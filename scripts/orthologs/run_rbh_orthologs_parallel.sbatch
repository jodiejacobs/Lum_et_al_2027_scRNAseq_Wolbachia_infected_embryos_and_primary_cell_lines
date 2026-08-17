#!/bin/bash
#SBATCH --job-name=dmel_dsim_rbh_par
#SBATCH --mail-user=jodiejacobs1021@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --output=rbh_orthologs_parallel_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64gb
#SBATCH --time=06:00:00
#SBATCH --partition=medium
# Optional: uncomment and set if you belong to multiple Slurm accounts
# #SBATCH --account=<GROUP>

# Parallel variant of run_rbh_orthologs.sbatch. Same pipeline, same
# node-local-scratch + timeout safeguards (see that script's header for
# the full backstory), but the dmel/dsim branches of each stage --
# copying inputs, gffread extraction, diamond makedb, and the two
# reciprocal blastp directions -- run concurrently as background jobs
# instead of one after another. They're independent until the final
# filter_rbh.py merge step, so there's no correctness reason to
# serialize them.
#
# Each diamond process gets half of $SLURM_CPUS_PER_TASK rather than
# all of it, since running one diamond job at a time with all 64
# threads doesn't actually use them efficiently on proteomes this small
# (~13-14k proteins/species) -- diamond's speedup flattens out well
# before 64 threads for a database this size, so two 32-thread jobs
# finishing side by side beats one 64-thread job with most of its
# threads idling on tiny per-thread chunks.
#
# NOTE: if the sequential version is also failing/hanging, parallelizing
# it will NOT fix that on its own -- it just does the same work with
# better resource utilization. If this version fails too, the .out log
# below will say which of the two concurrent branches failed and
# whether it was a timeout or a real error (each branch's full log is
# printed in full after both finish, prefixed with its label) -- that's
# the important part to send back for further debugging.
#
# Before running: same prerequisites as run_rbh_orthologs.sbatch --
# check WORKDIR, SCRIPT_DIR, DMEL_*/DSIM_* paths, SCRATCH_BASE, and the
# KALLISTO_ENV/DIAMOND_ENV mamba envs below.

set -euo pipefail

date; hostname; pwd
echo "=== memory ==="; free -h || true

WORKDIR=/private/groups/russelllab/jodie/scRNAseq/reference/orthologs
SCRIPT_DIR=/private/groups/russelllab/jodie/scRNAseq/Lum_et_al_2027_scRNAseq_Wolbachia_infected_embryos_and_primary_cell_lines/scripts/orthologs  # where gtf_id_map.py / filter_rbh.py live

DMEL_FASTA=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_melanogaster/dmel-all-chromosome-r6.68.fasta
DMEL_GTF=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_melanogaster/dmel-all-r6.68.gtf

DSIM_FASTA=/private/groups/russelllab/jodie/scRNAseq/reference/NCBI_ref/Dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna
DSIM_GTF=/private/groups/russelllab/jodie/scRNAseq/reference/NCBI_ref/Dsim/GCF_016746395.2_Prin_Dsim_3.1.gtf

KALLISTO_ENV=kallisto_bustools   # provides gffread
DIAMOND_ENV=diamond              # provides diamond

# Node-local scratch (see run_rbh_orthologs.sbatch's header for why).
SCRATCH_BASE="${TMPDIR:-/tmp}"
SCRATCH="$SCRATCH_BASE/dmel_dsim_rbh_par_${SLURM_JOB_ID:-$$}"
LOGDIR="$SCRATCH/logs"
mkdir -p "$LOGDIR"
echo "Using local scratch: $SCRATCH"
df -h "$SCRATCH" || true

mkdir -p "$WORKDIR"

# Half the allocated CPUs per concurrent diamond job (min 1).
THREADS_HALF=$(( ${SLURM_CPUS_PER_TASK:-2} / 2 ))
if [ "$THREADS_HALF" -lt 1 ]; then THREADS_HALF=1; fi
echo "Threads per parallel diamond job: $THREADS_HALF (of $SLURM_CPUS_PER_TASK total)"

# ---------------------------------------------------------------------
# bg_step <logfile> <timeout> <label> <command...>
#   Runs one command under `timeout`, with all its output (plus a
#   start/end/error banner) captured to <logfile>. Meant to be launched
#   with `&`; on failure it `exit`s the backgrounded subshell with the
#   command's exit code (124 for timeout), which `wait_all` below
#   picks up via `wait $pid`.
#
# launch <logfile> <timeout> <label> <command...>
#   Backgrounds a bg_step call and records its PID/logfile/label in the
#   PIDS/LOGS/LABELS arrays for wait_all to reconcile.
#
# wait_all
#   Waits for every job launched since the last wait_all, prints each
#   one's full log (so concurrent output doesn't interleave illegibly),
#   and aborts the script if any of them failed -- so a parallel stage
#   still fails loudly and specifically, same as the sequential version.
# ---------------------------------------------------------------------
declare -a PIDS=()
declare -a LOGS=()
declare -a LABELS=()

bg_step() {
    local logfile="$1" budget="$2" label="$3"
    shift 3
    {
        echo "--- [$label] starting ($(date))"
        if timeout "$budget" "$@"; then
            echo "--- [$label] done ($(date))"
        else
            local ec=$?
            if [ "$ec" -eq 124 ]; then
                echo "ERROR: [$label] timed out after $budget with no progress."
                echo "  Host: $(hostname)  Scratch: $SCRATCH"
                df -h "$SCRATCH" || true
            else
                echo "ERROR: [$label] failed with exit code $ec ($(date))"
            fi
            exit "$ec"
        fi
    } > "$logfile" 2>&1
}

launch() {
    local logfile="$1" budget="$2" label="$3"
    shift 3
    bg_step "$logfile" "$budget" "$label" "$@" &
    PIDS+=("$!")
    LOGS+=("$logfile")
    LABELS+=("$label")
}

wait_all() {
    local failed=0
    for i in "${!PIDS[@]}"; do
        local pid="${PIDS[$i]}" log="${LOGS[$i]}" label="${LABELS[$i]}"
        local ec=0
        wait "$pid" || ec=$?
        echo "===== [$label] log ====="
        cat "$log"
        if [ "$ec" -ne 0 ]; then
            echo "ERROR: [$label] failed (exit $ec)" >&2
            failed=1
        fi
    done
    PIDS=(); LOGS=(); LABELS=()
    if [ "$failed" -ne 0 ]; then
        echo "ERROR: one or more parallel steps above failed -- aborting." >&2
        exit 1
    fi
}

# ---- 0. copy inputs to local scratch, all 4 files concurrently ----
echo "=== copying inputs to local scratch (parallel) ==="
launch "$LOGDIR/copy_dmel_fasta.log" 15m "copy dmel fasta" cp "$DMEL_FASTA" "$SCRATCH/dmel_genome.fasta"
launch "$LOGDIR/copy_dmel_gtf.log"   5m  "copy dmel gtf"   cp "$DMEL_GTF"   "$SCRATCH/dmel.gtf"
launch "$LOGDIR/copy_dsim_fasta.log" 15m "copy dsim fasta" cp "$DSIM_FASTA" "$SCRATCH/dsim_genome.fasta"
launch "$LOGDIR/copy_dsim_gtf.log"   5m  "copy dsim gtf"   cp "$DSIM_GTF"   "$SCRATCH/dsim.gtf"
wait_all

cd "$SCRATCH"

source "$(dirname "$(dirname "$(which conda)")")/etc/profile.d/conda.sh"

# ---- 1. extract protein sequences from genome + GTF (gffread -y), dmel/dsim in parallel ----
conda activate "$KALLISTO_ENV"

# See run_rbh_orthologs.sbatch for why this filtering step exists
# (NCBI/Gnomon GTFs carry gene rows with a literal empty transcript_id
# "" that crashes gffread; harmless no-op on FlyBase's GTF).
grep -v 'transcript_id ""' dmel.gtf > dmel_filtered.gtf
grep -v 'transcript_id ""' dsim.gtf > dsim_filtered.gtf

echo "=== gffread extraction (parallel) ==="
launch "$LOGDIR/gffread_dmel.log" 20m "gffread dmel" gffread dmel_filtered.gtf -g dmel_genome.fasta -y dmel_proteins.raw.fasta
launch "$LOGDIR/gffread_dsim.log" 20m "gffread dsim" gffread dsim_filtered.gtf -g dsim_genome.fasta -y dsim_proteins.raw.fasta
wait_all

# gffread inserts non-standard placeholder characters (most commonly
# '.') for CDS regions it can't cleanly translate -- typically an
# incomplete/partial codon at the edge of a partial gene model, common
# in whole-genome NCBI/Gnomon annotations like dsim's. DIAMOND doesn't
# accept those and hard-fails ("Invalid character in sequence: '.'")
# on the whole file. Also re-wraps to a fixed line width as a defensive
# measure/diagnostic (see sanitize_protein_fasta.py's docstring for the
# 2026-08-16 dmel `diamond makedb` hang this was investigating).
echo "=== sanitizing protein FASTAs (parallel) ==="
launch "$LOGDIR/sanitize_dmel.log" 5m "sanitize dmel proteins" python "$SCRIPT_DIR/sanitize_protein_fasta.py" dmel_proteins.raw.fasta dmel_proteins.fasta
launch "$LOGDIR/sanitize_dsim.log" 5m "sanitize dsim proteins" python "$SCRIPT_DIR/sanitize_protein_fasta.py" dsim_proteins.raw.fasta dsim_proteins.fasta
wait_all

# Fail loudly rather than silently handing diamond an empty/tiny
# proteome -- see run_rbh_orthologs.sbatch for the incident this
# guards against.
n_dmel=$(grep -c '^>' dmel_proteins.fasta || true)
n_dsim=$(grep -c '^>' dsim_proteins.fasta || true)
echo "gffread extracted $n_dmel dmel proteins, $n_dsim dsim proteins"
if [ "$n_dmel" -lt 1000 ] || [ "$n_dsim" -lt 1000 ]; then
    echo "ERROR: expected >1000 proteins per species (D. melanogaster/" \
         "D. simulans each have ~13-14k protein-coding genes) -- got" \
         "$n_dmel dmel / $n_dsim dsim. Check the gffread logs above." >&2
    exit 1
fi

# ---- 2. {protein_id, transcript_id} -> gene_id maps, parsed straight
#         from each GTF (not from FASTA headers) -- fast, run in series ----
python "$SCRIPT_DIR/gtf_id_map.py" dmel.gtf > dmel_id_map.tsv
python "$SCRIPT_DIR/gtf_id_map.py" dsim.gtf > dsim_id_map.tsv

conda deactivate

# ---- 3. build diamond DBs, dmel/dsim in parallel ----
conda activate "$DIAMOND_ENV"

echo "=== diamond makedb (parallel, $THREADS_HALF threads each) ==="
launch "$LOGDIR/makedb_dmel.log" 30m "diamond makedb dmel" diamond makedb --in dmel_proteins.fasta -d dmel_db --threads "$THREADS_HALF"
launch "$LOGDIR/makedb_dsim.log" 30m "diamond makedb dsim" diamond makedb --in dsim_proteins.fasta -d dsim_db --threads "$THREADS_HALF"
wait_all

# ---- 4. reciprocal blastp, both directions in parallel ----
echo "=== diamond blastp (parallel, $THREADS_HALF threads each) ==="
launch "$LOGDIR/blastp_dsim_vs_dmel.log" 45m "diamond blastp dsim->dmel" \
    diamond blastp -q dsim_proteins.fasta -d dmel_db -o dsim_vs_dmel.tsv \
    --max-target-seqs 1 --evalue 1e-5 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --threads "$THREADS_HALF"

launch "$LOGDIR/blastp_dmel_vs_dsim.log" 45m "diamond blastp dmel->dsim" \
    diamond blastp -q dmel_proteins.fasta -d dsim_db -o dmel_vs_dsim.tsv \
    --max-target-seqs 1 --evalue 1e-5 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --threads "$THREADS_HALF"
wait_all

conda deactivate

# ---- 5. collapse to gene-level reciprocal best hits (fast, run in series) ----
# --name-a/--name-b: match the "Dsim"/"Dmel" column names config.yaml's
# ortholog_map and integrate_v2.py's load_ortholog_map() expect.
conda activate "$KALLISTO_ENV"   # just needs a stdlib-only python3

python "$SCRIPT_DIR/filter_rbh.py" \
    dsim_vs_dmel.tsv dmel_vs_dsim.tsv \
    dsim_id_map.tsv dmel_id_map.tsv \
    --name-a Dsim --name-b Dmel \
    > dmel_dsim_orthologs_rbh.tsv

conda deactivate

echo "=== copying results back to \$WORKDIR (NFS) ==="
cp dmel_dsim_orthologs_rbh.tsv dmel_proteins.fasta dsim_proteins.fasta \
   dmel_id_map.tsv dsim_id_map.tsv dsim_vs_dmel.tsv dmel_vs_dsim.tsv \
   "$WORKDIR/"

echo "Done. Genome-wide ortholog table: $WORKDIR/dmel_dsim_orthologs_rbh.tsv"
echo "Columns: Dsim(gene_id) Dmel(gene_id) pident evalue bitscore"

echo "Cleaning up local scratch: $SCRATCH"
rm -rf "$SCRATCH"

date
