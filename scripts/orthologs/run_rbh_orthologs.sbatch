#!/bin/bash
#SBATCH --job-name=dmel_dsim_rbh
#SBATCH --mail-user=jodiejacobs1021@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --output=rbh_orthologs_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64gb
#SBATCH --time=06:00:00
#SBATCH --partition=medium
# Optional: uncomment and set if you belong to multiple Slurm accounts
# #SBATCH --account=<GROUP>

# Genome-wide D. melanogaster <-> D. simulans ortholog table via
# reciprocal best BLAST hit (DIAMOND), since FlyBase no longer ships
# a bulk orthologs file (OrthoDB calls moved to a live API in 2022)
# and DIOPT doesn't cover D. simulans.
#
# Built directly from each species' genome FASTA + GTF (protein
# sequences extracted locally with gffread) rather than a pre-built
# FlyBase translation FASTA. This is required now that the D. simulans
# annotation comes from NCBI RefSeq (Prin_Dsim_3.1, GCF_016746395.2,
# Gnomon gene models with "LOC#####"/accession-style gene IDs) instead
# of FlyBase, which no longer maintains D. simulans sequence/annotation
# at all -- there is no FlyBase-style translation FASTA to download for
# it anymore. D. melanogaster still comes from FlyBase, but is handled
# the same way here for a consistent gffread-based extraction step on
# both sides.
#
# Requires: gtf_id_map.py and filter_rbh.py in the same directory.
#
# RUNTIME NOTE (learned the hard way on 2026-08-13): everything runs on
# node-local scratch ($SCRATCH_BASE below), not on the NFS-mounted
# $WORKDIR. DIAMOND's makedb/blastp hang at 0% CPU for hours when their
# input/output/temp files live on NFS -- a known upstream issue, not
# specific to this pipeline or cluster:
#   https://github.com/bbuchfink/diamond/issues/319
#   https://github.com/bbuchfink/diamond/issues/930
#   https://github.com/bbuchfink/diamond/issues/185
# Inputs are copied to local scratch once at the start, all gffread/
# diamond work happens there, and only the final results are copied
# back to $WORKDIR on NFS at the end. Every long-running external
# command is also wrapped in `timeout` so a *future* hang (NFS glitch,
# flaky node, whatever) fails within a bounded time with a clear
# message instead of silently eating the whole #SBATCH --time budget --
# that's what happened before: a hung `diamond makedb` sat at 0% CPU /
# 0 MB for 2.5+ hours with no indication anything was wrong until
# someone happened to check `seff`.
#
# Before running:
#   1. Set WORKDIR below (final NFS destination for results).
#   2. Set DMEL_FASTA/DMEL_GTF/DSIM_FASTA/DSIM_GTF to your genome FASTA +
#      GTF pairs (whole-genome/chromosome FASTA, NOT a translation/
#      protein FASTA -- proteins are extracted from these with gffread).
#   3. Verify SCRATCH_BASE below is actually node-local fast disk on
#      your cluster with enough free space (a few hundred MB is plenty
#      for these genome sizes) -- run `df -h "$TMPDIR"` (or /tmp) on a
#      compute node once to check. Override SCRATCH_BASE if your
#      cluster uses a different local-scratch convention (e.g.
#      /scratch/local, /localscratch).
#   4. mamba envs:
#        - KALLISTO_ENV below (this repo's existing kallisto_bustools
#          env) must have gffread on PATH.
#        - DIAMOND_ENV below, e.g.:
#            mamba create -n diamond -c bioconda -c conda-forge diamond
#   5. Adjust #SBATCH resources to your cluster's partition limits.

set -euo pipefail

# Log host, start time, and resources for troubleshooting
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

# Node-local scratch -- see RUNTIME NOTE above. $TMPDIR is typically set
# per-job to node-local storage by Slurm; falls back to /tmp otherwise.
SCRATCH_BASE="${TMPDIR:-/tmp}"
SCRATCH="$SCRATCH_BASE/dmel_dsim_rbh_${SLURM_JOB_ID:-$$}"
mkdir -p "$SCRATCH"
echo "Using local scratch: $SCRATCH"
df -h "$SCRATCH" || true

mkdir -p "$WORKDIR"

# run_step <timeout> <label> <command...>: runs a command under `timeout`,
# printing start/end timestamps, and on failure prints which step failed
# and whether it was a timeout (exit 124) vs a real error before exiting.
run_step() {
    local budget="$1" label="$2"
    shift 2
    echo "--- [$label] starting ($(date))"
    if timeout "$budget" "$@"; then
        echo "--- [$label] done ($(date))"
    else
        local ec=$?
        if [ "$ec" -eq 124 ]; then
            echo "ERROR: [$label] timed out after $budget with no progress." >&2
            echo "  Host: $(hostname)  Scratch: $SCRATCH" >&2
            df -h "$SCRATCH" >&2 || true
            echo "  If this keeps happening: check node health, local disk" \
                 "space above, and whether \$SCRATCH_BASE is genuinely" \
                 "node-local (not another network mount) on this cluster." >&2
        else
            echo "ERROR: [$label] failed with exit code $ec ($(date))" >&2
        fi
        exit "$ec"
    fi
}

echo "=== copying inputs to local scratch ==="
run_step 15m "copy dmel fasta" cp "$DMEL_FASTA" "$SCRATCH/dmel_genome.fasta"
run_step 5m  "copy dmel gtf"   cp "$DMEL_GTF"   "$SCRATCH/dmel.gtf"
run_step 15m "copy dsim fasta" cp "$DSIM_FASTA" "$SCRATCH/dsim_genome.fasta"
run_step 5m  "copy dsim gtf"   cp "$DSIM_GTF"   "$SCRATCH/dsim.gtf"

cd "$SCRATCH"

source "$(dirname "$(dirname "$(which conda)")")/etc/profile.d/conda.sh"

# ---- 1. extract protein sequences from genome + GTF (gffread -y) ----
conda activate "$KALLISTO_ENV"

# NCBI RefSeq/Gnomon GTFs (like the dsim GCF_016746395.2 one) carry
# gene-level rows with a literal, empty transcript_id "" attribute
# (rather than omitting transcript_id entirely, as FlyBase GTFs do).
# gffread's parser tries to use that as the record's ID, gets an empty
# string, and dies with "Error: no valid ID found for GFF record" --
# a known gffread/NCBI-GTF issue, not specific to this pipeline. Strip
# those rows before handing the GTF to gffread (harmless no-op on GTFs
# that don't have any, like FlyBase's, so applied to both defensively).
# gtf_id_map.py doesn't need this filtering -- it already skips rows
# with an empty transcript_id on its own.
grep -v 'transcript_id ""' dmel.gtf > dmel_filtered.gtf
grep -v 'transcript_id ""' dsim.gtf > dsim_filtered.gtf

run_step 20m "gffread dmel" gffread dmel_filtered.gtf -g dmel_genome.fasta -y dmel_proteins.fasta
run_step 20m "gffread dsim" gffread dsim_filtered.gtf -g dsim_genome.fasta -y dsim_proteins.fasta

# Fail loudly here rather than silently handing diamond an empty/tiny
# proteome -- cheap sanity check given how easy this is to miss (this
# is exactly what happened on 2026-08-13: a previous run's gffread call
# died partway through, left a 0-byte dsim_proteins.fasta, and a later
# run's -f-only existence check let it slip through unnoticed).
n_dmel=$(grep -c '^>' dmel_proteins.fasta || true)
n_dsim=$(grep -c '^>' dsim_proteins.fasta || true)
echo "gffread extracted $n_dmel dmel proteins, $n_dsim dsim proteins"
if [ "$n_dmel" -lt 1000 ] || [ "$n_dsim" -lt 1000 ]; then
    echo "ERROR: expected >1000 proteins per species (D. melanogaster/" \
         "D. simulans each have ~13-14k protein-coding genes) -- got" \
         "$n_dmel dmel / $n_dsim dsim. Check the gffread step above." >&2
    exit 1
fi

# ---- 2. {protein_id, transcript_id} -> gene_id maps, parsed straight
#         from each GTF (not from FASTA headers) ----
python "$SCRIPT_DIR/gtf_id_map.py" dmel.gtf > dmel_id_map.tsv
python "$SCRIPT_DIR/gtf_id_map.py" dsim.gtf > dsim_id_map.tsv

conda deactivate

# ---- 3. build diamond DBs ----
conda activate "$DIAMOND_ENV"

run_step 30m "diamond makedb dmel" diamond makedb --in dmel_proteins.fasta -d dmel_db --threads "$SLURM_CPUS_PER_TASK"
run_step 30m "diamond makedb dsim" diamond makedb --in dsim_proteins.fasta -d dsim_db --threads "$SLURM_CPUS_PER_TASK"

# ---- 4. reciprocal blastp, best hit per query only ----
run_step 45m "diamond blastp dsim->dmel" diamond blastp -q dsim_proteins.fasta -d dmel_db -o dsim_vs_dmel.tsv \
    --max-target-seqs 1 --evalue 1e-5 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --threads "$SLURM_CPUS_PER_TASK"

run_step 45m "diamond blastp dmel->dsim" diamond blastp -q dmel_proteins.fasta -d dsim_db -o dmel_vs_dsim.tsv \
    --max-target-seqs 1 --evalue 1e-5 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --threads "$SLURM_CPUS_PER_TASK"

conda deactivate

# ---- 5. collapse to gene-level reciprocal best hits ----
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
