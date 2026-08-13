#!/bin/bash
#SBATCH --job-name=dmel_dsim_rbh
#SBATCH --mail-user=jodiejacobs1021@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --output=rbh_orthologs_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
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
# Before running:
#   1. Set WORKDIR below.
#   2. Set DMEL_FASTA/DMEL_GTF/DSIM_FASTA/DSIM_GTF to your genome FASTA +
#      GTF pairs (whole-genome/chromosome FASTA, NOT a translation/
#      protein FASTA -- proteins are extracted from these with gffread).
#   3. mamba envs:
#        - KALLISTO_ENV below (this repo's existing kallisto_bustools
#          env) must have gffread on PATH.
#        - DIAMOND_ENV below, e.g.:
#            mamba create -n diamond -c bioconda -c conda-forge diamond
#   4. Adjust #SBATCH resources to your cluster's partition limits.

set -euo pipefail

# Log host, start time, and working directory for troubleshooting
date; hostname; pwd

WORKDIR=/private/groups/russelllab/jodie/scRNAseq/reference/orthologs
SCRIPT_DIR=/private/groups/russelllab/jodie/scRNAseq/Lum_et_al_2027_scRNAseq_Wolbachia_infected_embryos_and_primary_cell_lines/scripts/orthologs  # where gtf_id_map.py / filter_rbh.py live

DMEL_FASTA=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_melanogaster/dmel-all-chromosome-r6.68.fasta
DMEL_GTF=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_melanogaster/dmel-all-r6.68.gtf

DSIM_FASTA=/private/groups/russelllab/jodie/scRNAseq/reference/NCBI_ref/Dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna
DSIM_GTF=/private/groups/russelllab/jodie/scRNAseq/reference/NCBI_ref/Dsim/GCF_016746395.2_Prin_Dsim_3.1.gtf

KALLISTO_ENV=kallisto_bustools   # provides gffread
DIAMOND_ENV=diamond              # provides diamond

mkdir -p "$WORKDIR"
cd "$WORKDIR"

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
grep -v 'transcript_id ""' "$DMEL_GTF" > dmel_filtered.gtf
grep -v 'transcript_id ""' "$DSIM_GTF" > dsim_filtered.gtf

if [ ! -f dmel_proteins.fasta ]; then
    gffread dmel_filtered.gtf -g "$DMEL_FASTA" -y dmel_proteins.fasta
fi
if [ ! -f dsim_proteins.fasta ]; then
    gffread dsim_filtered.gtf -g "$DSIM_FASTA" -y dsim_proteins.fasta
fi

# ---- 2. transcript_id -> gene_id maps, parsed straight from each GTF ----
# (not from FASTA headers -- the NCBI GTF has no FlyBase-style
# parent=FBgn convention to parse, so gtf_id_map.py reads gene_id/
# transcript_id out of the GTF attribute column directly, which works
# for both the FlyBase dmel GTF and the NCBI/Gnomon dsim GTF.)
python "$SCRIPT_DIR/gtf_id_map.py" "$DMEL_GTF" > dmel_id_map.tsv
python "$SCRIPT_DIR/gtf_id_map.py" "$DSIM_GTF" > dsim_id_map.tsv

conda deactivate

# ---- 3. build diamond DBs ----
conda activate "$DIAMOND_ENV"

diamond makedb --in dmel_proteins.fasta -d dmel_db --threads "$SLURM_CPUS_PER_TASK"
diamond makedb --in dsim_proteins.fasta -d dsim_db --threads "$SLURM_CPUS_PER_TASK"

# ---- 4. reciprocal blastp, best hit per query only ----
diamond blastp -q dsim_proteins.fasta -d dmel_db -o dsim_vs_dmel.tsv \
    --max-target-seqs 1 --evalue 1e-5 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --threads "$SLURM_CPUS_PER_TASK"

diamond blastp -q dmel_proteins.fasta -d dsim_db -o dmel_vs_dsim.tsv \
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

echo "Done. Genome-wide ortholog table: $WORKDIR/dmel_dsim_orthologs_rbh.tsv"
echo "Columns: Dsim(gene_id) Dmel(gene_id) pident evalue bitscore"

date
