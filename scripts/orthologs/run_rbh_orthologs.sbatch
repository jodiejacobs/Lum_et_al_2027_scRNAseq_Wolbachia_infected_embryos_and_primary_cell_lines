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
# Requires: get_id_map.py and filter_rbh.py in the same directory.
#
# Before running:
#   1. Set WORKDIR below.
#   2. Set DMEL_PEP_URL / DSIM_PEP_URL to the FlyBase translation FASTA
#      URLs for the releases you're actually using (check
#      https://flybase.org/downloads/bulkdata for current dmel release,
#      and https://flybase.org/genomes/Drosophila_simulans/ for dsim
#      r2.01 — do NOT use the r1.0 GLEAN release, it has no standard
#      gene models).
#   3. mamba env with diamond installed, e.g.:
#        mamba create -n diamond -c bioconda -c conda-forge diamond
#   4. Adjust #SBATCH resources to your cluster's partition limits.

set -euo pipefail

# Log host, start time, and working directory for troubleshooting
date; hostname; pwd

WORKDIR=/private/groups/russelllab/jodie/scRNAseq/reference/orthologs
SCRIPT_DIR=/private/groups/russelllab/jodie/scRNAseq/reference/orthologs  # where get_id_map.py / filter_rbh.py live

DMEL_PEP_URL="https://ftp.flybase.net/releases/current/precomputed_files/genes/fasta/dmel-all-translation-current.fasta.gz"
DSIM_PEP_URL="https://ftp.flybase.net/genomes/Drosophila_simulans/dsim_r2.01_FB2015_01/fasta/dsim-all-translation-r2.01.fasta.gz"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

source "$(mamba info --base)/etc/profile.d/conda.sh"
mamba activate diamond

# ---- 1. download proteomes ----
if [ ! -f dmel_proteins.fasta ]; then
    wget -O dmel.fasta.gz "$DMEL_PEP_URL"
    gunzip -c dmel.fasta.gz > dmel_proteins.fasta
fi
if [ ! -f dsim_proteins.fasta ]; then
    wget -O dsim.fasta.gz "$DSIM_PEP_URL"
    gunzip -c dsim.fasta.gz > dsim_proteins.fasta
fi

# ---- 2. protein_id -> gene_id maps ----
python "$SCRIPT_DIR/get_id_map.py" dmel_proteins.fasta > dmel_id_map.tsv
python "$SCRIPT_DIR/get_id_map.py" dsim_proteins.fasta > dsim_id_map.tsv

# ---- 3. build diamond DBs ----
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

# ---- 5. collapse to gene-level reciprocal best hits ----
python "$SCRIPT_DIR/filter_rbh.py" \
    dsim_vs_dmel.tsv dmel_vs_dsim.tsv \
    dsim_id_map.tsv dmel_id_map.tsv \
    > dmel_dsim_orthologs_rbh.tsv

echo "Done. Genome-wide ortholog table: $WORKDIR/dmel_dsim_orthologs_rbh.tsv"
echo "Columns: geneA(dsim) geneB(dmel) pident evalue bitscore"

date
