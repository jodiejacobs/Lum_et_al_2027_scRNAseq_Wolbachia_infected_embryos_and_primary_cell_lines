#!/usr/bin/env bash
# Reciprocal-best-BLAST-hit (RBH) D. simulans <-> D. melanogaster protein
# orthologs -- a sequence-based fallback/cross-check for
# find_dmel_dsim_orthologs.py, which calls orthologs from FlyBase's own
# gene_symbol projection instead. Use this for:
#   - dsim genes find_dmel_dsim_orthologs.py left unmatched (no dmel gene
#     shares that symbol) -- pass its --unmatched-out file as RESTRICT_IDS
#   - independently verifying the symbol-based calls (leave RESTRICT_IDS
#     unset to run genome-wide instead)
#
# Needs blast (blastp/makeblastdb) on PATH, in addition to gffread. Neither
# of this repo's mamba envs currently has blast; add it to kallisto_bustools
# (which already has gffread) or make a dedicated env:
#   mamba install -n kallisto_bustools -c bioconda blast
#   # or:
#   mamba create -n blast -c bioconda -c conda-forge blast gffread
#
# Run on the server:
#   mamba activate kallisto_bustools   # (or your blast env)
#   bash run_reciprocal_blast_orthologs.sh

set -euo pipefail

# ---- EDIT THESE ----
DMEL_FASTA=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_melanogaster/dmel-all-aligned-r6.68.fasta
DMEL_GTF=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_melanogaster/dmel-all-r6.68.gtf
DSIM_FASTA=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_simulans/dsim-all-chromosome-r2.02.fasta
DSIM_GTF=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_simulans/dsim-all-r2.02.gtf
OUTDIR=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/dsim_dmel_orthologs
THREADS=16
SCRIPT_DIR=/private/groups/russelllab/jodie/scRNAseq/Lum_et_al_2027_scRNAseq_Wolbachia_infected_embryos_and_primary_cell_lines/snakemake_scripts/reference
# Optional: dsim gene_ids to restrict RBH output to (e.g. the
# --unmatched-out file from find_dmel_dsim_orthologs.py). Leave empty to
# run genome-wide.
RESTRICT_IDS=""
# ---------------------

mkdir -p "$OUTDIR"

echo "=== Extracting protein sequences (gffread -y, one per transcript with a CDS) ==="
gffread "$DMEL_GTF" -g "$DMEL_FASTA" -y "$OUTDIR/dmel_proteins.fa"
gffread "$DSIM_GTF" -g "$DSIM_FASTA" -y "$OUTDIR/dsim_proteins.fa"

echo "=== Building BLAST protein databases ==="
makeblastdb -in "$OUTDIR/dmel_proteins.fa" -dbtype prot -out "$OUTDIR/dmel_db"
makeblastdb -in "$OUTDIR/dsim_proteins.fa" -dbtype prot -out "$OUTDIR/dsim_db"

echo "=== dsim -> dmel BLASTP ==="
blastp -query "$OUTDIR/dsim_proteins.fa" -db "$OUTDIR/dmel_db" \
  -out "$OUTDIR/dsim_vs_dmel.tsv" -outfmt 6 -max_target_seqs 5 -evalue 1e-5 \
  -num_threads "$THREADS"

echo "=== dmel -> dsim BLASTP ==="
blastp -query "$OUTDIR/dmel_proteins.fa" -db "$OUTDIR/dsim_db" \
  -out "$OUTDIR/dmel_vs_dsim.tsv" -outfmt 6 -max_target_seqs 5 -evalue 1e-5 \
  -num_threads "$THREADS"

echo "=== Calling reciprocal best hits (gene-level) ==="
restrict_flag=()
if [[ -n "$RESTRICT_IDS" ]]; then
  restrict_flag=(--restrict-ids "$RESTRICT_IDS")
fi
python3 "$SCRIPT_DIR/reciprocal_best_hits.py" \
  --forward "$OUTDIR/dsim_vs_dmel.tsv" \
  --reverse "$OUTDIR/dmel_vs_dsim.tsv" \
  --dsim-gtf "$DSIM_GTF" --dmel-gtf "$DMEL_GTF" \
  -o "$OUTDIR/dsim_dmel_rbh_orthologs.tsv" \
  "${restrict_flag[@]}"

echo "Done -- $OUTDIR/dsim_dmel_rbh_orthologs.tsv"
echo "To fold these into the symbol-based table, use dsim_gene_id from"
echo "dsim_to_dmel_orthologs.tsv rows with match_status == unmatched, and"
echo "join in dmel_gene_id from dsim_dmel_rbh_orthologs.tsv (match_status ="
echo "'rbh_blast')."
