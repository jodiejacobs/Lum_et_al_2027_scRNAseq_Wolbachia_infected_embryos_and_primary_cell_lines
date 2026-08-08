#!/usr/bin/env bash
# Build kallisto|bustools transcriptomes for dmel, dsim, dwil with gene
# SYMBOLS as the t2g "gene" column -> becomes adata.var_names after
# `kb count` / bustools loading.
#
# WHY THIS EXISTS:
#   dmel-all-r6.68.gtf and dwil-all-r1.04.gtf are both standard FlyBase GTF
#   (gene_id "FBgn...", gene_symbol "...") -- already the same format,
#   nothing to convert. dsim-all-r1.0.gff is NOT the same format: it's a
#   legacy 2007 GLEANR alignment annotation (match/match_part features)
#   with no native dsim FBgn id, only a dmel_ortholog= tag. Rather than
#   hack that into a fake GTF, this pulls the current official FlyBase
#   dsim annotation (r2.02), which is native-FBgn, standard gene/mRNA/exon
#   GTF -- same shape as dmel/dwil. Genome fasta + GTF are pulled from the
#   same release (dsim_r2.02_FB2017_04) so coordinates are guaranteed to
#   match; FlyBase hasn't republished dsim gff/gtf since (the r2.02
#   assembly/annotation itself hasn't changed, only that later release
#   folders stopped mirroring gff/gtf on this FTP).
#
#   Every GTF then goes through swap_gene_id_to_symbol.py, which rewrites
#   gene_id to hold the gene_symbol value (stripping the "Dsim\"/"Dwil\"
#   species-prefix FlyBase puts on non-melanogaster symbols, so orthologs
#   share the same base symbol across species). kb ref keys t2g off
#   gene_id, so the transcriptome/index/t2g built from the swapped GTF has
#   symbols, not FBgn IDs, as the gene column.
#
#   NOTE: gene symbols are not guaranteed unique the way FBgn IDs are.
#   swap_gene_id_to_symbol.py prints a warning listing any symbol shared by
#   more than one gene -- those will collide in adata.var_names and need
#   `.var_names_make_unique()` downstream.
#
# Run on the server:
#   mamba activate kallisto_bustools   # kallisto + bustools + gffread + kb-python

set -euo pipefail

# ---- EDIT THESE ----
GENOME_ROOT=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes
SCRIPT_DIR=/private/groups/russelllab/jodie/scRNAseq/Lum_et_al_2027_scRNAseq_Wolbachia_infected_embryos_and_primary_cell_lines/snakemake_scripts/reference
# ---------------------

SWAP_SCRIPT="$SCRIPT_DIR/swap_gene_id_to_symbol.py"
RRNA_SCRIPT="$SCRIPT_DIR/find_rrna_genes.py"

DMEL_DIR="$GENOME_ROOT/Drosophila_melanogaster"
DSIM_DIR="$GENOME_ROOT/Drosophila_simulans"
DWIL_DIR="$GENOME_ROOT/Drosophila_willistoni"
mkdir -p "$DMEL_DIR" "$DSIM_DIR" "$DWIL_DIR"

### 1. dsim: replace legacy r1.0 GLEANR gff with current FlyBase r2.02 GTF ###
cd "$DSIM_DIR"
wget -N https://s3ftp.flybase.org/genomes/Drosophila_simulans/dsim_r2.02_FB2017_04/fasta/dsim-all-chromosome-r2.02.fasta.gz
wget -N https://s3ftp.flybase.org/genomes/Drosophila_simulans/dsim_r2.02_FB2017_04/gtf/dsim-all-r2.02.gtf.gz
gunzip -kf dsim-all-chromosome-r2.02.fasta.gz dsim-all-r2.02.gtf.gz
echo "dsim sanity check (should show native dsim FBgn, not dmel_ortholog=):"
head -3 dsim-all-r2.02.gtf

### 2. dwil: already FlyBase GTF. Current release is r1.05 -- optional,
###    your existing r1.04 is already in the right format.
# cd "$DWIL_DIR"
# wget -N https://s3ftp.flybase.org/genomes/Drosophila_willistoni/dwil_r1.05_FB2016_05/fasta/dwil-all-chromosome-r1.05.fasta.gz
# wget -N https://s3ftp.flybase.org/genomes/Drosophila_willistoni/dwil_r1.05_FB2016_05/gtf/dwil-all-r1.05.gtf.gz
# gunzip -kf dwil-all-chromosome-r1.05.fasta.gz dwil-all-r1.05.gtf.gz

### 3. Build kallisto|bustools transcriptome per species (symbol-keyed) ###
build_ref () {
  local name=$1 fasta=$2 gtf=$3 outdir=$4
  mkdir -p "$outdir"
  echo "=== Building $name transcriptome (gene symbols) ==="

  local symbol_gtf="$outdir/${name}.symbol.gtf"
  python3 "$SWAP_SCRIPT" "$gtf" "$symbol_gtf"

  gffread "$symbol_gtf" -g "$fasta" -w "$outdir/transcripts.fa"
  kb ref \
    -i "$outdir/index.idx" \
    -g "$outdir/t2g.txt" \
    -f1 "$outdir/transcripts.fa" \
    "$fasta" "$symbol_gtf"

  # rRNA genes for this genome, by symbol -- feeds the titer script instead
  # of a hardcoded gene list
  python3 "$RRNA_SCRIPT" "$symbol_gtf" > "$outdir/host_rrna_genes.txt"
}

build_ref dmel \
  "$DMEL_DIR/dmel-all-aligned-r6.68.fasta" \
  "$DMEL_DIR/dmel-all-r6.68.gtf" \
  "$DMEL_DIR/kallisto_ref"

build_ref dsim \
  "$DSIM_DIR/dsim-all-chromosome-r2.02.fasta" \
  "$DSIM_DIR/dsim-all-r2.02.gtf" \
  "$DSIM_DIR/kallisto_ref"

build_ref dwil \
  "$DWIL_DIR/dwil-all-chromosome-r1.04.fasta" \
  "$DWIL_DIR/dwil-all-r1.04.gtf" \
  "$DWIL_DIR/kallisto_ref"

### 4. Confirm the gene column is now symbols, and show the rRNA genes found ###
for sp in Drosophila_melanogaster Drosophila_simulans Drosophila_willistoni; do
  echo "--- $sp t2g.txt sample ---"
  head -3 "$GENOME_ROOT/$sp/kallisto_ref/t2g.txt"
  echo "--- $sp host_rrna_genes.txt ---"
  cat "$GENOME_ROOT/$sp/kallisto_ref/host_rrna_genes.txt"
done

# For a combined host+Wolbachia reference (e.g. Dmel_wMel), also run
# find_rrna_genes.py against the Wolbachia GTF on its own before
# concatenating it with the host GTF, e.g.:
#   python3 "$RRNA_SCRIPT" wMel.gtf > wMel_kallisto_ref/symbiont_rrna_genes.txt
# Wolbachia gene IDs are NCBI locus tags (e.g. GQX67_05945), untouched by
# swap_gene_id_to_symbol.py -- that script only rewrites host FlyBase GTFs.
