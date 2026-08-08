#!/usr/bin/env bash
# Build combined host + Wolbachia strain + 16S kallisto|bustools references,
# one per genome key used in samples.csv / config.yaml's genome_components
# (Dmel_wMel, Dmel_wRi_Riv84, Dmel_wWil, Dsim_wMel, Dsim_wRi_M23,
# Dsim_wRi_Riv84, Dwil_wWil).
#
# Each combined reference = that host's genome + symbol-keyed GTF (from
# build_dmel_dsim_dwil_transcriptomes.sh -- regenerated here if missing)
#                          + that Wolbachia strain's own genome + GTF
#                            (untouched, NCBI locus-tag gene_ids)
#                          + the shared 16S reference
#                            (highly_represented_16S_filtered.renamed.fasta /
#                            16S_fixed.gtf)
# concatenated and run through gffread + kb ref together, so a single
# kallisto index/t2g covers host transcripts (by symbol), symbiont
# transcripts (by locus tag), and 16S sequences at once.
#
# Also writes each Wolbachia strain's own rrna_genes.txt (via
# find_rrna_genes.py against that strain's *own* GTF, not the combined one)
# next to its source files -- config.yaml's symbiont_rrna_genes should point
# there.
#
# Run AFTER build_dmel_dsim_dwil_transcriptomes.sh.
# Run on the server:
#   mamba activate kallisto_bustools

set -euo pipefail

# ---- EDIT THESE ----
FLYBASE_ROOT=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes
REF_ROOT=/private/groups/russelllab/jodie/scRNAseq/reference
COMBINED_ROOT=/private/groups/russelllab/jodie/scRNAseq/Lum_et_al_2027_scRNAseq_Wolbachia_infected_embryos_and_primary_cell_lines/reference
SCRIPT_DIR=/private/groups/russelllab/jodie/scRNAseq/Lum_et_al_2027_scRNAseq_Wolbachia_infected_embryos_and_primary_cell_lines/snakemake_scripts/reference
# ---------------------

SWAP_SCRIPT="$SCRIPT_DIR/swap_gene_id_to_symbol.py"
RRNA_SCRIPT="$SCRIPT_DIR/find_rrna_genes.py"

SIXTEEN_S_FASTA="$REF_ROOT/highly_represented_16S/highly_represented_16S_filtered.renamed.fasta"
SIXTEEN_S_GTF="$REF_ROOT/16S_fixed.gtf"

declare -A HOST_FASTA=(
  [Dmel]="$FLYBASE_ROOT/Drosophila_melanogaster/dmel-all-aligned-r6.68.fasta"
  [Dsim]="$FLYBASE_ROOT/Drosophila_simulans/dsim-all-chromosome-r2.02.fasta"
  [Dwil]="$FLYBASE_ROOT/Drosophila_willistoni/dwil-all-chromosome-r1.04.fasta"
)
declare -A HOST_GTF=(
  [Dmel]="$FLYBASE_ROOT/Drosophila_melanogaster/dmel-all-r6.68.gtf"
  [Dsim]="$FLYBASE_ROOT/Drosophila_simulans/dsim-all-r2.02.gtf"
  [Dwil]="$FLYBASE_ROOT/Drosophila_willistoni/dwil-all-r1.04.gtf"
)
declare -A HOST_SYMBOL_GTF=(
  [Dmel]="$FLYBASE_ROOT/Drosophila_melanogaster/kallisto_ref/dmel.symbol.gtf"
  [Dsim]="$FLYBASE_ROOT/Drosophila_simulans/kallisto_ref/dsim.symbol.gtf"
  [Dwil]="$FLYBASE_ROOT/Drosophila_willistoni/kallisto_ref/dwil.symbol.gtf"
)

declare -A WOLBACHIA_FASTA=(
  [wMel]="$REF_ROOT/wMel_GCF_016584425.1/GCF_016584425.1_ASM1658442v1_genomic.fna"
  [wRi_Riv84]="$REF_ROOT/wRi_Riv84_GCF_000022285.1/GCF_000022285.1_ASM2228v1_genomic.fna"
  [wRi_M23]="$REF_ROOT/wRi_M23_GCA_979474595.1/20260402_wRi_M23_pilon.fasta"
  [wWil]="$REF_ROOT/wWil_GCF_040084705.1/wWil_GCF_040084705.fna"
)
declare -A WOLBACHIA_GTF=(
  [wMel]="$REF_ROOT/wMel_GCF_016584425.1/wMel_GCF_016584425.1.gtf"
  [wRi_Riv84]="$REF_ROOT/wRi_Riv84_GCF_000022285.1/wRi_Riv84_GCF_000022285.1.gtf"
  [wRi_M23]="$REF_ROOT/wRi_M23_GCA_979474595.1/wRi_M23_v3.gtf"
  [wWil]="$REF_ROOT/wWil_GCF_040084705.1/wWil_GCF_040084705.1.gtf"
)

# genome key (as used in samples.csv / config.yaml) -> host key, symbiont key
declare -A GENOME_HOST=(
  [Dmel_wMel]=Dmel      [Dmel_wRi_Riv84]=Dmel  [Dmel_wWil]=Dmel
  [Dsim_wMel]=Dsim      [Dsim_wRi_M23]=Dsim    [Dsim_wRi_Riv84]=Dsim
  [Dwil_wWil]=Dwil
)
declare -A GENOME_SYMBIONT=(
  [Dmel_wMel]=wMel      [Dmel_wRi_Riv84]=wRi_Riv84  [Dmel_wWil]=wWil
  [Dsim_wMel]=wMel      [Dsim_wRi_M23]=wRi_M23      [Dsim_wRi_Riv84]=wRi_Riv84
  [Dwil_wWil]=wWil
)
# output dirs -- match config.yaml's genome-key -> path entries exactly
declare -A GENOME_OUTDIR=(
  [Dmel_wMel]="$COMBINED_ROOT/Drosophila_melanogaster_wMel_combined_16S"
  [Dmel_wRi_Riv84]="$COMBINED_ROOT/Drosophila_melanogaster_wRi_Riv84_combined_16S"
  [Dmel_wWil]="$COMBINED_ROOT/Drosophila_melanogaster_wWil_combined_16S"
  [Dsim_wMel]="$COMBINED_ROOT/Drosophila_simulans_wMel_combined_16S"
  [Dsim_wRi_M23]="$COMBINED_ROOT/Drosophila_simulans_wRi_M23_combined_16S"
  [Dsim_wRi_Riv84]="$COMBINED_ROOT/Drosophila_simulans_wRi_Riv84_combined_16S"
  [Dwil_wWil]="$COMBINED_ROOT/Drosophila_willistoni_wWil_combined_16S"
)

# Concatenate files with a guaranteed newline between them, so a missing
# trailing newline in one source file can't glue its last line to the next
# file's first line.
concat_safe () {
  local out=$1; shift
  awk 'FNR==1 && NR!=1 {print ""} {print}' "$@" > "$out"
}

# Warn (don't fail) if any GTF seqname isn't present in the combined fasta --
# gffread silently drops features on unmatched seqnames, which would
# otherwise show up downstream as "fewer transcripts than expected" with no
# obvious cause.
check_seqnames () {
  local fasta=$1 gtf=$2
  local missing
  missing=$(comm -23 \
    <(grep -v '^$' "$gtf" | cut -f1 | sort -u) \
    <(grep '^>' "$fasta" | sed 's/^>//; s/[[:space:]].*//' | sort -u) \
    || true)
  if [[ -n "$missing" ]]; then
    echo "  WARNING: GTF seqname(s) not found in combined fasta -- gffread will silently skip these features:"
    echo "$missing" | sed 's/^/    /'
  fi
}

# One rrna_genes.txt per Wolbachia strain (against its own GTF, not the
# combined one), written next to its source files.
for symbiont in "${!WOLBACHIA_GTF[@]}"; do
  out="$(dirname "${WOLBACHIA_GTF[$symbiont]}")/rrna_genes.txt"
  echo "=== rRNA genes: $symbiont -> $out ==="
  python3 "$RRNA_SCRIPT" "${WOLBACHIA_GTF[$symbiont]}" > "$out"
done

build_combined () {
  local genome_key=$1
  local host=${GENOME_HOST[$genome_key]}
  local symbiont=${GENOME_SYMBIONT[$genome_key]}
  local outdir=${GENOME_OUTDIR[$genome_key]}

  echo "=== Building $genome_key ($host + $symbiont + 16S) -> $outdir ==="
  mkdir -p "$outdir"

  # Make sure the host's symbol-keyed GTF exists; regenerate if this is
  # being run before/without build_dmel_dsim_dwil_transcriptomes.sh.
  if [[ ! -s "${HOST_SYMBOL_GTF[$host]}" ]]; then
    echo "  ${HOST_SYMBOL_GTF[$host]} missing -- generating from ${HOST_GTF[$host]}"
    mkdir -p "$(dirname "${HOST_SYMBOL_GTF[$host]}")"
    python3 "$SWAP_SCRIPT" "${HOST_GTF[$host]}" "${HOST_SYMBOL_GTF[$host]}"
  fi

  concat_safe "$outdir/combined.fasta" \
    "${HOST_FASTA[$host]}" "${WOLBACHIA_FASTA[$symbiont]}" "$SIXTEEN_S_FASTA"
  concat_safe "$outdir/combined.gtf" \
    "${HOST_SYMBOL_GTF[$host]}" "${WOLBACHIA_GTF[$symbiont]}" "$SIXTEEN_S_GTF"

  check_seqnames "$outdir/combined.fasta" "$outdir/combined.gtf"

  gffread "$outdir/combined.gtf" -g "$outdir/combined.fasta" -w "$outdir/transcripts.fa"
  kb ref \
    -i "$outdir/index.idx" \
    -g "$outdir/t2g.txt" \
    -f1 "$outdir/transcripts.fa" \
    "$outdir/combined.fasta" "$outdir/combined.gtf"
}

for genome_key in "${!GENOME_OUTDIR[@]}"; do
  build_combined "$genome_key"
done

echo "=== Done. t2g.txt gene-column sample per genome: ==="
for genome_key in "${!GENOME_OUTDIR[@]}"; do
  echo "--- $genome_key ---"
  head -3 "${GENOME_OUTDIR[$genome_key]}/t2g.txt"
done
