# D. melanogaster – D. simulans ortholog mapping

Genome-wide table mapping *D. simulans* genes to their *D. melanogaster*
orthologs, built via reciprocal best BLAST hit (RBH) with DIAMOND.

## Why this exists

The legacy FlyBase `dsim-all-r1.0.gff` annotation includes a
`dmel_ortholog=` attribute, but that file is from the original ~2005
GLEAN-based *D. simulans* annotation and predates standard GFF3 gene
models (it's built from `match`/`match_part` alignment records, not
`gene`/`mRNA`/`exon`/`CDS`), so it's unusable with standard tools like
`gffread` and its ortholog calls are outdated.

The obvious replacement — FlyBase's precomputed
`dmel_orthologs_in_drosophila_species` bulk TSV — no longer exists.
FlyBase moved OrthoDB-derived orthology calls to a live per-gene API
in 2022 and dropped them from precomputed files, GFF, and JBrowse.
DIOPT, the other standard Drosophila ortholog tool, doesn't cover
*D. simulans* either. Querying the FlyBase API gene-by-gene isn't
practical at genome scale.

So this pipeline computes orthologs directly via reciprocal best BLAST
hit between the two species' proteomes, which is self-contained and
doesn't depend on any external ortholog database's availability.

## Method

1. Download translation (protein) FASTAs for both species from FlyBase.
   - *D. melanogaster*: current release (r6.68, FB2026_02)
   - *D. simulans*: r2.01 (FB2015_01) — the last FlyBase-curated
     annotation with standard gene models and FBgn IDs. (FlyBase no
     longer updates *D. simulans* sequence/annotation at all; NCBI's
     GNOMON pipeline now maintains it, but under NCBI gene IDs rather
     than FBgn, which don't cross-reference this pipeline's output —
     hence sticking with FlyBase r2.01 here.)
2. Parse each FASTA header to map protein ID → FBgn gene ID
   (`get_id_map.py`).
3. Build DIAMOND databases for both proteomes and run `blastp` in
   both directions (dsim→dmel and dmel→dsim), keeping only the single
   best hit per query protein.
4. Collapse protein-level hits to gene level and keep only gene pairs
   that are each other's best hit in both directions
   (`filter_rbh.py`).

RBH is stricter than FlyBase's old ortholog calls — it won't capture
many-to-many relationships from gene duplications/paralogs — so expect
somewhat fewer pairs than older `dmel_ortholog=` annotations gave,
that's expected.

## Files

| File | Purpose |
|---|---|
| `run_rbh_orthologs.sbatch` | SLURM pipeline: download, DIAMOND makedb, reciprocal blastp, RBH filtering |
| `get_id_map.py` | Extracts protein ID → FBgn gene ID from a FlyBase translation FASTA header |
| `filter_rbh.py` | Collapses two one-directional best-hit tables into gene-level reciprocal best hits |
| `README.md` | This file |

## Running it

Requires a `diamond` mamba env:

```bash
mamba create -n diamond -c bioconda -c conda-forge diamond
```

Compute nodes on this cluster (`phoenix-08.prism`, `medium` partition)
don't have outbound internet access, so the proteome FASTAs must be
downloaded from a login node first — the sbatch script skips the
download step automatically if the files already exist in `WORKDIR`:

```bash
mkdir -p /private/groups/russelllab/jodie/scRNAseq/reference/orthologs
cd /private/groups/russelllab/jodie/scRNAseq/reference/orthologs
wget -O dmel.fasta.gz "https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.68_FB2026_02/fasta/dmel-all-translation-r6.68.fasta.gz"
gunzip -c dmel.fasta.gz > dmel_proteins.fasta
wget -O dsim.fasta.gz "https://s3ftp.flybase.org/genomes/Drosophila_simulans/dsim_r2.01_FB2015_01/fasta/dsim-all-translation-r2.01.fasta.gz"
gunzip -c dsim.fasta.gz > dsim_proteins.fasta
```

Then submit the job:

```bash
sbatch /private/groups/russelllab/jodie/scRNAseq/Lum_et_al_2027_scRNAseq_Wolbachia_infected_embryos_and_primary_cell_lines/scripts/orthologs/run_rbh_orthologs.sbatch
```

Double-check `WORKDIR` and `SCRIPT_DIR` at the top of the sbatch
script point to the right paths on your system before submitting.

## Output

`dmel_dsim_orthologs_rbh.tsv` in `WORKDIR`
(`/private/groups/russelllab/jodie/scRNAseq/reference/orthologs/`):

| Column | Meaning |
|---|---|
| `geneA` | *D. simulans* FBgn |
| `geneB` | *D. melanogaster* FBgn |
| `pident` | % identity of the reciprocal best hit |
| `evalue` | BLAST e-value |
| `bitscore` | BLAST bitscore |

Intermediate files also left in `WORKDIR` (useful for debugging):

- `dmel_id_map.tsv` / `dsim_id_map.tsv` — protein ID → FBgn maps
- `dmel_db.dmnd` / `dsim_db.dmnd` — DIAMOND databases
- `dsim_vs_dmel.tsv` / `dmel_vs_dsim.tsv` — one-directional best-hit
  results, protein-level, before reciprocal filtering

## Result (2026-08-11 run)

- 13,794 *D. simulans* genes queried, 13,761 *D. melanogaster* genes
  queried
- **13,299 reciprocal best-hit ortholog pairs**, consistent with
  *D. melanogaster*'s ~13,900 protein-coding genes
