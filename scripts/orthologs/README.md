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

On top of that, FlyBase no longer maintains *D. simulans* sequence or
annotation at all — the current *D. simulans* reference used in this
project is NCBI RefSeq Prin_Dsim_3.1 (`GCF_016746395.2`, Gnomon gene
models, `LOC#####`/accession-style gene IDs), which shares no ID
namespace with FlyBase FBgn IDs and has no FlyBase-style translation
FASTA to download.

So this pipeline computes orthologs directly via reciprocal best BLAST
hit between the two species' proteomes — extracted locally from each
species' own genome FASTA + GTF — which is self-contained and doesn't
depend on any external ortholog database or a FlyBase-formatted
proteome being available for *D. simulans*.

## Method

1. Extract protein sequences directly from each species' genome FASTA +
   GTF with `gffread -y` (one protein per transcript with a CDS):
   - *D. melanogaster*: FlyBase r6.68 (`dmel-all-chromosome-r6.68.fasta`
     / `dmel-all-r6.68.gtf`)
   - *D. simulans*: NCBI RefSeq Prin_Dsim_3.1, `GCF_016746395.2`
     (`GCF_016746395.2_Prin_Dsim_3.1_genomic.fna` /
     `GCF_016746395.2_Prin_Dsim_3.1.gtf`)
2. Map each gffread protein FASTA header to its gene ID by reading
   `gene_id`/`transcript_id`/`protein_id` straight out of the GTF
   attribute column (`gtf_id_map.py`). This is generic GTF parsing, not
   FlyBase-header parsing, so it works the same way on the FlyBase dmel
   GTF and the NCBI/Gnomon dsim GTF.
   **Important:** `gffread -y` writes `protein_id` (not `transcript_id`)
   as the FASTA header whenever the GTF's CDS rows carry one -- true for
   both FlyBase (`protein_id "FBpp######"`) and NCBI/Gnomon
   (`protein_id "XP_######"`) GTFs. `gtf_id_map.py` maps both
   `protein_id` and `transcript_id` to `gene_id` so the lookup resolves
   either way; mapping transcript_id alone would let raw FBpp/XP_
   accessions leak through into the final `Dsim`/`Dmel` gene_id columns
   instead of FBgn/LOC gene IDs.
3. Build DIAMOND databases for both proteomes and run `blastp` in
   both directions (dsim→dmel and dmel→dsim), keeping only the single
   best hit per query protein.
4. Collapse protein-level hits to gene level and keep only gene pairs
   that are each other's best hit in both directions
   (`filter_rbh.py`), output with columns named `Dsim`/`Dmel` to match
   what `config.yaml`'s `ortholog_map` and the integration script
   (`integrate_v2.py`) expect.

RBH is stricter than FlyBase's old ortholog calls — it won't capture
many-to-many relationships from gene duplications/paralogs — so expect
somewhat fewer pairs than older `dmel_ortholog=` annotations gave,
that's expected. Also note: `Dsim` gene IDs in the output table are now
NCBI/Gnomon IDs (e.g. `LOC120284240`), not FBgn — make sure whatever
built the *D. simulans* kallisto\|bustools t2g uses the same
`GCF_016746395.2` GTF, so `adata.var_names` for Dsim samples land in the
same ID namespace as this table's `Dsim` column.

## Files

| File | Purpose |
|---|---|
| `run_rbh_orthologs.sbatch` | SLURM pipeline: gffread protein extraction, DIAMOND makedb, reciprocal blastp, RBH filtering |
| `gtf_id_map.py` | Extracts protein_id/transcript_id → gene_id from any GTF's attribute column (FlyBase or NCBI/Gnomon) |
| `get_id_map.py` | Legacy: extracts protein_id → FBgn from a FlyBase translation FASTA header (`parent=FBgn...`) — only usable if you have that kind of FASTA for both species |
| `filter_rbh.py` | Collapses two one-directional best-hit tables into gene-level reciprocal best hits |
| `README.md` | This file |

## Running it

Requires two mamba envs:

```bash
# gffread -- already present in this repo's kallisto_bustools env
# diamond:
mamba create -n diamond -c bioconda -c conda-forge diamond
```

Everything now runs from genome FASTA + GTF files already on disk (no
download step, so no login-node/compute-node internet dependency).
Double-check `WORKDIR`, `SCRIPT_DIR`, and the `DMEL_FASTA`/`DMEL_GTF`/
`DSIM_FASTA`/`DSIM_GTF`/env-name variables at the top of
`run_rbh_orthologs.sbatch` point to the right paths on your system, then
submit:

```bash
sbatch /private/groups/russelllab/jodie/scRNAseq/Lum_et_al_2027_scRNAseq_Wolbachia_infected_embryos_and_primary_cell_lines/scripts/orthologs/run_rbh_orthologs.sbatch
```

## Output

`dmel_dsim_orthologs_rbh.tsv` in `WORKDIR`
(`/private/groups/russelllab/jodie/scRNAseq/reference/orthologs/`):

| Column | Meaning |
|---|---|
| `Dsim` | *D. simulans* gene ID (NCBI/Gnomon, e.g. `LOC120284240`) |
| `Dmel` | *D. melanogaster* FlyBase gene ID (FBgn) |
| `pident` | % identity of the reciprocal best hit |
| `evalue` | BLAST e-value |
| `bitscore` | BLAST bitscore |

Intermediate files also left in `WORKDIR` (useful for debugging):

- `dmel_proteins.fasta` / `dsim_proteins.fasta` — gffread-extracted
  proteomes
- `dmel_id_map.tsv` / `dsim_id_map.tsv` — transcript_id → gene_id maps
- `dmel_db.dmnd` / `dsim_db.dmnd` — DIAMOND databases
- `dsim_vs_dmel.tsv` / `dmel_vs_dsim.tsv` — one-directional best-hit
  results, protein-level, before reciprocal filtering

## Result (2026-08-11 run, FlyBase r2.01 D. simulans)

The counts below are from the prior run against FlyBase's r2.01
*D. simulans* annotation, before the switch to NCBI RefSeq
Prin_Dsim_3.1. Re-run and update this section after the next run against
`GCF_016746395.2`.

- 13,794 *D. simulans* genes queried, 13,761 *D. melanogaster* genes
  queried
- **13,299 reciprocal best-hit ortholog pairs**, consistent with
  *D. melanogaster*'s ~13,900 protein-coding genes
