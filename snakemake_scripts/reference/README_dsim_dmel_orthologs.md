# D. simulans -> D. melanogaster ortholog mapping

Four scripts in this directory build a per-gene ortholog table between the
FlyBase D. simulans annotation (`dsim-all-r2.02.gtf`) and the D.
melanogaster reference annotation (`dmel-all-r6.68.gtf`). Run them in
order: 1 -> 2 -> 3 -> 4.

| # | Script | What it does | Needs |
|---|--------|---------------|-------|
| 1 | `find_dmel_dsim_orthologs.py` | Curated symbol-based calls, no new tools | your existing GTFs only |
| 2 | `run_reciprocal_blast_orthologs.sh` | Sequence-based calls (reciprocal best BLAST hit) | `gffread`, `blast` |
| 3 | `reciprocal_best_hits.py` | Called automatically by script 2 -- not usually run by hand | -- |
| 4 | `merge_ortholog_calls.py` | Combines 1 + 2 into one final table | output of 1 and 2 |

## Why two methods

FlyBase gives a D. simulans gene its D. melanogaster homolog's symbol
*only when it has confidently identified that homolog* (e.g.
`Dsim\Adh`, `Dsim\CG1704`). Every other gene keeps a `GD#####` systematic
ID -- Dsim's equivalent of Dmel's own `CG` numbers -- which carries no
ortholog information at all. In this annotation, most genes are
`GD`-numbered, so the symbol method (script 1) only resolves a minority.
Script 2 fills in the rest by BLASTing predicted proteins between species
and keeping reciprocal best hits (RBH): each gene's best hit in the other
species, and vice versa.

Script 1's calls are curated by FlyBase and preferred where they exist;
script 2's calls are a sequence-similarity inference and are used for
everything else. Script 4 merges the two with that precedence.

## 1. Symbol-based calls

```bash
mamba activate kallisto_bustools   # just needs python3

DMEL_GTF=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_melanogaster/dmel-all-r6.68.gtf
DSIM_GTF=/private/groups/russelllab/jodie/scRNAseq/reference/Flybase_genomes/Drosophila_simulans/dsim-all-r2.02.gtf

python3 snakemake_scripts/reference/find_dmel_dsim_orthologs.py \
  --dsim-gtf "$DSIM_GTF" --dmel-gtf "$DMEL_GTF" \
  -o dsim_to_dmel_orthologs.tsv \
  --unmatched-out dsim_unmatched.txt
```

Output `dsim_to_dmel_orthologs.tsv` columns: `dsim_gene_id`, `dsim_symbol`,
`dsim_symbol_clean`, `dmel_gene_id`, `dmel_symbol`, `match_status`.

`match_status` values:
- `matched` -- confident 1:1 call.
- `ambiguous_dsim_paralogs` -- two or more Dsim genes project to the same
  Dmel gene (recent duplication). Each still gets that Dmel gene as its
  call; flagged so you know it's not unique.
- `ambiguous_dmel_symbol` -- the (rare) case where the Dsim symbol matches
  more than one Dmel gene. `dmel_gene_id` lists all candidates,
  comma-separated.
- `unmatched` -- no Dmel gene shares that symbol (typically a
  `GD#####`-numbered gene). Written to `dsim_unmatched.txt` for step 2.

## 2. BLAST reciprocal-best-hit calls

Needs `blast` (`blastp`/`makeblastdb`), not currently in any mamba env in
this repo:

```bash
mamba install -n kallisto_bustools -c bioconda blast
# or: mamba create -n blast -c bioconda -c conda-forge blast gffread
```

Open `run_reciprocal_blast_orthologs.sh` and check the `EDIT THESE` block
at the top (genome fasta/GTF paths, `OUTDIR`, `THREADS`). Since most genes
need this method here (not just the leftovers), leave `RESTRICT_IDS=""` to
run genome-wide rather than pointing it at `dsim_unmatched.txt`:

```bash
mamba activate kallisto_bustools   # or your blast env
bash snakemake_scripts/reference/run_reciprocal_blast_orthologs.sh
```

This extracts protein sequences with `gffread -y`, BLASTs each proteome
against the other (`-evalue 1e-5 -max_target_seqs 5`), and calls
`reciprocal_best_hits.py` to collapse hits to the gene level and keep only
reciprocal best hits. Output: `$OUTDIR/dsim_dmel_rbh_orthologs.tsv`, columns
`dsim_gene_id`, `dmel_gene_id`, `pident`, `evalue`, `bitscore`.

Runtime: proteome-wide BLASTP for two ~13-14k gene genomes is minutes to
low tens of minutes at `THREADS=16`; request that many cores if running
under SLURM.

You normally don't run `reciprocal_best_hits.py` directly -- it's called
by the shell script -- but its `--restrict-ids` flag is there if you ever
want to rerun RBH-calling alone against a subset of Dsim genes without
re-BLASTing.

## 3. Merge into one final table

```bash
python3 snakemake_scripts/reference/merge_ortholog_calls.py \
  --symbol-table dsim_to_dmel_orthologs.tsv \
  --rbh-table /path/to/OUTDIR/dsim_dmel_rbh_orthologs.tsv \
  --dmel-gtf "$DMEL_GTF" \
  -o dsim_to_dmel_orthologs_final.tsv
```

Output `dsim_to_dmel_orthologs_final.tsv` columns: `dsim_gene_id`,
`dsim_symbol`, `dmel_gene_id`, `dmel_symbol`, `method`, `flag`, `pident`,
`evalue`, `bitscore` (BLAST stats blank for symbol-based rows).

`method`:
- `symbol` -- from step 1 (`matched` or `ambiguous_dsim_paralogs`).
- `blast_rbh` -- from step 2, used because step 1 left this gene
  `unmatched` or `ambiguous_dmel_symbol`.
- `none` -- neither method produced a call; `dmel_gene_id` is blank.

`flag` worth checking by hand:
- `ambiguous_dsim_paralogs` -- shared Dmel target, not unique.
- `blast_conflicts_symbol_candidates` -- step 1 had specific Dmel
  candidates from a shared symbol, but BLAST's reciprocal best hit picked
  a different gene entirely. Worth a manual look.
- `blast_confirms_symbol_candidate` -- BLAST's pick matches one of step
  1's candidates -- resolves the ambiguity.
- `unmatched` (with `method` = `none`) -- no ortholog call from either
  method; likely a Dsim-specific gene or one BLAST also couldn't place
  confidently (consider loosening `-evalue`/`-max_target_seqs` in step 2
  if you want to push harder on these).

The script prints a summary of how many genes fall into each
`(method, flag)` combination to stderr when it runs.

## Quick sanity checks

```bash
# overall match rate
awk -F'\t' 'NR>1{c[$5]++} END{for (m in c) print m, c[m]}' dsim_to_dmel_orthologs_final.tsv

# genes with no call at all
awk -F'\t' 'NR>1 && $3==""' dsim_to_dmel_orthologs_final.tsv | wc -l
```
