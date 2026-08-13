# Prompt for Claude Code (run on razzmatazz/emerald, in this repo)

## Context

This is a Snakemake pipeline (`Snakefile`, `config/config.yaml`) for scRNA-seq
of Wolbachia-infected Drosophila (Dmel/Dsim/Dwil) embryos and primary cell
lines, with several Wolbachia strains (wMel, wRi_Riv84, wRi_M23, wWil).
Samples are processed with kallisto|bustools against combined host+Wolbachia
+16S references, then QC-filtered and a per-cell `wolbachia_titer` is
calculated in `snakemake_scripts/filtering/kallisto_bustools_qcfilter_adata_no_pybiomart.py`.

I've been debugging rRNA gene identification for the Wolbachia titer
calculation in a chat session (not on this server) with no direct filesystem
access, so a lot of this was diagnosed by asking me to paste `bash` output
back and forth. That's slow -- please pick up from here with direct access to
`/private/groups/russelllab/jodie/scRNAseq/` and actually verify/finish this
end to end.

## What's already been fixed (already applied to the repo)

1. **`Snakefile`, `get_symbiont_16s_region()`** (~line 160-210): matched
   `locus_tag`/`gene_id` too strictly against `config['symbiont_16s_gene']`
   bare locus tags (e.g. `WRI_RS06005`). The Wolbachia GTFs' 16S entries come
   from `cmsearch`/Infernal and have `gene_id "gene-WRI_RS06005"` /
   `transcript_id "rna-WRI_RS06005"` (prefixed, no `locus_tag` attribute at
   all). Fixed to also match with an optional `gene-`/`rna-` prefix. `import re`
   was added.

2. **`snakemake_scripts/reference/find_rrna_genes.py`**: only parsed
   GTF-quoted attributes (`key "value";`). The strains' *source* GFF3 files
   (`genomic.gff` / `genomic.fixed.gff` next to each Wolbachia genome fasta)
   use unquoted GFF3 syntax (`key=value;key=value`) and are the ones that
   actually carry proper `rRNA` feature-typing + `gbkey=rRNA` + bare
   `locus_tag` -- the *derived* `.gtf` files (e.g.
   `wMel_GCF_016584425.1.gtf`) do NOT reliably carry this (see below).
   Fixed `parse_attrs()` to fall back to GFF3-style parsing when GTF-style
   finds nothing. Also changed the default `--id-attr` priority from
   `gene, gene_id, locus_tag` to `locus_tag, gene, gene_id`, since
   `locus_tag` is the one that's always present/bare in the GFF3s.

3. **New script `snakemake_scripts/reference/find_host_rrna_genes.py`**:
   the host (FlyBase) GTFs don't reliably expose `rRNA` feature type or
   `gbkey`/`gene_biotype` attributes the way NCBI GTFs do, so this instead
   matches FlyBase's own rRNA gene *symbol* convention
   (`^\d+(\.\d+)?SrRNA` e.g. `18SrRNA:CR41548`, `28SrRNA-Psi:CR40596`),
   reading `gene_symbol` directly off the raw FlyBase GTF (no need to run it
   through `swap_gene_id_to_symbol.py` first). Strips the `Dsim\`/`Dwil\`
   species prefix. Mitochondrial rRNA (`mt:lrRNA`/`mt:srRNA`) excluded by
   default (Wolbachia alters host mtDNA copy number, which would confound
   titer). **Confirmed working for Dmel** (`dmel-all-r6.68.gtf`).

4. **Dsim has no rRNA genes annotated at all** in its current FlyBase r2.02
   GTF (`dsim-all-r2.02.gtf`) -- confirmed via
   `grep -io 'gene_symbol "[^"]*rrna[^"]*"'` returning nothing, and no
   `SrRNA`-pattern symbols anywhere in the file. This is a real annotation
   gap (rDNA repeats are hard to assemble), not a script bug. **Still open**:
   I asked whether an alternative NCBI RefSeq Dsim assembly
   (`/private/groups/russelllab/jodie/scRNAseq/reference/kb_references_species_specific_id/Drosophila_simulans_GCF_016746395.1/Drosophila_simulans_GCF_016746395.1.gtf`)
   has rRNA annotated instead -- **please check this now** (feature-type
   counts, any `rRNA` feature rows, `gbkey`/`gene_biotype` values, grep for
   "ribosomal"). If it does, we may want to build a proper Dsim reference
   from that NCBI assembly instead of the old FlyBase one, or at minimum use
   its GTF just to get a host rRNA gene list.

5. **Wolbachia titer calculation changed to a stopgap** in
   `snakemake_scripts/filtering/kallisto_bustools_qcfilter_adata_no_pybiomart.py`,
   function `calculate_wolbachia_titer()`: it used to be a true ratio
   `symbiont / (symbiont + host)` rRNA reads, but with host rRNA gene lists
   missing for some species (Dsim) that ratio silently collapsed to ~1.0 for
   any cell with symbiont reads (misleading, not a real titer). **Per my
   request**, it's now just **raw total symbiont (Wolbachia) rRNA transcript
   counts per cell**, not normalized against host at all. This is
   intentional and temporary -- I said I'd deal with proper normalization
   later. `--host_rrna_genes` CLI arg is still accepted but currently
   unused (kept for when normalization is revisited). Don't "fix" this back
   to a ratio without checking with me first -- I know it's not a real titer
   right now.

## Important finding that still needs investigating: AGAT-mangled Wolbachia gene IDs

The derived Wolbachia GTFs (e.g.
`/private/groups/russelllab/jodie/scRNAseq/reference/wMel_GCF_016584425.1/wMel_GCF_016584425.1.gtf`)
were built with **AGAT** (source column shows `AGAT`/`tRNAscan-SE`, not
`cmsearch` for everything). Sample lines:

```
NZ_CP046925.1	AGAT	CDS	1	1383	.	+	0	transcript_id "gene-GQX67_RS00005"; gene_id "agat-gene-1"; gene_name "dnaA";
NZ_CP046925.1	AGAT	transcript	1	1383	.	+	.	transcript_id "gene-GQX67_RS00005"; gene_id "agat-gene-1"; gene_name "dnaA"
```

For **protein-coding genes**, AGAT rewrote `gene_id` to synthetic sequential
IDs (`agat-gene-1`, `agat-gene-2`, ...) completely disconnected from the
NCBI locus tag. Only the non-coding features (`cmsearch`/`tRNAscan-SE`
sourced -- rRNA, tRNA) kept `gene_id "gene-<locus_tag>"`.

**This means**: if the combined kb reference for any Wolbachia strain was
built from this derived `.gtf` (via `gffread` + `kb ref`, see
`snakemake_scripts/alignment/build_combined_host_wolbachia_references.sh`),
then `adata.var_names` for essentially every *protein-coding* Wolbachia gene
would be a meaningless `agat-gene-N` string -- not just an rRNA-list problem,
but potentially invalidating any per-gene Wolbachia differential
expression/pathway analysis downstream. **Please check whether this is
actually the case** in the transcriptomes currently in use (see next
section for exactly which files/paths to check), and flag it clearly if so.
If it is, we likely need to rebuild the Wolbachia GTF from
`genomic.fixed.gff` (or `genomic.gff`) directly via `gffread`/`agat` with
settings that preserve `locus_tag` as `gene_id`, rather than the current
derived `.gtf`.

## Also still open: reference path mismatch

`config.yaml`'s genome dict (`Dmel_wMel`, `Dsim_wRi_Riv84`, etc., ~line
14-20) and `build_combined_host_wolbachia_references.sh`'s `COMBINED_ROOT`
both point to
`/private/groups/russelllab/jodie/scRNAseq/Lum_et_al_2027_scRNAseq_Wolbachia_infected_embryos_and_primary_cell_lines/reference/...`.
But when I searched for existing `t2g.txt` files, I found them at different
locations instead:

```
/private/groups/russelllab/jodie/scRNAseq/reference/kb_references_species_specific_id/Drosophila_simulans_wRi_Riv84_combined_16S/t2g.txt
/private/groups/russelllab/jodie/scRNAseq/reference/kb_references_species_specific_id/Drosophila_simulans_wRi_M23_combined_16S/t2g.txt
/private/groups/russelllab/jodie/scRNAseq/reference/kb_references_species_specific_id/Drosophila_simulans_wMel_combined_16S/t2g.txt
/private/groups/russelllab/jodie/scRNAseq/reference/kb_references_species_specific_id/Drosophila_willistoni_wWil_combined_16S/t2g.txt
/private/groups/russelllab/jodie/scRNAseq/reference/Drosophila_melanogaster_wWil_combined_16S/t2g.txt
```

None of these match the `config.yaml`-expected path exactly. **Please
confirm**: which `t2g.txt`/`index.idx` are actually referenced by
`config.yaml`'s genome dict right now (i.e. does
`Dmel_wMel: ".../Lum_et_al_2027_.../reference/Drosophila_melanogaster_wMel_combined_16S"`
actually exist and get used, or is `config.yaml` stale / pointing somewhere
empty while the pipeline actually reads from
`kb_references_species_specific_id/`)? If `config.yaml` is wrong, fix the
paths to point at the real, currently-used references.

## What I need you to do, in order

1. **Check the alternative Dsim NCBI assembly** for rRNA annotation (see
   item 4 above). If it has rRNA, figure out the best way to get a
   `host_rrna_genes.txt` for Dsim from it (may require running
   `find_rrna_genes.py`, not `find_host_rrna_genes.py`, if it's NCBI-style
   with `gbkey`/`locus_tag` rather than FlyBase-style symbols -- check the
   attribute format first).

2. **Regenerate `rrna_genes.txt` for every Wolbachia strain** from each
   strain's `genomic.fixed.gff` (fall back to `genomic.gff` if that doesn't
   exist) using the fixed `find_rrna_genes.py`, for:
   - `wMel_GCF_016584425.1`
   - `wRi_Riv84_GCF_000022285.1`
   - `wRi_M23_GCA_979474595.1`
   - `wWil_GCF_040084705.1`

   I expected **4 rRNA genes** for at least one of these strains (possibly a
   duplicated rRNA operon) but only saw 3 (23S, 5S, 16S) in a `head -5` of
   wMel's GFF3 earlier -- please get the actual full count per strain and
   let me know if any have 4 (or more) and what they are.

3. **Resolve the AGAT `agat-gene-N` question** (see section above) -- check
   whether the Wolbachia references actually in use have this problem, and
   report back clearly. This might be a bigger issue than the rRNA lists.

4. **Confirm `adata.var_names` format** for Wolbachia genes in a real
   filtered h5ad (e.g. something under `results/h5ad_results/` or
   `results/filtered_h5ad/` for a wMel or wRi sample) -- specifically check
   whether the regenerated `rrna_genes.txt` locus tags (bare, e.g.
   `GQX67_RS05935`) actually appear in `adata.var_names`, or whether
   `var_names` carry a `gene-`/`rna-` prefix or (per the AGAT issue) are
   `agat-gene-N` instead. Use `--id-attr` on `find_rrna_genes.py` to emit
   whatever format actually matches if bare locus tags don't.

5. **Resolve the reference path mismatch** (see section above) and fix
   `config.yaml` if it's pointing at the wrong/stale reference directories.

6. Once 1-5 are sorted, **re-run `filter_h5ad` for a couple of test samples**
   (one wMel, one wRi strain, one Dmel host, one Dsim host if possible) and
   confirm in the log (`logs/filter/{sample_id}.log`) that
   `calculate_wolbachia_titer` reports non-zero "Symbiont rRNA genes
   present" counts, and note down what `wolbachia_titer` values look like
   (mean/median) so we can sanity check them against expectations.

## Conventions / environment

- mamba environments; `scanpy_env` etc. are defined in `config/config.yaml`.
  Activate with `source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh && conda activate <env>`
  (see any rule's `shell:` block in `Snakefile` for the pattern).
- SLURM cluster (`snakemake --executor slurm ...`), but for these ad hoc
  diagnostic checks just run things directly on the login/interactive node.
- Please don't re-run the full `build_dmel_dsim_dwil_transcriptomes.sh` --
  I already built the host transcriptomes manually and don't want them
  regenerated/overwritten.
- Show me exact commands + output as you go, and call out clearly anywhere
  you had to make a judgment call (e.g. picking an `--id-attr`, choosing
  which GFF3 file to use) so I can sanity check it.
