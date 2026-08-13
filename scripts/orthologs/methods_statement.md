## Ortholog identification (methods statement)

Orthologous gene pairs between *D. melanogaster* and *D. simulans*
were identified using a reciprocal best BLAST hit (RBH) approach.
Protein sequences were extracted from each species' reference genome
and annotation with gffread v[X.X.X] (Pertea & Pertea 2020): *D.
melanogaster* release r6.68 (FlyBase) and *D. simulans* Prin_Dsim_3.1
(NCBI RefSeq assembly GCF_016746395.2, Gnomon annotation — FlyBase no
longer maintains *D. simulans* sequence or annotation). The resulting
proteomes were searched against each other using DIAMOND blastp
v[X.X.X] (Buchfink et al. 2021) with an e-value threshold of 1e-5,
retaining a single best hit per query. Gene pairs that were each
other's best hit in both directions were retained as putative
one-to-one orthologs.

[Re-run and fill in: N reciprocal best-hit ortholog pairs from N
*D. simulans* and N *D. melanogaster* genes queried. The previous run
(13,299 pairs from 13,794 *D. simulans* / 13,761 *D. melanogaster*
genes) used FlyBase r2.01 for *D. simulans* and is no longer current
now that *D. simulans* is sourced from NCBI RefSeq.]

**Citation for DIAMOND:**
Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at
tree-of-life scale using DIAMOND. *Nat Methods*. 2021;18(4):366-368.

**Citation for gffread:**
Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. *F1000Res*.
2020;9:ISCB Comm J-304.

**Before using this in a manuscript:**
- Fill in the DIAMOND and gffread versions — run `diamond --version`
  and `gffread --version` in the relevant mamba envs and drop them in
  above.
- Fill in the updated gene/pair counts from the NCBI RefSeq
  (`GCF_016746395.2`) run.
- If reviewers want it, the full pipeline and scripts are in this
  repo at `scripts/orthologs/` (see `README.md`) for a methods/code
  availability statement.
