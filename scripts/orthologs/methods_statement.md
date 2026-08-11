## Ortholog identification (methods statement)

Orthologous gene pairs between *D. melanogaster* and *D. simulans*
were identified using a reciprocal best BLAST hit (RBH) approach.
Protein sequences were obtained from FlyBase (*D. melanogaster*
release r6.68, FB2026_02; *D. simulans* release r2.01, FB2015_01,
the most recent FlyBase-curated annotation with standard gene models
for this species) and searched against each other using DIAMOND
blastp v[X.X.X] (Buchfink et al. 2021) with an e-value threshold of
1e-5, retaining a single best hit per query. Gene pairs that were each
other's best hit in both directions were retained as putative
one-to-one orthologs. This yielded 13,299 reciprocal best-hit ortholog
pairs from 13,794 *D. simulans* and 13,761 *D. melanogaster* genes
queried.

**Citation for DIAMOND:**
Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at
tree-of-life scale using DIAMOND. *Nat Methods*. 2021;18(4):366-368.

**Before using this in a manuscript:**
- Fill in the DIAMOND version — run `diamond --version` in the
  `diamond` mamba env and drop it in above.
- If reviewers want it, the full pipeline and scripts are in this
  repo at `scripts/orthologs/` (see `README.md`) for a methods/code
  availability statement.
