#!/usr/bin/env python3
"""
species_candidates.py
=====================
Rule pseudotime_species_candidates. Rank candidate genes that could explain
why D. simulans establishes cell lines more readily than D. melanogaster, for
knockout in Dsim or transfer of the Dsim allele into Dmel.

Why not the species x stage interaction: the Dsim embryo samples differ from the
Dmel embryos (e.g. ecdysone-cascade genes), consistent with a stage mismatch,
so 'primary vs embryo' interactions are dominated by embryo differences.
Instead we rank CULTURE-EMERGENT species differences, using pseudobulk
log2CPM averaged per lineage (each lineage counts once):
  - embryos similar:            |mean Dsim embryo - mean Dmel embryo| <= --max_embryo_diff
  - primary cell lines differ:  every Dsim primary > every Dmel primary by >= --min_gap
                                (or the reverse for Dmel-high genes)
  - expressed:                  higher species' primaries >= --min_expr log2CPM
  - maintained_in_lines:        the same ordering holds in every cell line by >= 1,
                                including the uninfected Dsim lines, so the
                                difference is not kept up by current infection
Annotations: shared cell-line signature membership, and (if the
pseudotime_primary_proliferative output exists) the gene's log2FC in the
proliferative population vs the rest in each primary cell line.
Cross-species levels carry mapping bias (Dsim reads mapped to the Dsim genome,
then to 1:1 Dmel orthologs; gene models differ), so genes without a symbol,
non-coding genes and genes at 0 in every Dmel sample are flagged for manual
checking of the gene model and read coverage before any experiment.

Outputs: candidates_Dsim_high.csv, candidates_Dmel_high.csv, species_gap_all.csv,
candidates_heatmap.pdf
"""
import os
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from pt_utils import savefig


def lineage_means(cpm, S, sp, st):
    s = S[(S["species"] == sp) & (S["stage"] == st)]
    return pd.concat({l: cpm[g.index].mean(axis=1) for l, g in s.groupby("lineage")}, axis=1)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--de_dir", required=True)
    p.add_argument("--prolif_dir", default=None, help="pseudotime_primary_proliferative output")
    p.add_argument("--max_embryo_diff", type=float, default=1.0)
    p.add_argument("--min_gap", type=float, default=1.5)
    p.add_argument("--min_expr", type=float, default=4.0)
    p.add_argument("--top_n", type=int, default=40)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    D = args.de_dir

    cnt = pd.read_csv(os.path.join(D, "pseudobulk_counts.csv.gz"), index_col=0)
    S = pd.read_csv(os.path.join(D, "pseudobulk_samples.csv"), index_col=0)
    sym = pd.read_csv(os.path.join(D, "de_pooled_line_vs_primary.csv"), index_col=0)["symbol"]
    cpm = np.log2(cnt.div(cnt.sum()) * 1e6 + 1)
    E = {sp: lineage_means(cpm, S, sp, "embryo") for sp in ("Dmel", "Dsim")}
    P = {sp: lineage_means(cpm, S, sp, "primary") for sp in ("Dmel", "Dsim")}
    L = {sp: lineage_means(cpm, S, sp, "line") for sp in ("Dmel", "Dsim")}

    t = pd.DataFrame({"symbol": sym.reindex(cpm.index)}, index=cpm.index)
    t["embryo_diff_Dsim_minus_Dmel"] = E["Dsim"].mean(1) - E["Dmel"].mean(1)
    for tag, hi, lo in [("Dsim_high", "Dsim", "Dmel"), ("Dmel_high", "Dmel", "Dsim")]:
        t[f"primary_gap_{tag}"] = P[hi].min(1) - P[lo].max(1)
        t[f"line_gap_{tag}"] = L[hi].min(1) - L[lo].max(1)
    for sp in ("Dmel", "Dsim"):
        t[f"{sp}_embryo_mean"] = E[sp].mean(1)
        t[f"{sp}_primary_min"], t[f"{sp}_primary_max"] = P[sp].min(1), P[sp].max(1)
        t[f"{sp}_line_min"], t[f"{sp}_line_max"] = L[sp].min(1), L[sp].max(1)

    # shared cell-line signature membership
    f = [os.path.join(D, f"de_{sp}_line_vs_primary.csv") for sp in ("Dmel", "Dsim")]
    if all(os.path.exists(x) for x in f):
        a, b = (pd.read_csv(x, index_col=0) for x in f)
        j = a[["log2FoldChange", "padj"]].join(b[["log2FoldChange", "padj"]], lsuffix="_m",
                                                rsuffix="_s", how="inner").dropna()
        sig = j[(j.padj_m < 0.05) & (j.padj_s < 0.05)
                & (np.sign(j.log2FoldChange_m) == np.sign(j.log2FoldChange_s))]
        t["cellline_signature"] = np.where(t.index.isin(sig.index[sig.log2FoldChange_m > 0]), "up",
                                  np.where(t.index.isin(sig.index[sig.log2FoldChange_m < 0]), "down", ""))

    # proliferative-population markers per lineage
    if args.prolif_dir and os.path.isdir(args.prolif_dir):
        for mf in sorted(glob.glob(os.path.join(args.prolif_dir, "*", "markers.csv"))):
            lin = os.path.basename(os.path.dirname(mf))
            m = pd.read_csv(mf).set_index("names")
            t[f"prolif_lfc_{lin}"] = m["logfoldchanges"].reindex(t.index)
            t[f"prolif_padj_{lin}"] = m["pvals_adj"].reindex(t.index)

    no_symbol = t["symbol"].isna() | t["symbol"].astype(str).str.startswith("FBgn")
    ncrna = t["symbol"].astype(str).str.match(r"^(?:CR\d|lncRNA:|asRNA:|snoRNA|tRNA|mir-)")
    t["check_gene_model"] = no_symbol | ncrna | (P["Dmel"].max(1) < 0.5) & (E["Dmel"].mean(1) < 0.5) \
        | (P["Dsim"].max(1) < 0.5) & (E["Dsim"].mean(1) < 0.5)
    t.to_csv(os.path.join(args.out_dir, "species_gap_all.csv"))

    ok = t["embryo_diff_Dsim_minus_Dmel"].abs() <= args.max_embryo_diff
    out = {}
    for tag, hi in [("Dsim_high", "Dsim"), ("Dmel_high", "Dmel")]:
        c = t[ok & (t[f"primary_gap_{tag}"] >= args.min_gap)
              & (t[f"{hi}_primary_min"] >= args.min_expr)].copy()
        c["maintained_in_lines"] = c[f"line_gap_{tag}"] >= 1
        c = c.sort_values(["check_gene_model", "maintained_in_lines", f"primary_gap_{tag}"],
                          ascending=[True, False, False])
        c.to_csv(os.path.join(args.out_dir, f"candidates_{tag}.csv"))
        out[tag] = c
        clean = c[~c["check_gene_model"]]
        print(f"\n{tag}: {len(c)} genes ({int(c['maintained_in_lines'].sum())} maintained in all cell "
              f"lines; {int(c['check_gene_model'].sum())} flagged to check gene model)")
        cols = ["symbol", "embryo_diff_Dsim_minus_Dmel", f"primary_gap_{tag}", f"line_gap_{tag}",
                "maintained_in_lines"] + [x for x in ["cellline_signature"] if x in c] + \
               [x for x in c if x.startswith("prolif_lfc_")]
        print(clean[cols].head(args.top_n).round(2).to_string(index=False))

    # heatmap of the top clean candidates across all samples
    top = pd.concat([out[k][~out[k]["check_gene_model"]].head(args.top_n // 2) for k in out])
    if len(top):
        order = S.sort_values(["species", "lineage", "stage"]).index
        z = cpm.loc[top.index, order]
        z.index = top["symbol"].astype(str)
        z.columns = [f"{S.at[c, 'lineage']} | {S.at[c, 'stage']}" for c in order]
        fig, ax = plt.subplots(figsize=(0.35 * z.shape[1] + 4, 0.22 * z.shape[0] + 2))
        sns.heatmap(z, cmap="viridis", ax=ax, cbar_kws={"label": "log2CPM"})
        ax.tick_params(axis="y", labelsize=6); plt.xticks(rotation=60, ha="right", fontsize=7)
        ax.set_title("Culture-emergent species differences (top Dsim-high, then Dmel-high)", fontsize=9)
        savefig(fig, os.path.join(args.out_dir, "candidates_heatmap.pdf"))
    print("Done ->", args.out_dir)


if __name__ == "__main__":
    main()
