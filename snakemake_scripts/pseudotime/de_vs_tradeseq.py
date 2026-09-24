#!/usr/bin/env python3
"""
de_vs_tradeseq.py
=================
Rule pseudotime_de_concordance. Do the discrete pseudobulk DE (de_stages.py)
and the continuous tradeSeq transition tests (tradeseq_species.R) agree?

For every trajectory, its species' DE contrast is compared with that
trajectory's tradeSeq transition test:
  primary_vs_embryo  <->  tradeseq_embryo_to_primary.csv
  line_vs_primary    <->  tradeseq_primary_to_cellline.csv
Reported: Spearman rho of fold changes, overlap of significant genes, and sign
agreement among genes significant in both.
"""

import os
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from pt_utils import savefig

PAIRS = [("primary_vs_embryo", "embryo_to_primary"),
         ("line_vs_primary", "primary_to_cellline")]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--de_dir", required=True)
    p.add_argument("--groups", nargs="+", required=True, help="trajectory names")
    p.add_argument("--species", nargs="+", required=True, help="species of each trajectory")
    p.add_argument("--tradeseq_root", default="results/pseudotime")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    rows = []
    fig, axes = plt.subplots(len(args.groups), len(PAIRS),
                             figsize=(4.5 * len(PAIRS), 4 * len(args.groups)), squeeze=False)
    for i, (g, sp_) in enumerate(zip(args.groups, args.species)):
        for j, (de_name, ts_name) in enumerate(PAIRS):
            ax = axes[i, j]
            de_p = os.path.join(args.de_dir, f"de_{sp_}_{de_name}.csv")
            ts_p = os.path.join(args.tradeseq_root, g, "tradeseq", f"tradeseq_{ts_name}.csv")
            if not (os.path.exists(de_p) and os.path.exists(ts_p)):
                print(f"  missing input for {g} {de_name}; skipping")
                ax.axis("off")
                continue
            de = pd.read_csv(de_p, index_col=0)
            ts = pd.read_csv(ts_p).set_index("gene")
            lfc_col = [c for c in ts.columns if c.startswith("logFC")][0]
            u = de.index.intersection(ts.index)
            d = pd.DataFrame({"symbol": de.loc[u, "symbol"],
                              "de_log2FC": de.loc[u, "log2FoldChange"],
                              "de_padj": de.loc[u, "padj"],
                              # tradeSeq's own logFC scale (from its GAM link);
                              # only rank/sign are compared, so scale doesn't matter
                              "ts_logFC": ts.loc[u, lfc_col],
                              "ts_padj": ts.loc[u, "padj"]}).dropna(subset=["de_log2FC", "ts_logFC"])
            de_sig = d["de_padj"] < args.alpha
            ts_sig = d["ts_padj"] < args.alpha
            both = de_sig & ts_sig
            rho = spearmanr(d["de_log2FC"], d["ts_logFC"])[0]
            agree = (np.sign(d.loc[both, "de_log2FC"]) == np.sign(d.loc[both, "ts_logFC"])).mean()
            d.to_csv(os.path.join(args.out_dir, f"concordance_{g}_{de_name}.csv"))
            rows.append(dict(trajectory=g, species=sp_, contrast=de_name, n_genes=len(d),
                             de_sig=int(de_sig.sum()), tradeseq_sig=int(ts_sig.sum()),
                             both_sig=int(both.sum()),
                             jaccard=both.sum() / max(1, (de_sig | ts_sig).sum()),
                             logFC_spearman=rho, sign_agreement_both_sig=agree))
            ax.scatter(d.loc[~both, "de_log2FC"], d.loc[~both, "ts_logFC"], s=2, c="#d0d0d0")
            ax.scatter(d.loc[both, "de_log2FC"], d.loc[both, "ts_logFC"], s=3, c="#8e44ad")
            ax.axhline(0, c="grey", lw=0.5); ax.axvline(0, c="grey", lw=0.5)
            ax.set_xlabel(f"pseudobulk DE log2FC ({sp_})")
            ax.set_ylabel("tradeSeq logFC")
            ax.set_title(f"{g}: {de_name}\nrho = {rho:.2f}, sig in both = {int(both.sum())}",
                         fontsize=8)
    savefig(fig, os.path.join(args.out_dir, "de_vs_tradeseq_scatter.pdf"))
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(args.out_dir, "de_vs_tradeseq_summary.csv"), index=False)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
