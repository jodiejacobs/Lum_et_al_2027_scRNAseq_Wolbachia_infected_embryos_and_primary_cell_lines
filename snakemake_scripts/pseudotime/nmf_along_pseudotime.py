#!/usr/bin/env python3
"""
nmf_along_pseudotime.py
=======================
Step 5 (rule pseudotime_nmf). Per-species NMF gene programs and how their
usage changes along SCEPTIC pseudotime. Programs are fit independently per
species and matched across species later (compare_species.py) by
ortholog overlap of their top genes.

Same NMF setup as Jacobs et al. 2026 nmf_programs.py: sklearn NMF,
init='nndsvda', random_state=42, on log-normalised HVG expression.

Outputs (in --out_dir)
  nmf_top_genes_<species>.csv      : program, rank, gene, label, weight
  nmf_usage_<species>.csv          : per-cell usage (rows sum to 1)
  nmf_usage_along_pt_<species>.csv : mean usage per pseudotime bin x lineage
  nmf_pt_correlation_<species>.csv : Spearman usage ~ pseudotime per program
  nmf_usage_heatmap_<species>.pdf, nmf_usage_curves_<species>.pdf
"""

import os
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.decomposition import NMF
from statsmodels.stats.multitest import multipletests

from pt_utils import savefig as _savefig, load_gene_labels


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True, help="run_sceptic_stages.py output")
    p.add_argument("--species", required=True)
    p.add_argument("--flybase_annotation", default=None)
    p.add_argument("--ortholog_map", default=None)
    p.add_argument("--n_programs", type=int, default=12)
    p.add_argument("--n_top", type=int, default=100)
    p.add_argument("--n_bins", type=int, default=20)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    sp_ = args.species

    adata = sc.read_h5ad(args.h5ad)
    hvg = adata.var_names[adata.var["highly_variable"]]
    X = adata[:, hvg].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    print(f"NMF: {X.shape[0]} cells x {X.shape[1]} HVGs, k = {args.n_programs}")

    model = NMF(n_components=args.n_programs, init="nndsvda", random_state=42, max_iter=500)
    W = model.fit_transform(X)
    H = model.components_
    names = [f"{sp_}_P{i + 1:02d}" for i in range(args.n_programs)]

    labels = load_gene_labels(args.flybase_annotation, args.ortholog_map, sp_)
    rows = []
    for i, nm in enumerate(names):
        order = np.argsort(H[i])[::-1][:args.n_top]
        for r, j in enumerate(order):
            rows.append(dict(program=nm, rank=r + 1, gene=hvg[j],
                             label=labels.get(hvg[j], hvg[j]), weight=H[i, j]))
    top = pd.DataFrame(rows)
    top.to_csv(os.path.join(args.out_dir, f"nmf_top_genes_{sp_}.csv"), index=False)
    short = {nm: nm + " (" + ", ".join(top.loc[top.program == nm, "label"].head(4)) + ")"
             for nm in names}

    usage = pd.DataFrame(W / np.clip(W.sum(axis=1, keepdims=True), 1e-12, None),
                         index=adata.obs_names, columns=names)
    usage.to_csv(os.path.join(args.out_dir, f"nmf_usage_{sp_}.csv"))

    pt = adata.obs["sceptic_pseudotime"].values
    cor = []
    for nm in names:
        r, pv = spearmanr(usage[nm].values, pt)
        cor.append(dict(program=nm, top_genes=short[nm], rho=r, p=pv))
    cor = pd.DataFrame(cor)
    cor["padj"] = multipletests(cor["p"], method="fdr_bh")[1]
    cor.sort_values("rho").to_csv(
        os.path.join(args.out_dir, f"nmf_pt_correlation_{sp_}.csv"), index=False)

    edges = np.linspace(0, 2, args.n_bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    df = usage.copy()
    df["bin_center"] = centers[pd.cut(pt, bins=edges, labels=False, include_lowest=True)]
    df["lineage"] = adata.obs["lineage"].astype(str).values
    long = df.melt(id_vars=["lineage", "bin_center"], var_name="program", value_name="usage")
    summ = (long.groupby(["program", "lineage", "bin_center"])["usage"]
            .agg(["mean", "std", "count"]).reset_index())
    summ.to_csv(os.path.join(args.out_dir, f"nmf_usage_along_pt_{sp_}.csv"), index=False)

    # Heatmap: program x pseudotime bin (all lineages pooled), z-scored per program
    pooled = df.groupby("bin_center")[names].mean()
    pooled = pooled[df.groupby("bin_center").size() >= 20]
    z = (pooled - pooled.mean()) / pooled.std(ddof=0).replace(0, np.nan)
    order = cor.set_index("program").loc[names, "rho"].sort_values().index
    fig, ax = plt.subplots(figsize=(10, 0.4 * len(names) + 2))
    sns.heatmap(z[order].T.fillna(0), cmap="RdBu_r", center=0, ax=ax,
                yticklabels=[short[n] for n in order],
                xticklabels=[f"{c:.2f}" for c in z.index], cbar_kws={"label": "z (mean usage)"})
    ax.set_xlabel("SCEPTIC pseudotime bin centre")
    ax.set_title(f"{sp_}: NMF program usage along pseudotime")
    _savefig(fig, os.path.join(args.out_dir, f"nmf_usage_heatmap_{sp_}.pdf"))

    ncol = 4
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3 * nrow), sharex=True)
    for ax, nm in zip(axes.ravel(), names):
        s = summ[(summ.program == nm) & (summ["count"] >= 20)]
        for lin, d in s.groupby("lineage"):
            ax.plot(d["bin_center"], d["mean"], marker="o", ms=3, label=lin)
        ax.set_title(short[nm], fontsize=7)
        ax.axvline(1, ls="--", c="grey", lw=0.5)
    for ax in axes.ravel()[len(names):]:
        ax.axis("off")
    axes.ravel()[0].legend(fontsize=6)
    fig.supxlabel("SCEPTIC pseudotime (0 embryo, 1 primary, 2 cell line)")
    fig.supylabel("mean program usage")
    _savefig(fig, os.path.join(args.out_dir, f"nmf_usage_curves_{sp_}.pdf"))


if __name__ == "__main__":
    main()
