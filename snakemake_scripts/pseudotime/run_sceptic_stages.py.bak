#!/usr/bin/env python3
"""
run_sceptic_stages.py
=====================
Step 2 (rule pseudotime_sceptic). Supervised SCEPTIC pseudotime over the
three culture stages, per species, plus the readouts along pseudotime.

Adapted from Jacobs et al. 2026 run_sceptic.py. The 2026 version trained on
infection time points (D0..D56). Here the labels are the stage of
immortalization: embryo = 0, primary_cells = 1, cell_culture = 2. SCEPTIC
pseudotime is sum_k P(stage k) * k, so it runs on a 0-2 scale in both
species. That shared scale is what lets the joint tradeSeq model compare them.

SCEPTIC uses its own internal 3-fold KFold on cells, so every cell's
pseudotime is an out-of-fold prediction. Caveat: folds split cells, not
samples. With stage fully confounded with sequencing run, the classifier
can use batch signal, and pseudotime will pile up near 0/1/2. Read the
within-stage spread, and how the lineages order inside each stage, as the
informative part.

Readouts
--------
  confusion matrix, pseudotime by sample, UMAP, SCEPTIC vs DPT concordance,
  Wolbachia titer vs pseudotime, atlas identity + label-transfer confidence
  across pseudotime bins.

Also exports the stratified tradeSeq input (--tradeseq_cells_per_sample
cells per sample, genes detected in >= --tradeseq_min_frac of those cells).
"""

import os
import argparse
import warnings

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

from sceptic import run_sceptic_and_evaluate

from pt_utils import savefig as _savefig, load_gene_labels

warnings.filterwarnings("ignore", category=FutureWarning)

STAGE_ORDER = ["embryo", "primary_cells", "cell_culture"]
STAGE_COLORS = {"embryo": "#4C72B0", "primary_cells": "#DD8452", "cell_culture": "#C44E52"}


def _sp(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return np.nan, np.nan, int(m.sum())
    r, p = spearmanr(x[m], y[m])
    return r, p, int(m.sum())


def _condition_order(obs):
    return (obs[["condition", "stage_numeric", "lineage"]].drop_duplicates()
            .sort_values(["stage_numeric", "lineage", "condition"])["condition"].tolist())


# ─────────────────────────────────────────────────────────────────────────────
# SCEPTIC
# ─────────────────────────────────────────────────────────────────────────────

def run(adata, n_pcs, method):
    embed = adata.uns.get("pt_embedding", "X_pca")
    X = np.asarray(adata.obsm[embed][:, :n_pcs])
    labels = adata.obs["stage_numeric"].astype(int).values
    label_list = np.array(sorted(np.unique(labels)))
    print(f"SCEPTIC on {embed}[:, :{n_pcs}] | {len(labels)} cells | "
          f"label counts {pd.Series(labels).value_counts().sort_index().to_dict()}")
    cm, pred, pt, prob = run_sceptic_and_evaluate(X, labels, label_list=label_list,
                                                   method=method)
    stages = [STAGE_ORDER[i] for i in label_list]
    adata.obs["sceptic_pseudotime"] = pt
    adata.obs["sceptic_pred_stage"] = pd.Categorical(
        [stages[int(i)] for i in pred], categories=stages)
    for j, s in enumerate(stages):
        adata.obs[f"sceptic_prob_{s}"] = prob[:, j]
    return pd.DataFrame(cm, index=stages, columns=stages)


# ─────────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────────

def plot_confusion(cm, fig_dir, species):
    cm.to_csv(os.path.join(fig_dir, f"confusion_matrix_{species}.csv"))
    norm = cm.div(cm.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(4.5, 4))
    sns.heatmap(norm, annot=cm.astype(int), fmt="d", cmap="Blues", vmin=0, vmax=1, ax=ax,
                cbar_kws={"label": "row fraction"})
    ax.set_xlabel("Predicted stage"); ax.set_ylabel("True stage")
    ax.set_title(f"{species}: SCEPTIC out-of-fold predictions")
    _savefig(fig, os.path.join(fig_dir, f"confusion_matrix_{species}.pdf"))


def plot_pt_by_condition(obs, fig_dir, species):
    order = _condition_order(obs)
    fig, ax = plt.subplots(figsize=(max(6, len(order) * 0.8), 4.5))
    sns.violinplot(data=obs, x="condition", y="sceptic_pseudotime", order=order,
                   hue="sample_type", hue_order=STAGE_ORDER, palette=STAGE_COLORS,
                   dodge=False, cut=0, inner="quartile", ax=ax)
    ax.set_ylim(-0.05, 2.05); ax.set_ylabel("SCEPTIC pseudotime (0-2)")
    ax.set_title(f"{species}: pseudotime by sample"); plt.xticks(rotation=45, ha="right")
    _savefig(fig, os.path.join(fig_dir, f"pseudotime_by_condition_{species}.pdf"))


def plot_dpt_concordance(obs, fig_dir, species, stats):
    if "dpt_pseudotime" not in obs:
        return
    fig, ax = plt.subplots(figsize=(5.5, 5))
    for st in STAGE_ORDER:
        m = obs["sample_type"] == st
        ax.scatter(obs.loc[m, "dpt_pseudotime"], obs.loc[m, "sceptic_pseudotime"],
                   s=2, alpha=0.3, c=STAGE_COLORS[st], label=st)
    r, p, n = _sp(obs["dpt_pseudotime"].values, obs["sceptic_pseudotime"].values)
    stats.append(dict(test="sceptic_vs_dpt", subset="all", rho=r, p=p, n=n))
    for st in STAGE_ORDER:
        m = (obs["sample_type"] == st).values
        r_s, p_s, n_s = _sp(obs["dpt_pseudotime"].values[m], obs["sceptic_pseudotime"].values[m])
        stats.append(dict(test="sceptic_vs_dpt", subset=st, rho=r_s, p=p_s, n=n_s))
    ax.set_xlabel("DPT (unsupervised)"); ax.set_ylabel("SCEPTIC (supervised)")
    ax.set_title(f"{species}: Spearman rho = {r:.2f}"); ax.legend(markerscale=4)
    _savefig(fig, os.path.join(fig_dir, f"sceptic_vs_dpt_{species}.pdf"))


def plot_titer(obs, fig_dir, species, stats, n_bins):
    if "wolbachia_titer" not in obs or obs["wolbachia_titer"].notna().sum() < 10:
        print("  No titer values -- skipping titer plots")
        return
    df = obs[["sceptic_pseudotime", "wolbachia_titer", "condition", "sample_type",
              "lineage"]].dropna(subset=["wolbachia_titer"])
    r, p, n = _sp(df["sceptic_pseudotime"].values, df["wolbachia_titer"].values)
    stats.append(dict(test="titer_vs_pseudotime", subset="all", rho=r, p=p, n=n))
    for key, col in [("stage", "sample_type"), ("condition", "condition")]:
        for g, sub in df.groupby(col, observed=True):
            r_g, p_g, n_g = _sp(sub["sceptic_pseudotime"].values, sub["wolbachia_titer"].values)
            stats.append(dict(test="titer_vs_pseudotime", subset=f"{key}:{g}",
                              rho=r_g, p=p_g, n=n_g))

    order = _condition_order(obs)
    pal = dict(zip(order, sns.color_palette("tab20", len(order))))
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.scatterplot(data=df, x="sceptic_pseudotime", y="wolbachia_titer", hue="condition",
                    hue_order=order, palette=pal, s=4, alpha=0.4, linewidth=0, ax=ax)
    ax.set_title(f"{species}: titer vs pseudotime (rho = {r:.2f})")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, markerscale=3)
    _savefig(fig, os.path.join(fig_dir, f"titer_vs_pseudotime_{species}.pdf"))

    df["pt_bin"] = pd.cut(df["sceptic_pseudotime"], bins=np.linspace(0, 2, n_bins + 1),
                          include_lowest=True)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    sns.boxplot(data=df, x="pt_bin", y="wolbachia_titer", hue="lineage",
                showfliers=False, ax=ax)
    ax.set_xlabel("SCEPTIC pseudotime bin"); plt.xticks(rotation=60, ha="right", fontsize=7)
    ax.set_title(f"{species}: titer across pseudotime bins")
    _savefig(fig, os.path.join(fig_dir, f"titer_by_pseudotime_bin_{species}.pdf"))


def plot_identity(obs, fig_dir, species, n_bins, conf_thr, top_n=12, min_cells=20):
    """Atlas identity composition + confidence across pseudotime bins
    (replaces the 2026 leiden stacked-area/enrichment plots)."""
    for label in ["atlas_annotation", "atlas_germ_layer"]:
        conf = f"{label}_confidence"
        if label not in obs:
            continue
        df = obs[["sceptic_pseudotime", label, "lineage"] +
                 ([conf] if conf in obs else [])].copy()
        df[label] = df[label].astype(str)
        if conf in df:
            df.loc[df[conf] < conf_thr, label] = "low_confidence"
        top = df[label].value_counts().index[:top_n]
        df[label] = np.where(df[label].isin(top), df[label], "other")
        edges = np.linspace(0, 2, n_bins + 1)
        centers = (edges[:-1] + edges[1:]) / 2
        df["bin"] = pd.cut(df["sceptic_pseudotime"], bins=edges, labels=False,
                           include_lowest=True)
        frac = (df.groupby(["bin", label]).size().unstack(fill_value=0)
                .reindex(range(n_bins), fill_value=0))
        n_per_bin = frac.sum(axis=1)
        frac = frac.div(n_per_bin.replace(0, np.nan), axis=0)
        out = frac.copy(); out.insert(0, "n_cells", n_per_bin); out.insert(0, "bin_center", centers)
        out.to_csv(os.path.join(fig_dir, f"identity_along_pseudotime_{label}_{species}.csv"))

        # stacked bars; bins with < min_cells are left empty rather than
        # interpolated (SCEPTIC pseudotime is often clumped near 0/1/2)
        fig, ax = plt.subplots(figsize=(9, 4.5))
        f = frac[n_per_bin >= min_cells].fillna(0)
        bottom = np.zeros(len(f))
        for col, color in zip(f.columns, sns.color_palette("tab20", f.shape[1])):
            ax.bar(centers[f.index], f[col].values, width=2 / n_bins * 0.95,
                   bottom=bottom, color=color, label=col)
            bottom += f[col].values
        ax.set_xlim(0, 2); ax.set_ylim(0, 1)
        ax.set_xlabel("SCEPTIC pseudotime"); ax.set_ylabel("Fraction of cells")
        ax.set_title(f"{species}: {label} along pseudotime")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        _savefig(fig, os.path.join(fig_dir, f"identity_along_pseudotime_{label}_{species}.pdf"))

        if conf in obs:
            d = obs[["sceptic_pseudotime", conf, "lineage"]].copy()
            d["bin_center"] = centers[pd.cut(d["sceptic_pseudotime"], bins=edges,
                                             labels=False, include_lowest=True)]
            fig, ax = plt.subplots(figsize=(7, 4))
            sns.lineplot(data=d, x="bin_center", y=conf, hue="lineage",
                         errorbar="sd", marker="o", ax=ax)
            ax.set_xlabel("SCEPTIC pseudotime"); ax.set_ylabel(f"{conf}")
            ax.set_title(f"{species}: atlas label-transfer confidence along pseudotime")
            _savefig(fig, os.path.join(fig_dir, f"confidence_along_pseudotime_{label}_{species}.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# tradeSeq export
# ─────────────────────────────────────────────────────────────────────────────

def export_tradeseq(adata, out_dir, n_per_sample, min_frac, labels, seed=42):
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(seed)
    idx = []
    for _, pos in adata.obs.groupby("source_file", observed=True).indices.items():
        idx.extend(rng.choice(pos, size=min(n_per_sample, len(pos)), replace=False))
    sub = adata[np.sort(idx)]
    counts = sub.layers["counts"]
    counts = counts.tocsc() if sp.issparse(counts) else sp.csc_matrix(counts)
    frac = np.asarray((counts > 0).mean(axis=0)).ravel()
    keep = frac >= min_frac
    counts = counts[:, keep]
    genes = sub.var_names[keep]
    print(f"tradeSeq input: {sub.n_obs} cells "
          f"({n_per_sample}/sample cap), {keep.sum()} genes detected in >= {min_frac:.0%}")

    mat = sp.csr_matrix(counts.T)
    mat.data = np.round(mat.data).astype(np.int64)
    scipy.io.mmwrite(os.path.join(out_dir, "counts_genesXcells.mtx"), mat)
    pd.DataFrame({"gene": genes, "label": [labels.get(g, g) for g in genes]}).to_csv(
        os.path.join(out_dir, "genes.tsv"), sep="\t", index=False)
    sub.obs[["sceptic_pseudotime", "condition", "lineage", "sample_type",
             "source_file"]].to_csv(os.path.join(out_dir, "cells.csv"))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True, help="prepare_species.py output")
    p.add_argument("--species", required=True)
    p.add_argument("--n_pcs", type=int, default=30)
    p.add_argument("--method", default="xgboost", choices=["xgboost", "svm"])
    p.add_argument("--n_bins", type=int, default=20)
    p.add_argument("--conf_threshold", type=float, default=0.5)
    p.add_argument("--tradeseq_cells_per_sample", type=int, default=1000)
    p.add_argument("--tradeseq_min_frac", type=float, default=0.05)
    p.add_argument("--flybase_annotation", default=None)
    p.add_argument("--ortholog_map", default=None)
    p.add_argument("--out_h5ad", required=True)
    p.add_argument("--out_obs_csv", required=True)
    p.add_argument("--tradeseq_dir", required=True)
    p.add_argument("--fig_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.fig_dir, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    cm = run(adata, args.n_pcs, args.method)
    obs = adata.obs
    print("\nMedian pseudotime per sample:")
    print(obs.groupby(["sample_type", "condition"], observed=True)["sceptic_pseudotime"]
          .median().to_string())

    stats = []
    acc = np.trace(cm.values) / cm.values.sum()
    stats.append(dict(test="sceptic_accuracy", subset="all", rho=acc, p=np.nan,
                      n=int(cm.values.sum())))
    r, pv, n = _sp(obs["stage_numeric"].values.astype(float), obs["sceptic_pseudotime"].values)
    stats.append(dict(test="sceptic_vs_stage", subset="all", rho=r, p=pv, n=n))

    plot_confusion(cm, args.fig_dir, args.species)
    plot_pt_by_condition(obs, args.fig_dir, args.species)
    plot_dpt_concordance(obs, args.fig_dir, args.species, stats)
    plot_titer(obs, args.fig_dir, args.species, stats, args.n_bins)
    plot_identity(obs, args.fig_dir, args.species, args.n_bins, args.conf_threshold)
    sc.settings.figdir = args.fig_dir
    sc.pl.umap(adata, color=["sceptic_pseudotime", "sceptic_pred_stage"], show=False,
               save=f"_sceptic_{args.species}.pdf")

    pd.DataFrame(stats).to_csv(os.path.join(args.fig_dir, f"sceptic_stats_{args.species}.csv"),
                               index=False)
    cols = ["species", "lineage", "stage_numeric", "sceptic_pseudotime", "sceptic_pred_stage",
            "dpt_pseudotime"] + [c for c in obs if c.startswith("sceptic_prob_")]
    obs[[c for c in cols if c in obs]].to_csv(args.out_obs_csv)

    labels = load_gene_labels(args.flybase_annotation, args.ortholog_map, args.species)
    export_tradeseq(adata, args.tradeseq_dir, args.tradeseq_cells_per_sample,
                    args.tradeseq_min_frac, labels)
    adata.write(args.out_h5ad)
    print(f"\nWrote {args.out_h5ad}")


if __name__ == "__main__":
    main()
