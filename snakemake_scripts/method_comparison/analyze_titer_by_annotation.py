"""
analyze_titer_by_annotation.py
===============================
integrate_v2.py's analyze_titer_by_cluster(), ported to group by a
transferred cell-type label (e.g. atlas_annotation from
integrate_via_atlas_projection.py / annotate_with_flysta3d_ingest.py)
instead of obs['leiden'].

Why a groupby column instead of Leiden here specifically
----------------------------------------------------------
Leiden clusters on the atlas-projected object would be clusters of the
ATLAS's own developmental biology, re-discovered on an embedding that was
never fit to be sensitive to your titer/condition axis in the first place
(see integrate_via_atlas_projection.py's docstring) -- not a meaningful
grouping for "does titer differ across cell types". atlas_<label> IS a
meaningful grouping: it's your best call on what cell type each cell
actually is, independent of how well or badly any embedding separates it.
This script answers "does Wolbachia titer differ by cell-type identity",
which is the question you actually want, directly off that column, no
clustering step needed at all.

Every plot function below takes --groupby so you're not locked to
atlas_annotation specifically -- atlas_tissue, atlas_germ_layer, or
(after map_cellline_to_embryo.py) embryo_annotation all work the same way.

Run with:
    mamba activate scanpy
    python snakemake_scripts/method_comparison/analyze_titer_by_annotation.py \\
        --adata results/integrated/integrated_atlas.h5ad \\
        --groupby atlas_annotation \\
        --condition_col condition \\
        --fig_dir results/integrated/figures_atlas \\
        --sample wolbachia_infection
"""

import os
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import kruskal

from annotate_with_flysta3d import _savefig


def _group_colors(adata, groupby):
    groups = sorted(adata.obs[groupby].dropna().astype(str).unique())
    cmap = plt.colormaps["tab20"]
    return groups, [cmap(i % 20) for i in range(len(groups))]


def analyze_titer_by_annotation(adata, fig_dir, sample, groupby="atlas_annotation",
                                 condition_col="condition", origin_col=None,
                                 titer_col="wolbachia_titer"):
    """
    How does Wolbachia titer vary across transferred cell-type categories?

    Produces (mirrors analyze_titer_by_cluster's plot set, grouped by
    `groupby` instead of 'leiden'):
      - Boxplot + strip: titer per category
      - Violin: titer per category
      - Bar: % infected cells (titer > 0) per category
      - Heatmap: mean titer per category x condition (if condition_col present)
      - Boxplot: titer per category, split by origin (if origin_col/inferrable)
      - Kruskal-Wallis test across categories
    """
    if titer_col not in adata.obs.columns:
        print(f"  Skipping titer analysis: no '{titer_col}' column")
        return None
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"'{groupby}' not found in adata.obs -- available columns: "
            f"{list(adata.obs.columns)}"
        )

    os.makedirs(fig_dir, exist_ok=True)

    print("\n" + "=" * 60)
    print(f"WOLBACHIA TITER BY {groupby.upper()}")
    print("=" * 60)

    adata = adata.copy()
    adata.obs[groupby] = adata.obs[groupby].astype(str)

    if origin_col is None or origin_col not in adata.obs.columns:
        if "source_file" in adata.obs.columns:
            adata.obs["_origin"] = np.where(
                adata.obs["source_file"].astype(str).str.contains("embryo", case=False),
                "embryo", "cell_line",
            )
            origin_col = "_origin"
        else:
            origin_col = None

    groups, colors = _group_colors(adata, groupby)

    keep_cols = [groupby, titer_col]
    if condition_col in adata.obs.columns:
        keep_cols.append(condition_col)
    if origin_col:
        keep_cols.append(origin_col)
    obs = adata.obs[keep_cols].copy()
    obs_titer = obs.dropna(subset=[titer_col])

    # ── Kruskal-Wallis ────────────────────────────────────────────────
    value_groups = [obs_titer[obs_titer[groupby] == g][titer_col].values
                     for g in groups if (obs_titer[groupby] == g).sum() > 0]
    if len(value_groups) < 2:
        print(f"  Only {len(value_groups)} non-empty {groupby} group(s) with "
              f"'{titer_col}' -- skipping Kruskal-Wallis")
        h, p_kw = np.nan, np.nan
    else:
        h, p_kw = kruskal(*value_groups)
        print(f"\nKruskal-Wallis: H={h:.2f}  p={p_kw:.2e}")
        print(f"Titer {'SIGNIFICANTLY' if p_kw < 0.05 else 'does NOT significantly'} "
              f"differ across {groupby}")

    stats = obs_titer.groupby(groupby)[titer_col].agg(["mean", "median", "std", "count"])
    print(f"\nTiter summary per {groupby}:")
    print(stats.to_string())
    stats.to_csv(os.path.join(fig_dir, f"titer_stats_by_{groupby}_{sample}.csv"))

    title_stat = f"Kruskal-Wallis H={h:.2f}  p={p_kw:.2e}" if np.isfinite(h) else ""

    # ── Plot 1: boxplot + strip ──────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(10, len(groups) * 0.8), 6))
    sns.stripplot(data=obs_titer, x=groupby, y=titer_col, order=groups,
                  color="black", alpha=0.15, size=1.5, ax=ax)
    bp = ax.boxplot(
        [obs_titer[obs_titer[groupby] == g][titer_col].values for g in groups],
        positions=range(len(groups)), widths=0.55, patch_artist=True,
        whiskerprops=dict(alpha=0.7), capprops=dict(alpha=0.7),
        medianprops=dict(color="black", linewidth=2),
        flierprops=dict(markersize=1),
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=45, ha="right")
    ax.set_xlabel(groupby)
    ax.set_ylabel(titer_col)
    ax.set_title(f"Wolbachia titer by {groupby} -- {sample}\n{title_stat}")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_boxplot_by_{groupby}_{sample}.pdf"))

    # ── Plot 2: violin ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(10, len(groups) * 0.8), 6))
    sc.pl.violin(adata, titer_col, groupby=groupby, ax=ax, show=False, rotation=45)
    ax.set_title(f"Wolbachia titer distribution by {groupby} -- {sample}\n{title_stat}")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_violin_by_{groupby}_{sample}.pdf"))

    # ── Plot 3: % infected per category ──────────────────────────────
    pct_infected = obs_titer.groupby(groupby)[titer_col].apply(
        lambda x: (x > 0).sum() / len(x) * 100
    ).reindex(groups)
    fig, ax = plt.subplots(figsize=(max(8, len(groups) * 0.8), 5))
    ax.bar(range(len(pct_infected)), pct_infected.values, color=colors,
           alpha=0.85, edgecolor="black", linewidth=0.4)
    ax.set_xticks(range(len(pct_infected)))
    ax.set_xticklabels(pct_infected.index, rotation=45, ha="right")
    ax.set_xlabel(groupby)
    ax.set_ylabel("% Wolbachia-positive cells")
    ax.set_title(f"Infection rate by {groupby} -- {sample}")
    ax.grid(axis="y", alpha=0.3)
    for i, v in enumerate(pct_infected.values):
        if np.isfinite(v):
            ax.text(i, v + 0.5, f"{v:.1f}%", ha="center", va="bottom", fontsize=7)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"infection_pct_by_{groupby}_{sample}.pdf"))

    # ── Plot 4: mean titer heatmap -- category x condition ───────────
    if condition_col in obs_titer.columns:
        mean_titer = obs_titer.pivot_table(
            values=titer_col, index=groupby, columns=condition_col, aggfunc="mean",
            observed=True)
        fig, ax = plt.subplots(figsize=(max(10, mean_titer.shape[1] * 1.2),
                                        max(5, mean_titer.shape[0] * 0.5)))
        sns.heatmap(mean_titer, cmap="viridis", ax=ax,
                    cbar_kws={"label": f"Mean {titer_col}"}, linewidths=0.3)
        ax.set_xlabel(condition_col)
        ax.set_ylabel(groupby)
        ax.set_title(f"Mean Wolbachia titer -- {groupby} x {condition_col} -- {sample}")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        _savefig(fig, os.path.join(
            fig_dir, f"titer_heatmap_{groupby}_{condition_col}_{sample}.pdf"))
    else:
        print(f"  Skipping category x condition heatmap: no '{condition_col}' column")

    # ── Plot 5: titer by category x origin (infected cells only) ─────
    infected = obs_titer[obs_titer[titer_col] > 0]
    if origin_col and len(infected) > 0:
        fig, ax = plt.subplots(figsize=(12, 5))
        sns.boxplot(data=infected, x=groupby, y=titer_col, hue=origin_col,
                    ax=ax, flierprops=dict(markersize=1), order=groups)
        ax.set_xlabel(groupby)
        ax.set_ylabel(f"{titer_col} (infected cells only)")
        ax.set_title(f"Titer by {groupby}, by {origin_col} -- {sample}")
        ax.legend(title=origin_col, bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        _savefig(fig, os.path.join(
            fig_dir, f"titer_by_{groupby}_origin_{sample}.pdf"))

    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adata", required=True,
                         help="e.g. results/integrated/integrated_atlas.h5ad "
                              "(output of integrate_via_atlas_projection.py)")
    parser.add_argument("--groupby", default="atlas_annotation",
                         help="Cell-type column to group by, e.g. "
                              "atlas_annotation, atlas_tissue, "
                              "atlas_germ_layer, embryo_annotation")
    parser.add_argument("--condition_col", default="condition")
    parser.add_argument("--origin_col", default=None,
                         help="obs column distinguishing embryo vs. cell "
                              "line. Omit to infer from 'embryo' appearing "
                              "in source_file, matching the Snakefile's own "
                              "EMBRYO_SAMPLE_IDS convention.")
    parser.add_argument("--titer_col", default="wolbachia_titer")
    parser.add_argument("--sample", default="integrated")
    parser.add_argument("--fig_dir", required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.adata)
    print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes from {args.adata}")

    analyze_titer_by_annotation(
        adata, args.fig_dir, args.sample,
        groupby=args.groupby, condition_col=args.condition_col,
        origin_col=args.origin_col, titer_col=args.titer_col,
    )


if __name__ == "__main__":
    main()
