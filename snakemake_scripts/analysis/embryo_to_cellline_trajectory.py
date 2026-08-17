"""
embryo_to_cellline_trajectory.py
==================================
Exploratory analysis of results/integrated/integrated.h5ad (rule integrate's
output): how do the cultured primary cell lines relate to the embryonic
tissues they were derived from?

Cell line establishment isn't a developmental continuum (it's dissociation +
selection + culture adaptation), so this deliberately does NOT compute a
single pseudotime. Instead it runs several independent, complementary
readouts of "how embryo-like is this cell line, and which tissue does it
resemble":

  1. Tissue-of-origin composition  -- per-condition breakdown of the KNN-
     transferred cell_type_<label> calls (confidence-filtered).
  2. Label diversity / heterogeneity collapse -- Shannon entropy of cell
     type calls per condition (whole embryos should be far more diverse
     than an established line).
  3. Label-transfer confidence as a drift metric -- per condition, split by
     embryo vs. cell line.
  4. Pseudobulk correlation -- condition-level and tissue-level pseudobulk
     profiles (full gene set, log1p), hierarchically clustered, so you can
     see which conditions/tissues actually resemble each other
     transcriptome-wide, independent of the KNN calls.
  5. Marker-module scoring -- marker genes derived from the embryo cells'
     own tissue calls (rank_genes_groups), then scored (score_genes) in
     every cell/condition. A second, independent check on tissue identity
     that also reveals markers being RETAINED vs. LOST during culture.
  6. Cell cycle / proliferation shift -- phase composition and S/G2M scores
     per condition; cell lines should show a higher cycling fraction than
     the embryo tissue they came from.
  7. Wolbachia infection effects -- titer vs. confidence/composition/cell
     cycle within cell line cells (correlation + infected-vs-uninfected
     comparisons), since some lines have paired infected/uninfected
     versions.
  8. Host species check (Dmel vs. Dsim) -- since Dsim samples only reach
     the embryo atlas via ortholog remapping, a quick QC/biology check on
     whether their tissue calls are as confident/consistent as Dmel's.
  9. Leiden cluster composition -- cluster x embryo/cell-line and cluster x
     cell type, complementing integrate_v2.py's own condition-enrichment
     analysis with a cell-type-aware view.
  10. UMAP overview -- one plot per key variable for a quick visual pass.

Every plot has a matching CSV of the underlying numbers written next to it.
Nothing here is meant to be the final analysis -- it's meant to be looked at
and used to decide what's actually worth pursuing further.

Run with:
    mamba activate scanpy
    python snakemake_scripts/analysis/embryo_to_cellline_trajectory.py \\
        --input results/integrated/integrated.h5ad \\
        --fig_dir results/trajectory_analysis
"""

import os
import argparse

import numpy as np
import pandas as pd
import scipy.sparse
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _savefig(fig, path):
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"   Saved: {path}")


def _condition_order(obs):
    """Embryo conditions first (sorted), then cell line conditions (sorted)
    -- used to keep every bar/box plot in a consistent, biologically
    grouped order with embryos and cell lines visually separated."""
    df = obs[["condition", "is_embryo"]].drop_duplicates()
    embryo   = sorted(df.loc[df["is_embryo"].astype(bool), "condition"])
    cellline = sorted(df.loc[~df["is_embryo"].astype(bool), "condition"])
    return embryo + cellline


def _embryo_cellline_divider(ax, order, obs):
    """Draw a vertical dashed line between the embryo and cell line groups
    on a categorical x-axis using the order returned by _condition_order."""
    df = obs[["condition", "is_embryo"]].drop_duplicates().set_index("condition")
    n_embryo = int(df.loc[order, "is_embryo"].astype(bool).sum())
    if 0 < n_embryo < len(order):
        ax.axvline(n_embryo - 0.5, color="black", lw=1.2, ls="--", alpha=0.6)


_MISSING_LABEL_STRINGS = {"NA", "nan", "None", "none", ""}


def _clean_label_series(series):
    """cell_type_<label> columns are object-dtype and go through
    integrate_v2.py's _sanitize_obs() before being written to h5ad, which
    stringifies missing values to the literal "NA" (h5ad can't store mixed
    str/NaN object columns) rather than leaving real NaN. A plain
    .dropna()/.notna() on these columns after loading the integrated object
    back in would therefore silently miss all of the "missing" cells and
    let a spurious "NA" category leak into every composition/diversity/
    marker plot as if it were a real cell type. This restores real NaN for
    any of the stringified missing-value spellings so downstream
    .dropna()/.notna() calls behave as intended."""
    s = series.astype(str).str.strip()
    return s.mask(s.isin(_MISSING_LABEL_STRINGS), np.nan)


def detect_label_cols(adata):
    """cell_type_<label> columns written by integrate_v2.py's unification
    step (coalesced atlas_<label> / embryo_<label>), skipping the
    _confidence variants."""
    cols = sorted({
        c for c in adata.obs.columns
        if c.startswith("cell_type_") and not c.endswith("_confidence")
    })
    if not cols:
        print("  WARNING: no cell_type_<label> columns found in adata.obs -- "
              "did rule integrate run its label-unification step? Falling "
              "back to raw atlas_<label>/embryo_<label> columns if present.")
        cols = sorted({
            c[len("atlas_"):] for c in adata.obs.columns
            if c.startswith("atlas_") and not c.endswith("_confidence")
        })
        cols = [f"atlas_{c}" for c in cols]
    return cols


def get_raw_full_adata(adata):
    """Full-gene, log1p-normalised (pre-scale, pre-HVG-subset) AnnData with
    the complete obs metadata -- built from adata.raw (integrate_v2.py's
    preprocess() sets adata.raw = adata right after normalize_total+log1p,
    before HVG subsetting and scaling). Needed for pseudobulk / marker
    scoring, since adata.X itself is HVG-subset and z-scored."""
    if adata.raw is None:
        raise ValueError(
            "adata.raw is missing -- pseudobulk/marker analyses need the "
            "full-gene log1p-normalised matrix that integrate_v2.py's "
            "preprocess() stores there before HVG subsetting."
        )
    full = adata.raw.to_adata()
    full.obs = adata.obs.copy()
    return full


def pseudobulk_by_group(adata_full, groupby, min_cells=20):
    """Mean log1p expression per group (adata_full.X), full gene set.
    Loops over groups rather than building one giant per-cell DataFrame, to
    keep memory bounded regardless of gene count."""
    X = adata_full.X
    if scipy.sparse.issparse(X):
        X = X.tocsr()
    groups = adata_full.obs[groupby].astype(str)
    vc = groups.value_counts()
    keep = vc[vc >= min_cells].index.tolist()
    dropped = vc[vc < min_cells]
    if len(dropped):
        print(f"   pseudobulk[{groupby}]: dropping {len(dropped)} group(s) "
              f"with < {min_cells} cells: {dropped.to_dict()}")

    rows, kept = [], []
    for g in keep:
        mask = (groups == g).values
        sub = X[mask]
        mean = np.asarray(sub.mean(axis=0)).ravel() if scipy.sparse.issparse(sub) else sub.mean(axis=0)
        rows.append(mean)
        kept.append(g)
    return pd.DataFrame(rows, index=kept, columns=adata_full.var_names)


# ─────────────────────────────────────────────────────────────────────────────
# 1. Tissue-of-origin composition
# ─────────────────────────────────────────────────────────────────────────────

def analyze_composition(adata, label_cols, fig_dir, conf_threshold):
    print("\n" + "=" * 70)
    print("1. TISSUE-OF-ORIGIN COMPOSITION")
    print("=" * 70)
    order = _condition_order(adata.obs)

    for col in label_cols:
        conf_col = f"{col}_confidence"
        if conf_col in adata.obs.columns:
            mask = adata.obs[conf_col].astype(float) >= conf_threshold
            print(f"  {col}: {mask.sum()}/{adata.n_obs} cells pass "
                  f"confidence >= {conf_threshold}")
        else:
            mask = pd.Series(True, index=adata.obs_names)

        sub = adata.obs.loc[mask, ["condition", col]].copy()
        sub[col] = _clean_label_series(sub[col])
        sub = sub.dropna()
        if sub.empty:
            print(f"  {col}: no cells left after confidence filtering -- skipping")
            continue

        ct = pd.crosstab(sub["condition"], sub[col], normalize="index") * 100
        ct = ct.reindex([c for c in order if c in ct.index])
        ct.to_csv(os.path.join(fig_dir, f"composition_{col}.csv"))

        n_types = ct.shape[1]
        cmap = matplotlib.colormaps["tab20"].resampled(max(n_types, 1))
        fig, ax = plt.subplots(figsize=(max(10, len(ct) * 0.7), 6))
        ct.plot(kind="bar", stacked=True, ax=ax,
                color=[cmap(i) for i in range(n_types)],
                edgecolor="black", linewidth=0.3)
        _embryo_cellline_divider(ax, [c for c in order if c in ct.index], adata.obs)
        ax.set_xlabel("Condition")
        ax.set_ylabel(f"% of cells (confidence >= {conf_threshold})")
        ax.set_title(f"Cell type composition by condition -- {col}")
        ax.legend(title=col, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
        plt.xticks(rotation=45, ha="right")
        _savefig(fig, os.path.join(fig_dir, f"composition_{col}.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# 2. Label diversity / heterogeneity collapse
# ─────────────────────────────────────────────────────────────────────────────

def analyze_diversity(adata, label_cols, fig_dir, conf_threshold):
    print("\n" + "=" * 70)
    print("2. LABEL DIVERSITY (Shannon entropy)")
    print("=" * 70)
    order = _condition_order(adata.obs)
    rows = []
    for col in label_cols:
        conf_col = f"{col}_confidence"
        mask = (adata.obs[conf_col].astype(float) >= conf_threshold
                if conf_col in adata.obs.columns else pd.Series(True, index=adata.obs_names))
        sub = adata.obs.loc[mask].copy()
        sub[col] = _clean_label_series(sub[col])
        for cond, g in sub.groupby("condition"):
            vc = g[col].dropna().value_counts(normalize=True)
            if vc.empty:
                continue
            H = float(scipy.stats.entropy(vc.values, base=2))
            rows.append(dict(label_col=col, condition=cond,
                              is_embryo=bool(g["is_embryo"].iloc[0]),
                              shannon_entropy=H, n_types=len(vc), n_cells=len(g)))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(fig_dir, "diversity_shannon_entropy.csv"), index=False)

    for col in label_cols:
        sub = df[df["label_col"] == col].set_index("condition")
        if sub.empty:
            continue
        present = [c for c in order if c in sub.index]
        sub = sub.reindex(present).dropna(subset=["shannon_entropy"])
        colors = ["#4CAF50" if e else "#FF7043" for e in sub["is_embryo"]]
        fig, ax = plt.subplots(figsize=(max(10, len(sub) * 0.6), 5))
        ax.bar(range(len(sub)), sub["shannon_entropy"].values, color=colors,
               edgecolor="black", linewidth=0.4)
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(sub.index, rotation=45, ha="right")
        ax.set_ylabel("Shannon entropy (bits)")
        ax.set_title(f"Cell type diversity by condition -- {col}\n"
                      "(green=embryo, orange=cell line)")
        _savefig(fig, os.path.join(fig_dir, f"diversity_{col}.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# 3. Confidence as a drift metric
# ─────────────────────────────────────────────────────────────────────────────

def analyze_confidence(adata, label_cols, fig_dir):
    print("\n" + "=" * 70)
    print("3. LABEL-TRANSFER CONFIDENCE BY CONDITION")
    print("=" * 70)
    order = _condition_order(adata.obs)
    for col in label_cols:
        conf_col = f"{col}_confidence"
        if conf_col not in adata.obs.columns:
            continue
        obs = adata.obs[["condition", "is_embryo", conf_col]].dropna()
        stats = obs.groupby("condition")[conf_col].agg(["mean", "median", "std", "count"])
        stats.to_csv(os.path.join(fig_dir, f"confidence_summary_{col}.csv"))

        present_order = [c for c in order if c in obs["condition"].unique()]
        fig, ax = plt.subplots(figsize=(max(10, len(present_order) * 0.6), 5))
        sns.boxplot(data=obs, x="condition", y=conf_col, hue="is_embryo",
                    order=present_order, ax=ax, showfliers=False)
        _embryo_cellline_divider(ax, present_order, adata.obs)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel(f"KNN confidence ({col})")
        ax.set_title(f"Label transfer confidence by condition -- {col}")
        plt.xticks(rotation=45, ha="right")
        _savefig(fig, os.path.join(fig_dir, f"confidence_{col}.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# 4. Pseudobulk correlation
# ─────────────────────────────────────────────────────────────────────────────

def analyze_pseudobulk(adata_full, label_cols, fig_dir, min_cells):
    print("\n" + "=" * 70)
    print("4. PSEUDOBULK CORRELATION")
    print("=" * 70)

    # 4a. condition x condition correlation, hierarchically clustered
    pb_cond = pseudobulk_by_group(adata_full, "condition", min_cells)
    if pb_cond.shape[0] >= 2:
        corr = pb_cond.T.corr(method="spearman")
        corr.to_csv(os.path.join(fig_dir, "pseudobulk_condition_correlation.csv"))
        g = sns.clustermap(corr, cmap="RdBu_r", center=0, figsize=(10, 10),
                            cbar_kws={"label": "Spearman rho"})
        g.savefig(os.path.join(fig_dir, "pseudobulk_condition_correlation.pdf"),
                  bbox_inches="tight", dpi=150)
        plt.close(g.fig)
        print(f"   Saved: {os.path.join(fig_dir, 'pseudobulk_condition_correlation.pdf')}")

    # 4b. condition pseudobulk vs. embryo-tissue-level pseudobulk (embryo cells only)
    embryo_mask = adata_full.obs["is_embryo"].astype(bool).values
    for col in label_cols:
        if col not in adata_full.obs.columns:
            continue
        embryo_adata = adata_full[embryo_mask].copy()
        embryo_adata.obs[col] = _clean_label_series(embryo_adata.obs[col])
        embryo_adata = embryo_adata[embryo_adata.obs[col].notna()]
        pb_tissue = pseudobulk_by_group(embryo_adata, col, min_cells)
        if pb_tissue.shape[0] < 2 or pb_cond.shape[0] < 1:
            print(f"   {col}: not enough groups for tissue pseudobulk comparison -- skipping")
            continue

        shared = pb_tissue.columns.intersection(pb_cond.columns)
        corr = pd.DataFrame(index=pb_cond.index, columns=pb_tissue.index, dtype=float)
        for cond in pb_cond.index:
            for tissue in pb_tissue.index:
                rho, _ = scipy.stats.spearmanr(pb_cond.loc[cond, shared], pb_tissue.loc[tissue, shared])
                corr.loc[cond, tissue] = rho
        corr = corr.astype(float)
        corr.to_csv(os.path.join(fig_dir, f"pseudobulk_vs_tissue_{col}.csv"))

        order = [c for c in _condition_order(adata_full.obs) if c in corr.index]
        corr = corr.reindex(order)
        fig, ax = plt.subplots(figsize=(max(8, corr.shape[1] * 0.8),
                                         max(6, corr.shape[0] * 0.4)))
        sns.heatmap(corr, cmap="viridis", annot=(corr.size <= 150), fmt=".2f",
                    linewidths=0.3, ax=ax, cbar_kws={"label": "Spearman rho"})
        ax.set_xlabel(f"Embryo tissue ({col}, embryo cells only)")
        ax.set_ylabel("Condition (pseudobulk)")
        ax.set_title(f"Which embryo tissue does each condition's pseudobulk "
                      f"resemble? -- {col}")
        _savefig(fig, os.path.join(fig_dir, f"pseudobulk_vs_tissue_{col}.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# 5. Marker-module scoring (identity retention / loss)
# ─────────────────────────────────────────────────────────────────────────────

def analyze_markers(adata_full, label_cols, fig_dir, min_cells, top_n):
    print("\n" + "=" * 70)
    print("5. MARKER-MODULE SCORING (identity retention)")
    print("=" * 70)

    for col in label_cols:
        if col not in adata_full.obs.columns:
            continue
        embryo_mask = adata_full.obs["is_embryo"].astype(bool).values
        sub = adata_full[embryo_mask].copy()
        sub.obs[col] = _clean_label_series(sub.obs[col])
        vc = sub.obs[col].value_counts()
        keep = vc[vc >= min_cells].index
        if len(keep) < 2:
            print(f"   {col}: fewer than 2 embryo tissue groups with >= "
                  f"{min_cells} cells -- skipping marker DE")
            continue
        sub = sub[sub.obs[col].isin(keep)].copy()
        sub.obs[col] = sub.obs[col].astype("category")

        print(f"   {col}: computing markers (wilcoxon) for {len(keep)} "
              f"embryo tissue groups ...")
        sc.tl.rank_genes_groups(sub, groupby=col, method="wilcoxon", use_raw=False)

        marker_sets, marker_rows = {}, []
        for group in sub.obs[col].cat.categories:
            df = sc.get.rank_genes_groups_df(sub, group=group)
            df = df[df["logfoldchanges"] > 0].head(top_n)
            if df.empty:
                continue
            marker_sets[group] = df["names"].tolist()
            df = df.copy()
            df.insert(0, "tissue", group)
            marker_rows.append(df)

        if not marker_rows:
            print(f"   {col}: no positive markers found -- skipping")
            continue
        pd.concat(marker_rows, ignore_index=True).to_csv(
            os.path.join(fig_dir, f"markers_{col}_top{top_n}.csv"), index=False)

        score_cols = []
        for tissue, genes in marker_sets.items():
            present = [g for g in genes if g in adata_full.var_names]
            if len(present) < 5:
                continue
            score_name = f"_score_{col}_{tissue}"
            sc.tl.score_genes(adata_full, gene_list=present, score_name=score_name, use_raw=False)
            score_cols.append(score_name)

        if not score_cols:
            continue
        mean_scores = adata_full.obs.groupby("condition")[score_cols].mean()
        mean_scores.columns = [c[len(f"_score_{col}_"):] for c in mean_scores.columns]
        order = [c for c in _condition_order(adata_full.obs) if c in mean_scores.index]
        mean_scores = mean_scores.reindex(order)
        mean_scores.to_csv(os.path.join(fig_dir, f"marker_scores_{col}.csv"))

        fig, ax = plt.subplots(figsize=(max(8, mean_scores.shape[1] * 0.8),
                                         max(6, mean_scores.shape[0] * 0.4)))
        sns.heatmap(mean_scores.astype(float), cmap="RdBu_r", center=0,
                    linewidths=0.3, ax=ax,
                    cbar_kws={"label": "Mean marker module score"})
        ax.set_xlabel(f"Marker module (top {top_n} genes per embryo tissue, {col})")
        ax.set_ylabel("Condition")
        ax.set_title(f"Marker retention/loss by condition -- {col}")
        _savefig(fig, os.path.join(fig_dir, f"marker_scores_{col}.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# 6. Cell cycle / proliferation shift
# ─────────────────────────────────────────────────────────────────────────────

def analyze_cell_cycle(adata, fig_dir):
    print("\n" + "=" * 70)
    print("6. CELL CYCLE / PROLIFERATION SHIFT")
    print("=" * 70)
    if "phase" not in adata.obs.columns:
        print("   No 'phase' column -- skipping (preprocess() may not have "
              "found enough cell cycle marker genes)")
        return

    order = _condition_order(adata.obs)
    present_order = [c for c in order if c in adata.obs["condition"].unique()]

    ct = pd.crosstab(adata.obs["condition"], adata.obs["phase"], normalize="index") * 100
    ct = ct.reindex(present_order)
    ct.to_csv(os.path.join(fig_dir, "cellcycle_phase_composition.csv"))

    fig, ax = plt.subplots(figsize=(max(10, len(ct) * 0.6), 6))
    ct.plot(kind="bar", stacked=True, ax=ax, edgecolor="black", linewidth=0.3,
            colormap="Set2")
    _embryo_cellline_divider(ax, present_order, adata.obs)
    ax.set_ylabel("% of cells")
    ax.set_title("Cell cycle phase composition by condition")
    ax.legend(title="Phase", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    _savefig(fig, os.path.join(fig_dir, "cellcycle_phase_composition.pdf"))

    adata.obs["_cycling"] = ~adata.obs["phase"].astype(str).str.upper().isin(["G1", "G0/G1"])
    cyc = (adata.obs.groupby("condition")["_cycling"].mean() * 100).reindex(present_order)
    cyc.to_csv(os.path.join(fig_dir, "cellcycle_cycling_fraction.csv"))
    fig, ax = plt.subplots(figsize=(max(10, len(cyc) * 0.6), 5))
    colors = ["#4CAF50" if adata.obs.loc[adata.obs["condition"] == c, "is_embryo"].iloc[0]
              else "#FF7043" for c in cyc.index]
    ax.bar(range(len(cyc)), cyc.values, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_xticks(range(len(cyc)))
    ax.set_xticklabels(cyc.index, rotation=45, ha="right")
    ax.set_ylabel("% cycling cells (S/G2M)")
    ax.set_title("Cycling fraction by condition (green=embryo, orange=cell line)")
    _savefig(fig, os.path.join(fig_dir, "cellcycle_cycling_fraction.pdf"))

    for score_col in ["S_score", "G2M_score"]:
        if score_col not in adata.obs.columns:
            continue
        fig, ax = plt.subplots(figsize=(max(10, len(present_order) * 0.6), 5))
        sns.boxplot(data=adata.obs, x="condition", y=score_col, hue="is_embryo",
                    order=present_order, ax=ax, showfliers=False)
        _embryo_cellline_divider(ax, present_order, adata.obs)
        ax.set_title(f"{score_col} by condition")
        plt.xticks(rotation=45, ha="right")
        _savefig(fig, os.path.join(fig_dir, f"cellcycle_{score_col}.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# 7. Wolbachia infection effects (cell line cells only)
# ─────────────────────────────────────────────────────────────────────────────

def analyze_wolbachia(adata, label_cols, fig_dir):
    print("\n" + "=" * 70)
    print("7. WOLBACHIA INFECTION EFFECTS")
    print("=" * 70)
    if "wolbachia_titer" not in adata.obs.columns:
        print("   No 'wolbachia_titer' column -- skipping")
        return

    cellline = adata.obs.loc[~adata.obs["is_embryo"].astype(bool)].copy()
    cellline["_infected"] = cellline["wolbachia_titer"].astype(float) > 0

    # 7a. titer vs. confidence, per label
    for col in label_cols:
        conf_col = f"{col}_confidence"
        if conf_col not in cellline.columns:
            continue
        sub = cellline[["wolbachia_titer", conf_col]].dropna()
        sub = sub[np.isfinite(sub["wolbachia_titer"])]
        if len(sub) < 20:
            continue
        rho, p = scipy.stats.spearmanr(sub["wolbachia_titer"], sub[conf_col])
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(sub["wolbachia_titer"], sub[conf_col], s=3, alpha=0.25, rasterized=True)
        ax.set_xlabel("Wolbachia titer")
        ax.set_ylabel(f"KNN confidence ({col})")
        ax.set_title(f"Titer vs. {col} confidence\nSpearman rho={rho:.3f}, p={p:.2e}")
        _savefig(fig, os.path.join(fig_dir, f"wolbachia_titer_vs_confidence_{col}.pdf"))
        pd.DataFrame({"spearman_rho": [rho], "p_value": [p], "n_cells": [len(sub)]}).to_csv(
            os.path.join(fig_dir, f"wolbachia_titer_vs_confidence_{col}.csv"), index=False)

    # 7b. infected vs. uninfected: composition per condition
    for col in label_cols:
        if col not in cellline.columns:
            continue
        ct = pd.crosstab([cellline["condition"], cellline["_infected"]],
                          cellline[col], normalize="index") * 100
        ct.to_csv(os.path.join(fig_dir, f"wolbachia_composition_by_infection_{col}.csv"))

    # 7c. infected vs. uninfected: cycling fraction per condition
    if "phase" in cellline.columns:
        cellline["_cycling"] = ~cellline["phase"].astype(str).str.upper().isin(["G1", "G0/G1"])
        rate = (cellline.groupby(["condition", "_infected"])["_cycling"].mean() * 100).unstack()
        rate.to_csv(os.path.join(fig_dir, "wolbachia_cycling_fraction_by_infection.csv"))
        if rate.shape[1] >= 1:
            fig, ax = plt.subplots(figsize=(max(8, len(rate) * 0.8), 5))
            rate.plot(kind="bar", ax=ax, edgecolor="black", linewidth=0.4)
            ax.set_ylabel("% cycling cells (S/G2M)")
            ax.set_title("Cycling fraction: infected vs. uninfected, by condition")
            ax.legend(title="Infected")
            plt.xticks(rotation=45, ha="right")
            _savefig(fig, os.path.join(fig_dir, "wolbachia_cycling_fraction_by_infection.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# 8. Host species check (Dmel vs. Dsim)
# ─────────────────────────────────────────────────────────────────────────────

def analyze_species(adata, label_cols, fig_dir):
    print("\n" + "=" * 70)
    print("8. HOST SPECIES CHECK (Dmel vs. Dsim)")
    print("=" * 70)
    adata.obs["_species"] = np.where(
        adata.obs["condition"].astype(str).str.contains("dsim", case=False), "Dsim", "Dmel")
    cellline = adata.obs.loc[~adata.obs["is_embryo"].astype(bool)].copy()

    for col in label_cols:
        conf_col = f"{col}_confidence"
        if conf_col not in cellline.columns:
            continue
        sub = cellline[["_species", conf_col]].dropna()
        if sub["_species"].nunique() < 2:
            continue
        fig, ax = plt.subplots(figsize=(5, 5))
        sns.boxplot(data=sub, x="_species", y=conf_col, ax=ax, showfliers=False)
        ax.set_ylim(0, 1.05)
        ax.set_title(f"Label-transfer confidence: Dmel vs. Dsim cell lines -- {col}")
        _savefig(fig, os.path.join(fig_dir, f"species_confidence_{col}.pdf"))
        sub.groupby("_species")[conf_col].agg(["mean", "median", "std", "count"]).to_csv(
            os.path.join(fig_dir, f"species_confidence_{col}.csv"))


# ─────────────────────────────────────────────────────────────────────────────
# 9. Leiden cluster composition
# ─────────────────────────────────────────────────────────────────────────────

def analyze_cluster_composition(adata, label_cols, fig_dir):
    print("\n" + "=" * 70)
    print("9. LEIDEN CLUSTER COMPOSITION")
    print("=" * 70)
    if "leiden" not in adata.obs.columns:
        print("   No 'leiden' column -- skipping")
        return

    ct = pd.crosstab(adata.obs["leiden"], adata.obs["is_embryo"], normalize="index") * 100
    ct.to_csv(os.path.join(fig_dir, "cluster_composition_is_embryo.csv"))
    fig, ax = plt.subplots(figsize=(max(8, len(ct) * 0.5), 5))
    ct.plot(kind="bar", stacked=True, ax=ax, edgecolor="black", linewidth=0.3,
            color=["#FF7043", "#4CAF50"])
    ax.set_xlabel("Leiden cluster")
    ax.set_ylabel("% of cells")
    ax.set_title("Embryo vs. cell line composition per cluster")
    ax.legend(title="is_embryo", bbox_to_anchor=(1.02, 1), loc="upper left")
    _savefig(fig, os.path.join(fig_dir, "cluster_composition_is_embryo.pdf"))

    for col in label_cols:
        if col not in adata.obs.columns:
            continue
        ct2 = pd.crosstab(adata.obs["leiden"], adata.obs[col], normalize="index") * 100
        ct2.to_csv(os.path.join(fig_dir, f"cluster_composition_{col}.csv"))
        fig, ax = plt.subplots(figsize=(max(8, ct2.shape[1] * 0.6),
                                         max(6, ct2.shape[0] * 0.4)))
        sns.heatmap(ct2, cmap="viridis", linewidths=0.3, ax=ax,
                    cbar_kws={"label": "% of cluster"})
        ax.set_xlabel(col)
        ax.set_ylabel("Leiden cluster")
        ax.set_title(f"Cluster composition -- {col}")
        _savefig(fig, os.path.join(fig_dir, f"cluster_composition_{col}.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# 10. UMAP overview
# ─────────────────────────────────────────────────────────────────────────────

def plot_umap_overview(adata, label_cols, fig_dir):
    print("\n" + "=" * 70)
    print("10. UMAP OVERVIEW")
    print("=" * 70)
    if "X_umap" not in adata.obsm:
        print("   No X_umap embedding -- skipping")
        return
    sc.settings.figdir = fig_dir

    quick_vars = [v for v in ["is_embryo", "leiden", "method", "wolbachia_titer", "phase"]
                  if v in adata.obs.columns]
    if quick_vars:
        sc.pl.umap(adata, color=quick_vars, ncols=3, save="_overview_summary.pdf", show=False)

    for col in ["condition"] + label_cols:
        if col in adata.obs.columns:
            sc.pl.umap(adata, color=col, save=f"_overview_{col}.pdf", show=False,
                       legend_fontsize=6)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Exploratory embryo -> cell line trajectory/identity "
                     "analysis on the integrated object."
    )
    parser.add_argument("--input", default="results/integrated/integrated.h5ad")
    parser.add_argument("--fig_dir", default="results/trajectory_analysis")
    parser.add_argument("--conf_threshold", type=float, default=0.5,
                         help="Minimum KNN confidence to count a cell type "
                              "call in composition/diversity analyses.")
    parser.add_argument("--min_cells", type=int, default=20,
                         help="Minimum cells per group for pseudobulk/marker analyses.")
    parser.add_argument("--top_n_markers", type=int, default=50)
    parser.add_argument("--skip_markers", action="store_true",
                         help="Skip section 5 (marker DE + scoring) -- the "
                              "slowest step.")
    parser.add_argument("--skip_pseudobulk", action="store_true")
    parser.add_argument("--skip_umap", action="store_true")

    args = parser.parse_args()
    os.makedirs(args.fig_dir, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    for required in ("condition", "is_embryo"):
        if required not in adata.obs.columns:
            raise ValueError(
                f"adata.obs is missing '{required}' -- is {args.input} really "
                "the output of integrate_v2.py's add_metadata()?"
            )

    label_cols = detect_label_cols(adata)
    print(f"  Detected cell type label columns: {label_cols}")
    if not label_cols:
        print("  WARNING: no cell type label columns found at all -- "
              "composition/diversity/marker analyses will be skipped.")

    n_embryo   = int(adata.obs["is_embryo"].astype(bool).sum())
    n_cellline = adata.n_obs - n_embryo
    print(f"  {n_embryo:,} embryo cells, {n_cellline:,} cell line cells")
    print(f"  Conditions: {_condition_order(adata.obs)}")

    if label_cols:
        analyze_composition(adata, label_cols, args.fig_dir, args.conf_threshold)
        analyze_diversity(adata, label_cols, args.fig_dir, args.conf_threshold)
        analyze_confidence(adata, label_cols, args.fig_dir)

    analyze_cell_cycle(adata, args.fig_dir)
    analyze_wolbachia(adata, label_cols, args.fig_dir)
    analyze_species(adata, label_cols, args.fig_dir)
    analyze_cluster_composition(adata, label_cols, args.fig_dir)

    if not args.skip_pseudobulk or not args.skip_markers:
        adata_full = get_raw_full_adata(adata)
        if not args.skip_pseudobulk and label_cols:
            analyze_pseudobulk(adata_full, label_cols, args.fig_dir, args.min_cells)
        if not args.skip_markers and label_cols:
            analyze_markers(adata_full, label_cols, args.fig_dir, args.min_cells,
                             args.top_n_markers)

    if not args.skip_umap:
        plot_umap_overview(adata, label_cols, args.fig_dir)

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"All plots + CSVs -> {args.fig_dir}/")
    n_pdf = len([f for f in os.listdir(args.fig_dir) if f.endswith(".pdf")])
    n_csv = len([f for f in os.listdir(args.fig_dir) if f.endswith(".csv")])
    print(f"{n_pdf} PDFs, {n_csv} CSVs written.")


if __name__ == "__main__":
    main()
