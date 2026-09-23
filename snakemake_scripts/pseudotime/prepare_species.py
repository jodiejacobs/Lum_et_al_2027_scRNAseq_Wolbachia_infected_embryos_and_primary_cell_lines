#!/usr/bin/env python3
"""
prepare_species.py
==================
Step 1 of the embryo -> primary cells -> immortalized cell line pseudotime
analysis (rule pseudotime_prepare). Run once per host species.

Why a new object instead of integrated.h5ad's embedding
-------------------------------------------------------
integrated.h5ad's X_pca_atlas / X_umap_atlas come from a frozen projection
onto the Flysta3D-v2 atlas. That embedding is fit on embryonic developmental
biology and is not built to resolve culture adaptation, which is the axis this
analysis is after. So each species gets its own PCA fit on its own cells,
using native host gene IDs (Dsim keeps its NCBI IDs, so no genes are lost to
ortholog remapping). integrated.h5ad supplies only per-cell metadata: the cell
whitelist, atlas labels + confidences, sample_type, titer.

Steps
-----
  1. Load raw host-gene counts (adata.raw) from each filtered_h5ad sample of
     this species; drop Wolbachia + 16S features.
  2. Keep cells present in integrated.h5ad and copy its obs columns across.
  3. Add species / lineage / stage_numeric (embryo=0, primary_cells=1,
     cell_culture=2).
  4. Root filter: keep embryo cells only if their atlas_annotation (at
     confidence >= --conf_threshold) makes up >= --root_min_frac of this
     species' primary + cell_culture cells. All non-embryo cells are kept.
  5. normalize_total + log1p -> HVG -> scale -> PCA (optional Harmony) ->
     neighbors -> UMAP -> leiden -> diffusion map -> DPT (rooted in embryo
     cells) -> PAGA (by condition and by leiden).

Outputs
-------
  <out_h5ad>                      : trajectory cells only; X = log1p-normalised
                                    host genes, layers['counts'] = raw counts
  <out_cells_csv>                 : EVERY cell of this species with
                                    pt_in_trajectory flag (used by merge step)
  <fig_dir>/root_filter_*.csv     : which embryo tissues were kept / dropped
  <fig_dir>/*.pdf                 : UMAPs, diffusion map, DPT, PAGA

Pattern follows Jacobs et al. 2026 (run_sceptic.py / integrate_v2.py), but
built for three discrete culture stages instead of an infection time course.
"""

import os
import re
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from pt_utils import savefig as _savefig

STAGE_ORDER = ["embryo", "primary_cells", "cell_culture"]
STAGE_NUM = {s: i for i, s in enumerate(STAGE_ORDER)}
STAGE_COLORS = {"embryo": "#4C72B0", "primary_cells": "#DD8452", "cell_culture": "#C44E52"}

# obs columns copied from integrated.h5ad (whatever subset exists)
OBS_KEEP = ["condition", "replicate", "method", "source_file", "sample_type",
            "wolbachia_titer", "n_counts", "n_genes", "percent_mito",
            "atlas_annotation", "atlas_annotation_confidence",
            "atlas_tissue", "atlas_tissue_confidence",
            "atlas_germ_layer", "atlas_germ_layer_confidence"]


def gtf_gene_ids(path):
    """Set of gene_id values in a GTF ('gene-' prefix stripped)."""
    ids = set()
    pat = re.compile(r'gene_id "([^"]+)"')
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            m = pat.search(line)
            if m:
                ids.add(re.sub(r"^gene-", "", m.group(1)))
    return ids


def host_gene_mask(var_names, host_ids, symbiont_ids):
    """True for host genes. Uses host GTF membership when it covers the
    var_names well; otherwise falls back to excluding symbiont + 16S IDs."""
    v = pd.Index(var_names)
    in_host = v.isin(host_ids) if host_ids else np.zeros(len(v), bool)
    frac = in_host.mean()
    print(f"  {in_host.sum()}/{len(v)} features match host GTF gene_ids ({frac:.1%})")
    if frac >= 0.5:
        return np.asarray(in_host)
    print("  WARNING: <50% match host GTF -- falling back to excluding "
          "symbiont GTF gene_ids and 16S_* features")
    return np.asarray(~(v.isin(symbiont_ids) | v.str.startswith("16S_")))


def load_species_counts(paths, whitelist, host_ids, symbiont_ids):
    adatas = []
    for p in paths:
        base = os.path.splitext(os.path.basename(p))[0]
        print(f"\nLoading {base}")
        a = sc.read_h5ad(p)
        if a.raw is None:
            raise ValueError(f"{p} has no .raw (expected raw counts)")
        X = a.raw.X
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        b = ad.AnnData(X=X.astype(np.float32), var=pd.DataFrame(index=a.raw.var_names))
        b.obs_names = [f"{base}__{bc}" for bc in a.obs_names]
        b = b[:, host_gene_mask(b.var_names, host_ids, symbiont_ids)].copy()
        keep = b.obs_names.isin(whitelist)
        print(f"  {keep.sum()}/{b.n_obs} cells in integrated.h5ad; {b.n_vars} host genes")
        adatas.append(b[keep].copy())
    out = ad.concat(adatas, join="outer", fill_value=0)
    out.X = sp.csr_matrix(out.X)
    return out


def root_filter(obs, label_col, conf_col, conf_thr, min_frac, fig_dir, species):
    """Embryo cells kept only if their (confident) tissue call is a type the
    cultures actually contain."""
    is_emb = obs["sample_type"] == "embryo"
    conf_ok = obs[conf_col] >= conf_thr if conf_col in obs else pd.Series(True, index=obs.index)
    cult = obs.loc[~is_emb & conf_ok, label_col].astype(str)
    frac_cult = cult.value_counts(normalize=True)
    keep_types = set(frac_cult[frac_cult >= min_frac].index)

    emb = obs.loc[is_emb, label_col].astype(str)
    summary = pd.DataFrame({
        "frac_of_culture_cells": frac_cult,
        "n_embryo_cells": emb.value_counts(),
        "n_embryo_confident": obs.loc[is_emb & conf_ok, label_col].astype(str).value_counts(),
    }).fillna(0)
    summary["kept_as_root_type"] = summary.index.isin(keep_types)
    summary = summary.sort_values("frac_of_culture_cells", ascending=False)
    summary.to_csv(os.path.join(fig_dir, f"root_filter_tissues_{species}.csv"))

    in_traj = ~is_emb | (is_emb & conf_ok & obs[label_col].astype(str).isin(keep_types))
    print(f"\nRoot filter ({label_col}, conf>={conf_thr}, >= {min_frac:.1%} of culture cells):")
    print(f"  kept types: {sorted(keep_types)}")
    print(f"  embryo cells kept: {int((in_traj & is_emb).sum())}/{int(is_emb.sum())}")
    return in_traj


def choose_root(adata):
    """Embryo cell at the embryo-side extreme of diffusion component 1."""
    dc1 = adata.obsm["X_diffmap"][:, 1]
    emb = (adata.obs["sample_type"] == "embryo").values
    if np.median(dc1[emb]) <= np.median(dc1[~emb]):
        idx = np.where(emb)[0][np.argmin(dc1[emb])]
    else:
        idx = np.where(emb)[0][np.argmax(dc1[emb])]
    return int(idx)


def plot_overview(adata, fig_dir, species):
    sc.settings.figdir = fig_dir
    for col in ["sample_type", "condition", "lineage", "atlas_annotation",
                "dpt_pseudotime", "wolbachia_titer", "leiden"]:
        if col in adata.obs:
            sc.pl.umap(adata, color=col, show=False, title=f"{species}: {col}",
                       save=f"_pt_{species}_{col}.pdf")

    fig, ax = plt.subplots(figsize=(6, 5))
    dm = adata.obsm["X_diffmap"]
    for st in STAGE_ORDER:
        m = (adata.obs["sample_type"] == st).values
        ax.scatter(dm[m, 1], dm[m, 2], s=2, alpha=0.4, c=STAGE_COLORS[st], label=st)
    r = adata.uns["iroot"]
    ax.scatter(dm[r, 1], dm[r, 2], s=80, marker="*", c="black", label="DPT root")
    ax.set_xlabel("DC1"); ax.set_ylabel("DC2"); ax.legend(markerscale=4)
    ax.set_title(f"{species}: diffusion map")
    _savefig(fig, os.path.join(fig_dir, f"diffmap_{species}.pdf"))

    order = (adata.obs[["condition", "stage_numeric"]].drop_duplicates()
             .sort_values(["stage_numeric", "condition"])["condition"].tolist())
    fig, ax = plt.subplots(figsize=(max(6, len(order) * 0.8), 4.5))
    sns.violinplot(data=adata.obs, x="condition", y="dpt_pseudotime", order=order,
                   hue="sample_type", hue_order=STAGE_ORDER, palette=STAGE_COLORS,
                   dodge=False, cut=0, ax=ax)
    ax.set_title(f"{species}: DPT by sample"); plt.xticks(rotation=45, ha="right")
    _savefig(fig, os.path.join(fig_dir, f"dpt_by_condition_{species}.pdf"))


def run_paga(adata, groups, fig_dir, species):
    sc.tl.paga(adata, groups=groups)
    con = adata.uns["paga"]["connectivities"].toarray()
    cats = adata.obs[groups].cat.categories
    pd.DataFrame(con, index=cats, columns=cats).to_csv(
        os.path.join(fig_dir, f"paga_connectivities_{groups}_{species}.csv"))
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.paga(adata, threshold=0.05, ax=ax, show=False, fontsize=7,
               title=f"{species}: PAGA by {groups}")
    _savefig(fig, os.path.join(fig_dir, f"paga_{groups}_{species}.pdf"))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--integrated", required=True)
    p.add_argument("--filtered", required=True, nargs="+")
    p.add_argument("--species", required=True)
    p.add_argument("--lineages", required=True, help="TSV: condition<TAB>lineage")
    p.add_argument("--host_gtf", default=None)
    p.add_argument("--symbiont_gtfs", nargs="*", default=[])
    p.add_argument("--root_label_col", default="atlas_annotation")
    p.add_argument("--conf_threshold", type=float, default=0.5)
    p.add_argument("--root_min_frac", type=float, default=0.01)
    p.add_argument("--n_top_genes", type=int, default=2000)
    p.add_argument("--n_pcs", type=int, default=30)
    p.add_argument("--n_neighbors", type=int, default=30)
    p.add_argument("--leiden_res", type=float, default=1.0)
    p.add_argument("--harmony_key", default=None,
                   help="obs column for Harmony; omit for no batch correction")
    p.add_argument("--out_h5ad", required=True)
    p.add_argument("--out_cells_csv", required=True)
    p.add_argument("--fig_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.fig_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_h5ad) or ".", exist_ok=True)

    print(f"=== Pseudotime prep: {args.species} ===")
    integ = ad.read_h5ad(args.integrated, backed="r")
    iobs = integ.obs[[c for c in OBS_KEEP if c in integ.obs.columns]].copy()
    integ.file.close()

    host_ids = gtf_gene_ids(args.host_gtf) if args.host_gtf else set()
    sym_ids = set().union(*[gtf_gene_ids(g) for g in args.symbiont_gtfs]) if args.symbiont_gtfs else set()

    adata = load_species_counts(args.filtered, set(iobs.index), host_ids, sym_ids)
    adata.obs = iobs.loc[adata.obs_names].copy()

    lin = pd.read_csv(args.lineages, sep="\t", index_col=0).iloc[:, 0].to_dict()
    adata.obs["species"] = args.species
    adata.obs["lineage"] = adata.obs["condition"].map(lin).fillna("unassigned")
    adata.obs["sample_type"] = adata.obs["sample_type"].astype(str)
    unknown = set(adata.obs["sample_type"]) - set(STAGE_ORDER)
    if unknown:
        raise ValueError(f"Unknown sample_type values: {unknown}")
    adata.obs["stage_numeric"] = adata.obs["sample_type"].map(STAGE_NUM).astype(int)
    print("\nCells per condition:")
    print(adata.obs.groupby(["sample_type", "lineage", "condition"], observed=True)
          .size().to_string())

    in_traj = root_filter(adata.obs, args.root_label_col,
                          f"{args.root_label_col}_confidence", args.conf_threshold,
                          args.root_min_frac, args.fig_dir, args.species)
    adata.obs["pt_in_trajectory"] = in_traj.values
    adata.obs[["species", "lineage", "condition", "sample_type", "stage_numeric",
               "pt_in_trajectory"]].to_csv(args.out_cells_csv)
    adata = adata[adata.obs["pt_in_trajectory"]].copy()
    n_stage = adata.obs["sample_type"].value_counts()
    if n_stage.get("embryo", 0) == 0 or n_stage.size < 2:
        raise ValueError(
            f"Need embryo cells plus at least one later stage after the root filter; "
            f"got {n_stage.to_dict()}. Check --filtered includes embryo samples, or "
            "lower --root_min_frac / --conf_threshold.")

    # ── Preprocess ───────────────────────────────────────────────────────────
    sc.pp.filter_genes(adata, min_cells=10)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, flavor="seurat")
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=max(50, args.n_pcs))
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.uns["pca"] = hvg.uns["pca"]
    embed = "X_pca"
    if args.harmony_key:
        sc.external.pp.harmony_integrate(adata, args.harmony_key, basis="X_pca",
                                         adjusted_basis="X_pca_harmony")
        embed = "X_pca_harmony"
    adata.uns["pt_embedding"] = embed
    print(f"\nEmbedding for neighbors/SCEPTIC: {embed}")

    sc.pp.neighbors(adata, use_rep=embed, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.leiden_res)
    sc.tl.diffmap(adata, n_comps=15)
    adata.uns["iroot"] = choose_root(adata)
    sc.tl.dpt(adata)
    print(f"DPT root: {adata.obs_names[adata.uns['iroot']]}")
    n_inf = int((~np.isfinite(adata.obs["dpt_pseudotime"])).sum())
    if n_inf:
        print(f"  WARNING: {n_inf} cells have infinite DPT (graph component not "
              "connected to the root); they are excluded from DPT comparisons")

    for col in ["condition", "lineage", "sample_type"]:
        adata.obs[col] = pd.Categorical(adata.obs[col])
    adata.obs["sample_type"] = adata.obs["sample_type"].cat.reorder_categories(
        [s for s in STAGE_ORDER if s in adata.obs["sample_type"].cat.categories])

    plot_overview(adata, args.fig_dir, args.species)
    run_paga(adata, "condition", args.fig_dir, args.species)
    run_paga(adata, "leiden", args.fig_dir, args.species)

    adata.write(args.out_h5ad)
    print(f"\nWrote {args.out_h5ad} ({adata.n_obs} cells x {adata.n_vars} genes)")


if __name__ == "__main__":
    main()
