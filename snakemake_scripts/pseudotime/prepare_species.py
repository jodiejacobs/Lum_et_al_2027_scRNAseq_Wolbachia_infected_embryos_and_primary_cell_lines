#!/usr/bin/env python3
"""
prepare_species.py
==================
Step 1 of the embryo -> primary cells -> immortalized cell line pseudotime
analysis (rule pseudotime_prepare). Run once per trajectory (one matched
embryo -> primary cells -> cell line lineage; --group), or per species.

Why a new object instead of integrated.h5ad's embedding
-------------------------------------------------------
integrated.h5ad's X_pca_atlas / X_umap_atlas come from a frozen projection
onto the Flysta3D-v2 atlas. That embedding is fit on embryonic developmental
biology and is not built to resolve culture adaptation, which is the axis this
analysis is after. So each species gets its own PCA fit on its own cells,
using Dmel FlyBase gene IDs. Species passed --ortholog_map (Dsim) are remapped
to 1:1 RBH Dmel orthologs first, then filtered against the Dmel GTF, so the
annotation matches integrated.h5ad and Dmel. integrated.h5ad supplies only
per-cell metadata: the cell
whitelist, atlas labels + confidences, sample_type, titer.

Steps
-----
  1. Load raw counts (adata.raw) from each filtered_h5ad sample of this
     species; optionally remap to 1:1 Dmel orthologs; keep host genes only
     (drops Wolbachia + 16S features).
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

from pt_utils import savefig as _savefig, load_orthologs, remap_to_dmel, load_flybase_symbols

# Genes left out of the embedding (HVG/PCA) because they mostly track library
# prep and dissociation rather than cell state: cytosolic + mito ribosomal
# proteins, mitochondrial genome, heat shock, and fly immediate-early /
# dissociation-stress genes. Matched on FlyBase symbols. They stay in the
# object (and in tradeSeq); override with --exclude_gene_regex.
DEFAULT_EXCLUDE_REGEX = (r"^(?:Rp[LS]\d|RpLP|mRp[LS]|mt:|Hsp\d|Hsc70|Hsromega)"
                         r"|^(?:kay|Jra|puc|Hr38|sr|Ets21C|Thor)$")

STAGE_ORDER = ["embryo", "primary_cells", "cell_culture"]
STAGE_NUM = {s: i for i, s in enumerate(STAGE_ORDER)}
STAGE_COLORS = {"embryo": "#4C72B0", "primary_cells": "#DD8452", "cell_culture": "#C44E52"}

# obs columns copied from integrated.h5ad (whatever subset exists)
OBS_KEEP = ["condition", "replicate", "method", "source_file", "sample_type",
            "wolbachia_titer", "n_counts", "n_genes", "percent_mito", "doublet_score",
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


def load_species_counts(paths, whitelist, host_ids, symbiont_ids, to_dmel=None):
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
        if to_dmel:
            b = remap_to_dmel(b, to_dmel, label=base)
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


def plot_qc(adata, fig_dir, species):
    """Is the stage separation partly technical? QC metrics on the UMAP and a
    per-sample QC table."""
    sc.settings.figdir = fig_dir
    cols = [c for c in ["n_counts", "n_genes", "percent_mito", "doublet_score"]
            if c in adata.obs]
    if not cols:
        return
    sc.pl.umap(adata, color=cols, show=False, ncols=2, cmap="viridis",
               save=f"_pt_{species}_qc.pdf")
    tab = (adata.obs.groupby(["sample_type", "condition"], observed=True)[cols]
           .median().round(3))
    tab.insert(0, "n_cells", adata.obs.groupby(["sample_type", "condition"],
                                               observed=True).size())
    tab.to_csv(os.path.join(fig_dir, f"qc_by_sample_{species}.csv"))
    print("\nMedian QC per sample:\n" + tab.to_string())


def knn_stage_mixing(adata, fig_dir, species):
    """Fraction of each cell's kNN neighbours from each stage, averaged per
    sample. If primary cells never neighbour embryo or cell-line cells, the
    stages are disconnected and no continuous trajectory exists in this
    embedding."""
    conn = adata.obsp["distances"].tocsr()
    stage = adata.obs["sample_type"].astype(str).values
    stages = [s for s in STAGE_ORDER if s in set(stage)]
    rows, cols_ = conn.nonzero()
    nb = pd.DataFrame({"cell": rows, "nb_stage": stage[cols_]})
    frac = (nb.groupby(["cell", "nb_stage"]).size().unstack(fill_value=0)
            .reindex(columns=stages, fill_value=0))
    frac = frac.div(frac.sum(axis=1), axis=0)
    frac["condition"] = adata.obs["condition"].astype(str).values[frac.index]
    frac["sample_type"] = stage[frac.index]
    out = frac.groupby(["sample_type", "condition"])[stages].mean().round(4)
    out.columns = [f"frac_neighbours_{c}" for c in stages]
    out.to_csv(os.path.join(fig_dir, f"knn_stage_mixing_{species}.csv"))
    print("\nkNN neighbours by stage (mean fraction per sample):\n" + out.to_string())


def run_paga(adata, groups, fig_dir, species):
    try:
        sc.tl.paga(adata, groups=groups)
    except ValueError as e:
        # igraph fails when no kNN edges connect the groups at all
        print(f"  WARNING: PAGA by {groups} failed ({e}); groups are disconnected "
              "in the kNN graph -- skipping")
        return
    con = adata.uns["paga"]["connectivities"].toarray()
    cats = adata.obs[groups].cat.categories
    pd.DataFrame(con, index=cats, columns=cats).to_csv(
        os.path.join(fig_dir, f"paga_connectivities_{groups}_{species}.csv"))
    # connectivity heatmap instead of sc.pl.paga's graph drawing, which breaks
    # with newer scipy sparse arrays (TypeError: sparse array length is ambiguous)
    n = len(cats)
    fig, ax = plt.subplots(figsize=(0.5 * n + 3, 0.45 * n + 2.5))
    sns.heatmap(pd.DataFrame(con, index=cats, columns=cats), cmap="viridis",
                vmin=0, vmax=1, annot=n <= 12, fmt=".2f", annot_kws={"size": 7},
                square=True, cbar_kws={"label": "PAGA connectivity"}, ax=ax)
    ax.set_title(f"{species}: PAGA connectivity by {groups}")
    plt.xticks(rotation=45, ha="right")
    _savefig(fig, os.path.join(fig_dir, f"paga_{groups}_{species}.pdf"))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--integrated", required=True)
    p.add_argument("--filtered", required=True, nargs="+")
    p.add_argument("--species", required=True, help="host species (Dmel / Dsim)")
    p.add_argument("--group", default=None,
                   help="trajectory name used for files/titles (e.g. a lineage); default = species")
    p.add_argument("--lineages", required=True, help="TSV: condition<TAB>lineage")
    p.add_argument("--host_gtf", default=None,
                   help="GTF of the gene space used (Dmel GTF when --ortholog_map is set)")
    p.add_argument("--ortholog_map", default=None,
                   help="RBH table (Dsim, Dmel columns); remaps this species to Dmel FBgn")
    p.add_argument("--symbiont_gtfs", nargs="*", default=[])
    p.add_argument("--root_label_col", default="atlas_annotation")
    p.add_argument("--conf_threshold", type=float, default=0.5)
    p.add_argument("--root_min_frac", type=float, default=0.01)
    p.add_argument("--n_top_genes", type=int, default=2000)
    p.add_argument("--n_pcs", type=int, default=30)
    p.add_argument("--n_neighbors", type=int, default=30)
    p.add_argument("--leiden_res", type=float, default=1.0)
    p.add_argument("--hvg_batch_key", default="source_file",
                   help="select HVGs within each value of this obs column (per sample), "
                        "so genes that only differ BETWEEN samples/stages don't drive "
                        "the embedding; '' = pooled HVGs")
    p.add_argument("--exclude_gene_regex", default=DEFAULT_EXCLUDE_REGEX,
                   help="FlyBase-symbol regex of genes kept out of HVG/PCA; '' = none")
    p.add_argument("--flybase_annotation", default=None,
                   help="fbgn_annotation_ID*.tsv.gz, for symbol-based gene exclusion")
    p.add_argument("--harmony_key", default=None,
                   help="obs column for Harmony; omit for no batch correction")
    p.add_argument("--out_h5ad", required=True)
    p.add_argument("--out_cells_csv", required=True)
    p.add_argument("--fig_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.fig_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_h5ad) or ".", exist_ok=True)

    name = args.group or args.species
    print(f"=== Pseudotime prep: {name} ({args.species}) ===")
    integ = ad.read_h5ad(args.integrated, backed="r")
    iobs = integ.obs[[c for c in OBS_KEEP if c in integ.obs.columns]].copy()
    integ.file.close()

    host_ids = gtf_gene_ids(args.host_gtf) if args.host_gtf else set()
    # symbiont GTFs are only used as a fallback (host GTF match < 50%), so a
    # missing one is a warning, not an error
    sym_ids = set()
    for g in args.symbiont_gtfs:
        if os.path.exists(g):
            sym_ids |= gtf_gene_ids(g)
        else:
            print(f"  WARNING: symbiont GTF not found, skipping: {g}")

    to_dmel = load_orthologs(args.ortholog_map) if args.ortholog_map else None
    adata = load_species_counts(args.filtered, set(iobs.index), host_ids, sym_ids, to_dmel)
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
                          args.root_min_frac, args.fig_dir, name)
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
    # genes excluded from the embedding (see DEFAULT_EXCLUDE_REGEX)
    sym = load_flybase_symbols(args.flybase_annotation)
    adata.var["symbol"] = [s if isinstance(s, str) else g
                           for g, s in ((g, sym.get(g)) for g in adata.var_names)]
    excl = (adata.var["symbol"].str.contains(args.exclude_gene_regex, regex=True, na=False)
            if args.exclude_gene_regex else pd.Series(False, index=adata.var_names))
    adata.var["pt_excluded"] = excl.values
    print(f"\nExcluded from HVG/PCA: {int(excl.sum())} genes "
          f"(e.g. {', '.join(adata.var.loc[excl, 'symbol'].head(8))})")
    cand = adata[:, ~adata.var["pt_excluded"]].copy()
    bk = args.hvg_batch_key if args.hvg_batch_key else None
    sc.pp.highly_variable_genes(cand, n_top_genes=args.n_top_genes, flavor="seurat",
                                batch_key=bk)
    print(f"HVGs: {int(cand.var['highly_variable'].sum())} "
          f"({'within each ' + bk if bk else 'pooled'})")
    adata.var["highly_variable"] = adata.var_names.isin(
        cand.var_names[cand.var["highly_variable"]])
    adata.uns["pt_hvg_batch_key"] = bk or ""
    del cand
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

    plot_overview(adata, args.fig_dir, name)
    plot_qc(adata, args.fig_dir, name)
    knn_stage_mixing(adata, args.fig_dir, name)
    run_paga(adata, "condition", args.fig_dir, name)
    run_paga(adata, "leiden", args.fig_dir, name)

    adata.write(args.out_h5ad)
    print(f"\nWrote {args.out_h5ad} ({adata.n_obs} cells x {adata.n_vars} genes)")


if __name__ == "__main__":
    main()
