#!/usr/bin/env python3
"""
primary_proliferative.py
========================
Rule pseudotime_primary_proliferative. Cluster each primary cell line on its
own, find the proliferative population (the candidate founders of the cell
line), and describe it against the rest of the primary cell line.

Per primary cell line (one per lineage, from prepared_<lineage>.h5ad):
  1. Re-embed the primary-cell-line cells alone: within-sample HVGs (stress /
     ribosomal / mito genes excluded, as in prepare_species.py), PCA, kNN,
     Leiden (--resolution), UMAP. Score proliferation (cell_states.MODULES)
     and S / G2M phase (scanpy score_genes_cell_cycle with fly genes).
  2. Proliferative population, two definitions:
       cluster : Leiden clusters whose mean proliferation z-score (within this
                 primary cell line) is >= --cluster_min_z
       cell    : cells with proliferation z >= --cell_min_z (cross-check;
                 matches cell_states.py's 'proliferating' flag)
     The cluster definition is the primary one; if no cluster passes, the
     cell definition is used and the table says so.
  3. What separates them (proliferative vs rest of the same primary cell line;
     markers.csv flags cell-cycle genes, which are expected by construction):
       qc_comparison.csv      : UMIs, genes detected, % mito, doublet score,
                                medians + Mann-Whitney p
       markers.csv            : Wilcoxon (scanpy), all genes, with symbols
       state_mix.csv          : marker-based cell states (cell_states.py modules)
       similarity_to_line.csv : Spearman rho of the population's pseudobulk
                                with the cell line's pseudobulk, vs the rest
  4. Wolbachia load (wolbachia_titer from integrated.h5ad):
       wolbachia.csv : fraction infected (titer > 0), median titer, Mann-Whitney
                       p, and an OLS coefficient for proliferative status with
                       log UMI depth as covariate (titer depends on depth)
  5. Embryonic origin, two methods:
       origin_atlas.csv : Flysta3D atlas_annotation (confidence >= --conf) of
                          proliferative vs rest (chi-square)
       origin_embryo.csv: kNN transfer from the lineage's own embryo cells in a
                          PCA fit on the embryo, using only genes that do NOT
                          change between embryo and primary cell line
                          (|pooled log2FC| < --stable_lfc or padj >= 0.05), so
                          the culture injury/immune program cannot drive the
                          match; label = majority of k embryo neighbours,
                          confidence = vote fraction
  6. GSEA (gsea/): prerank on the Wilcoxon z-scores per lineage and combined
     across lineages (Stouffer), each also without cell-cycle genes;
     gsea_<key>.csv, gsea_nes_heatmap.pdf, combined_ranking.csv
Across lineages (out_dir root):
  summary.csv, consistent_markers.csv (genes up in the proliferative
  population in >= n-1 lineages), wolbachia_summary.csv, origin_summary.csv,
  and summary figures.
"""
import os
import re
import argparse
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, chi2_contingency, spearmanr
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm

from pt_utils import savefig
from cell_states import MODULES, IDENTITY, score_modules

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

EXCLUDE = (r"^(?:Rp[LS]\d|RpLP|mRp[LS]|mt:|Hsp\d|Hsc70|Hsromega)"
           r"|^(?:kay|Jra|puc|Hr38|sr|Ets21C|Thor)$")
S_GENES = ["PCNA", "Mcm2", "Mcm3", "Mcm5", "Mcm6", "Mcm7", "Mcm10", "RnrL", "RnrS", "dup",
           "Orc1", "Cdc6", "E2f1", "CycE", "DNApol-alpha50", "DNApol-alpha73", "Fen1", "RPA2"]
G2M_GENES = ["CycB", "CycA", "Cdk1", "polo", "aurA", "aurB", "stg", "pav", "feo", "Incenp",
             "BubR1", "cmet", "mad2", "Det", "tum", "sti", "Klp61F", "asp", "Map205", "CycB3"]


def sym_map(adata):
    return pd.Series(adata.var["symbol"].astype(str).values, index=adata.var_names)


def ids_for(adata, symbols):
    s = sym_map(adata)
    return [g for g in adata.var_names[s.isin(symbols).values]]


def zscore(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / (x.std() or 1.0)


def mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def embed_primary(prim, resolution, seed):
    sym = sym_map(prim)
    keep = ~sym.str.contains(EXCLUDE, regex=True).values
    X = prim.copy()
    sc.pp.highly_variable_genes(X, n_top_genes=2000, flavor="seurat")
    X.var.loc[~keep, "highly_variable"] = False
    try:
        sc.pp.pca(X, n_comps=30, mask_var="highly_variable", random_state=seed)
    except TypeError:  # scanpy < 1.10
        sc.pp.pca(X, n_comps=30, use_highly_variable=True, random_state=seed)
    sc.pp.neighbors(X, n_neighbors=15, random_state=seed)
    sc.tl.leiden(X, resolution=resolution, random_state=seed, key_added="pl_leiden")
    sc.tl.umap(X, random_state=seed)
    return X


def pseudobulk(adata, mask):
    C = adata.layers["counts"] if "counts" in adata.layers else adata.X
    v = np.asarray(C[mask].sum(axis=0)).ravel()
    return np.log1p(1e6 * v / max(v.sum(), 1))


def embryo_origin(adata, prim_mask, emb_mask, stable_genes, k, conf, seed):
    """kNN label transfer from the lineage's embryo cells, stage-stable genes only."""
    lab = adata.obs["atlas_annotation"].astype(str)
    if "atlas_annotation_confidence" in adata.obs:
        ok = adata.obs["atlas_annotation_confidence"].astype(float) >= conf
    else:
        ok = pd.Series(True, index=adata.obs_names)
    ref = adata[emb_mask & ok.values]
    if ref.n_obs < 50:
        return None
    genes = [g for g in stable_genes if g in adata.var_names]
    R = ref[:, genes].copy()
    sc.pp.highly_variable_genes(R, n_top_genes=min(1500, len(genes)), flavor="seurat")
    hv = R.var_names[R.var["highly_variable"]]
    def dense(a):
        x = a[:, hv].X
        return x.toarray() if sp.issparse(x) else np.asarray(x)
    Xr = dense(ref)
    mu, sd = Xr.mean(0), Xr.std(0) + 1e-6
    Zr = (Xr - mu) / sd
    U, S, Vt = np.linalg.svd(Zr - Zr.mean(0), full_matrices=False)
    P = Vt[:30].T
    Er = Zr @ P
    Eq = ((dense(adata[prim_mask]) - mu) / sd) @ P
    nn = NearestNeighbors(n_neighbors=k).fit(Er)
    _, idx = nn.kneighbors(Eq)
    labs = ref.obs["atlas_annotation"].astype(str).values[idx]
    best, frac = [], []
    for row in labs:
        v, c = np.unique(row, return_counts=True)
        best.append(v[c.argmax()]); frac.append(c.max() / k)
    return pd.DataFrame({"embryo_origin": best, "embryo_origin_conf": frac},
                        index=adata.obs_names[prim_mask]), len(hv)


def compare_labels(labels, prol, min_frac=0.0):
    t = pd.crosstab(labels, prol)
    t.columns = ["rest" if not c else "proliferative" for c in t.columns]
    frac = t / t.sum()
    try:
        p = chi2_contingency(t.values)[1] if t.shape[0] > 1 and t.shape[1] == 2 else np.nan
    except ValueError:
        p = np.nan
    out = pd.concat([t.add_prefix("n_"), frac.add_prefix("frac_")], axis=1)
    if {"frac_proliferative", "frac_rest"} <= set(out):
        out["log2_enrichment"] = np.log2((out["frac_proliferative"] + 0.005) / (out["frac_rest"] + 0.005))
    return out.sort_values("n_proliferative" if "n_proliferative" in out else out.columns[0],
                           ascending=False), p


def analyze(path, args, de_dir):
    adata = sc.read_h5ad(path)
    if "symbol" not in adata.var:
        adata.var["symbol"] = adata.var_names
    name = os.path.basename(path).removesuffix(".h5ad").removeprefix("prepared_")
    species = str(adata.obs["species"].iloc[0])
    out = os.path.join(args.out_dir, name)
    os.makedirs(out, exist_ok=True)
    st = adata.obs["sample_type"].astype(str)
    prim_mask = (st == "primary_cells").values
    if prim_mask.sum() < 100:
        print(f"  {name}: <100 primary cells, skipped"); return None
    print(f"\n=== {name}: {prim_mask.sum()} primary-cell-line cells ===")

    prim = adata[prim_mask].copy()
    score_modules(prim, {k: MODULES[k] for k in ["proliferation"] + IDENTITY
                         + ["immune_AMP", "injury_JAK_JNK", "apoptosis"] if k in MODULES})
    s_ids, g_ids = ids_for(prim, S_GENES), ids_for(prim, G2M_GENES)
    if len(s_ids) >= 3 and len(g_ids) >= 3:
        sc.tl.score_genes_cell_cycle(prim, s_genes=s_ids, g2m_genes=g_ids, random_state=0)
    prim = embed_primary(prim, args.resolution, args.seed)
    prim.obs["prolif_z"] = zscore(prim.obs["score_proliferation"])

    # proliferative population
    cz = prim.obs.groupby("pl_leiden", observed=True)["prolif_z"].mean().sort_values(ascending=False)
    clusters = list(cz.index[cz >= args.cluster_min_z])
    cell_def = (prim.obs["prolif_z"] >= args.cell_min_z).values
    if clusters:
        prol = prim.obs["pl_leiden"].isin(clusters).values
        definition = "cluster"
    else:
        prol, definition = cell_def, "cell (no cluster passed)"
    prim.obs["proliferative"] = prol
    print(f"  clusters by mean proliferation z: {cz.round(2).to_dict()}")
    print(f"  proliferative ({definition}): {prol.sum()} cells "
          f"({100 * prol.mean():.1f}%); cell-level definition: {cell_def.sum()} cells; "
          f"overlap {int((prol & cell_def).sum())}")
    cz.rename("mean_prolif_z").to_csv(os.path.join(out, "cluster_proliferation.csv"))

    # 3a. QC
    obs = prim.obs
    C = prim.layers["counts"] if "counts" in prim.layers else prim.X
    obs["genes_detected"] = np.asarray((C > 0).sum(axis=1)).ravel()
    obs["umis"] = np.asarray(C.sum(axis=1)).ravel()
    qrows = []
    for c in ["umis", "genes_detected", "percent_mito", "doublet_score", "prolif_z",
              "S_score", "G2M_score"]:
        if c in obs:
            a, b = obs.loc[prol, c].astype(float), obs.loc[~prol, c].astype(float)
            ratio = (a.median() / b.median() if c in ("umis", "genes_detected", "percent_mito")
                     and b.median() else np.nan)
            qrows.append(dict(metric=c, median_proliferative=a.median(), median_rest=b.median(),
                              difference=a.median() - b.median(), ratio=ratio, p=mwu(a, b)))
    if "phase" in obs:
        ph = pd.crosstab(obs["phase"], prol, normalize="columns")
        ph.columns = ["rest" if not c else "proliferative" for c in ph.columns]
        ph.to_csv(os.path.join(out, "phase.csv"))
    qc = pd.DataFrame(qrows)
    qc.to_csv(os.path.join(out, "qc_comparison.csv"), index=False)
    print(qc.round(3).to_string(index=False))

    # 3b. markers
    prim.obs["group"] = np.where(prol, "proliferative", "rest")
    sc.tl.rank_genes_groups(prim, "group", groups=["proliferative"], reference="rest",
                            method="wilcoxon", use_raw=False, pts=True)
    mk = sc.get.rank_genes_groups_df(prim, "proliferative")
    mk["symbol"] = mk["names"].map(sym_map(prim))
    cc = set(S_GENES) | set(G2M_GENES) | set(MODULES.get("proliferation", []))
    mk["cell_cycle_gene"] = mk["symbol"].isin(cc) | mk["symbol"].str.contains(EXCLUDE, regex=True)
    mk.to_csv(os.path.join(out, "markers.csv"), index=False)
    up = mk[(mk["pvals_adj"] < 0.05) & (mk["logfoldchanges"] > 0.5)]
    print(f"  markers: {len(up)} up (padj<0.05, log2FC>0.5); top non-cell-cycle: "
          + ", ".join(up.loc[~up["cell_cycle_gene"], "symbol"].head(25).astype(str)))

    # 3c. state mix
    zc = [f"score_{m}" for m in IDENTITY if f"score_{m}" in obs]
    if zc:
        z = (obs[zc] - obs[zc].mean()) / obs[zc].std(ddof=0)
        state = z.idxmax(axis=1).str.replace("score_", "", regex=False)
        state[z.max(axis=1) <= 0.5] = "unassigned"
        sm_, _ = compare_labels(state, pd.Series(prol, index=obs.index))
        sm_.to_csv(os.path.join(out, "state_mix.csv"))

    # 3d. similarity to the cell line
    line_mask = (st == "cell_culture").values
    sim = {}
    if line_mask.sum() >= 50:
        hv = prim.var_names[prim.var["highly_variable"]]
        gi = adata.var_names.get_indexer(hv)
        full_prim = adata[prim_mask]
        pl = pseudobulk(adata[line_mask], np.ones(line_mask.sum(), bool))[gi]
        pp = pseudobulk(full_prim, prol)[gi]
        pr = pseudobulk(full_prim, ~prol)[gi]
        sim = dict(rho_proliferative_line=spearmanr(pp, pl).correlation,
                   rho_rest_line=spearmanr(pr, pl).correlation)
        # null: random sets of the same size
        rng = np.random.default_rng(args.seed)
        null = [spearmanr(pseudobulk(full_prim, np.isin(np.arange(prim.n_obs),
                                                        rng.choice(prim.n_obs, prol.sum(), replace=False)))[gi],
                          pl).correlation for _ in range(args.n_perm)]
        sim["rho_random_mean"] = float(np.mean(null))
        sim["p_perm"] = (1 + sum(n >= sim["rho_proliferative_line"] for n in null)) / (1 + args.n_perm)
        pd.DataFrame([sim]).to_csv(os.path.join(out, "similarity_to_line.csv"), index=False)
        print(f"  similarity to cell line: {sim}")

    # 4. Wolbachia
    wol = {}
    if "wolbachia_titer" in obs:
        t = obs["wolbachia_titer"].astype(float)
        a, b = t[prol], t[~prol]
        wol = dict(frac_infected_proliferative=float((a > 0).mean()), frac_infected_rest=float((b > 0).mean()),
                   median_titer_proliferative=float(a.median()), median_titer_rest=float(b.median()),
                   p_mwu=mwu(a, b))
        d = pd.DataFrame({"titer": t, "prol": prol.astype(float),
                          "log_umis": np.log10(obs["umis"].clip(lower=1))}).dropna()
        if len(d) > 20:
            fit = sm.OLS(d["titer"], sm.add_constant(d[["prol", "log_umis"]])).fit(cov_type="HC3")
            wol.update(ols_coef_proliferative=float(fit.params["prol"]),
                       ols_p_proliferative=float(fit.pvalues["prol"]))
        pd.DataFrame([wol]).to_csv(os.path.join(out, "wolbachia.csv"), index=False)
        print(f"  Wolbachia: {wol}")

    # 5. origin
    orig = {}
    if "atlas_annotation" in obs:
        lab = obs["atlas_annotation"].astype(str)
        if "atlas_annotation_confidence" in obs:
            lab = lab.where(obs["atlas_annotation_confidence"].astype(float) >= args.conf, "low_confidence")
        oa, pa = compare_labels(lab, pd.Series(prol, index=obs.index))
        oa.to_csv(os.path.join(out, "origin_atlas.csv"))
        orig["atlas_top_proliferative"] = oa.index[0]
        orig["atlas_chi2_p"] = pa
    de_p = os.path.join(de_dir, f"de_{species}_primary_vs_embryo.csv")
    emb_mask = (st == "embryo").values
    if os.path.exists(de_p) and emb_mask.sum() >= 50 and "atlas_annotation" in adata.obs:
        de = pd.read_csv(de_p, index_col=0)
        changed = de.index[(de["padj"] < 0.05) & (de["log2FoldChange"].abs() >= args.stable_lfc)]
        sym = sym_map(adata)
        stable = [g for g in adata.var_names
                  if g not in set(changed) and not re.search(EXCLUDE, sym[g])]
        res = embryo_origin(adata, prim_mask, emb_mask, stable, args.k, args.conf, args.seed)
        if res is not None:
            eo, n_hv = res
            obs = obs.join(eo)
            lab = eo["embryo_origin"].where(eo["embryo_origin_conf"] >= 0.5, "ambiguous")
            oe, pe = compare_labels(lab, pd.Series(prol, index=obs.index))
            oe.to_csv(os.path.join(out, "origin_embryo.csv"))
            orig.update(embryo_top_proliferative=oe.index[0], embryo_chi2_p=pe,
                        n_stable_genes=len(stable), n_stable_hvg=n_hv)
            print(f"  embryo-origin transfer on {len(stable)} stage-stable genes "
                  f"({n_hv} HVGs):\n" + oe.head(8).round(3).to_string())
    obs.to_csv(os.path.join(out, "cells.csv.gz"))

    # figures
    prim.obs = obs.reindex(prim.obs_names).combine_first(prim.obs)
    cols = [c for c in ["pl_leiden", "prolif_z", "proliferative", "phase", "wolbachia_titer",
                        "atlas_annotation", "embryo_origin"] if c in prim.obs]
    prim.obs["proliferative"] = prim.obs["proliferative"].astype(str)
    fig = sc.pl.umap(prim, color=cols, ncols=4, show=False, return_fig=True, size=8)
    savefig(fig, os.path.join(out, "umap_primary.pdf"))
    if "wolbachia_titer" in obs:
        fig, ax = plt.subplots(figsize=(3.5, 3.5))
        sns.violinplot(x=np.where(prol, "proliferative", "rest"), y=obs["wolbachia_titer"].astype(float),
                       ax=ax, cut=0, inner="quartile")
        ax.set_ylabel("Wolbachia titer"); ax.set_title(name, fontsize=9)
        savefig(fig, os.path.join(out, "wolbachia_violin.pdf"))

    row = dict(lineage=name, species=species, n_primary=int(prim_mask.sum()),
               definition=definition, clusters=";".join(map(str, clusters)),
               n_proliferative=int(prol.sum()), frac_proliferative=float(prol.mean()),
               n_cell_definition=int(cell_def.sum()), n_markers_up=len(up))
    row.update({f"qc_{r['metric']}_" + ("ratio" if r["metric"] in ("umis", "genes_detected", "percent_mito")
                else "diff"): (r["ratio"] if r["metric"] in ("umis", "genes_detected", "percent_mito")
                               else r["difference"]) for r in qrows})
    row["n_markers_up_non_cell_cycle"] = int((~up["cell_cycle_gene"]).sum())
    row.update(sim); row.update(wol); row.update(orig)
    return row, mk.assign(lineage=name)


def run_gsea(m, args):
    """GSEA prerank on the Wilcoxon z-scores (proliferative vs rest): per
    lineage and combined across lineages (Stouffer: sum of z / sqrt(n)), each
    with and without cell-cycle genes (which are expected by construction)."""
    from de_stages import load_gene_sets, gsea_all
    gene_sets = load_gene_sets(args.gmt, args.gene_set_libraries, args.gene_set_organism)
    if not gene_sets:
        print("  GSEA skipped: no gene sets loaded"); return
    gdir = os.path.join(args.out_dir, "gsea")
    os.makedirs(gdir, exist_ok=True)
    m = m.dropna(subset=["symbol", "scores"])
    z = m.pivot_table(index="symbol", columns="lineage", values="scores", aggfunc="first")
    cc = m.groupby("symbol")["cell_cycle_gene"].first().reindex(z.index).fillna(False).astype(bool)
    n = z.notna().sum(axis=1)
    comb = (z.sum(axis=1) / np.sqrt(n))[n >= max(2, z.shape[1] - 1)]
    pd.DataFrame({"stouffer_z": comb, "n_lineages": n.reindex(comb.index),
                  "cell_cycle_gene": cc.reindex(comb.index)}).sort_values(
        "stouffer_z", ascending=False).to_csv(os.path.join(gdir, "combined_ranking.csv"))
    results = {}
    for lin in z.columns:
        s_ = z[lin].dropna()
        results[f"{lin}"] = pd.DataFrame({"symbol": s_.index, "stat": s_.values})
        s2 = s_[~cc.reindex(s_.index).values]
        results[f"{lin}_no_cell_cycle"] = pd.DataFrame({"symbol": s2.index, "stat": s2.values})
    results["combined"] = pd.DataFrame({"symbol": comb.index, "stat": comb.values})
    c2 = comb[~cc.reindex(comb.index).values]
    results["combined_no_cell_cycle"] = pd.DataFrame({"symbol": c2.index, "stat": c2.values})
    gsea_all(results, gene_sets, gdir, args.gsea_permutations)
    for key in ["combined", "combined_no_cell_cycle"]:
        f = os.path.join(gdir, f"gsea_{key}.csv")
        if os.path.exists(f):
            g = pd.read_csv(f)
            g = g[g["FDR q-val"] < 0.05].sort_values("NES", ascending=False)
            print(f"\nGSEA {key}: {len(g)} terms FDR < 0.05; top up:\n"
                  + g.head(12)[["Term", "NES", "FDR q-val"]].round(3).to_string(index=False)
                  + "\ntop down:\n"
                  + g.tail(8)[["Term", "NES", "FDR q-val"]].round(3).to_string(index=False))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ads", nargs="+", required=True)
    p.add_argument("--de_dir", required=True)
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--cluster_min_z", type=float, default=1.0)
    p.add_argument("--cell_min_z", type=float, default=1.0)
    p.add_argument("--stable_lfc", type=float, default=1.0)
    p.add_argument("--conf", type=float, default=0.5)
    p.add_argument("--k", type=int, default=15)
    p.add_argument("--n_perm", type=int, default=200)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--gene_set_libraries", nargs="*", default=["GO_Biological_Process_2018"])
    p.add_argument("--gene_set_organism", default="Fly")
    p.add_argument("--gmt", default=None)
    p.add_argument("--gsea_permutations", type=int, default=1000)
    p.add_argument("--skip_gsea", action="store_true")
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    rows, marks = [], []
    for h in args.h5ads:
        r = analyze(h, args, args.de_dir)
        if r:
            rows.append(r[0]); marks.append(r[1])
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(args.out_dir, "summary.csv"), index=False)
    print("\n" + summ.T.to_string())
    if marks:
        m = pd.concat(marks)
        up = m[(m["pvals_adj"] < 0.05) & (m["logfoldchanges"] > 0.5)]
        cons = up.groupby("symbol")["lineage"].nunique().sort_values(ascending=False)
        is_cc = m.groupby("symbol")["cell_cycle_gene"].first()
        need = max(2, len(rows) - 1)
        cons = cons[cons >= need].rename("n_lineages_up").to_frame()
        cons["cell_cycle_gene"] = is_cc.reindex(cons.index).values
        lfc = m.pivot_table(index="symbol", columns="lineage", values="logfoldchanges", aggfunc="first")
        cons.join(lfc).to_csv(os.path.join(args.out_dir, "consistent_markers.csv"))
        print(f"\n{len(cons)} genes up in the proliferative population in >= {need} lineages "
              f"({int((~cons['cell_cycle_gene'].astype(bool)).sum())} not cell-cycle):\n"
              + ", ".join(cons.index[~cons["cell_cycle_gene"].astype(bool)][:60]))
    wcols = ["lineage"] + [c for c in summ if c.startswith(("frac_infected", "median_titer", "p_mwu", "ols_"))]
    summ[wcols].to_csv(os.path.join(args.out_dir, "wolbachia_summary.csv"), index=False)
    ocols = ["lineage"] + [c for c in summ if c.startswith(("atlas_", "embryo_", "n_stable"))]
    summ[ocols].to_csv(os.path.join(args.out_dir, "origin_summary.csv"), index=False)
    if marks and not args.skip_gsea:
        run_gsea(pd.concat(marks), args)
    print("Done ->", args.out_dir)


if __name__ == "__main__":
    main()
