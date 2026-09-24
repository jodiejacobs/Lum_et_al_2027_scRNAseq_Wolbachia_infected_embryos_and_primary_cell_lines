#!/usr/bin/env python3
"""
de_stages.py
============
Rule pseudotime_de. Discrete, pseudobulk DE counterpart to the pseudotime
analysis: embryo -> primary cells -> cell line, with each sample (library)
as the unit and the matched lineages as replicates. Uses PyDESeq2.

Input: the prepare_species.py objects (prepared_<trajectory>.h5ad), so the
cells are the same trajectory cells as the pseudotime (embryo cells limited
to atlas types present in culture) and every trajectory is in Dmel FlyBase
gene space (Dsim remapped to 1:1 orthologs).

Pseudobulk: raw counts (layers['counts']) summed per sample; samples with
< --min_cells cells are dropped.

Models (stage levels: embryo, primary, line)
  per species   ~ lineage + stage          (2 lineages per species)
  pooled        ~ lineage + stage          (all 4 lineages; genes detected in
                                            both species)
  interaction   ~ lineage + stage + dsim_primary + dsim_line
                (dsim_* = 1 for Dsim samples at that stage, else 0; their
                coefficients are the species x stage interaction: how much
                the Dsim stage change differs from the Dmel one)

Contrasts: primary_vs_embryo, line_vs_primary, line_vs_embryo (per species,
pooled) and the matching interaction terms.

Also
  - per-lineage consistency: CPM log2FC of each transition within every
    lineage; for pooled-significant genes, how many lineages agree in sign
    (+ sign test), heatmap of the top genes
  - preranked GSEA on every contrast (rank = Wald stat)

Caveat: embryo and primary stages have one library per lineage, and stage is
confounded with sequencing run, so stage effects include batch. n is small
(2 lineages per species); treat per-species results as lower-powered than
the pooled model.
"""

import os
import argparse
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binomtest

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from pt_utils import savefig, load_flybase_symbols

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

STAGE_MAP = {"embryo": "embryo", "primary_cells": "primary", "cell_culture": "line"}
CONTRASTS = [("primary_vs_embryo", "primary", "embryo"),
             ("line_vs_primary", "line", "primary"),
             ("line_vs_embryo", "line", "embryo")]


# ─────────────────────────────────────────────────────────────────────────────
# Pseudobulk
# ─────────────────────────────────────────────────────────────────────────────

def pseudobulk(paths, min_cells):
    counts, meta, genes_by_species = [], [], {}
    for p in paths:
        a = ad.read_h5ad(p)
        X = a.layers["counts"]
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        sp_ = str(a.obs["species"].iloc[0])
        genes_by_species.setdefault(sp_, set()).update(a.var_names)
        for s, idx in a.obs.groupby("source_file", observed=True).indices.items():
            o = a.obs.iloc[idx[0]]
            counts.append(pd.Series(np.asarray(X[idx].sum(axis=0)).ravel(),
                                    index=a.var_names, name=s))
            meta.append(dict(sample=s, condition=o["condition"], lineage=o["lineage"],
                             species=sp_, stage=STAGE_MAP[str(o["sample_type"])],
                             n_cells=len(idx)))
        print(f"  {os.path.basename(p)}: {a.obs['source_file'].nunique()} samples")
    cm = pd.concat(counts, axis=1).fillna(0).T.round().astype(int)
    md = pd.DataFrame(meta).set_index("sample")
    drop = md.index[md["n_cells"] < min_cells]
    if len(drop):
        print(f"  WARNING: dropping samples with < {min_cells} cells: {list(drop)}")
    md = md.drop(drop)
    return cm.loc[md.index], md, genes_by_species


def filter_genes(cm, min_count=10, min_samples=2):
    return cm.loc[:, (cm >= min_count).sum(axis=0) >= min_samples]


# ─────────────────────────────────────────────────────────────────────────────
# DESeq2
# ─────────────────────────────────────────────────────────────────────────────

def fit(cm, md, design, n_cpus):
    md = md.copy()
    for c in ["lineage", "stage", "species"]:
        md[c] = md[c].astype(str)
    dds = DeseqDataSet(counts=cm, metadata=md, design=design,
                       refit_cooks=True, n_cpus=n_cpus, quiet=True)
    dds.deseq2()
    return dds


def run_contrast(dds, contrast, sym, n_cpus):
    ds = DeseqStats(dds, contrast=contrast, n_cpus=n_cpus, quiet=True)
    ds.summary()
    res = ds.results_df.copy()
    res.insert(0, "symbol", [sym.get(g, g) for g in res.index])
    return res.sort_values("pvalue")


def coef_vector(dds, col):
    cols = list(dds.obsm["design_matrix"].columns)
    v = np.zeros(len(cols))
    for c, w in col.items():
        v[cols.index(c)] = w
    return v


def volcano(res, title, path, alpha):
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    r = res.dropna(subset=["padj"])
    sig = r["padj"] < alpha
    ax.scatter(r.loc[~sig, "log2FoldChange"], -np.log10(r.loc[~sig, "pvalue"]), s=3,
               c="#bdc3c7")
    ax.scatter(r.loc[sig, "log2FoldChange"], -np.log10(r.loc[sig, "pvalue"]), s=4,
               c=np.where(r.loc[sig, "log2FoldChange"] > 0, "#c0392b", "#2980b9"))
    for g, row in r[sig].nsmallest(15, "pvalue").iterrows():
        ax.annotate(row["symbol"], (row["log2FoldChange"], -np.log10(row["pvalue"])),
                    fontsize=6)
    ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10 p")
    ax.set_title(f"{title}\n{int(sig.sum())} genes padj < {alpha}", fontsize=9)
    savefig(fig, path)


# ─────────────────────────────────────────────────────────────────────────────
# Per-lineage consistency
# ─────────────────────────────────────────────────────────────────────────────

def lineage_consistency(cm, md, pooled_res, sym, out, alpha, top_n=50):
    cpm = cm.div(cm.sum(axis=1), axis=0) * 1e6
    stage_mean = np.log2(cpm + 1).groupby([md["lineage"], md["stage"]]).mean()
    rows = []
    for name, b, a in CONTRASTS:
        lfc = {}
        for lin in md["lineage"].unique():
            if (lin, a) in stage_mean.index and (lin, b) in stage_mean.index:
                lfc[lin] = stage_mean.loc[(lin, b)] - stage_mean.loc[(lin, a)]
        lfc = pd.DataFrame(lfc)
        res = pooled_res[name]
        sig = res.index[(res["padj"] < alpha)].intersection(lfc.index)
        pooled_sign = np.sign(res.loc[sig, "log2FoldChange"])
        agree = (np.sign(lfc.loc[sig]).eq(pooled_sign, axis=0)).sum(axis=1)
        n_lin = lfc.shape[1]
        tab = lfc.loc[sig].copy()
        tab.columns = [f"log2FC_{c}" for c in tab.columns]
        tab.insert(0, "symbol", [sym.get(g, g) for g in sig])
        tab["pooled_log2FC"] = res.loc[sig, "log2FoldChange"]
        tab["pooled_padj"] = res.loc[sig, "padj"]
        tab["n_lineages_agree"] = agree
        tab["sign_test_p"] = [binomtest(int(k), n_lin, 0.5, alternative="greater").pvalue
                              for k in agree]
        tab.sort_values(["n_lineages_agree", "pooled_padj"], ascending=[False, True]).to_csv(
            os.path.join(out, f"lineage_consistency_{name}.csv"))
        rows.append(dict(contrast=name, n_pooled_sig=len(sig), n_lineages=n_lin,
                         all_agree=int((agree == n_lin).sum()),
                         frac_all_agree=float((agree == n_lin).mean()) if len(sig) else np.nan))

        top = res.loc[sig].nsmallest(top_n, "padj").index
        if len(top) >= 5:
            h = lfc.loc[top]
            h.index = [sym.get(g, g) for g in top]
            fig, ax = plt.subplots(figsize=(0.8 * n_lin + 3, 0.18 * len(top) + 1.5))
            lim = np.nanpercentile(np.abs(h.values), 98)
            sns.heatmap(h, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, ax=ax,
                        yticklabels=True, cbar_kws={"label": "log2FC (CPM)"})
            ax.tick_params(axis="y", labelsize=5)
            ax.set_title(f"{name}: top {len(top)} pooled DE genes per lineage", fontsize=9)
            savefig(fig, os.path.join(out, f"lineage_consistency_{name}.pdf"))
    pd.DataFrame(rows).to_csv(os.path.join(out, "lineage_consistency_summary.csv"), index=False)


# ─────────────────────────────────────────────────────────────────────────────
# GSEA
# ─────────────────────────────────────────────────────────────────────────────

def load_gene_sets(gmt, libs, organism):
    import gseapy as gp
    if gmt:
        return {os.path.basename(gmt): gp.read_gmt(gmt)}
    out = {}
    for lib in libs:
        try:
            out[lib] = gp.get_library(name=lib, organism=organism)
        except Exception as e:
            print(f"  WARNING: could not fetch {lib} ({e}); pass --gmt for offline use")
    return out


def gsea_all(results, gene_sets, out, permutations, top_terms=25):
    import gseapy as gp
    nes = {}
    for key, res in results.items():
        rnk = res.dropna(subset=["stat"]).set_index("symbol")["stat"]
        rnk = rnk[~rnk.index.duplicated()].sort_values(ascending=False)
        frames = []
        for lib, gs in gene_sets.items():
            r = gp.prerank(rnk=rnk, gene_sets=gs, min_size=10, max_size=500,
                           permutation_num=permutations, seed=42, threads=4,
                           outdir=None, verbose=False).res2d
            r["library"] = lib
            frames.append(r)
        r = pd.concat(frames)
        r["NES"] = pd.to_numeric(r["NES"], errors="coerce")
        r["FDR q-val"] = pd.to_numeric(r["FDR q-val"], errors="coerce")
        r.sort_values("FDR q-val").to_csv(os.path.join(out, f"gsea_{key}.csv"), index=False)
        nes[key] = r.set_index("Term")["NES"].groupby(level=0).first()
        sig = r.loc[r["FDR q-val"] < 0.25, "Term"]
        nes[key + "__sig"] = pd.Series(True, index=sig.unique())
    keys = list(results)
    mat = pd.DataFrame({k: nes[k] for k in keys})
    sig_any = pd.concat([nes[k + "__sig"] for k in keys]).index.unique()
    mat = mat.loc[mat.index.intersection(sig_any)]
    if len(mat):
        top = mat.abs().max(axis=1).nlargest(top_terms * 2).index
        m = mat.loc[top].fillna(0)
        fig, ax = plt.subplots(figsize=(0.7 * len(keys) + 6, 0.22 * len(top) + 2))
        sns.heatmap(m, cmap="RdBu_r", center=0, ax=ax, yticklabels=[t[:60] for t in m.index],
                    cbar_kws={"label": "NES"})
        ax.tick_params(axis="y", labelsize=6); plt.xticks(rotation=45, ha="right")
        ax.set_title("GSEA NES across DE contrasts (terms FDR < 0.25 in any)", fontsize=9)
        savefig(fig, os.path.join(out, "gsea_nes_heatmap.pdf"))
        mat.to_csv(os.path.join(out, "gsea_nes_matrix.csv"))


# ─────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ads", nargs="+", required=True, help="prepared_<trajectory>.h5ad files")
    p.add_argument("--flybase_annotation", default=None)
    p.add_argument("--min_cells", type=int, default=30)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--n_cpus", type=int, default=4)
    p.add_argument("--gene_set_libraries", nargs="*", default=["GO_Biological_Process_2018"])
    p.add_argument("--gene_set_organism", default="Fly")
    p.add_argument("--gmt", default=None)
    p.add_argument("--gsea_permutations", type=int, default=1000)
    p.add_argument("--skip_gsea", action="store_true")
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    out = args.out_dir
    os.makedirs(out, exist_ok=True)
    sym = load_flybase_symbols(args.flybase_annotation)

    print("Building pseudobulk")
    cm, md, genes_by_species = pseudobulk(args.h5ads, args.min_cells)
    md.to_csv(os.path.join(out, "pseudobulk_samples.csv"))
    cm.T.to_csv(os.path.join(out, "pseudobulk_counts.csv.gz"))
    print(md.to_string())

    results, summary = {}, []

    def record(key, res):
        results[key] = res
        res.to_csv(os.path.join(out, f"de_{key}.csv"))
        volcano(res, key, os.path.join(out, f"volcano_{key}.pdf"), args.alpha)
        s = res["padj"] < args.alpha
        summary.append(dict(contrast=key, n_tested=int(res["padj"].notna().sum()),
                            n_sig=int(s.sum()),
                            n_up=int((s & (res["log2FoldChange"] > 0)).sum()),
                            n_down=int((s & (res["log2FoldChange"] < 0)).sum())))
        print(f"  {key}: {summary[-1]['n_sig']} genes padj < {args.alpha}")

    # per species
    for sp_ in sorted(md["species"].unique()):
        m = md[md["species"] == sp_]
        c = filter_genes(cm.loc[m.index])
        design = "~lineage + stage" if m["lineage"].nunique() > 1 else "~stage"
        print(f"\n[{sp_}] {len(m)} samples x {c.shape[1]} genes, design {design}")
        dds = fit(c, m, design, args.n_cpus)
        for name, b, a in CONTRASTS:
            record(f"{sp_}_{name}", run_contrast(dds, ["stage", b, a], sym, args.n_cpus))

    # pooled + interaction (genes present in every species)
    if md["species"].nunique() > 1:
        shared = sorted(set.intersection(*genes_by_species.values()))
        c = filter_genes(cm[shared])
        print(f"\n[pooled] {len(md)} samples x {c.shape[1]} shared genes, design ~lineage + stage")
        dds = fit(c, md, "~lineage + stage", args.n_cpus)
        pooled = {}
        for name, b, a in CONTRASTS:
            pooled[name] = run_contrast(dds, ["stage", b, a], sym, args.n_cpus)
            record(f"pooled_{name}", pooled[name])

        other = [s for s in sorted(md["species"].unique()) if s != "Dmel"][0]
        mi = md.copy()
        mi["dsim_primary"] = ((mi["species"] == other) & (mi["stage"] == "primary")).astype(float)
        mi["dsim_line"] = ((mi["species"] == other) & (mi["stage"] == "line")).astype(float)
        design = "~lineage + stage + dsim_primary + dsim_line"
        print(f"\n[interaction] {design} ({other} vs Dmel)")
        dds = fit(c, mi, design, args.n_cpus)
        cols = list(dds.obsm["design_matrix"].columns)
        print(f"  design columns: {cols}")
        ip = [x for x in cols if "dsim_primary" in x][0]
        il = [x for x in cols if "dsim_line" in x][0]
        for name, vec in [("primary_vs_embryo", {ip: 1}),
                          ("line_vs_primary", {il: 1, ip: -1}),
                          ("line_vs_embryo", {il: 1})]:
            record(f"interaction_{other}_minus_Dmel_{name}",
                   run_contrast(dds, coef_vector(dds, vec), sym, args.n_cpus))

        print("\nPer-lineage consistency")
        lineage_consistency(cm[c.columns], md, pooled, sym, out, args.alpha)

    pd.DataFrame(summary).to_csv(os.path.join(out, "de_summary.csv"), index=False)
    print("\n" + pd.DataFrame(summary).to_string(index=False))

    if not args.skip_gsea:
        gs = load_gene_sets(args.gmt, args.gene_set_libraries, args.gene_set_organism)
        if gs:
            print("\nGSEA")
            gsea_all(results, gs, out, args.gsea_permutations)
    print("\nDone ->", out)


if __name__ == "__main__":
    main()
