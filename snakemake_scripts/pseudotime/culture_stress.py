#!/usr/bin/env python3
"""
culture_stress.py
=================
Rule pseudotime_culture_stress. Could Wolbachia stabilize early primary
culture by limiting the injury / immune response to dissociation, by buffering
oxidative stress, or by provisioning nutrients (iron, heme)?

Gene programs are fixed in advance (MODULES_STRESS below; genes not found are
listed in modules_used.csv). OXPHOS is added from the regex used elsewhere.

Predictions of a stabilizing effect, stated before running:
  P1 culture level  : lineages whose primary culture carries more Wolbachia show
                      a smaller embryo -> primary rise in injury_JAK_JNK,
                      immune_AMP, redox_NRF2 and nutrient_starvation.
  P2 cell level     : within a culture, cells with more Wolbachia express LESS
                      of these programs (if the effect is cell-autonomous).
                      A population-wide effect predicts no within-culture
                      association, so P2 can only support, not refute, P1.
  P3 infected line  : Dsim6B-wMel expresses less of these programs than the
                      doxycycline-cured Dsim6B (same origin).

Outputs (out_dir):
  modules_used.csv          genes found per module
  shift_embryo_to_primary.csv  per lineage x module: median score in embryo and
                            primary, standardized shift (Cohen's d) with a
                            bootstrap 95% CI over cells, and the pseudobulk
                            module log2 fold change (primary vs embryo)
  culture_level.csv         per module: Spearman rho across the 4 lineages
                            between primary Wolbachia load (and load retained
                            from the embryo) and the shift; exact permutation p
                            over all 24 orderings (smallest attainable 0.083).
                            Species, strain and load are confounded: descriptive.
  within_culture_modules.csv / within_culture_genes.csv
                            per lineage x condition (median Wolbachia UMIs >=
                            --min_wol): partial Spearman of log1p(Wolbachia UMIs)
                            with each module score / gene | log host UMIs; BH
  infected_vs_cured.csv     per module: mean log2FC Dsim6B-wMel vs Dsim6B, genes
                            up / down (padj < 0.05), Mann-Whitney of module gene
                            Wald stats vs all other genes
  culture_stress.pdf        shift vs load per module; within-culture rho heatmap
"""
import os
import re
import itertools
import math
import argparse
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, spearmanr, mannwhitneyu, t as tdist
from statsmodels.stats.multitest import multipletests

from pt_utils import savefig
from cell_states import MODULES, score_modules
from primary_proliferative import OXPHOS_REGEX

warnings.filterwarnings("ignore", category=FutureWarning)
STAGES = {"embryo": "embryo", "primary_cells": "primary", "cell_culture": "line"}
MODULES_STRESS = {
    "injury_JAK_JNK": MODULES["injury_JAK_JNK"],
    "immune_AMP": MODULES["immune_AMP"],
    "redox_NRF2": ["cnc", "Keap1", "Gclc", "Gclm", "GstD1", "GstD2", "GstD5", "GstE1", "GstE6",
                   "GstE7", "Mgst1", "Jafrac1", "Jafrac2", "Prx3", "Prx5", "Trxr-1", "Sod1",
                   "Sod2", "Cat", "Cyp6a2", "Cyp6a8"],
    "iron_heme": ["Fer1HCH", "Fer2LCH", "Fer3HCH", "Tsf1", "Tsf2", "Tsf3", "Mvl", "Irp-1A",
                  "Irp-1B", "Alas", "Pbgs", "Ppox", "FeCH", "Coprox", "Ho"],
    "nutrient_starvation": ["Thor", "InR", "Pepck", "Pepck1", "bmm", "Lip3", "sug"],
}
STRESS = ["injury_JAK_JNK", "immune_AMP", "redox_NRF2", "nutrient_starvation"]


def partial_spearman(y, x, z):
    r = lambda v: pd.Series(v).rank().values
    Z = np.column_stack([np.ones(len(y)), r(z)])
    res = lambda v: v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
    a, b = res(r(y)), res(r(x))
    if a.std() == 0 or b.std() == 0:
        return np.nan, np.nan
    rho = float(np.corrcoef(a, b)[0, 1])
    df = len(y) - 3
    tt = rho * np.sqrt(df / max(1 - rho ** 2, 1e-12))
    return rho, 2 * tdist.sf(abs(tt), df)


def cohens_d(a, b, n_boot, rng):
    def d(x, y):
        s = np.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / (len(x) + len(y) - 2))
        return (y.mean() - x.mean()) / s if s > 0 else np.nan
    est = d(a, b)
    boots = [d(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(n_boot)]
    return est, np.nanpercentile(boots, 2.5), np.nanpercentile(boots, 97.5)


def exact_spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    obs = spearmanr(x, y).correlation
    null = [spearmanr(x, np.array(p)).correlation for p in itertools.permutations(y)]
    return obs, float(np.mean(np.abs(null) >= abs(obs) - 1e-12))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad_dir", required=True)
    p.add_argument("--lineages", nargs="+", required=True)
    p.add_argument("--infection_dir", required=True)
    p.add_argument("--de_dir", required=True, help="results/pseudotime/de (pseudobulk counts)")
    p.add_argument("--infection_de", required=True, help="de_Dsim6B-wMel_vs_Dsim6B.csv")
    p.add_argument("--min_wol", type=float, default=5)
    p.add_argument("--n_boot", type=int, default=500)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out_dir", required=True)
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    pd.set_option("display.width", 250)

    wol = pd.read_csv(os.path.join(a.infection_dir, "cells_infection_immune.csv.gz"), index_col=0,
                      usecols=lambda c: c in ("Unnamed: 0", "lineage", "sample_type", "wolbachia_umis",
                                              "wolbachia_gene_frac"))
    shifts, within_m, within_g, used_rows = [], [], [], []
    for lin in a.lineages:
        print(f"\n=== {lin}")
        A = sc.read_h5ad(os.path.join(a.h5ad_dir, lin, f"prepared_{lin}.h5ad"))
        A.obs["stage"] = A.obs["sample_type"].astype(str).map(STAGES)
        A = A[A.obs["stage"].notna()].copy()
        sym = A.var["symbol"].astype(str)
        mods = dict(MODULES_STRESS)
        mods["oxphos"] = [s for s in sym if re.search(OXPHOS_REGEX, s)]
        used = score_modules(A, mods)
        for m, g in used.items():
            used_rows.append(dict(lineage=lin, module=m, n_genes=len(g),
                                  genes=",".join(g) if m != "oxphos" else f"{len(g)} OXPHOS genes"))
        obs = A.obs[["stage", "condition"]].astype(str).copy()
        for m in used:
            obs[m] = A.obs[f"score_{m}"].values
        C = A.layers["counts"] if "counts" in A.layers else None
        obs["host_umis"] = np.asarray(C.sum(1)).ravel() if C is not None else np.nan
        panel = sorted({g for m in MODULES_STRESS.values() for g in m} & set(sym))
        s2i = {s: i for i, s in reversed(list(enumerate(sym)))}
        X = A.X
        for g in panel:
            x = X[:, s2i[g]]
            obs[f"g_{g}"] = np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()

        # shift from embryo to primary (cell scores)
        e, pr = obs[obs.stage == "embryo"], obs[obs.stage == "primary"]
        for m in used:
            d, lo, hi = cohens_d(e[m].values, pr[m].values, a.n_boot, rng)
            shifts.append(dict(lineage=lin, module=m, median_embryo=e[m].median(),
                               median_primary=pr[m].median(), cohens_d=d, d_ci_low=lo, d_ci_high=hi))

        # within-culture: load vs module / gene
        o = obs.join(wol[["wolbachia_umis"]], how="inner")
        o["load"] = np.log1p(o["wolbachia_umis"])
        o["log_host"] = np.log(o["host_umis"].clip(lower=1))
        for (st, cond), dd in o.groupby(["stage", "condition"]):
            if dd["wolbachia_umis"].median() < a.min_wol:
                continue
            for m in used:
                r_, p_ = partial_spearman(dd[m].values, dd["load"].values, dd["log_host"].values)
                within_m.append(dict(lineage=lin, stage=st, condition=cond, module=m,
                                     n=len(dd), rho=r_, p=p_))
            for g in panel:
                v = dd[f"g_{g}"].values
                r_, p_ = partial_spearman(v, dd["load"].values, dd["log_host"].values)
                within_g.append(dict(lineage=lin, stage=st, condition=cond, gene=g,
                                     frac_expressing=float((v > 0).mean()), rho=r_, p=p_))
        del A

    used_df = pd.DataFrame(used_rows)
    used_df.to_csv(os.path.join(a.out_dir, "modules_used.csv"), index=False)

    # pseudobulk module log2FC, primary vs embryo, per lineage
    cnt = pd.read_csv(os.path.join(a.de_dir, "pseudobulk_counts.csv.gz"), index_col=0)
    S = pd.read_csv(os.path.join(a.de_dir, "pseudobulk_samples.csv"), index_col=0)
    symtab = pd.read_csv(os.path.join(a.de_dir, "de_pooled_line_vs_primary.csv"), index_col=0)["symbol"]
    cpm = np.log2(cnt.div(cnt.sum()) * 1e6 + 1)
    sh = pd.DataFrame(shifts)
    pb = []
    for lin in a.lineages:
        s = S[S["lineage"] == lin]
        e = cpm[s.index[s["stage"] == "embryo"]].mean(1)
        pr = cpm[s.index[s["stage"] == "primary"]].mean(1)
        for m, genes in MODULES_STRESS.items():
            ids = symtab.index[symtab.isin(genes)].intersection(cpm.index)
            ids = ids[(e[ids] >= 1) | (pr[ids] >= 1)]
            pb.append(dict(lineage=lin, module=m, pseudobulk_n_genes=len(ids),
                           pseudobulk_mean_log2FC=float((pr[ids] - e[ids]).mean()) if len(ids) else np.nan))
    sh = sh.merge(pd.DataFrame(pb), on=["lineage", "module"], how="left")

    # culture-level load
    wm = pd.DataFrame(within_m)
    load = wol.groupby(["lineage", "sample_type"])["wolbachia_gene_frac"].median().unstack()
    load = pd.DataFrame({"primary_load_pct": 100 * load["primary_cells"],
                         "load_retained": load["primary_cells"] / load["embryo"]})
    sh = sh.join(load, on="lineage")
    sh.to_csv(os.path.join(a.out_dir, "shift_embryo_to_primary.csv"), index=False)
    print("\nEmbryo -> primary shift per module:\n" + sh.round(3).to_string(index=False))

    cl = []
    for m, d in sh.groupby("module"):
        d = d.dropna(subset=["cohens_d"])
        if len(d) < 3:
            continue
        for xcol in ["primary_load_pct", "load_retained"]:
            for ycol in ["cohens_d", "pseudobulk_mean_log2FC"]:
                dd = d.dropna(subset=[ycol])
                if len(dd) < 3:
                    continue
                rho, pe = exact_spearman(dd[xcol], dd[ycol])
                cl.append(dict(module=m, load_measure=xcol, shift_measure=ycol, n_lineages=len(dd),
                               rho=rho, p_exact=pe,
                               min_attainable_p=2 / math.factorial(len(dd)) if len(dd) > 1 else np.nan))
    cl = pd.DataFrame(cl)
    cl.to_csv(os.path.join(a.out_dir, "culture_level.csv"), index=False)
    print("\nCulture level (descriptive; species, strain and load confounded):\n"
          + cl.round(3).to_string(index=False))

    if len(wm):
        wm["padj"] = multipletests(wm["p"].fillna(1), method="fdr_bh")[1]
    wg = pd.DataFrame(within_g)
    if len(wg):
        wg["padj"] = multipletests(wg["p"].fillna(1), method="fdr_bh")[1]
    wm.to_csv(os.path.join(a.out_dir, "within_culture_modules.csv"), index=False)
    wg.to_csv(os.path.join(a.out_dir, "within_culture_genes.csv"), index=False)
    print("\nWithin culture, load vs module | log host UMIs:\n"
          + wm.pivot_table(index=["lineage", "condition"], columns="module", values="rho").round(3).to_string())

    # infected vs cured line
    de = pd.read_csv(a.infection_de, index_col=0)
    rows = []
    for m, genes in MODULES_STRESS.items():
        inm = de["symbol"].isin(genes)
        d = de[inm]
        if len(d) < 3:
            continue
        u = mannwhitneyu(d["stat"].dropna(), de.loc[~inm, "stat"].dropna(), alternative="two-sided")
        rows.append(dict(module=m, n_genes=len(d), mean_log2FC=d["log2FoldChange"].mean(),
                         n_up=int(((d.padj < 0.05) & (d.log2FoldChange > 0)).sum()),
                         n_down=int(((d.padj < 0.05) & (d.log2FoldChange < 0)).sum()),
                         p_mwu_vs_other_genes=u.pvalue,
                         genes_up=",".join(d.loc[(d.padj < 0.05) & (d.log2FoldChange > 0), "symbol"]),
                         genes_down=",".join(d.loc[(d.padj < 0.05) & (d.log2FoldChange < 0), "symbol"])))
    iv = pd.DataFrame(rows)
    iv.to_csv(os.path.join(a.out_dir, "infected_vs_cured.csv"), index=False)
    print("\nDsim6B-wMel vs Dsim6B:\n" + iv.drop(columns=["genes_up", "genes_down"]).round(4).to_string(index=False))

    # figure
    mods = [m for m in MODULES_STRESS if m in set(sh["module"])] + (["oxphos"] if "oxphos" in set(sh["module"]) else [])
    fig, axes = plt.subplots(2, len(mods), figsize=(3.2 * len(mods), 6.4), squeeze=False)
    for j, m in enumerate(mods):
        ax = axes[0, j]
        d = sh[sh["module"] == m]
        ax.errorbar(d["primary_load_pct"], d["cohens_d"],
                    yerr=[d["cohens_d"] - d["d_ci_low"], d["d_ci_high"] - d["cohens_d"]], fmt="o", color="k")
        for _, r in d.iterrows():
            ax.annotate(r["lineage"], (r["primary_load_pct"], r["cohens_d"]), fontsize=6,
                        xytext=(3, 3), textcoords="offset points")
        ax.set_xscale("log"); ax.axhline(0, color="grey", lw=0.5)
        ax.set_title(m, fontsize=9); ax.set_xlabel("primary Wolbachia (% UMIs)", fontsize=8)
        if j == 0:
            ax.set_ylabel("embryo -> primary shift (d)", fontsize=8)
        ax2 = axes[1, j]
        w = wm[wm["module"] == m]
        ax2.barh(w["lineage"] + " | " + w["stage"], w["rho"], color="grey")
        ax2.axvline(0, color="k", lw=0.5); ax2.tick_params(labelsize=6)
        ax2.set_xlabel("within-culture rho (load)", fontsize=8)
    fig.tight_layout()
    savefig(fig, os.path.join(a.out_dir, "culture_stress.pdf"))
    print("Done ->", a.out_dir)


if __name__ == "__main__":
    main()
