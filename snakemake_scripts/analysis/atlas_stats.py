#!/usr/bin/env python3
"""
atlas_stats.py
==============
Rule atlas_stats. Statistical tests for the atlas-projection claims, computed
from the tables written by embryo_to_cellline_trajectory.py (condition level,
so conditions are the replicates; no per-cell pseudoreplication).

  1. diversity_tests.csv       : Shannon entropy, embryos vs cell lines,
                                 exact two-sided Mann-Whitney per label column
  2. correlation_permutation.csv: mean within-group minus between-group
                                 pseudobulk Spearman rho; p from all
                                 C(n, n_embryo) relabelings of embryo/cell line
     parent_rank.csv           : each cell line's rho with every embryo and the
                                 rank of its own parental embryo (--parents)
  3. tissue_specificity.csv    : per condition, top tissue rho, gap to the
                                 second tissue; exact Mann-Whitney embryo vs line
  4. titer_by_condition.csv    : (with --integrated) Wolbachia titer among
                                 titer > 0 cells per condition: n, median, IQR,
                                 fraction >= 0.9; exact Mann-Whitney on
                                 condition medians, infected cell lines vs embryos
"""
import os
import argparse
import itertools

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    r = mannwhitneyu(a, b, alternative="two-sided", method="exact")
    return dict(n_embryo=len(a), n_line=len(b), mean_embryo=a.mean(), mean_line=b.mean(),
                U=r.statistic, p_exact=r.pvalue, min_possible_p=2 / len(list(
                    itertools.combinations(range(len(a) + len(b)), len(a)))))


def diversity(traj, out):
    d = pd.read_csv(os.path.join(traj, "diversity_shannon_entropy.csv"))
    d["is_embryo"] = d["is_embryo"].astype(str).str.lower() == "true"
    rows = []
    for col, g in d.groupby("label_col"):
        rows.append(dict(label_col=col, **mwu(g.loc[g.is_embryo, "shannon_entropy"],
                                              g.loc[~g.is_embryo, "shannon_entropy"])))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(out, "diversity_tests.csv"), index=False)
    print("\nShannon entropy, embryo vs cell line:\n" + res.round(4).to_string(index=False))
    return d


def correlation(traj, out, embryos, parents):
    c = pd.read_csv(os.path.join(traj, "pseudobulk_condition_correlation.csv"), index_col=0)
    conds = list(c.index)
    iu = np.triu_indices(len(conds), 1)
    vals = c.values[iu]

    def stat(emb):
        e = np.array([x in emb for x in conds])
        same = e[iu[0]] == e[iu[1]]
        within_e = vals[same & e[iu[0]]].mean()
        within_l = vals[same & ~e[iu[0]]].mean()
        return vals[same].mean() - vals[~same].mean(), within_e, within_l, vals[~same].mean()

    obs, we, wl, bt = stat(set(embryos))
    k = sum(x in embryos for x in conds)
    null = np.array([stat(set(s))[0] for s in itertools.combinations(conds, k)])
    res = pd.DataFrame([dict(mean_within_embryo=we, mean_within_line=wl, mean_between=bt,
                             within_minus_between=obs, n_relabelings=len(null),
                             p_perm=(null >= obs).mean())])
    res.to_csv(os.path.join(out, "correlation_permutation.csv"), index=False)
    print("\nPseudobulk correlation, within vs between:\n" + res.round(4).to_string(index=False))
    rows = []
    for line, parent in parents.items():
        if line not in c.index or parent not in c.columns:
            continue
        r = c.loc[line, embryos].sort_values(ascending=False)
        rows.append(dict(cell_line=line, parent=parent, rho_parent=r[parent],
                         parent_rank=int(list(r.index).index(parent)) + 1, n_embryos=len(r),
                         best_embryo=r.index[0], rho_best=r.iloc[0]))
    pr = pd.DataFrame(rows)
    pr.to_csv(os.path.join(out, "parent_rank.csv"), index=False)
    print("\nRank of parental embryo:\n" + pr.round(3).to_string(index=False))


def tissue(traj, out, embryos, label="cell_type_tissue"):
    t = pd.read_csv(os.path.join(traj, f"pseudobulk_vs_tissue_{label}.csv"), index_col=0)
    rows = []
    for cond, r in t.iterrows():
        s = r.dropna().sort_values(ascending=False)
        rows.append(dict(condition=cond, is_embryo=cond in embryos, top_tissue=s.index[0],
                         rho_top=s.iloc[0], rho_second=s.iloc[1], gap=s.iloc[0] - s.iloc[1],
                         rho_range=s.iloc[0] - s.iloc[-1]))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(out, "tissue_specificity.csv"), index=False)
    tests = pd.DataFrame([dict(metric=m, **mwu(d.loc[d.is_embryo, m], d.loc[~d.is_embryo, m]))
                          for m in ["rho_top", "gap", "rho_range"]])
    tests.to_csv(os.path.join(out, "tissue_specificity_tests.csv"), index=False)
    print("\nTissue specificity:\n" + d.round(3).to_string(index=False))
    print(tests.round(4).to_string(index=False))


def titer(integrated, out, embryos, infected_lines):
    import anndata as ad
    a = ad.read_h5ad(integrated, backed="r")
    obs = a.obs[["condition", "wolbachia_titer"]].copy()
    a.file.close()
    obs["condition"] = obs["condition"].astype(str)
    pos = obs[obs["wolbachia_titer"] > 0]
    g = pos.groupby("condition")["wolbachia_titer"]
    d = pd.DataFrame({"n_cells": obs.groupby("condition").size(), "n_titer_pos": g.size(),
                      "median": g.median(), "q25": g.quantile(0.25), "q75": g.quantile(0.75),
                      "frac_ge_0.9": g.apply(lambda x: (x >= 0.9).mean())})
    d["group"] = np.where(d.index.isin(embryos), "embryo",
                          np.where(d.index.isin(infected_lines), "infected_line", "uninfected_line"))
    d.to_csv(os.path.join(out, "titer_by_condition.csv"))
    tests = pd.DataFrame([dict(metric=m, **mwu(d.loc[d.group == "embryo", m],
                                               d.loc[d.group == "infected_line", m]))
                          for m in ["median", "frac_ge_0.9"]])
    tests.to_csv(os.path.join(out, "titer_tests.csv"), index=False)
    print("\nTiter among titer > 0 cells:\n" + d.round(3).to_string())
    print(tests.round(4).to_string(index=False))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--traj_dir", required=True, help="embryo_to_cellline_trajectory.py output")
    p.add_argument("--integrated", default=None, help="integrated.h5ad (for titer tests)")
    p.add_argument("--parents", nargs="+", default=[], help="cell_line:parental_embryo")
    p.add_argument("--infected_lines", nargs="*", default=[])
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    d = diversity(args.traj_dir, args.out_dir)
    embryos = sorted(d.loc[d.is_embryo, "condition"].unique())
    parents = dict(x.split(":", 1) for x in args.parents)
    correlation(args.traj_dir, args.out_dir, embryos, parents)
    tissue(args.traj_dir, args.out_dir, embryos)
    if args.integrated:
        titer(args.integrated, args.out_dir, embryos, args.infected_lines)
    print("Done ->", args.out_dir)


if __name__ == "__main__":
    main()
