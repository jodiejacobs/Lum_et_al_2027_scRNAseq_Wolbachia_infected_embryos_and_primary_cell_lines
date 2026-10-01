#!/usr/bin/env python3
"""
wolbachia_target_states.py
==========================
Rule pseudotime_wolbachia_states. Does Wolbachia load go with the cell states
that matter for immortalization?

Design: the three states are defined FIRST, from host genes only, and written
to disk before any Wolbachia data are read. Wolbachia load is joined
afterwards and tested with a fixed set of tests (no re-cutting).

States (per lineage; embryo, primary cell line and cell line cells are scored
together, so every call is on one shared scale):

  1. dividing          S-phase and G2/M gene-set scores, each compared with
                       --n_null random gene sets of the same size drawn from
                       the same expression bins (per-cell z). dividing =
                       S or G2/M z >= --div_z. Absolute call: a cell is called
                       dividing only if its cell-cycle genes stand out from
                       expression-matched random genes, not because it is in
                       the top tail of its culture.
  2. line_logit        out-of-fold logit from a classifier (L2 logistic
                       regression, 5-fold CV) trained on cell line vs primary
                       cell line cells of the same lineage. Features exclude
                       cell-cycle genes and genes correlated with the
                       cell-cycle score (|r| > --cc_corr), ribosomal /
                       mitochondrial / heat-shock / immediate-early genes,
                       antimicrobial-peptide and injury genes, and
                       Wolbachia-responsive genes (Dsim6B-wMel vs Dsim6B DE,
                       padj < 0.05 and |log2FC| >= --wol_de_lfc). Counts are
                       downsampled to --downsample UMIs per cell first, so
                       depth cannot separate the classes. line_like =
                       line_logit > 0 (classified as cell line).
  3. death_resistance  mean z of death inhibitors (Diap1, Diap2, Buffy, Bruce)
                       minus mean z of H99 pro-apoptotic genes (rpr, hid,
                       grim, skl); H99_silent = no UMI from any H99 gene.

Validation (states_validation.csv, classifier_auc.csv), before any test:
  per lineage x stage: % dividing, % of dividing / non-dividing cells
  expressing >= 2 core cell-cycle markers, agreement with cell_states.py's
  'proliferating' flag; classifier AUC; % line_like; death-resistance medians.

Tests (per lineage x condition, wherever median Wolbachia UMIs >= --min_wol):
  load = log1p(Wolbachia UMIs); covariate = log host UMIs (no ratio)
  tests_continuous.csv : partial Spearman of load with dividing z (max of S,
                         G2/M), line_logit (not in cell lines), and
                         death_resistance | log host UMIs
  tests_binary.csv     : OLS log1p(Wolbachia UMIs) ~ group + log host UMIs for
                         dividing and line_like (>= --min_cells per group);
                         fold = exp(coef)
  tests_primary_combined.csv : primary cell lines only, Stouffer across
                         lineages; BH across all tests in each file.
"""
import os
import argparse
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, t as tdist
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from pt_utils import savefig
from cell_states import MODULES
from primary_proliferative import S_GENES, G2M_GENES, EXCLUDE

warnings.filterwarnings("ignore", category=FutureWarning)
PRO = ["rpr", "hid", "grim", "skl"]
ANTI = ["Diap1", "Diap2", "Buffy", "Bruce"]
CORE = ["PCNA", "Mcm2", "Mcm5", "CycA", "CycB", "Cdk1", "polo", "stg", "aurB"]
STAGES = {"embryo": "embryo", "primary_cells": "primary", "cell_culture": "line"}


def dense(x):
    return np.asarray(x.todense()) if sp.issparse(x) else np.asarray(x)


def sym2pos(adata):
    d = {}
    for i, s in enumerate(adata.var["symbol"].astype(str).values):
        d.setdefault(s, i)
    return d


def null_z(X, idx, gene_means, n_null, rng, n_bins=25):
    """Per-cell z of a gene-set mean vs bin-matched random gene sets."""
    bins = pd.Series(gene_means).rank(method="first")
    bins = pd.qcut(bins, n_bins, labels=False).values
    pool = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    cols = [np.asarray(idx)] + [np.array([rng.choice(pool[bins[g]]) for g in idx])
                                for _ in range(n_null)]
    W = sp.lil_matrix((X.shape[1], len(cols)))
    for j, c in enumerate(cols):
        for g in c:
            W[g, j] += 1.0 / len(c)
    S = dense(X @ W.tocsr())
    obs, null = S[:, 0], S[:, 1:]
    return (obs - null.mean(1)) / np.clip(null.std(1), 1e-6, None)


def gene_corr(X, y):
    """Pearson r of every gene (columns of sparse X) with vector y."""
    n = X.shape[0]
    yc = (y - y.mean()) / (y.std() + 1e-12)
    mu = np.asarray(X.mean(0)).ravel()
    sq = np.asarray(X.multiply(X).mean(0)).ravel() if sp.issparse(X) else (X ** 2).mean(0)
    sd = np.sqrt(np.clip(sq - mu ** 2, 1e-12, None))
    cov = np.asarray(X.T @ yc).ravel() / n
    return cov / sd


def partial_spearman(y, x, z):
    r = lambda v: pd.Series(v).rank().values
    Z = np.column_stack([np.ones(len(y)), r(z)])
    res = lambda v: v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
    a, b = res(r(y)), res(r(x))
    rho = float(np.corrcoef(a, b)[0, 1])
    df = len(y) - 3
    tt = rho * np.sqrt(df / max(1 - rho ** 2, 1e-12))
    p = 2 * tdist.sf(abs(tt), df)
    return rho, p, np.sign(rho) * norm.isf(p / 2)


def define_states(path, lin, wol_genes, a, rng):
    """Host-gene-only state calls for every cell of one lineage."""
    A = sc.read_h5ad(path)
    A.obs["stage"] = A.obs["sample_type"].astype(str).map(STAGES)
    A = A[A.obs["stage"].notna()].copy()
    X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
    s2p = sym2pos(A)
    sym = A.var["symbol"].astype(str)
    obs = pd.DataFrame(index=A.obs_names)
    obs["lineage"], obs["stage"] = lin, A.obs["stage"].values
    obs["condition"] = A.obs["condition"].astype(str).values
    C = A.layers["counts"] if "counts" in A.layers else None
    obs["host_umis"] = np.asarray(C.sum(1)).ravel() if C is not None else np.nan

    # 1. dividing (absolute, null-calibrated)
    gm = np.asarray(X.mean(0)).ravel()
    for name, genes in [("S", S_GENES), ("G2M", G2M_GENES)]:
        idx = [s2p[g] for g in genes if g in s2p]
        obs[f"z_{name}"] = null_z(X, idx, gm, a.n_null, rng)
    obs["z_cycle"] = obs[["z_S", "z_G2M"]].max(1)
    obs["dividing"] = obs["z_cycle"] >= a.div_z
    core = [s2p[g] for g in CORE if g in s2p]
    obs["n_core_expressed"] = np.asarray((X[:, core] > 0).sum(1)).ravel()

    # 3. death resistance
    def zmean(genes):
        idx = [s2p[g] for g in genes if g in s2p]
        if not idx:
            return np.full(A.n_obs, np.nan), []
        M = dense(X[:, idx])
        M = (M - M.mean(0)) / np.clip(M.std(0), 1e-6, None)
        return M.mean(1), [g for g in genes if g in s2p]
    pro, pro_used = zmean(PRO)
    anti, anti_used = zmean(ANTI)
    obs["pro_apoptotic"], obs["anti_apoptotic"] = pro, anti
    obs["death_resistance"] = anti - pro
    if C is not None and pro_used:
        obs["H99_silent"] = np.asarray(C[:, [s2p[g] for g in pro_used]].sum(1)).ravel() == 0

    # 2. line-likeness classifier (primary vs cell line, out-of-fold)
    sub = obs["stage"].isin(["primary", "line"]).values
    y = (obs.loc[sub, "stage"] == "line").values.astype(int)
    auc, n_feat, coefs = np.nan, 0, None
    obs["line_logit"] = np.nan
    if C is not None and y.sum() >= a.min_cells and (1 - y).sum() >= a.min_cells:
        B = ad.AnnData(X=sp.csr_matrix(C[sub]).astype(np.float32), var=A.var[["symbol"]].copy())
        sc.pp.downsample_counts(B, counts_per_cell=a.downsample, random_state=0)
        sc.pp.normalize_total(B, target_sum=1e4)
        sc.pp.log1p(B)
        Xb = B.X.tocsr()
        cc_genes = set(S_GENES) | set(G2M_GENES) | set(MODULES.get("proliferation", []))
        drop_mod = set(MODULES.get("immune_AMP", [])) | set(MODULES.get("injury_JAK_JNK", []))
        r_cc = gene_corr(Xb, obs.loc[sub, "z_cycle"].values)
        keep = (~sym.str.contains(EXCLUDE, regex=True, na=False).values
                & ~sym.isin(cc_genes | drop_mod | set(PRO) | set(ANTI)).values
                & ~A.var_names.isin(wol_genes)
                & (np.abs(np.nan_to_num(r_cc)) <= a.cc_corr)
                & (np.asarray((Xb > 0).mean(0)).ravel() >= 0.02))
        Xk = Xb[:, keep]
        mu = np.asarray(Xk.mean(0)).ravel()
        var = np.asarray(Xk.multiply(Xk).mean(0)).ravel() - mu ** 2
        top = np.argsort(var)[::-1][:a.n_genes]
        D = dense(Xk[:, top])
        D = (D - D.mean(0)) / np.clip(D.std(0), 1e-6, None)
        clf = LogisticRegression(C=a.C, class_weight="balanced", max_iter=3000)
        cv = StratifiedKFold(5, shuffle=True, random_state=0)
        oof = cross_val_predict(clf, D, y, cv=cv, method="decision_function")
        auc = roc_auc_score(y, oof)
        obs.loc[sub, "line_logit"] = oof
        n_feat = len(top)
        fit = clf.fit(D, y)
        g = A.var_names[keep][top]
        coefs = pd.DataFrame({"symbol": sym.values[keep][top], "coef": fit.coef_[0]}, index=g) \
            .sort_values("coef")
    obs["line_like"] = obs["line_logit"] > 0
    return obs, dict(lineage=lin, classifier_auc=auc, n_features=n_feat,
                     n_wolbachia_responsive_excluded=int(A.var_names.isin(wol_genes).sum()),
                     pro_genes=",".join(pro_used), anti_genes=",".join(anti_used)), coefs


def validation(obs, cs_dir):
    flag = None
    if cs_dir:
        f = os.path.join(cs_dir, f"states_{obs['lineage'].iloc[0]}.csv.gz")
        if os.path.exists(f):
            flag = pd.read_csv(f, index_col=0, usecols=lambda c: c in ("Unnamed: 0", "proliferating"))
            flag = (flag["proliferating"].astype(str) == "True").reindex(obs.index)
    rows = []
    for (lin, st, cond), d in obs.groupby(["lineage", "stage", "condition"]):
        dv = d["dividing"]
        r = dict(lineage=lin, stage=st, condition=cond, n_cells=len(d),
                 pct_dividing=100 * dv.mean(), n_dividing=int(dv.sum()),
                 pct_core2_in_dividing=100 * (d.loc[dv, "n_core_expressed"] >= 2).mean() if dv.any() else np.nan,
                 pct_core2_in_other=100 * (d.loc[~dv, "n_core_expressed"] >= 2).mean(),
                 pct_line_like=100 * d["line_like"].mean() if d["line_logit"].notna().any() else np.nan,
                 median_line_logit=d["line_logit"].median(),
                 median_death_resistance=d["death_resistance"].median(),
                 pct_H99_silent=100 * d["H99_silent"].mean() if "H99_silent" in d else np.nan)
        if flag is not None:
            f = flag.reindex(d.index)
            r["pct_flag_proliferating"] = 100 * f.mean()
            r["pct_flagged_called_dividing"] = 100 * d.loc[f.fillna(False).astype(bool), "dividing"].mean()
        rows.append(r)
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad_dir", required=True, help="results/pseudotime")
    p.add_argument("--lineages", nargs="+", required=True)
    p.add_argument("--infection_dir", required=True)
    p.add_argument("--infection_de", required=True, help="de_Dsim6B-wMel_vs_Dsim6B.csv")
    p.add_argument("--cell_states_dir", default=None)
    p.add_argument("--n_null", type=int, default=100)
    p.add_argument("--div_z", type=float, default=3.0)
    p.add_argument("--downsample", type=int, default=2000)
    p.add_argument("--n_genes", type=int, default=2000)
    p.add_argument("--cc_corr", type=float, default=0.2)
    p.add_argument("--wol_de_lfc", type=float, default=1.0)
    p.add_argument("--C", type=float, default=0.1)
    p.add_argument("--min_cells", type=int, default=50)
    p.add_argument("--min_wol", type=float, default=5)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out_dir", required=True)
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    de = pd.read_csv(a.infection_de, index_col=0)
    wol_genes = set(de.index[(de["padj"] < 0.05) & (de["log2FoldChange"].abs() >= a.wol_de_lfc)])
    print(f"Wolbachia-responsive genes excluded from the classifier: {len(wol_genes)}")

    # ---- Stage A: states from host genes only (no Wolbachia data read yet)
    all_obs, cls, val = [], [], []
    for lin in a.lineages:
        path = os.path.join(a.h5ad_dir, lin, f"prepared_{lin}.h5ad")
        print(f"\n=== {lin}: defining states")
        obs, info, coefs = define_states(path, lin, wol_genes, a, rng)
        obs.to_csv(os.path.join(a.out_dir, f"states_{lin}.csv.gz"))
        if coefs is not None:
            coefs.to_csv(os.path.join(a.out_dir, f"classifier_coefs_{lin}.csv"))
        cls.append(info)
        val.append(validation(obs, a.cell_states_dir))
        all_obs.append(obs)
        print(f"  classifier AUC {info['classifier_auc']:.3f} on {info['n_features']} genes")
    val = pd.concat(val)
    val.to_csv(os.path.join(a.out_dir, "states_validation.csv"), index=False)
    pd.DataFrame(cls).to_csv(os.path.join(a.out_dir, "classifier_auc.csv"), index=False)
    pd.set_option("display.width", 250)
    print("\nState validation (written before any Wolbachia data were read):\n"
          + val.round(2).to_string(index=False))

    # ---- Stage B: join Wolbachia load and run the fixed tests
    obs = pd.concat(all_obs)
    wol = pd.read_csv(os.path.join(a.infection_dir, "cells_infection_immune.csv.gz"), index_col=0,
                      usecols=lambda c: c in ("Unnamed: 0", "wolbachia_umis"))
    obs = obs.join(wol, how="inner")
    obs["load"] = np.log1p(obs["wolbachia_umis"])
    obs["log_host"] = np.log(obs["host_umis"].clip(lower=1))
    cont, binr = [], []
    for (lin, st, cond), d in obs.groupby(["lineage", "stage", "condition"]):
        if d["wolbachia_umis"].median() < a.min_wol:
            continue
        for state in ["z_cycle", "line_logit", "death_resistance"]:
            if state == "line_logit" and st != "primary":
                continue
            v = d[state].values.astype(float)
            ok = ~np.isnan(v)
            if ok.sum() < a.min_cells:
                continue
            rho, pv, z = partial_spearman(v[ok], d["load"].values[ok], d["log_host"].values[ok])
            cont.append(dict(lineage=lin, stage=st, condition=cond, state=state, n=int(ok.sum()),
                             rho=rho, p=pv, z=z))
        for grp in ["dividing", "line_like"]:
            if grp == "line_like" and st != "primary":
                continue
            g = d[grp].astype(bool)
            if g.sum() < a.min_cells or (~g).sum() < a.min_cells:
                binr.append(dict(lineage=lin, stage=st, condition=cond, group=grp,
                                 n_in=int(g.sum()), n_out=int((~g).sum()), note="too few cells"))
                continue
            Xd = sm.add_constant(np.column_stack([g.astype(float), d["log_host"]]))
            f = sm.OLS(d["load"].values, Xd).fit()
            binr.append(dict(lineage=lin, stage=st, condition=cond, group=grp, n_in=int(g.sum()),
                             n_out=int((~g).sum()),
                             median_wol_umis_in=d.loc[g, "wolbachia_umis"].median(),
                             median_wol_umis_out=d.loc[~g, "wolbachia_umis"].median(),
                             fold_in_vs_out=float(np.exp(f.params[1])),
                             ci_low=float(np.exp(f.conf_int()[1][0])), ci_high=float(np.exp(f.conf_int()[1][1])),
                             p=float(f.pvalues[1])))
    cont, binr = pd.DataFrame(cont), pd.DataFrame(binr)
    if len(cont):
        cont["padj"] = multipletests(cont["p"], method="fdr_bh")[1]
    if "p" in binr and binr["p"].notna().any():
        m = binr["p"].notna()
        binr.loc[m, "padj"] = multipletests(binr.loc[m, "p"], method="fdr_bh")[1]
    cont.to_csv(os.path.join(a.out_dir, "tests_continuous.csv"), index=False)
    binr.to_csv(os.path.join(a.out_dir, "tests_binary.csv"), index=False)
    pc = cont[cont["stage"] == "primary"]
    comb = pc.groupby("state").agg(n_lineages=("lineage", "size"),
                                   n_positive=("rho", lambda v: int((v > 0).sum())),
                                   rho_min=("rho", "min"), rho_max=("rho", "max"),
                                   z_stouffer=("z", lambda v: v.sum() / np.sqrt(len(v)))).reset_index()
    comb["p_stouffer"] = 2 * norm.sf(comb["z_stouffer"].abs())
    comb.to_csv(os.path.join(a.out_dir, "tests_primary_combined.csv"), index=False)
    print("\nContinuous (load vs state | log host UMIs):\n" + cont.round(4).to_string(index=False))
    print("\nBinary (Wolbachia UMIs in group vs rest | log host UMIs):\n" + binr.round(4).to_string(index=False))
    print("\nPrimary cell lines, combined across lineages:\n" + comb.round(4).to_string(index=False))

    # ---- figure: primary cell lines, load vs each state
    prim = obs[obs["stage"] == "primary"]
    lins = sorted(prim["lineage"].unique())
    fig, axes = plt.subplots(3, len(lins), figsize=(3.4 * len(lins), 9), squeeze=False)
    for j, lin in enumerate(lins):
        d = prim[prim["lineage"] == lin]
        resid = d["load"] - np.polyval(np.polyfit(d["log_host"], d["load"], 1), d["log_host"])
        for i, state in enumerate(["z_cycle", "line_logit", "death_resistance"]):
            ax = axes[i, j]
            ax.hexbin(d[state], resid, gridsize=35, mincnt=1, cmap="Greys", bins="log")
            ax.set_title(lin if i == 0 else "", fontsize=9)
            ax.set_xlabel(state, fontsize=8)
            if j == 0:
                ax.set_ylabel("Wolbachia load\n(residual on log host UMIs)", fontsize=8)
    fig.tight_layout()
    savefig(fig, os.path.join(a.out_dir, "load_vs_states_primary.pdf"))
    print("Done ->", a.out_dir)


if __name__ == "__main__":
    main()
