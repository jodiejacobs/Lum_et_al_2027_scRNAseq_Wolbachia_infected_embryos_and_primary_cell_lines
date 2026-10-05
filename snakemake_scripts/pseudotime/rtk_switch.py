#!/usr/bin/env python3
"""
rtk_switch.py
=============
Rule pseudotime_rtk_switch. Cells lose Egfr and gain Pvr from embryo ->
primary culture -> cell line. What is the switch, and does it track Wolbachia?

Per cell (each lineage's prepared h5ad, all stages):
  receptor_class  Egfr_only / Pvr_only / both / neither (UMI detected > 0)
  switch_index    z(log expr Pvr) - z(log expr Egfr), z across the lineage
  modules         PVR_ligand (Pvf1-3), EGFR_ligand (spi, Krn, vn, grk),
                  EGFR_processing (rho, ru, S), MAPK targets, hemocyte and
                  mesoderm identity (gene lists below; genes found -> modules_used.csv)

Questions (Q) and outputs (out_dir):
  Q1 switch        switch_by_sample.csv     class fractions, mean expression and
                                            pseudobulk CPM of the panel per sample
                   switch_stage_tests.csv   per lineage, embryo->primary and
                                            primary->line: Cohen's d of switch_index
                                            (bootstrap CI), Fisher on Pvr_only
                                            fraction, pseudobulk log2 Pvr/Egfr shift
  Q2 ligands       autocrine.csv            per sample: is a ligand detected more
                                            often in cells with its receptor?
                                            Logistic regression ligand ~ receptor +
                                            log host UMIs (depth-adjusted odds ratio;
                                            NaN under separation) and raw Haldane OR
                   lr_by_state.csv          sender state x receiver state mean-product
                                            score per ligand-receptor pair,
                                            permutation p (state labels shuffled,
                                            --n_perm, default 2000)
  Q3 identity      identity.csv             identity module z in Pvr_only vs Egfr_only
                                            and Pvr+ vs Pvr- cells per sample (MWU)
                   identity_states.csv      receptor class x cell_state (cell_states.py)
                                            and x atlas_germ_layer
  Q4 proliferative proliferative.csv        primary cultures: Pvr+/Egfr+ detection
                                            in the proliferative population vs rest
                                            (depth-adjusted OR) and switch_index MWU;
                                            all samples: partial Spearman switch_index
                                            vs proliferation score | log host UMIs
  Q5 MAPK          mapk.csv                 MAPK-target score by receptor class;
                                            partial Spearman of Egfr, Pvr and aos
                                            (EGFR-biased target) with MAPK score
                                            (aos vs MAPK_no_aos, to avoid circularity)
  Q6 species       species_switch.csv       per lineage: Egfr retained into primary
                                            (CPM ratio) and log2 Pvr/Egfr per stage
                   species_de.csv           panel genes from the stage DE and
                                            Dsim - Dmel interaction tables
  Q7 Wolbachia     wolbachia_within.csv     per sample (median Wolbachia UMIs >=
                                            --min_wol): partial Spearman of
                                            log1p(Wolbachia UMIs) with switch_index,
                                            Egfr, Pvr, ligand modules, MAPK | log host
                                            UMIs; BH across all tests
                   wolbachia_by_class.csv   Wolbachia UMIs per receptor class; MWU of
                                            host-UMI-adjusted load, Pvr_only vs Egfr_only
                   wolbachia_quintiles.csv  within-sample load quintile -> Egfr+, Pvr+,
                                            switch_index
                   wolbachia_culture_level.csv  per lineage: primary load and load
                                            retained vs Egfr loss and switch shift;
                                            exact Spearman, n = 4 (descriptive)
                   wolbachia_infected_vs_cured.csv  panel genes, Dsim6B-wMel vs Dsim6B
  rtk_switch.pdf

Caveats: primary cultures are not ancestors of the cell lines, so a stage
shift cannot separate selection of Pvr+ cells from receptor switching within
cells. Detection depends on depth, so detection tests adjust for log host UMIs.
Lineage, host species and Wolbachia strain are confounded (n = 4).
"""
import os
import argparse
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from pt_utils import savefig
from cell_states import MODULES, score_modules, gene_expr
from culture_stress import partial_spearman, cohens_d, exact_spearman

warnings.filterwarnings("ignore")
STAGES = {"embryo": "embryo", "primary_cells": "primary", "cell_culture": "line"}
STAGE_ORDER = ["embryo", "primary", "line"]
RECEPTORS = ["Egfr", "Pvr"]
LIGANDS = {"Pvr": ["Pvf1", "Pvf2", "Pvf3"], "Egfr": ["spi", "Krn", "vn", "grk"]}
PANEL = RECEPTORS + LIGANDS["Pvr"] + LIGANDS["Egfr"] + ["rho", "ru", "S", "aos"]
MODULES_RTK = {
    "PVR_ligand": LIGANDS["Pvr"],
    "EGFR_ligand": LIGANDS["Egfr"],
    "EGFR_processing": ["rho", "ru", "S"],
    "MAPK_targets": MODULES["MAPK_targets"],
    "MAPK_no_aos": [g for g in MODULES["MAPK_targets"] if g != "aos"],
    "proliferation": MODULES["proliferation"],
    "plasmatocyte": MODULES["plasmatocyte"],
    "crystal_cell": MODULES["crystal_cell"],
    "lamellocyte": MODULES["lamellocyte"],
    "hemocyte_TF": ["srp", "gcm", "gcm2", "ush"],          # Pvr excluded (circular)
    "mesoderm": ["twi", "Mef2", "htl", "zfh1", "tin", "bap"],
    "muscle": MODULES["muscle"],
    "neural": MODULES["neural"],
    "epidermis": MODULES["epidermis"],
    "fat_body": MODULES["fat_body"],
}
IDENTITY = ["plasmatocyte", "crystal_cell", "lamellocyte", "hemocyte_TF", "mesoderm",
            "muscle", "neural", "epidermis", "fat_body"]
LR_PAIRS = [(l, r) for r in RECEPTORS for l in LIGANDS[r]]
WOL_READOUTS = ["switch_index", "expr_Egfr", "expr_Pvr", "PVR_ligand", "EGFR_ligand",
                "EGFR_processing", "MAPK_targets"]


def zs(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std() if x.std() > 0 else x * 0


def logit_or(y, x, cov, min_n=10):
    """Odds ratio of binary y for binary x, adjusted for cov."""
    y, x = np.asarray(y, int), np.asarray(x, int)
    if min(x.sum(), (1 - x).sum(), y.sum(), (1 - y).sum()) < min_n:
        return np.nan, np.nan
    X = sm.add_constant(np.column_stack([x, cov]))
    try:
        f = sm.Logit(y, X).fit(disp=0, maxiter=200)
    except Exception:
        return np.nan, np.nan
    # quasi-separation gives unstable, huge coefficients: report NaN instead
    if not f.mle_retvals.get("converged", False) or abs(f.params[1]) > 10:
        return np.nan, np.nan
    return float(np.exp(f.params[1])), float(f.pvalues[1])


def haldane_or(y, x):
    """Unadjusted odds ratio with 0.5 added to every cell."""
    y, x = np.asarray(y, bool), np.asarray(x, bool)
    a_, b_ = (y & x).sum() + .5, (~y & x).sum() + .5
    c_, d_ = (y & ~x).sum() + .5, (~y & ~x).sum() + .5
    return float(a_ * d_ / (b_ * c_))


def mwu_p(a, b, min_n=10):
    a, b = np.asarray(a)[~np.isnan(a)], np.asarray(b)[~np.isnan(b)]
    if len(a) < min_n or len(b) < min_n:
        return np.nan
    return mannwhitneyu(a, b, alternative="two-sided").pvalue


def lr_scores(X, labels, states, n_perm, rng):
    """Mean-product ligand-receptor scores, sender x receiver, with permutation p.
    X: cells x 2 (ligand, receptor)."""
    M = (np.asarray(labels)[:, None] == np.asarray(states)[None, :]).astype(float)
    M /= M.sum(0)                                   # cells x states, column means
    m = M.T @ X                                     # states x 2
    obs = np.outer(m[:, 0], m[:, 1])
    ge = np.zeros_like(obs)
    for _ in range(n_perm):                         # shuffling cells = shuffling labels
        mp = M.T @ X[rng.permutation(len(X))]
        ge += np.outer(mp[:, 0], mp[:, 1]) >= obs
    return obs, (ge + 1) / (n_perm + 1)


def load_lineage(lin, a):
    A = sc.read_h5ad(os.path.join(a.h5ad_dir, lin, f"prepared_{lin}.h5ad"))
    A.obs["stage"] = A.obs["sample_type"].astype(str).map(STAGES)
    A = A[A.obs["stage"].notna()].copy()
    used = score_modules(A, MODULES_RTK)
    keep = [c for c in ["stage", "condition", "source_file", "species", "atlas_germ_layer"]
            if c in A.obs]
    obs = A.obs[keep].astype(str).copy()
    obs["lineage"] = lin
    for m in used:
        obs[m] = A.obs[f"score_{m}"].values
    obs = obs.join(gene_expr(A, PANEL))
    C = A.layers["counts"] if "counts" in A.layers else None
    obs["host_umis"] = (np.asarray(C.sum(1)).ravel() if C is not None
                        else A.obs["n_counts"].astype(float).values)
    obs["log_host"] = np.log(obs["host_umis"].clip(lower=1))
    for g in RECEPTORS:
        if f"expr_{g}" not in obs:
            obs[f"expr_{g}"] = 0.0
    e, p = obs["expr_Egfr"] > 0, obs["expr_Pvr"] > 0
    obs["receptor_class"] = np.select([e & ~p, p & ~e, e & p], ["Egfr_only", "Pvr_only", "both"],
                                      "neither")
    obs["switch_index"] = zs(obs["expr_Pvr"]) - zs(obs["expr_Egfr"])
    for r in RECEPTORS:
        cols = [f"expr_{l}" for l in LIGANDS[r] if f"expr_{l}" in obs]
        obs[f"det_{r}_ligand"] = (obs[cols] > 0).any(axis=1) if cols else False
    procs = [f"expr_{g}" for g in ["rho", "ru"] if f"expr_{g}" in obs]
    obs["det_EGFR_processing"] = (obs[procs] > 0).any(axis=1) if procs else False

    # cell states, proliferative calls, Wolbachia
    st = os.path.join(a.states_dir, f"states_{lin}.csv.gz")
    if os.path.exists(st):
        obs = obs.join(pd.read_csv(st, index_col=0, usecols=lambda c: c in ("Unnamed: 0", "cell_state")))
    pf = os.path.join(a.prolif_dir, lin, "cells.csv.gz")
    if os.path.exists(pf):
        obs = obs.join(pd.read_csv(pf, index_col=0,
                                   usecols=lambda c: c in ("Unnamed: 0", "proliferative", "prolif_z")))
    lab = {}
    LR = []
    for smp, d in obs.groupby("source_file"):
        if "cell_state" not in d:
            break
        lab_ = d["cell_state"].astype(str).values
        states = [s for s, n in pd.Series(lab_).value_counts().items() if n >= a.min_cells]
        mask = np.isin(lab_, states)
        if len(states) < 2:
            continue
        for l, r in LR_PAIRS:
            if f"expr_{l}" not in d or f"expr_{r}" not in d:
                continue
            X = d.loc[mask, [f"expr_{l}", f"expr_{r}"]].values
            sc_, p_ = lr_scores(X, lab_[mask], states, a.n_perm, a.rng)
            for i, s in enumerate(states):
                for j, t in enumerate(states):
                    LR.append(dict(lineage=lin, sample=smp, stage=d["stage"].iloc[0],
                                   condition=d["condition"].iloc[0], ligand=l, receptor=r,
                                   sender=s, receiver=t, score=sc_[i, j], p_perm=p_[i, j]))
    del A
    return obs, used, LR


def pseudobulk_panel(de_dir):
    cnt = pd.read_csv(os.path.join(de_dir, "pseudobulk_counts.csv.gz"), index_col=0)
    S = pd.read_csv(os.path.join(de_dir, "pseudobulk_samples.csv"), index_col=0)
    sym = pd.read_csv(os.path.join(de_dir, "de_pooled_line_vs_primary.csv"), index_col=0)["symbol"]
    cpm = cnt.div(cnt.sum()) * 1e6
    ids = sym[sym.isin(PANEL)]
    pb = cpm.loc[ids.index.intersection(cpm.index)].rename(index=ids).T
    pb.columns = [f"cpm_{c}" for c in pb.columns]
    pb["log2_Pvr_over_Egfr"] = np.log2((pb["cpm_Pvr"] + 1) / (pb["cpm_Egfr"] + 1))
    return S.join(pb)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad_dir", required=True)
    p.add_argument("--lineages", nargs="+", required=True)
    p.add_argument("--states_dir", required=True)
    p.add_argument("--prolif_dir", required=True)
    p.add_argument("--infection_dir", required=True)
    p.add_argument("--de_dir", required=True)
    p.add_argument("--infection_de", required=True)
    p.add_argument("--min_wol", type=float, default=5)
    p.add_argument("--min_cells", type=int, default=20)
    p.add_argument("--n_perm", type=int, default=2000)
    p.add_argument("--n_boot", type=int, default=500)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out_dir", required=True)
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    a.rng = np.random.default_rng(a.seed)
    pd.set_option("display.width", 250)
    out = lambda f: os.path.join(a.out_dir, f)

    wol = pd.read_csv(os.path.join(a.infection_dir, "cells_infection_immune.csv.gz"), index_col=0,
                      usecols=lambda c: c in ("Unnamed: 0", "wolbachia_umis", "wolbachia_gene_frac"))
    cells, used_rows, LR = [], [], []
    for lin in a.lineages:
        print(f"\n=== {lin}")
        o, used, lr = load_lineage(lin, a)
        cells.append(o)
        LR += lr
        used_rows += [dict(lineage=lin, module=m, n_genes=len(g), genes=",".join(g))
                      for m, g in used.items()]
    obs = pd.concat(cells).join(wol, how="left")
    obs["wol_load"] = np.log1p(obs["wolbachia_umis"])
    pd.DataFrame(used_rows).to_csv(out("modules_used.csv"), index=False)
    obs.to_csv(out("cells_rtk.csv.gz"))
    pb = pseudobulk_panel(a.de_dir)

    # ---- Q1 switch per sample and per stage transition
    rows = []
    for smp, d in obs.groupby("source_file"):
        r = dict(sample=smp, lineage=d["lineage"].iloc[0], species=d.get("species", pd.Series([""])).iloc[0],
                 stage=d["stage"].iloc[0], condition=d["condition"].iloc[0], n_cells=len(d),
                 median_switch_index=d["switch_index"].median(),
                 mean_expr_Egfr=d["expr_Egfr"].mean(), mean_expr_Pvr=d["expr_Pvr"].mean())
        r.update({f"frac_{k}": v for k, v in
                  d["receptor_class"].value_counts(normalize=True).items()})
        rows.append(r)
    bys = pd.DataFrame(rows).fillna({c: 0 for c in
                                     ["frac_Egfr_only", "frac_Pvr_only", "frac_both", "frac_neither"]})
    bys = bys.merge(pb.drop(columns=["condition", "lineage", "species", "stage"], errors="ignore"),
                    left_on="sample", right_index=True, how="left")
    order = {s: i for i, s in enumerate(STAGE_ORDER)}
    bys = bys.sort_values(["lineage", "stage"], key=lambda c: c.map(order) if c.name == "stage" else c)
    bys.to_csv(out("switch_by_sample.csv"), index=False)
    print("\nQ1 switch per sample:\n" + bys[["sample", "stage", "frac_Egfr_only", "frac_Pvr_only",
          "frac_both", "frac_neither", "cpm_Egfr", "cpm_Pvr", "log2_Pvr_over_Egfr"]].round(3).to_string(index=False))

    tr = []
    for lin, d in obs.groupby("lineage"):
        for s0, s1 in [("embryo", "primary"), ("primary", "line")]:
            x, y = d[d.stage == s0], d[d.stage == s1]
            if not len(x) or not len(y):
                continue
            dd, lo, hi = cohens_d(x["switch_index"].values, y["switch_index"].values, a.n_boot, a.rng)
            t = [[(y.receptor_class == "Pvr_only").sum(), (y.receptor_class != "Pvr_only").sum()],
                 [(x.receptor_class == "Pvr_only").sum(), (x.receptor_class != "Pvr_only").sum()]]
            orr, fp = fisher_exact(t)
            b = pb[pb["lineage"] == lin]
            tr.append(dict(lineage=lin, transition=f"{s0}->{s1}", cohens_d_switch=dd, d_ci_low=lo,
                           d_ci_high=hi, frac_Pvr_only_from=(x.receptor_class == "Pvr_only").mean(),
                           frac_Pvr_only_to=(y.receptor_class == "Pvr_only").mean(),
                           frac_Egfr_pos_from=(x.expr_Egfr > 0).mean(),
                           frac_Egfr_pos_to=(y.expr_Egfr > 0).mean(), OR_Pvr_only=orr, fisher_p=fp,
                           pseudobulk_log2ratio_shift=b.loc[b.stage == s1, "log2_Pvr_over_Egfr"].mean()
                           - b.loc[b.stage == s0, "log2_Pvr_over_Egfr"].mean()))
    pd.DataFrame(tr).to_csv(out("switch_stage_tests.csv"), index=False)
    print("\nQ1 stage transitions:\n" + pd.DataFrame(tr).round(3).to_string(index=False))

    # ---- Q2 ligands: autocrine co-detection and sender/receiver
    au = []
    for smp, d in obs.groupby("source_file"):
        base = dict(sample=smp, lineage=d["lineage"].iloc[0], stage=d["stage"].iloc[0],
                    condition=d["condition"].iloc[0])
        for r, ycol in [("Pvr", "det_Pvr_ligand"), ("Egfr", "det_Egfr_ligand"),
                        ("Egfr", "det_EGFR_processing")]:
            rec = d[f"expr_{r}"] > 0
            orr, pv = logit_or(d[ycol], rec, d["log_host"].values)
            au.append(dict(base, receptor=r, readout=ycol, n_receptor_pos=int(rec.sum()),
                           frac_readout_in_rec_pos=d.loc[rec, ycol].mean() if rec.any() else np.nan,
                           frac_readout_in_rec_neg=d.loc[~rec, ycol].mean() if (~rec).any() else np.nan,
                           OR_raw_haldane=haldane_or(d[ycol], rec), OR_depth_adj=orr, p=pv))
        for l in LIGANDS["Pvr"]:
            if f"expr_{l}" in d:
                orr, pv = logit_or(d[f"expr_{l}"] > 0, d["expr_Pvr"] > 0, d["log_host"].values)
                au.append(dict(base, receptor="Pvr", readout=f"det_{l}", OR_depth_adj=orr, p=pv,
                               OR_raw_haldane=haldane_or(d[f"expr_{l}"] > 0, d["expr_Pvr"] > 0),
                               frac_readout_in_rec_pos=(d.loc[d.expr_Pvr > 0, f"expr_{l}"] > 0).mean()))
    au = pd.DataFrame(au)
    au["padj"] = multipletests(au["p"].fillna(1), method="fdr_bh")[1]
    au.to_csv(out("autocrine.csv"), index=False)
    print("\nQ2 autocrine (ligand detected | receptor+, depth-adjusted OR):\n"
          + au.pivot_table(index=["lineage", "stage", "sample"], columns="readout",
                           values="OR_depth_adj").round(2).to_string())
    lr = pd.DataFrame(LR)
    if len(lr):
        lr["padj"] = multipletests(lr["p_perm"], method="fdr_bh")[1]
        lr.to_csv(out("lr_by_state.csv"), index=False)

    # ---- Q3 identity of Pvr+ cells
    idn, ist = [], []
    idmods = [m for m in IDENTITY if m in obs]
    for smp, d in obs.groupby("source_file"):
        z = d[idmods].apply(zs)
        for lbl, m1, m0 in [("Pvr_only_vs_Egfr_only", d.receptor_class == "Pvr_only",
                             d.receptor_class == "Egfr_only"),
                            ("Pvr_pos_vs_neg", d.expr_Pvr > 0, d.expr_Pvr == 0)]:
            for m in idmods:
                idn.append(dict(sample=smp, lineage=d["lineage"].iloc[0], stage=d["stage"].iloc[0],
                                contrast=lbl, module=m, n1=int(m1.sum()), n0=int(m0.sum()),
                                mean_z_1=z.loc[m1, m].mean(), mean_z_0=z.loc[m0, m].mean(),
                                diff=z.loc[m1, m].mean() - z.loc[m0, m].mean(),
                                p=mwu_p(z.loc[m1, m].values, z.loc[m0, m].values)))
        for col in ["cell_state", "atlas_germ_layer"]:
            if col in d:
                ct = pd.crosstab(d["receptor_class"], d[col], normalize="index")
                for cls, row in ct.iterrows():
                    for k, v in row.items():
                        ist.append(dict(sample=smp, lineage=d["lineage"].iloc[0], stage=d["stage"].iloc[0],
                                        annotation=col, receptor_class=cls, label=k, frac=v))
    idn = pd.DataFrame(idn)
    idn["padj"] = multipletests(idn["p"].fillna(1), method="fdr_bh")[1]
    idn.to_csv(out("identity.csv"), index=False)
    pd.DataFrame(ist).to_csv(out("identity_states.csv"), index=False)
    print("\nQ3 identity, Pvr_only - Egfr_only (module z):\n"
          + idn[idn.contrast == "Pvr_only_vs_Egfr_only"].pivot_table(
              index=["lineage", "stage"], columns="module", values="diff").round(2).to_string())

    # ---- Q4 proliferative population
    pr = []
    for smp, d in obs.groupby("source_file"):
        r = dict(sample=smp, lineage=d["lineage"].iloc[0], stage=d["stage"].iloc[0])
        if "proliferation" in d:
            r["rho_switch_vs_prolif"], r["p_switch_vs_prolif"] = partial_spearman(
                d["switch_index"].values, d["proliferation"].values, d["log_host"].values)
        if "proliferative" in d and d["proliferative"].notna().any():
            pm = d["proliferative"].fillna(False).astype(bool)
            r.update(n_proliferative=int(pm.sum()),
                     frac_Pvr_pos_prolif=(d.loc[pm, "expr_Pvr"] > 0).mean(),
                     frac_Pvr_pos_rest=(d.loc[~pm, "expr_Pvr"] > 0).mean(),
                     frac_Egfr_pos_prolif=(d.loc[pm, "expr_Egfr"] > 0).mean(),
                     frac_Egfr_pos_rest=(d.loc[~pm, "expr_Egfr"] > 0).mean(),
                     switch_prolif_minus_rest=d.loc[pm, "switch_index"].median()
                     - d.loc[~pm, "switch_index"].median(),
                     p_switch_mwu=mwu_p(d.loc[pm, "switch_index"].values, d.loc[~pm, "switch_index"].values))
            for g in RECEPTORS:
                r[f"OR_{g}_pos_in_prolif"], r[f"p_OR_{g}"] = logit_or(
                    d[f"expr_{g}"] > 0, pm, d["log_host"].values)
        pr.append(r)
    pr = pd.DataFrame(pr)
    pr.to_csv(out("proliferative.csv"), index=False)
    print("\nQ4 proliferative:\n" + pr.round(3).to_string(index=False))

    # ---- Q5 MAPK readout
    mk = []
    for smp, d in obs.groupby("source_file"):
        if "MAPK_targets" not in d:
            continue
        r = dict(sample=smp, lineage=d["lineage"].iloc[0], stage=d["stage"].iloc[0])
        for cls, dd in d.groupby("receptor_class"):
            r[f"MAPK_mean_{cls}"] = dd["MAPK_targets"].mean()
            r[f"n_{cls}"] = len(dd)
        r["p_MAPK_Pvr_only_vs_neither"] = mwu_p(d.loc[d.receptor_class == "Pvr_only", "MAPK_targets"].values,
                                                d.loc[d.receptor_class == "neither", "MAPK_targets"].values)
        for g in ["expr_Egfr", "expr_Pvr", "expr_aos"]:
            # aos is in MAPK_targets, so test it against the score without aos
            score = "MAPK_no_aos" if g == "expr_aos" and "MAPK_no_aos" in d else "MAPK_targets"
            if g in d:
                r[f"rho_{g}_MAPK"], r[f"p_{g}_MAPK"] = partial_spearman(
                    d[g].values, d[score].values, d["log_host"].values)
        mk.append(r)
    pd.DataFrame(mk).to_csv(out("mapk.csv"), index=False)
    print("\nQ5 MAPK:\n" + pd.DataFrame(mk).round(3).to_string(index=False))

    # ---- Q6 species
    sw = []
    for lin, b in pb.groupby("lineage"):
        g = b.groupby("stage")[["cpm_Egfr", "cpm_Pvr", "log2_Pvr_over_Egfr"]].mean()
        r = dict(lineage=lin, species=b["species"].iloc[0])
        for s in g.index:
            r[f"log2_Pvr_over_Egfr_{s}"] = g.loc[s, "log2_Pvr_over_Egfr"]
            r[f"cpm_Egfr_{s}"] = g.loc[s, "cpm_Egfr"]
        if {"embryo", "primary"} <= set(g.index):
            r["Egfr_retained_primary"] = g.loc["primary", "cpm_Egfr"] / max(g.loc["embryo", "cpm_Egfr"], 1e-9)
            r["Pvr_fold_primary"] = g.loc["primary", "cpm_Pvr"] / max(g.loc["embryo", "cpm_Pvr"], 1e-9)
        sw.append(r)
    sw = pd.DataFrame(sw)
    sw.to_csv(out("species_switch.csv"), index=False)
    print("\nQ6 species:\n" + sw.round(3).to_string(index=False))
    sde = []
    for f in sorted(os.listdir(a.de_dir)):
        if f.startswith("de_") and f.endswith(".csv"):
            t = pd.read_csv(os.path.join(a.de_dir, f), index_col=0)
            if "symbol" in t:
                t = t[t["symbol"].isin(PANEL)].assign(table=f[3:-4])
                sde.append(t)
    if sde:
        sde = pd.concat(sde)
        sde.to_csv(out("species_de.csv"))
        print("\nQ6 Dsim - Dmel interaction, panel genes (log2FC):\n" + sde[sde.table.str.startswith(
            "interaction")].pivot_table(index="symbol", columns="table",
                                        values="log2FoldChange").round(2).to_string())

    # ---- Q7 Wolbachia UMIs vs the switch
    ww, wc, wq = [], [], []
    for smp, d in obs.groupby("source_file"):
        d = d.dropna(subset=["wolbachia_umis"])
        if not len(d) or d["wolbachia_umis"].median() < a.min_wol:
            continue
        base = dict(sample=smp, lineage=d["lineage"].iloc[0], stage=d["stage"].iloc[0],
                    condition=d["condition"].iloc[0], n=len(d),
                    median_wol_umis=d["wolbachia_umis"].median())
        for y in WOL_READOUTS:
            if y in d:
                rho, pv = partial_spearman(d[y].values, d["wol_load"].values, d["log_host"].values)
                ww.append(dict(base, readout=y, rho=rho, p=pv))
        # host-depth-adjusted load (residual of log1p Wolbachia UMIs on log host UMIs)
        Z = np.column_stack([np.ones(len(d)), d["log_host"].values])
        resid = d["wol_load"].values - Z @ np.linalg.lstsq(Z, d["wol_load"].values, rcond=None)[0]
        d = d.assign(wol_resid=resid)
        r = dict(base)
        for cls, dd in d.groupby("receptor_class"):
            r[f"median_wol_umis_{cls}"] = dd["wolbachia_umis"].median()
            r[f"mean_wol_resid_{cls}"] = dd["wol_resid"].mean()
            r[f"n_{cls}"] = len(dd)
        r["p_resid_Pvr_only_vs_Egfr_only"] = mwu_p(d.loc[d.receptor_class == "Pvr_only", "wol_resid"].values,
                                                   d.loc[d.receptor_class == "Egfr_only", "wol_resid"].values)
        r["p_resid_Pvr_pos_vs_neg"] = mwu_p(d.loc[d.expr_Pvr > 0, "wol_resid"].values,
                                            d.loc[d.expr_Pvr == 0, "wol_resid"].values)
        wc.append(r)
        q = pd.qcut(d["wol_resid"].rank(method="first"), 5, labels=[1, 2, 3, 4, 5])
        for k, dd in d.groupby(q):
            wq.append(dict(sample=smp, lineage=base["lineage"], stage=base["stage"], quintile=int(k),
                           median_wol_umis=dd["wolbachia_umis"].median(),
                           frac_Egfr_pos=(dd.expr_Egfr > 0).mean(), frac_Pvr_pos=(dd.expr_Pvr > 0).mean(),
                           mean_switch_index=dd["switch_index"].mean()))
    ww = pd.DataFrame(ww)
    if len(ww):
        ww["padj"] = multipletests(ww["p"].fillna(1), method="fdr_bh")[1]
    ww.to_csv(out("wolbachia_within.csv"), index=False)
    pd.DataFrame(wc).to_csv(out("wolbachia_by_class.csv"), index=False)
    wq = pd.DataFrame(wq)
    wq.to_csv(out("wolbachia_quintiles.csv"), index=False)
    if len(ww):
        print("\nQ7 within sample, log1p(Wolbachia UMIs) vs readout | log host UMIs (rho):\n"
              + ww.pivot_table(index=["lineage", "stage", "sample"], columns="readout",
                               values="rho").round(3).to_string())
    print("\nQ7 Wolbachia by receptor class:\n" + pd.DataFrame(wc).round(3).to_string(index=False))

    w2 = obs.dropna(subset=["wolbachia_gene_frac"]).groupby(["lineage", "stage"])["wolbachia_gene_frac"].median().unstack()
    cl = sw.set_index("lineage").join(pd.DataFrame({
        "primary_load_pct": 100 * w2.get("primary"),
        "load_retained": w2.get("primary") / w2.get("embryo")}))
    cl["switch_shift_primary"] = cl.get("log2_Pvr_over_Egfr_primary") - cl.get("log2_Pvr_over_Egfr_embryo")
    clr = []
    for x in ["primary_load_pct", "load_retained"]:
        for y in ["Egfr_retained_primary", "switch_shift_primary"]:
            dd = cl[[x, y]].dropna()
            if len(dd) >= 3:
                rho, pe = exact_spearman(dd[x], dd[y])
                clr.append(dict(load_measure=x, switch_measure=y, n_lineages=len(dd), rho=rho, p_exact=pe))
    cl.to_csv(out("wolbachia_culture_lineages.csv"))
    pd.DataFrame(clr).to_csv(out("wolbachia_culture_level.csv"), index=False)
    print("\nQ7 culture level (n = 4, descriptive):\n" + cl.round(3).to_string()
          + "\n" + pd.DataFrame(clr).round(3).to_string(index=False))

    de = pd.read_csv(a.infection_de, index_col=0)
    iv = de[de["symbol"].isin(PANEL)][["symbol", "baseMean", "log2FoldChange", "padj"]]
    iv.to_csv(out("wolbachia_infected_vs_cured.csv"))
    print("\nQ7 Dsim6B-wMel vs Dsim6B, panel genes:\n" + iv.round(4).to_string())

    # ---- figure
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    b = bys.set_index("sample")
    cols = {"Egfr_only": "#3B6FB6", "both": "#8C6BB1", "Pvr_only": "#C0392B", "neither": "#D0D0D0"}
    bottom = np.zeros(len(b))
    for k, c in cols.items():
        v = b.get(f"frac_{k}", pd.Series(0, index=b.index)).values
        ax[0, 0].bar(range(len(b)), v, bottom=bottom, color=c, label=k)
        bottom += v
    ax[0, 0].set_xticks(range(len(b)))
    ax[0, 0].set_xticklabels([f"{s} ({st})" for s, st in zip(b.index.str.replace("_pipseq", ""), b["stage"])],
                             rotation=90, fontsize=6)
    ax[0, 0].set_ylabel("fraction of cells")
    ax[0, 0].legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    ax[0, 0].set_title("Receptor detection per sample", fontsize=10)

    for lin, d in pb.groupby("lineage"):
        g = d.groupby("stage")["log2_Pvr_over_Egfr"].mean().reindex(STAGE_ORDER)
        ax[0, 1].plot(range(3), g.values, "o-", label=lin,
                      ls="-" if d["species"].iloc[0] == "Dmel" else "--")
    ax[0, 1].set_xticks(range(3)); ax[0, 1].set_xticklabels(STAGE_ORDER)
    ax[0, 1].axhline(0, color="grey", lw=0.5); ax[0, 1].legend(fontsize=7, frameon=False)
    ax[0, 1].set_ylabel("pseudobulk log2 (Pvr+1)/(Egfr+1) CPM")
    ax[0, 1].set_title("Switch by lineage (dashed = Dsim)", fontsize=10)

    if len(ww):
        h = ww.pivot_table(index="sample", columns="readout", values="rho").reindex(columns=WOL_READOUTS)
        h.index = h.index.str.replace("_pipseq", "")
        im = ax[1, 0].imshow(h.values, cmap="RdBu_r", vmin=-0.3, vmax=0.3, aspect="auto")
        ax[1, 0].set_xticks(range(h.shape[1])); ax[1, 0].set_xticklabels(h.columns, rotation=45, ha="right", fontsize=7)
        ax[1, 0].set_yticks(range(h.shape[0])); ax[1, 0].set_yticklabels(h.index, fontsize=6)
        sig = ww.pivot_table(index="sample", columns="readout", values="padj").reindex(columns=WOL_READOUTS)
        for i in range(h.shape[0]):
            for j in range(h.shape[1]):
                if sig.values[i, j] < 0.05:
                    ax[1, 0].text(j, i, "*", ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax[1, 0], label="partial rho (| log host UMIs)")
        ax[1, 0].set_title("Within sample: Wolbachia UMIs vs readout (* BH < 0.05)", fontsize=10)

    if len(wq):
        for smp, d in wq.groupby("sample"):
            ax[1, 1].plot(d["quintile"], d["mean_switch_index"], "o-", ms=3,
                          label=f"{smp.replace('_pipseq', '')} ({d['stage'].iloc[0]})")
        ax[1, 1].set_xlabel("Wolbachia load quintile (host-depth adjusted)")
        ax[1, 1].set_ylabel("mean switch index (Pvr - Egfr)")
        ax[1, 1].legend(fontsize=5, frameon=False)
        ax[1, 1].set_title("Switch index by within-sample Wolbachia load", fontsize=10)
    fig.tight_layout()
    savefig(fig, out("rtk_switch.pdf"))
    print("Done ->", a.out_dir)


if __name__ == "__main__":
    main()
