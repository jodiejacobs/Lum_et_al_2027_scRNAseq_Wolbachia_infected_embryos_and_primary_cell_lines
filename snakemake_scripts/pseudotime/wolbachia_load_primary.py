#!/usr/bin/env python3
"""
wolbachia_load_primary.py
=========================
Rule pseudotime_wolbachia_load. Does more Wolbachia in primary culture go with
a more cell-line-like, proliferative, or death-resistant state?

Load = per-cell Wolbachia gene UMI fraction (Wolbachia gene UMIs / all UMIs,
from infection_by_stage.py). The rRNA-based wolbachia_titer saturates near 1
in D. simulans, so it is reported but not tested.

  1. load_by_sample.csv : per lineage x stage, median per-cell load (IQR),
     % cells with >= --min_umis Wolbachia UMIs, primary/embryo load ratio,
     and the primary cell line's % proliferating and median cell-line
     signature. Lineages are confounded with host species AND Wolbachia strain
     (all Dsim carry wRi, all Dmel wMel), n = 4: descriptive only.
  2. within_culture.csv : within each primary cell line (cells as units),
     partial Spearman rho between load and each outcome, controlling for
     log host UMIs, then also for marker-based cell state. BH across all
     tests; combined across lineages by Stouffer z (equal weights).
     Two load measures: 'frac' (Wolbachia UMIs / all UMIs; Wolbachia / host
     UMIs gives identical ranks, so identical rho) and 'counts' (Wolbachia UMIs,
     with log host UMIs as the covariate, so no ratio). With --h5ad_dir, also
     per-gene log-normalized expression of --genes (default: apoptosis and
     death-inhibitor genes) from prepared_<lineage>.h5ad, to test whether
     heavily infected cells express more pro-apoptotic genes.
  3. quintiles.csv / quintiles.pdf : outcomes by within-culture load quintile.
Caveat: load is a ratio to host UMIs, so cells with less host transcription
(e.g. quiescent cells) score higher; the host-UMI covariate reduces but does
not remove this.
"""
import os
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, t as tdist
from statsmodels.stats.multitest import multipletests

from pt_utils import savefig

OUTCOMES = ["z_proliferation", "cellline_signature", "sig_up", "sig_down", "score_apoptosis",
            "score_injury_JAK_JNK", "score_immune_AMP", "score_oxphos"]


def partial_spearman(y, x, Z):
    r = lambda v: pd.Series(v).rank().values
    Zr = np.column_stack([np.ones(len(y))] + [Z[:, j] for j in range(Z.shape[1])])
    res = lambda v: v - Zr @ np.linalg.lstsq(Zr, v, rcond=None)[0]
    ey, ex = res(r(y)), res(r(x))
    rho = np.corrcoef(ey, ex)[0, 1]
    df = len(y) - 2 - (Zr.shape[1] - 1)
    tt = rho * np.sqrt(df / (1 - rho ** 2))
    return rho, 2 * tdist.sf(abs(tt), df), np.sign(rho) * norm.isf(tdist.sf(abs(tt), df))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--infection_dir", required=True)
    p.add_argument("--prolif_dir", required=True)
    p.add_argument("--min_umis", type=int, default=5)
    p.add_argument("--h5ad_dir", default=None, help="results/pseudotime (prepared_<lineage>.h5ad)")
    p.add_argument("--genes", nargs="*", default=["rpr", "hid", "grim", "skl", "Buffy", "Debcl", "Diap1",
                                                  "Diap2", "Dronc", "Drice", "Dcp-1", "p53", "puc"])
    p.add_argument("--out_dir", required=True)
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    cells = pd.read_csv(os.path.join(a.infection_dir, "cells_infection_immune.csv.gz"), index_col=0)
    cells["host_umis"] = (cells["total_umis"] - cells["wolbachia_umis"]).clip(lower=1)
    cells["load"] = cells["wolbachia_gene_frac"]

    # 1. between lineages (descriptive)
    g = cells.groupby(["lineage", "species", "sample_type"])
    by = pd.DataFrame({"n_cells": g.size(), "median_load": g["load"].median(),
                       "q25_load": g["load"].quantile(0.25), "q75_load": g["load"].quantile(0.75),
                       "frac_cells_ge_min_umis": g["wolbachia_umis"].apply(lambda v: (v >= a.min_umis).mean()),
                       "median_rRNA_titer": g["wolbachia_titer"].median(),
                       "frac_proliferating": g["proliferating"].mean()}).reset_index()
    emb = by[by.sample_type == "embryo"].set_index("lineage")["median_load"]
    by["load_ratio_vs_embryo"] = by["median_load"] / by["lineage"].map(emb)

    # 2-3. within each primary cell line
    rows, qrows = [], []
    prim = cells[cells.sample_type == "primary_cells"]
    for lin, d in prim.groupby("lineage"):
        f = os.path.join(a.prolif_dir, lin, "cells.csv.gz")
        if os.path.exists(f):
            extra = pd.read_csv(f, index_col=0)
            d = d.join(extra[[c for c in ["cellline_signature", "sig_up", "sig_down", "score_oxphos"]
                              if c in extra]], how="left")
        genes = []
        if a.h5ad_dir:
            import anndata as ad
            from cell_states import gene_expr
            h = ad.read_h5ad(os.path.join(a.h5ad_dir, lin, f"prepared_{lin}.h5ad"))
            ge = gene_expr(h[h.obs_names.isin(d.index)], a.genes)
            d = d.join(ge, how="left")
            genes = list(ge.columns)
            del h
        ld = np.log(d["host_umis"].values)[:, None]
        st = pd.get_dummies(d["cell_state"], drop_first=True).values.astype(float)
        for measure, x in [("frac", d["load"].values), ("counts", d["wolbachia_umis"].values)]:
            for y in [o for o in OUTCOMES if o in d] + genes:
                v = d[y].astype(float).values
                ok = ~np.isnan(v)
                r1, p1, z1 = partial_spearman(v[ok], x[ok], ld[ok])
                r2, p2, z2 = partial_spearman(v[ok], x[ok], np.column_stack([ld, st])[ok])
                rows.append(dict(lineage=lin, species=d["species"].iloc[0], load_measure=measure,
                                 outcome=y, n=int(ok.sum()), frac_expressing=float((v[ok] > 0).mean()),
                                 rho_host_umis=r1, p_host_umis=p1, z_host_umis=z1,
                                 rho_host_umis_state=r2, p_host_umis_state=p2, z_host_umis_state=z2))
        q = pd.qcut(d["load"].rank(method="first"), 5, labels=[1, 2, 3, 4, 5])
        qq = d.groupby(q, observed=True).agg(median_load=("load", "median"),
                                              pct_proliferating=("proliferating", lambda v: 100 * v.mean()),
                                              **{f"median_{o}": (o, "median") for o in OUTCOMES if o in d})
        qrows.append(qq.reset_index(names="load_quintile").assign(lineage=lin))
        by.loc[(by.lineage == lin) & (by.sample_type == "primary_cells"), "median_cellline_signature"] = \
            d["cellline_signature"].median() if "cellline_signature" in d else np.nan
    by.to_csv(os.path.join(a.out_dir, "load_by_sample.csv"), index=False)
    print("Wolbachia load by lineage and stage (descriptive; species, strain and lineage confounded):\n"
          + by.round(4).to_string(index=False))

    w = pd.DataFrame(rows)
    for k in ["host_umis", "host_umis_state"]:
        w[f"padj_{k}"] = multipletests(w[f"p_{k}"], method="fdr_bh")[1]
    comb = w.groupby(["load_measure", "outcome"]).agg(
        n_lineages=("lineage", "size"),
        n_positive=("rho_host_umis_state", lambda v: int((v > 0).sum())),
        rho_min=("rho_host_umis_state", "min"), rho_max=("rho_host_umis_state", "max"),
        z_stouffer=("z_host_umis_state", lambda v: v.sum() / np.sqrt(len(v)))).reset_index()
    comb["p_stouffer"] = 2 * norm.sf(comb["z_stouffer"].abs())
    comb["padj_stouffer"] = multipletests(comb["p_stouffer"], method="fdr_bh")[1]
    w.to_csv(os.path.join(a.out_dir, "within_culture.csv"), index=False)
    comb.to_csv(os.path.join(a.out_dir, "within_culture_combined.csv"), index=False)
    print("\nWithin-culture partial Spearman (load vs outcome | log host UMIs [+ cell state]):\n"
          + w.drop(columns=[c for c in w if c.startswith("z_")]).round(4).to_string(index=False))
    print("\nCombined across lineages:\n" + comb.round(4).to_string(index=False))

    qt = pd.concat(qrows)
    qt.to_csv(os.path.join(a.out_dir, "quintiles.csv"), index=False)
    show = [c for c in ["pct_proliferating", "median_cellline_signature", "median_score_apoptosis",
                        "median_score_immune_AMP"] if c in qt]
    fig, axes = plt.subplots(1, len(show), figsize=(3.2 * len(show), 3))
    for ax, c in zip(np.atleast_1d(axes), show):
        for lin, d in qt.groupby("lineage"):
            ax.plot(d["load_quintile"].astype(int), d[c], marker="o", label=lin)
        ax.set_xlabel("Wolbachia load quintile\n(within primary cell line)"); ax.set_title(c, fontsize=9)
    np.atleast_1d(axes)[0].legend(fontsize=7)
    fig.tight_layout()
    savefig(fig, os.path.join(a.out_dir, "quintiles.pdf"))
    print("Done ->", a.out_dir)


if __name__ == "__main__":
    main()
