#!/usr/bin/env python3
"""
infection_by_stage.py
=====================
Rule pseudotime_infection. Why is the antimicrobial-peptide (AMP) / immune
program so strong in primary culture? Separates three candidate causes:
dissociation injury (no bacterial signal needed), Wolbachia load, or
bacterial contamination of the cultures.

From each sample's raw counts (filtered_h5ad .raw: host + Wolbachia + 16S
features from the combined kallisto reference):
  * host genes      = gene_ids in that species' host GTF
  * 16S features    = "16S_*" (the highly-represented 16S panel); labelled from
                      the 16S GTF (seqname / name attributes) and split into
                      Wolbachia vs. other bacteria by label
  * Wolbachia genes = everything else (the strain's genome)

Outputs
  sample_summary.csv            per sample: UMIs, Wolbachia gene fraction,
                                Wolbachia-16S and other-bacteria-16S per million
                                UMIs, cells with any other-bacteria 16S read
  sixteenS_by_sample.csv/.pdf   per-16S-feature counts per million, top features
  sixteenS_labels.csv           16S feature -> label used (check this!)
  wolbachia_by_stage.pdf        per-cell Wolbachia fraction + rRNA titer by stage
  bacteria_by_stage.pdf         other-bacteria 16S per million UMIs by sample
  immune_vs_bacteria.csv/.pdf   within each stage x trajectory: AMP module score
                                vs. per-cell bacterial / Wolbachia signal
                                (needs cell_states.py tables)
Caveat: bacterial reads in droplets can be ambient; the per-cell association
with the AMP score is the more informative readout than totals.
"""

import os
import re
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, mannwhitneyu

from pt_utils import savefig

STAGES = ["embryo", "primary_cells", "cell_culture"]


def gtf_gene_ids(path):
    ids, pat = set(), re.compile(r'gene_id "([^"]+)"')
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                m = pat.search(line)
                if m:
                    ids.add(re.sub(r"^gene-", "", m.group(1)))
    return ids


def sixteen_s_labels(path):
    """16S gene_id -> readable label (first informative attribute, else seqname)."""
    lab = {}
    if not path or not os.path.exists(path):
        return lab
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            gid = attrs.get("gene_id")
            if not gid or gid in lab:
                continue
            cand = [attrs.get(k) for k in ("gene_name", "Name", "product", "organism",
                                           "note", "transcript_id")]
            cand = [c for c in cand if c and c != gid]
            lab[gid] = cand[0] if cand else f[0]
    return lab


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--filtered", nargs="+", required=True, help="filtered_h5ad files")
    p.add_argument("--sample_species", nargs="+", required=True,
                   help="species for each --filtered file (same order)")
    p.add_argument("--host_gtf_dmel", required=True)
    p.add_argument("--host_gtf_dsim", required=True)
    p.add_argument("--condition_sample_type", required=True)
    p.add_argument("--lineages", required=True)
    p.add_argument("--sixteen_s_gtf", default=None)
    p.add_argument("--states", nargs="*", default=[], help="cell_states.py states_*.csv.gz")
    p.add_argument("--infection_status", default=None,
                   help="TSV: condition<TAB>expected (infected/uninfected); overrides names")
    p.add_argument("--min_wolbachia_frac", type=float, default=1e-3,
                   help="sample-level Wolbachia gene UMI fraction called 'infected'")
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    out = args.out_dir
    os.makedirs(out, exist_ok=True)

    host_ids = {"Dmel": gtf_gene_ids(args.host_gtf_dmel), "Dsim": gtf_gene_ids(args.host_gtf_dsim)}
    stype = pd.read_csv(args.condition_sample_type, sep="\t", index_col=0).iloc[:, 0].to_dict()
    lin = pd.read_csv(args.lineages, sep="\t", index_col=0).iloc[:, 0].to_dict()
    labels16 = sixteen_s_labels(args.sixteen_s_gtf)

    rows, per16, cells = [], {}, []
    for path, species in zip(args.filtered, args.sample_species):
        base = os.path.splitext(os.path.basename(path))[0]
        cond = re.sub(r"-\d+_(pipseq|10x)$", "", base)
        a = ad.read_h5ad(path)
        X = a.raw.X if a.raw is not None else a.X
        X = X.tocsc() if sp.issparse(X) else sp.csc_matrix(X)
        v = pd.Index(a.raw.var_names if a.raw is not None else a.var_names)
        is16 = v.str.startswith("16S_")
        ishost = v.isin(host_ids[species]) & ~is16
        iswol = ~ishost & ~is16
        lab = pd.Series([labels16.get(g, g) for g in v[is16]], index=v[is16])
        wol16 = lab.str.contains("wolbachia", case=False).values
        tot = np.asarray(X.sum(axis=1)).ravel()
        c_host = np.asarray(X[:, ishost].sum(axis=1)).ravel()
        c_wol = np.asarray(X[:, iswol].sum(axis=1)).ravel()
        X16 = X[:, is16]
        c16_wol = np.asarray(X16[:, wol16].sum(axis=1)).ravel() if wol16.any() else np.zeros(len(tot))
        c16_oth = np.asarray(X16[:, ~wol16].sum(axis=1)).ravel()
        per16[base] = pd.Series(np.asarray(X16.sum(axis=0)).ravel(), index=v[is16])
        T = tot.sum()
        rows.append(dict(sample=base, condition=cond, species=species,
                         sample_type=stype.get(cond, "unknown"), lineage=lin.get(cond, "none"),
                         n_cells=a.n_obs, total_umis=int(T),
                         host_frac=c_host.sum() / T, wolbachia_gene_frac=c_wol.sum() / T,
                         wolbachia_16S_per_M=1e6 * c16_wol.sum() / T,
                         other_bacteria_16S_per_M=1e6 * c16_oth.sum() / T,
                         frac_cells_any_other_16S=float((c16_oth > 0).mean()),
                         frac_cells_any_wolbachia=float((c_wol > 0).mean())))
        cells.append(pd.DataFrame({
            "wolbachia_gene_frac": c_wol / np.maximum(tot, 1),
            "wolbachia_umis": c_wol, "other_16S_umis": c16_oth, "wolbachia_16S_umis": c16_wol,
            "total_umis": tot}, index=[f"{base}__{bc}" for bc in a.obs_names]))
        print(f"  {base}: {a.n_obs} cells; host {c_host.sum()/T:.1%}, "
              f"Wolbachia genes {c_wol.sum()/T:.2%}, other-16S {1e6*c16_oth.sum()/T:.1f}/M")

    summ = pd.DataFrame(rows)
    summ["sample_type"] = pd.Categorical(summ["sample_type"], categories=STAGES + ["unknown"])
    summ = summ.sort_values(["lineage", "sample_type", "sample"])
    expected = {}
    if args.infection_status and os.path.getsize(args.infection_status) > 0:
        expected = pd.read_csv(args.infection_status, sep="\t", index_col=0).iloc[:, 0].to_dict()
    summ["expected_infection"] = summ["condition"].map(expected).fillna("unspecified")
    summ["observed_infection"] = np.where(summ["wolbachia_gene_frac"] >= args.min_wolbachia_frac,
                                          "infected", "uninfected")
    summ["status_mismatch"] = ((summ["expected_infection"] != "unspecified")
                               & (summ["expected_infection"] != summ["observed_infection"]))
    for _, r in summ[summ["status_mismatch"]].iterrows():
        print(f"WARNING: {r['sample']} expected {r['expected_infection']}, "
              f"observed {r['observed_infection']} (Wolbachia genes {r['wolbachia_gene_frac']:.2%})")
    summ.to_csv(os.path.join(out, "sample_summary.csv"), index=False)
    print("\n" + summ.drop(columns=["condition"]).round(4).to_string(index=False))

    pd.Series({g: labels16.get(g, "(no label found)") for g in
               pd.concat(per16, axis=1).index}, name="label").to_csv(
        os.path.join(out, "sixteenS_labels.csv"))
    m16 = pd.concat(per16, axis=1).fillna(0)
    tot = summ.set_index("sample")["total_umis"]
    m16 = m16.div(tot[m16.columns], axis=1) * 1e6
    m16.index = [f"{labels16.get(g, g)} ({g})" for g in m16.index]
    m16.to_csv(os.path.join(out, "sixteenS_by_sample.csv"))
    order = summ["sample"].tolist()
    top = m16.max(axis=1).nlargest(25).index
    if len(top):
        fig, ax = plt.subplots(figsize=(0.45 * len(order) + 5, 0.3 * len(top) + 2))
        sns.heatmap(np.log10(m16.loc[top, order] + 1), cmap="magma", ax=ax,
                    cbar_kws={"label": "log10(16S UMIs per M + 1)"})
        ax.tick_params(axis="x", labelsize=7, rotation=90); ax.tick_params(axis="y", labelsize=6)
        ax.set_title("16S features by sample (top 25)", fontsize=9)
        savefig(fig, os.path.join(out, "sixteenS_by_sample.pdf"))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, col, t in [(axes[0], "other_bacteria_16S_per_M", "Other-bacteria 16S per M UMIs"),
                       (axes[1], "wolbachia_gene_frac", "Wolbachia gene UMI fraction")]:
        sns.barplot(data=summ, x="lineage", y=col, hue="sample_type", hue_order=STAGES,
                    ax=ax, errorbar=None)
        ax.set_title(t, fontsize=9); ax.tick_params(axis="x", rotation=30)
    savefig(fig, os.path.join(out, "bacteria_by_stage.pdf"))

    cell = pd.concat(cells)
    if args.states:
        st = pd.concat([pd.read_csv(f, index_col=0) for f in args.states])
        cell = st.join(cell, how="inner")
        cell["sample_type"] = pd.Categorical(cell["sample_type"], categories=STAGES)
        cell.to_csv(os.path.join(out, "cells_infection_immune.csv.gz"))

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        sns.boxplot(data=cell, x="trajectory", y="wolbachia_gene_frac", hue="sample_type",
                    hue_order=STAGES, showfliers=False, ax=axes[0])
        axes[0].set_title("Per-cell Wolbachia gene UMI fraction", fontsize=9)
        if "wolbachia_titer" in cell:
            sns.boxplot(data=cell, x="trajectory", y="wolbachia_titer", hue="sample_type",
                        hue_order=STAGES, showfliers=False, ax=axes[1])
            axes[1].set_title("rRNA-based Wolbachia titer (filter step)", fontsize=9)
        for ax in axes:
            ax.tick_params(axis="x", rotation=30)
        savefig(fig, os.path.join(out, "wolbachia_by_stage.pdf"))

        if "score_immune_AMP" in cell:
            res = []
            for (tr, stg), d in cell.groupby(["trajectory", "sample_type"], observed=True):
                if len(d) < 20:
                    continue
                r1 = spearmanr(d["score_immune_AMP"], d["wolbachia_gene_frac"])
                any16 = d["other_16S_umis"] > 0
                row = dict(trajectory=tr, sample_type=stg, n_cells=len(d),
                           rho_AMP_vs_wolbachia_frac=r1.correlation, p_wolbachia=r1.pvalue,
                           frac_cells_other_16S=float(any16.mean()))
                if 5 <= any16.sum() <= len(d) - 5:
                    u = mannwhitneyu(d.loc[any16, "score_immune_AMP"], d.loc[~any16, "score_immune_AMP"])
                    row.update(AMP_with_16S=d.loc[any16, "score_immune_AMP"].mean(),
                               AMP_without_16S=d.loc[~any16, "score_immune_AMP"].mean(),
                               p_AMP_16S=u.pvalue)
                res.append(row)
            res = pd.DataFrame(res)
            res.to_csv(os.path.join(out, "immune_vs_bacteria.csv"), index=False)
            print("\nAMP score vs. Wolbachia / bacterial signal:\n" + res.round(4).to_string(index=False))

            prim = cell[cell["sample_type"] == "primary_cells"]
            if len(prim):
                fig, axes = plt.subplots(1, 2, figsize=(11, 4))
                sns.scatterplot(data=prim, x="wolbachia_gene_frac", y="score_immune_AMP",
                                hue="trajectory", s=4, alpha=0.4, linewidth=0, ax=axes[0])
                axes[0].set_xscale("symlog", linthresh=1e-3)
                axes[0].set_title("Primary cells: AMP score vs Wolbachia fraction", fontsize=9)
                prim = prim.assign(other_16S=np.where(prim["other_16S_umis"] > 0, "yes", "no"))
                sns.boxplot(data=prim, x="trajectory", y="score_immune_AMP", hue="other_16S",
                            showfliers=False, ax=axes[1])
                axes[1].set_title("Primary cells: AMP score by other-bacteria 16S read", fontsize=9)
                axes[1].tick_params(axis="x", rotation=30)
                savefig(fig, os.path.join(out, "immune_vs_bacteria.pdf"))
    print("\nDone ->", out)


if __name__ == "__main__":
    main()
