#!/usr/bin/env python3
"""
composition_by_stage.py
=======================
Rule pseudotime_composition. Which embryonic cell types carry through into
culture? For each lineage (pseudotime_lineages), atlas cell-type composition
(Flysta3D label transfer from integrated.h5ad, confidence-filtered) at each
stage: whole embryo, primary cells, cell line.

Outputs (per label column, e.g. atlas_annotation / atlas_germ_layer)
  composition_<label>.csv          : fraction of confident cells per type,
                                     lineage x stage
  composition_<label>.pdf          : stacked bars, one panel per lineage
  enrichment_<label>.csv / .pdf    : log2(fraction in stage / fraction in
                                     embryo) per lineage, pseudocount 0.005
  confidence_by_stage.csv          : median transfer confidence per sample
"""

import os
import argparse

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from pt_utils import savefig

STAGES = ["embryo", "primary_cells", "cell_culture"]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--integrated", required=True)
    p.add_argument("--lineages", required=True, help="TSV: condition<TAB>lineage")
    p.add_argument("--label_cols", nargs="+", default=["atlas_annotation", "atlas_germ_layer"])
    p.add_argument("--conf_threshold", type=float, default=0.5)
    p.add_argument("--top_n", type=int, default=15)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    a = ad.read_h5ad(args.integrated, backed="r")
    obs = a.obs.copy()
    a.file.close()
    lin = pd.read_csv(args.lineages, sep="\t", index_col=0).iloc[:, 0].to_dict()
    obs["lineage"] = obs["condition"].astype(str).map(lin)
    obs = obs.dropna(subset=["lineage"])
    obs["sample_type"] = pd.Categorical(obs["sample_type"].astype(str), categories=STAGES)
    lineages = sorted(obs["lineage"].unique())

    conf_cols = [f"{c}_confidence" for c in args.label_cols if f"{c}_confidence" in obs]
    if conf_cols:
        (obs.groupby(["lineage", "sample_type", "condition"], observed=True)[conf_cols]
         .median().round(3).to_csv(os.path.join(args.out_dir, "confidence_by_stage.csv")))

    for label in [c for c in args.label_cols if c in obs]:
        conf = f"{label}_confidence"
        d = obs if conf not in obs else obs[obs[conf] >= args.conf_threshold]
        comp = (d.groupby(["lineage", "sample_type", label], observed=True).size()
                .unstack(fill_value=0))
        comp = comp.div(comp.sum(axis=1), axis=0)
        comp.to_csv(os.path.join(args.out_dir, f"composition_{label}.csv"))

        top = comp.mean().nlargest(args.top_n).index
        plot = comp[top].copy()
        plot["other"] = 1 - plot.sum(axis=1)
        colors = sns.color_palette("tab20", plot.shape[1])
        fig, axes = plt.subplots(1, len(lineages), figsize=(3 * len(lineages) + 3, 4.5),
                                 sharey=True, squeeze=False)
        for ax, l in zip(axes[0], lineages):
            sub = plot.loc[l] if l in plot.index.get_level_values(0) else None
            if sub is None:
                ax.axis("off"); continue
            sub.plot(kind="bar", stacked=True, ax=ax, color=colors, legend=False, width=0.8)
            ax.set_title(l, fontsize=9); ax.set_xlabel("")
            ax.tick_params(axis="x", rotation=30)
        axes[0, 0].set_ylabel(f"fraction of cells (conf >= {args.conf_threshold})")
        axes[0, -1].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        fig.suptitle(f"{label} composition by stage", fontsize=10)
        savefig(fig, os.path.join(args.out_dir, f"composition_{label}.pdf"))

        # enrichment relative to the lineage's own embryo
        rows = {}
        for l in lineages:
            if (l, "embryo") not in comp.index:
                continue
            emb = comp.loc[(l, "embryo")]
            for st in ["primary_cells", "cell_culture"]:
                if (l, st) in comp.index:
                    rows[(l, st)] = np.log2((comp.loc[(l, st)] + 0.005) / (emb + 0.005))
        if rows:
            enr = pd.DataFrame(rows).T
            enr.index.names = ["lineage", "sample_type"]
            enr.to_csv(os.path.join(args.out_dir, f"enrichment_{label}.csv"))
            e = enr[top]
            fig, ax = plt.subplots(figsize=(0.45 * len(top) + 3, 0.45 * len(e) + 2))
            lim = np.nanmax(np.abs(e.values))
            sns.heatmap(e, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, ax=ax,
                        annot=True, fmt=".1f", annot_kws={"size": 6},
                        yticklabels=[f"{a} | {b}" for a, b in e.index],
                        cbar_kws={"label": "log2(stage / embryo)"})
            plt.xticks(rotation=45, ha="right")
            ax.set_title(f"{label}: enrichment vs. the lineage's embryo", fontsize=9)
            savefig(fig, os.path.join(args.out_dir, f"enrichment_{label}.pdf"))
        print(f"  {label}: {comp.shape[1]} types")
    print("Done ->", args.out_dir)


if __name__ == "__main__":
    main()
