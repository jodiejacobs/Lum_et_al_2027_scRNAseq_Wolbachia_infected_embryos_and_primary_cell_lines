#!/usr/bin/env python3
"""
infection_line_de.py
====================
Rule pseudotime_infection_de. Pseudobulk DE between an infected cell line and
its uninfected parent line (default Dsim6B-wMel vs Dsim6B), using the
per-sample pseudobulk counts written by de_stages.py (replicate samples of
each line are the DESeq2 replicates).

Outputs (per pair, prefix <infected>_vs_<uninfected>)
  de_<pair>.csv          : PyDESeq2 results (design ~infection)
  volcano_<pair>.pdf
  gsea_<pair>.csv        : GSEA prerank on the Wald statistic
  vs_line_program_<pair>.csv / .pdf : does infection push expression along
                           the line-vs-primary axis? Spearman of the infection
                           log2FC with the species' line-vs-primary log2FC,
                           on all genes and on line-vs-primary DE genes
  panel_<pair>.csv       : infection log2FC for immune/AMP and pathway genes
"""
import os
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from pydeseq2.dds import DeseqDataSet

from pt_utils import savefig
from de_stages import filter_genes, run_contrast, volcano, load_gene_sets, gsea_all

PANEL = ["Drs", "DptA", "DptB", "AttA", "AttC", "CecA1", "Mtk", "Def", "IM1", "PGRP-LC",
         "PGRP-LB", "PGRP-SC2", "Rel", "dl", "Dif", "imd", "Tl", "upd3", "Socs36E", "TotA",
         "Pvr", "Pvf2", "sty", "pnt", "kek1", "Mkp3", "stg", "CycB", "hid", "rpr", "grim",
         "Diap1", "Hsp70Aa", "Hsp83", "Ubx", "abd-A", "Antp"]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--de_dir", required=True, help="de_stages.py output directory")
    p.add_argument("--pairs", nargs="+", default=["Dsim6B-wMel:Dsim6B"],
                   help="infected_condition:uninfected_condition")
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

    cnt = pd.read_csv(os.path.join(args.de_dir, "pseudobulk_counts.csv.gz"), index_col=0)
    S = pd.read_csv(os.path.join(args.de_dir, "pseudobulk_samples.csv"), index_col=0)
    sym = pd.read_csv(os.path.join(args.de_dir, "de_pooled_line_vs_primary.csv"),
                      index_col=0)["symbol"].to_dict()
    gene_sets = {} if args.skip_gsea else load_gene_sets(args.gmt, args.gene_set_libraries,
                                                         args.gene_set_organism)
    results = {}
    for pair in args.pairs:
        inf, uninf = pair.split(":")
        md = S[S["condition"].isin([inf, uninf])].copy()
        if (md["condition"] == inf).sum() < 2 or (md["condition"] == uninf).sum() < 2:
            print(f"  {pair}: need >= 2 samples per condition -- skipped")
            continue
        md["infection"] = np.where(md["condition"] == inf, "infected", "uninfected")
        species = md["species"].iloc[0]
        cm = filter_genes(cnt[md.index].T.astype(int))
        dds = DeseqDataSet(counts=cm, metadata=md[["infection"]], design="~infection",
                           refit_cooks=True, n_cpus=args.n_cpus, quiet=True)
        dds.deseq2()
        res = run_contrast(dds, ["infection", "infected", "uninfected"], sym, args.n_cpus)
        key = f"{inf}_vs_{uninf}"
        res.to_csv(os.path.join(out, f"de_{key}.csv"))
        volcano(res, f"{inf} vs {uninf}", os.path.join(out, f"volcano_{key}.pdf"), args.alpha)
        sig = res["padj"] < args.alpha
        print(f"  {key}: {sig.sum()} genes padj < {args.alpha} "
              f"({(sig & (res['log2FoldChange'] > 0)).sum()} up, "
              f"{(sig & (res['log2FoldChange'] < 0)).sum()} down in infected)")
        results[key] = res

        # position along the line-vs-primary axis
        lvp_p = os.path.join(args.de_dir, f"de_{species}_line_vs_primary.csv")
        if os.path.exists(lvp_p):
            lvp = pd.read_csv(lvp_p, index_col=0)
            j = res[["log2FoldChange"]].join(lvp[["log2FoldChange", "padj"]], rsuffix="_lvp",
                                              how="inner").dropna(subset=["log2FoldChange",
                                                                          "log2FoldChange_lvp"])
            rows = []
            for lab, d in [("all_genes", j), ("lvp_DE_genes", j[j["padj"] < 0.05])]:
                r = spearmanr(d["log2FoldChange"], d["log2FoldChange_lvp"])
                rows.append(dict(gene_set=lab, n=len(d), rho=r.correlation, p=r.pvalue))
            pd.DataFrame(rows).to_csv(os.path.join(out, f"vs_line_program_{key}.csv"), index=False)
            print(pd.DataFrame(rows).round(4).to_string(index=False))
            fig, ax = plt.subplots(figsize=(4, 4))
            d = j[j["padj"] < 0.05]
            ax.scatter(d["log2FoldChange_lvp"], d["log2FoldChange"], s=3, alpha=0.4)
            ax.axhline(0, c="k", lw=0.5); ax.axvline(0, c="k", lw=0.5)
            ax.set_xlabel(f"{species} line vs primary log2FC")
            ax.set_ylabel(f"{inf} vs {uninf} log2FC")
            ax.set_title(f"rho = {rows[1]['rho']:.2f} (line-vs-primary DE genes)", fontsize=9)
            savefig(fig, os.path.join(out, f"vs_line_program_{key}.pdf"))

        panel = res[res["symbol"].isin(PANEL)][["symbol", "baseMean", "log2FoldChange", "padj"]]
        panel.to_csv(os.path.join(out, f"panel_{key}.csv"))
        print(panel.round(3).to_string())
    if results and gene_sets:
        gsea_all(results, gene_sets, out, args.gsea_permutations)
    print("Done ->", out)


if __name__ == "__main__":
    main()
