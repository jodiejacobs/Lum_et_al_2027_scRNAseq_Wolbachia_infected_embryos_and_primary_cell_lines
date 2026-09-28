#!/usr/bin/env python3
"""
silencing_h3k27me3.py
=====================
Rule pseudotime_silencing. Are the genes the cell lines switch off
Polycomb (H3K27me3) targets? Gene categories come from the pseudobulk stage
DE (de_stages.py); H3K27me3 domains come from user-supplied BED files
(dm6 / FlyBase r6 coordinates; "chr" prefixes are ignored), e.g. embryo and
S2/Kc/BG3 modENCODE H3K27me3 peaks from the ENCODE portal.

Gene categories (genes expressed in primary culture: log2CPM >= --min_expr in
every lineage's primary sample; the Hox positive control uses all genes):
  off_all_lines   : log2CPM < 1 in every lineage's line
  down_consistent : line-vs-primary log2FC < -1 in all lineages, pooled padj < 0.05
                    (and not off_all_lines)
  unchanged       : |pooled log2FC| < 0.5
  up_consistent   : log2FC > 0.5 in all lineages, pooled padj < 0.05
A gene is "marked" if >= --min_cov of its body is covered by the BED.

Outputs (per BED label)
  gene_marks.csv                  : per gene: category, length, coverage per BED
  h3k27me3_enrichment.csv         : fraction marked per category; Fisher test vs
                                    unchanged; logistic OR for off/down adjusted
                                    for log gene length and primary expression
  h3k27me3_by_lfc.pdf             : fraction marked by line-vs-primary log2FC bin
Positive control printed to the log: fraction of Hox genes marked (should be
high for embryo H3K27me3; if ~0, check the genome build of the BED).
"""
import os
import re
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import fisher_exact

from pt_utils import savefig

HOX = ["lab", "pb", "Dfd", "Scr", "Antp", "Ubx", "abd-A", "Abd-B"]


def gene_spans(gtf):
    rows = []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            m = re.search(r'gene_id "([^"]+)"', f[8])
            if m:
                rows.append((m.group(1), f[0].removeprefix("chr"), int(f[3]) - 1, int(f[4]), f[2]))
    d = pd.DataFrame(rows, columns=["gene_id", "chrom", "start", "end", "feature"])
    if (d["feature"] == "gene").any():
        d = d[d["feature"] == "gene"]
    return d.groupby("gene_id").agg(chrom=("chrom", "first"), start=("start", "min"),
                                    end=("end", "max"))


def merged_bed(path):
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], comment="#",
                    names=["chrom", "start", "end"])
    b = b[~b["chrom"].astype(str).str.startswith(("track", "browser"))]
    b["chrom"] = b["chrom"].astype(str).str.removeprefix("chr")
    b[["start", "end"]] = b[["start", "end"]].astype(int)
    out = {}
    for c, d in b.sort_values(["chrom", "start"]).groupby("chrom"):
        s, e = d["start"].values, d["end"].values
        ms, me = [s[0]], [e[0]]
        for a, z in zip(s[1:], e[1:]):
            if a <= me[-1]:
                me[-1] = max(me[-1], z)
            else:
                ms.append(a); me.append(z)
        out[c] = (np.array(ms), np.array(me))
    return out


def covered_bp(start, end, iv):
    """bp of [start, end) covered by merged, sorted intervals iv=(s, e)."""
    s, e = iv
    cum = np.concatenate([[0], np.cumsum(e - s)])

    def upto(x):  # covered bp in (-inf, x)
        i = np.searchsorted(s, x, side="right")
        full = cum[np.maximum(i - 1, 0)]
        part = np.where(i > 0, np.minimum(x, e[np.maximum(i - 1, 0)]) - s[np.maximum(i - 1, 0)], 0)
        return np.where(i > 0, full + np.clip(part, 0, None), 0)
    return upto(end) - upto(start)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--de_dir", required=True)
    p.add_argument("--gtf", required=True, help="Dmel GTF (DE gene IDs are Dmel FBgn)")
    p.add_argument("--beds", nargs="+", required=True, help="label=path.bed")
    p.add_argument("--min_expr", type=float, default=2.0)
    p.add_argument("--min_cov", type=float, default=0.5)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    D = args.de_dir

    cnt = pd.read_csv(os.path.join(D, "pseudobulk_counts.csv.gz"), index_col=0)
    S = pd.read_csv(os.path.join(D, "pseudobulk_samples.csv"), index_col=0)
    cpm = np.log2(cnt.div(cnt.sum()) * 1e6 + 1)
    g = cpm.T.groupby([S["lineage"], S["stage"]]).mean().T
    pri, lin = g.xs("primary", axis=1, level=1), g.xs("line", axis=1, level=1)
    pooled = pd.read_csv(os.path.join(D, "de_pooled_line_vs_primary.csv"), index_col=0)
    cons = pd.read_csv(os.path.join(D, "lineage_consistency_line_vs_primary.csv"), index_col=0)
    fc = cons[[c for c in cons if c.startswith("log2FC_")]]

    t = pd.DataFrame({"symbol": pooled["symbol"], "lfc": pooled["log2FoldChange"],
                      "padj": pooled["padj"]})
    t = t.join(pri.min(axis=1).rename("pri_min")).join(lin.max(axis=1).rename("line_max"))
    t["expressed"] = t["pri_min"] >= args.min_expr
    down = fc.index[(fc < -1).all(axis=1)]
    up = fc.index[(fc > 0.5).all(axis=1)]
    t["category"] = "not_expressed_in_primary"
    t.loc[t["expressed"], "category"] = "other"
    e = t["expressed"]
    t.loc[e & (t["lfc"].abs() < 0.5), "category"] = "unchanged"
    t.loc[e & t.index.isin(up) & (t["padj"] < 0.05), "category"] = "up_consistent"
    t.loc[e & t.index.isin(down) & (t["padj"] < 0.05), "category"] = "down_consistent"
    t.loc[t["expressed"] & (t["line_max"] < 1), "category"] = "off_all_lines"
    print(t["category"].value_counts().to_string())

    spans = gene_spans(args.gtf)
    t = t.join(spans, how="inner")
    t["length"] = t["end"] - t["start"]
    rows = []
    for spec in args.beds:
        label, path = spec.split("=", 1)
        iv = merged_bed(path)
        cov = np.zeros(len(t))
        for c, idx in t.groupby("chrom").indices.items():
            if c in iv:
                sub = t.iloc[idx]
                cov[idx] = covered_bp(sub["start"].values, sub["end"].values, iv[c]) / sub["length"].values
        t[f"cov_{label}"] = cov
        mk = cov >= args.min_cov
        t[f"marked_{label}"] = mk
        hox = t["symbol"].isin(HOX)
        print(f"\n[{label}] {len(iv)} chromosomes; Hox genes marked: "
              f"{mk[hox.values].sum()}/{hox.sum()} (positive control)")
        ref = t["category"] == "unchanged"
        X = pd.DataFrame({"log_len": np.log10(t["length"].clip(lower=100)),
                          "pri_expr": t["pri_min"]})
        for cat in ["off_all_lines", "down_consistent", "up_consistent", "unchanged"]:
            sel = t["category"] == cat
            r = dict(bed=label, category=cat, n=int(sel.sum()),
                     frac_marked=float(mk[sel.values].mean()) if sel.any() else np.nan)
            if cat != "unchanged" and sel.sum() >= 5:
                tab = [[mk[sel.values].sum(), (~mk[sel.values]).sum()],
                       [mk[ref.values].sum(), (~mk[ref.values]).sum()]]
                r["fisher_OR_vs_unchanged"], r["fisher_p"] = fisher_exact(tab)
                keep = (sel | ref).values
                Xm = sm.add_constant(X[keep].assign(cat=sel[keep].astype(float)))
                try:
                    fit = sm.Logit(mk[keep].astype(float), Xm).fit(disp=0)
                    r["logit_OR_adj"] = float(np.exp(fit.params["cat"]))
                    r["logit_p"] = float(fit.pvalues["cat"])
                except Exception as e:
                    print(f"  logit failed for {cat}: {e}")
            rows.append(r)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(args.out_dir, "h3k27me3_enrichment.csv"), index=False)
    print("\n" + res.round(4).to_string(index=False))
    t.drop(columns=["start", "end"]).to_csv(os.path.join(args.out_dir, "gene_marks.csv"))

    t = t[t["expressed"]]
    bins = pd.cut(t["lfc"], [-np.inf, -6, -4, -2, -1, -0.5, 0.5, 1, np.inf])
    labels = [s.split("=", 1)[0] for s in args.beds]
    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    for lab in labels:
        f = t.groupby(bins, observed=True)[f"marked_{lab}"].mean()
        ax.plot(range(len(f)), f.values, marker="o", label=lab)
    ax.set_xticks(range(len(f)))
    ax.set_xticklabels([str(i) for i in f.index], rotation=30, fontsize=7)
    ax.set_xlabel("line vs primary log2FC (pooled)")
    ax.set_ylabel(f"fraction of genes marked\n(>= {args.min_cov:.0%} of body)")
    ax.set_title("H3K27me3 marking by expression change in the cell lines", fontsize=9)
    ax.legend(fontsize=7)
    savefig(fig, os.path.join(args.out_dir, "h3k27me3_by_lfc.pdf"))
    print("Done ->", args.out_dir)


if __name__ == "__main__":
    main()
