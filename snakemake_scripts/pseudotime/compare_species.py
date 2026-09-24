#!/usr/bin/env python3
"""
compare_species.py
==================
Step 6 (rule pseudotime_compare). Compares two trajectories A and B
(--name_a / --name_b): two lineages of one species, or one D. melanogaster
and one D. simulans lineage. Everything runs in Dmel gene space (Dsim genes
are already remapped to 1:1 RBH orthologs in prepare_species.py; native
Dsim IDs are still mapped through --ortholog_map if present).

  1. Ortholog-level gene overlap: pseudotime-associated genes (per-trajectory
     tradeSeq associationTest) shared vs species-specific (Fisher test over
     the universe of genes tested in both), Wald statistic
     rank correlation, and logFC direction concordance for each stage
     transition (embryo -> primary, primary -> cell line).
  2. Shared GO/pathway dynamics: preranked GSEA per trajectory per transition
     (rank = sign(logFC) * waldStat), then A NES vs B NES.
  3. Joint tradeSeq: conditionTest + shape similarity summarised into gene
     classes (conserved dynamics / diverged shape / level offset only).
  4. NMF programs: Jaccard overlap of top genes between A and B
     programs, Hungarian matching, and usage-along-pseudotime curves for
     matched pairs.
"""

import os
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact, spearmanr
from scipy.optimize import linear_sum_assignment

from pt_utils import savefig, load_orthologs, load_flybase_symbols

warnings.filterwarnings("ignore", category=FutureWarning)
TRANSITIONS = ["embryo_to_primary", "primary_to_cellline"]
A, B = "Dmel", "Dsim"   # trajectory names; set from --name_a / --name_b


def read_ts(d, name, d2m=None):
    path = os.path.join(d, name)
    if not os.path.exists(path):
        print(f"  missing: {path}")
        return None
    df = pd.read_csv(path)
    if d2m is not None:
        df["gene"] = df["gene"].map(lambda g: d2m.get(g, g))
        df = df.dropna(subset=["gene"])
    return df.set_index("gene")


def lfc_col(df):
    return [c for c in df.columns if c.startswith("logFC")][0]


# ─────────────────────────────────────────────────────────────────────────────
# 1. Gene overlap
# ─────────────────────────────────────────────────────────────────────────────

def gene_overlap(mel, sim, sym, out, alpha):
    universe = mel.index.intersection(sim.index)
    m_sig = set(mel.index[mel["padj"] < alpha]) & set(universe)
    s_sig = set(sim.index[sim["padj"] < alpha]) & set(universe)
    both = m_sig & s_sig
    table = [[len(both), len(m_sig - s_sig)],
             [len(s_sig - m_sig), len(universe) - len(m_sig | s_sig)]]
    odds, p = fisher_exact(table, alternative="greater")
    rho, rp = spearmanr(mel.loc[universe, "waldStat"], sim.loc[universe, "waldStat"])
    summary = {"universe": len(universe), f"{A}_sig": len(m_sig), f"{B}_sig": len(s_sig)}
    summary.update(shared_sig=len(both), jaccard=len(both) / max(1, len(m_sig | s_sig)),
                   fisher_odds=odds, fisher_p=p, waldStat_spearman=rho, waldStat_p=rp)
    print("  associationTest overlap:", summary)

    cls = pd.DataFrame(index=universe)
    cls["symbol"] = [sym.get(g, g) for g in universe]
    cls[f"{A}_waldStat"] = mel.loc[universe, "waldStat"]
    cls[f"{B}_waldStat"] = sim.loc[universe, "waldStat"]
    cls[f"{A}_padj"] = mel.loc[universe, "padj"]
    cls[f"{B}_padj"] = sim.loc[universe, "padj"]
    cls["class"] = np.select(
        [cls.index.isin(both), cls.index.isin(m_sig), cls.index.isin(s_sig)],
        ["shared_dynamic", f"{A}_only", f"{B}_only"], "not_dynamic")
    cls.sort_values(["class", f"{A}_waldStat"], ascending=[True, False]).to_csv(
        os.path.join(out, "association_gene_classes.csv"))

    fig, ax = plt.subplots(figsize=(5.5, 5))
    pal = {"shared_dynamic": "#8e44ad", f"{A}_only": "#2980b9",
           f"{B}_only": "#c0392b", "not_dynamic": "#d0d0d0"}
    for c in ["not_dynamic", f"{A}_only", f"{B}_only", "shared_dynamic"]:
        d = cls[cls["class"] == c]
        ax.scatter(np.log10(d[f"{A}_waldStat"] + 1), np.log10(d[f"{B}_waldStat"] + 1),
                   s=4, alpha=0.6, c=pal[c], label=f"{c} ({len(d)})")
    ax.set_xlabel(f"{A} log10(Wald + 1)"); ax.set_ylabel(f"{B} log10(Wald + 1)")
    ax.set_title(f"associationTest, 1:1 orthologs (rho = {rho:.2f})")
    ax.legend(fontsize=7, markerscale=3)
    savefig(fig, os.path.join(out, "association_overlap_scatter.pdf"))
    return summary


def transition_concordance(mel, sim, sym, name, out, alpha):
    u = mel.index.intersection(sim.index)
    lm, ls = lfc_col(mel), lfc_col(sim)
    df = pd.DataFrame({"symbol": [sym.get(g, g) for g in u],
                       f"{A}_logFC": mel.loc[u, lm], f"{B}_logFC": sim.loc[u, ls],
                       f"{A}_padj": mel.loc[u, "padj"], f"{B}_padj": sim.loc[u, "padj"]},
                      index=u)

    def call(lfc, padj):
        return np.where(padj >= alpha, "ns", np.where(lfc > 0, "up", "down"))
    df[f"{A}_call"] = call(df[f"{A}_logFC"], df[f"{A}_padj"])
    df[f"{B}_call"] = call(df[f"{B}_logFC"], df[f"{B}_padj"])
    df.to_csv(os.path.join(out, f"transition_{name}_concordance.csv"))
    ct = pd.crosstab(df[f"{A}_call"], df[f"{B}_call"])
    ct.to_csv(os.path.join(out, f"transition_{name}_crosstab.csv"))
    both = df[(df[f"{A}_call"] != "ns") & (df[f"{B}_call"] != "ns")]
    conc = (both[f"{A}_call"] == both[f"{B}_call"]).mean() if len(both) else np.nan
    rho = spearmanr(df[f"{A}_logFC"], df[f"{B}_logFC"])[0]

    fig, ax = plt.subplots(figsize=(5.5, 5))
    sig_both = (df[f"{A}_call"] != "ns") & (df[f"{B}_call"] != "ns")
    ax.scatter(df.loc[~sig_both, f"{A}_logFC"], df.loc[~sig_both, f"{B}_logFC"],
               s=3, c="#d0d0d0", alpha=0.5)
    ax.scatter(df.loc[sig_both, f"{A}_logFC"], df.loc[sig_both, f"{B}_logFC"],
               s=5, c="#8e44ad", alpha=0.7, label=f"sig in both ({sig_both.sum()})")
    ax.axhline(0, c="grey", lw=0.5); ax.axvline(0, c="grey", lw=0.5)
    ax.set_xlabel(f"{A} logFC"); ax.set_ylabel(f"{B} logFC"); ax.legend(fontsize=7)
    ax.set_title(f"{name.replace('_', ' ')}: rho = {rho:.2f}, "
                 f"direction concordance = {conc:.0%}")
    savefig(fig, os.path.join(out, f"transition_{name}_logFC_scatter.pdf"))
    return dict(transition=name, logFC_spearman=rho, n_sig_both=int(sig_both.sum()),
                direction_concordance=conc)


# ─────────────────────────────────────────────────────────────────────────────
# 2. GSEA
# ─────────────────────────────────────────────────────────────────────────────

def load_gene_sets(gmt, libs, organism):
    import gseapy as gp
    if gmt:
        return {os.path.basename(gmt): gp.read_gmt(gmt)}
    out = {}
    for lib in libs:
        try:
            out[lib] = gp.get_library(name=lib, organism=organism)
        except Exception as e:  # compute nodes are often offline
            print(f"  WARNING: could not fetch {lib} ({e}); pass --gmt for offline use")
    return out


def gsea(mel, sim, sym, gene_sets, out, name, permutations):
    import gseapy as gp
    res = {}
    for sp_, df in [(f"{A}", mel), (f"{B}", sim)]:
        rnk = (np.sign(df[lfc_col(df)]) * df["waldStat"]).rename(index=lambda g: sym.get(g, g))
        rnk = rnk[~rnk.index.duplicated()].dropna().sort_values(ascending=False)
        for lib, gs in gene_sets.items():
            r = gp.prerank(rnk=rnk, gene_sets=gs, min_size=10, max_size=500,
                           permutation_num=permutations, seed=42, threads=4,
                           outdir=None, verbose=False).res2d
            r["library"] = lib
            res.setdefault(sp_, []).append(r)
    if not res:
        return
    merged = None
    for sp_, frames in res.items():
        d = pd.concat(frames)[["library", "Term", "NES", "FDR q-val", "Lead_genes"]]
        d = d.rename(columns={"NES": f"{sp_}_NES", "FDR q-val": f"{sp_}_FDR",
                              "Lead_genes": f"{sp_}_lead_genes"})
        merged = d if merged is None else merged.merge(d, on=["library", "Term"], how="outer")
    for c in [f"{A}_NES", f"{B}_NES", f"{A}_FDR", f"{B}_FDR"]:
        merged[c] = pd.to_numeric(merged[c], errors="coerce")
    sig = (merged[f"{A}_FDR"] < 0.25) | (merged[f"{B}_FDR"] < 0.25)
    merged["pattern"] = np.select(
        [(merged[f"{A}_FDR"] < 0.25) & (merged[f"{B}_FDR"] < 0.25) &
         (np.sign(merged[f"{A}_NES"]) == np.sign(merged[f"{B}_NES"])),
         (merged[f"{A}_FDR"] < 0.25) & (merged[f"{B}_FDR"] < 0.25),
         merged[f"{A}_FDR"] < 0.25, merged[f"{B}_FDR"] < 0.25],
        ["shared_same_direction", "shared_opposite_direction", f"{A}_only", f"{B}_only"], "ns")
    merged.sort_values(f"{A}_NES").to_csv(os.path.join(out, f"gsea_{name}_species.csv"),
                                          index=False)

    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.scatter(merged.loc[~sig, f"{A}_NES"], merged.loc[~sig, f"{B}_NES"], s=4, c="#d0d0d0")
    ax.scatter(merged.loc[sig, f"{A}_NES"], merged.loc[sig, f"{B}_NES"], s=8, c="#8e44ad")
    lab = merged[merged["pattern"].str.startswith("shared")].copy()
    lab["m"] = lab[[f"{A}_NES", f"{B}_NES"]].abs().min(axis=1)
    for _, r in lab.nlargest(12, "m").iterrows():
        ax.annotate(r["Term"][:45], (r[f"{A}_NES"], r[f"{B}_NES"]), fontsize=5)
    ax.axhline(0, c="grey", lw=0.5); ax.axvline(0, c="grey", lw=0.5)
    ax.set_xlabel(f"{A} NES"); ax.set_ylabel(f"{B} NES")
    ax.set_title(f"GSEA, {name.replace('_', ' ')} (purple: FDR < 0.25 in either)")
    savefig(fig, os.path.join(out, f"gsea_{name}_nes_scatter.pdf"))


# ─────────────────────────────────────────────────────────────────────────────
# 3. Joint tradeSeq
# ─────────────────────────────────────────────────────────────────────────────

def joint_summary(joint_dir, classes, out, shape_thr):
    path = os.path.join(joint_dir, "tradeseq_shape_similarity.csv")
    if not os.path.exists(path):
        print(f"  missing: {path}")
        return
    sh = pd.read_csv(path).set_index("gene")
    sh = sh.join(classes[["class"]], how="left")
    sh["joint_class"] = np.select(
        [~sh["condition_sig"].astype(bool),
         sh["shape_r"] >= shape_thr,
         sh["shape_r"] < shape_thr],
        ["same_trajectory", "same_shape_level_offset", "diverged_shape"], "NA")
    sh.to_csv(os.path.join(out, "joint_gene_classes.csv"))
    ct = pd.crosstab(sh["joint_class"], sh["class"].fillna("not_tested"))
    ct.to_csv(os.path.join(out, "joint_vs_association_crosstab.csv"))
    print("  joint classes:\n" + ct.to_string())


# ─────────────────────────────────────────────────────────────────────────────
# 4. NMF matching
# ─────────────────────────────────────────────────────────────────────────────

def nmf_matching(nmf_mel, nmf_sim, d2m, sym, out, n_top):
    tm = pd.read_csv(os.path.join(nmf_mel, f"nmf_top_genes_{A}.csv"))
    ts = pd.read_csv(os.path.join(nmf_sim, f"nmf_top_genes_{B}.csv"))
    ts["gene"] = ts["gene"].map(lambda g: d2m.get(g, g))
    tm, ts = tm[tm["rank"] <= n_top], ts[ts["rank"] <= n_top].dropna(subset=["gene"])
    pm = {p: set(d.gene) for p, d in tm.groupby("program")}
    ps = {p: set(d.gene) for p, d in ts.groupby("program")}
    J = pd.DataFrame([[len(pm[a] & ps[b]) / len(pm[a] | ps[b]) for b in ps] for a in pm],
                     index=list(pm), columns=list(ps))
    J.to_csv(os.path.join(out, "nmf_program_jaccard.csv"))
    r, c = linear_sum_assignment(-J.values)
    pairs = pd.DataFrame({f"{A}_program": J.index[r], f"{B}_program": J.columns[c],
                          "jaccard": J.values[r, c]}).sort_values("jaccard", ascending=False)
    pairs["shared_top_genes"] = [
        ", ".join(sorted({sym.get(g, g) for g in pm[a] & ps[b]})[:15])
        for a, b in zip(pairs[f"{A}_program"], pairs[f"{B}_program"])]
    pairs.to_csv(os.path.join(out, "nmf_program_matches.csv"), index=False)

    fig, ax = plt.subplots(figsize=(0.5 * J.shape[1] + 3, 0.45 * J.shape[0] + 2))
    sns.heatmap(J, cmap="viridis", annot=True, fmt=".2f", annot_kws={"size": 6}, ax=ax)
    ax.set_title(f"NMF program overlap (Jaccard of top {n_top} genes, 1:1 orthologs)")
    savefig(fig, os.path.join(out, "nmf_program_jaccard.pdf"))

    um = pd.read_csv(os.path.join(nmf_mel, f"nmf_usage_along_pt_{A}.csv"))
    us = pd.read_csv(os.path.join(nmf_sim, f"nmf_usage_along_pt_{B}.csv"))
    pooled = lambda u: (u.assign(w=u["mean"] * u["count"]).groupby(["program", "bin_center"])
                        [["w", "count"]].sum().assign(mean=lambda d: d.w / d["count"])
                        .reset_index().query("count >= 20"))
    um, us = pooled(um), pooled(us)
    n = len(pairs)
    ncol = 4
    fig, axes = plt.subplots(int(np.ceil(n / ncol)), ncol, figsize=(4 * ncol, 3 * np.ceil(n / ncol)),
                             squeeze=False)
    for ax, (_, row) in zip(axes.ravel(), pairs.iterrows()):
        a = um[um.program == row[f"{A}_program"]]; b = us[us.program == row[f"{B}_program"]]
        ax.plot(a.bin_center, a["mean"], "o-", ms=3, c="#2980b9", label=row[f"{A}_program"])
        ax.plot(b.bin_center, b["mean"], "o-", ms=3, c="#c0392b", label=row[f"{B}_program"])
        ax.axvline(1, ls="--", c="grey", lw=0.5)
        ax.set_title(f"J = {row.jaccard:.2f}", fontsize=8); ax.legend(fontsize=6)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.supxlabel("SCEPTIC pseudotime"); fig.supylabel("mean usage")
    savefig(fig, os.path.join(out, "nmf_matched_program_curves.pdf"))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a_tradeseq", required=True)
    p.add_argument("--b_tradeseq", required=True)
    p.add_argument("--joint_tradeseq", required=True)
    p.add_argument("--a_nmf", required=True)
    p.add_argument("--b_nmf", required=True)
    p.add_argument("--ortholog_map", required=True)
    p.add_argument("--flybase_annotation", default=None)
    p.add_argument("--gene_set_libraries", nargs="*", default=["GO_Biological_Process_2018"])
    p.add_argument("--gene_set_organism", default="Fly")
    p.add_argument("--gmt", default=None, help="local GMT (Dmel symbols); overrides libraries")
    p.add_argument("--gsea_permutations", type=int, default=1000)
    p.add_argument("--skip_gsea", action="store_true")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--shape_r_threshold", type=float, default=0.7)
    p.add_argument("--nmf_top", type=int, default=50)
    p.add_argument("--name_a", default="Dmel", help="name of the first trajectory (reference)")
    p.add_argument("--name_b", default="Dsim", help="name of the second trajectory")
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    global A, B
    A, B = args.name_a, args.name_b
    out = args.out_dir
    os.makedirs(out, exist_ok=True)

    d2m = load_orthologs(args.ortholog_map)
    sym = load_flybase_symbols(args.flybase_annotation)

    print("\n[1] Ortholog-level gene overlap")
    mel = read_ts(args.a_tradeseq, "tradeseq_association.csv")
    sim = read_ts(args.b_tradeseq, "tradeseq_association.csv", d2m)
    summary = gene_overlap(mel, sim, sym, out, args.alpha)
    classes = pd.read_csv(os.path.join(out, "association_gene_classes.csv"), index_col=0)
    trans = []
    tr_tables = {}
    for t in TRANSITIONS:
        m = read_ts(args.a_tradeseq, f"tradeseq_{t}.csv")
        s = read_ts(args.b_tradeseq, f"tradeseq_{t}.csv", d2m)
        if m is None or s is None:
            continue
        tr_tables[t] = (m, s)
        trans.append(transition_concordance(m, s, sym, t, out, args.alpha))
    pd.DataFrame([summary]).to_csv(os.path.join(out, "association_overlap_summary.csv"), index=False)
    pd.DataFrame(trans).to_csv(os.path.join(out, "transition_concordance_summary.csv"), index=False)

    if not args.skip_gsea:
        print("\n[2] GSEA per transition")
        gene_sets = load_gene_sets(args.gmt, args.gene_set_libraries, args.gene_set_organism)
        if gene_sets:
            for t, (m, s) in tr_tables.items():
                gsea(m, s, sym, gene_sets, out, t, args.gsea_permutations)

    print("\n[3] Joint tradeSeq summary")
    joint_summary(args.joint_tradeseq, classes, out, args.shape_r_threshold)

    print("\n[4] NMF program matching")
    nmf_matching(args.a_nmf, args.b_nmf, d2m, sym, out, args.nmf_top)
    print("\nDone ->", out)


if __name__ == "__main__":
    main()
