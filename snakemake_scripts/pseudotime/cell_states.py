#!/usr/bin/env python3
"""
cell_states.py
==============
Rule pseudotime_cell_states. Marker-based cell states, replacing the
Flysta3D embryo-atlas labels for cultured cells. The atlas calls most
cell-line cells "CNS" although they have lost neural/neuroblast markers, and
most primary cells "fat body", driven by the injury/immune secretome.

Part 1 -- module scores and states (every trajectory cell)
  * sc.tl.score_genes for each marker module (MODULES below; override with
    --markers_tsv: columns module, gene) on the log-normalised expression.
  * Module scores are z-scored within the trajectory. cell_state = the
    identity module with the highest z (if > --state_min_z), else
    "unassigned". Two independent flags: proliferating (proliferation z >
    --flag_z) and immune_active (AMP z > --flag_z).
  Outputs: states_<trajectory>.csv.gz (per cell), state_composition.csv/.pdf
  (state x stage x lineage), module_scores_by_sample.csv, module_scores.pdf,
  marker_dotplot_<trajectory>.pdf.

Part 2 -- are there cell-line precursors in the primary culture?
  * line_similarity = Pearson r(cell, cell-line centroid) - r(cell, primary
    centroid), on the trajectory's HVGs. > 0 means the cell is closer to the
    cell line than to its own culture's average.
  * precursor candidates = primary cells with line_similarity > 0 (or, if
    fewer than --min_candidates, the top --top_frac of primary cells).
  * Rank test (precursor_rank_test.csv/.pdf): Spearman rho between each
    candidate set's per-gene shift (vs. the rest of the primary culture) and
    the pseudobulk line-vs-primary log2FC, on non-HVG genes only (candidates
    are picked on HVGs, so this avoids circularity), against random
    primary-cell sets of equal size. Candidate sets: most line-similar, most
    proliferative, least AMP/immune-active (same size each).
  * Candidates vs other primary cells: module scores (Wilcoxon), state mix,
    SCEPTIC P(cell line), exploratory Wilcoxon marker genes, and overlap of
    those markers with the pseudobulk line-vs-primary up genes (Fisher test).
  Outputs: line_similarity_by_stage.pdf, precursors_summary.csv,
  precursor_modules.csv, precursor_markers_<trajectory>.csv.
  Caveat: cell-level tests within one library; exploratory.
"""

import os
import argparse
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests

from pt_utils import savefig

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

STAGES = ["embryo", "primary_cells", "cell_culture"]

# FlyBase symbols. identity modules compete for cell_state; the others are flags/readouts
MODULES = {
    # identity
    "plasmatocyte":      ["Hml", "Pxn", "NimC1", "eater", "crq", "He", "Ppn", "NimB4"],
    "crystal_cell":      ["lz", "PPO1", "PPO2", "peb", "Bc"],
    "lamellocyte":       ["atilla", "PPO3", "mys", "Itgbn", "cher"],
    "fat_body":          ["Lsp2", "Lsp1alpha", "Lsp1beta", "Lsp1gamma", "Fbp1", "Fbp2",
                          "Yp1", "Yp2", "Yp3", "Tsf1"],
    "neural":            ["elav", "nSyb", "brp", "fne", "Syt1", "nrv3", "Appl", "cpo"],
    "neuroblast":        ["dpn", "mira", "wor", "ase", "insc", "klu", "dap"],
    "muscle":            ["Mhc", "Act57B", "Mef2", "Tm1", "Mlc1", "Mlc2", "up", "Prm"],
    "epidermis":         ["grh", "shg", "crb", "Cht2", "Gasp", "ttv"],
    "germline":          ["vas", "nos", "aub", "osk", "pgc", "bam"],
    # flags / readouts
    "proliferation":     ["PCNA", "stg", "CycB", "CycA", "polo", "Mcm2", "Mcm3", "Mcm5",
                          "Mcm6", "Mcm7", "aurB", "Cdk1", "RnrL", "RnrS", "dup"],
    "immune_AMP":        ["Drs", "DptA", "DptB", "AttA", "AttB", "AttC", "CecA1", "CecA2",
                          "CecC", "Mtk", "Def", "IM1", "IM2", "IM3", "Dro"],
    "injury_JAK_JNK":    ["upd3", "upd2", "Socs36E", "Mmp1", "puc", "TotA", "TotM"],
    "apoptosis":         ["hid", "rpr", "grim", "Dronc", "p53", "skl"],
    "hemocyte_core":     ["srp", "Pvr", "gcm", "gcm2"],
}
IDENTITY = ["plasmatocyte", "crystal_cell", "lamellocyte", "fat_body", "neural",
            "neuroblast", "muscle", "epidermis", "germline"]


def load_modules(path):
    if not path:
        return MODULES
    t = pd.read_csv(path, sep="\t")
    return {m: list(d["gene"]) for m, d in t.groupby("module")}


def score_modules(adata, modules, min_genes=3):
    sym = adata.var["symbol"].astype(str) if "symbol" in adata.var else pd.Series(
        adata.var_names, index=adata.var_names)
    s2id = pd.Series(adata.var_names, index=sym.values)
    s2id = s2id[~s2id.index.duplicated()]
    used = {}
    for m, genes in modules.items():
        ids = [s2id[g] for g in genes if g in s2id.index]
        if len(ids) < min_genes:
            print(f"  module {m}: only {len(ids)} genes found -- skipped")
            continue
        sc.tl.score_genes(adata, ids, score_name=f"score_{m}", random_state=0, use_raw=False)
        used[m] = [g for g in genes if g in s2id.index]
    return used


def call_states(obs, used, state_min_z, flag_z):
    cols = [f"score_{m}" for m in used]
    z = (obs[cols] - obs[cols].mean()) / obs[cols].std(ddof=0).replace(0, np.nan)
    z.columns = [c.replace("score_", "z_") for c in cols]
    ident = [f"z_{m}" for m in IDENTITY if m in used]
    obs = pd.concat([obs, z], axis=1)
    if ident:
        best = z[ident].idxmax(axis=1).str.replace("z_", "", regex=False)
        obs["cell_state"] = np.where(z[ident].max(axis=1) > state_min_z, best, "unassigned")
    else:
        print("  WARNING: no identity module scored; cell_state = unassigned")
        obs["cell_state"] = "unassigned"
    if "z_proliferation" in z:
        obs["proliferating"] = z["z_proliferation"] > flag_z
    if "z_immune_AMP" in z:
        obs["immune_active"] = z["z_immune_AMP"] > flag_z
    return obs


def line_similarity(adata):
    genes = adata.var["highly_variable"].values if "highly_variable" in adata.var else \
        np.ones(adata.n_vars, bool)
    X = adata.X[:, genes]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    st = adata.obs["sample_type"].astype(str).values
    cent = {s: X[st == s].mean(axis=0) for s in ["primary_cells", "cell_culture"] if (st == s).any()}
    if len(cent) < 2:
        return None

    def rowcorr(M, v):
        Mc = M - M.mean(axis=1, keepdims=True)
        vc = v - v.mean()
        return (Mc @ vc) / (np.linalg.norm(Mc, axis=1) * np.linalg.norm(vc) + 1e-12)
    return rowcorr(X, cent["cell_culture"]) - rowcorr(X, cent["primary_cells"])


def rank_test(adata, prim_idx, masks, de, n_perm, seed=0):
    """Does a candidate set's expression shift (vs. the rest of the primary
    culture) point the same way as the pseudobulk line-vs-primary change,
    transcriptome-wide? Spearman rho between per-gene mean log-expression
    difference (candidates - other primary cells) and the DE log2FC, on
    EVALUATION genes only: genes tested in the DE that are NOT HVGs, because
    line_similarity (which picks one candidate set) is computed on the HVGs.
    Null: the same statistic for n_perm random primary-cell sets of equal size."""
    hv = adata.var["highly_variable"].values if "highly_variable" in adata.var else \
        np.zeros(adata.n_vars, bool)
    genes = adata.var_names[~hv].intersection(de.index[de["log2FoldChange"].notna()])
    if len(genes) < 100:
        return []
    cols = adata.var_names.get_indexer(genes)
    X = adata[prim_idx].X
    X = X.tocsc()[:, cols] if sp.issparse(X) else sp.csc_matrix(np.asarray(X)[:, cols])
    X = X.tocsr()
    tot = np.asarray(X.sum(axis=0)).ravel()
    n = X.shape[0]
    lfc = de.loc[genes, "log2FoldChange"].values
    rng = np.random.default_rng(seed)

    def stat(mask):
        k = mask.sum()
        s_in = np.asarray(X[mask].sum(axis=0)).ravel()
        diff = s_in / k - (tot - s_in) / (n - k)
        return spearmanr(diff, lfc).correlation

    rows = []
    for name, mask in masks.items():
        mask = np.asarray(mask, bool)
        k = int(mask.sum())
        if k < 10 or n - k < 10:
            continue
        obs_rho = stat(mask)
        null = np.array([stat(np.isin(np.arange(n), rng.choice(n, k, replace=False)))
                         for _ in range(n_perm)])
        rows.append(dict(candidate_set=name, n_cells=k, n_eval_genes=len(genes),
                         rho=obs_rho, null_mean=null.mean(), null_sd=null.std(ddof=1),
                         z=(obs_rho - null.mean()) / (null.std(ddof=1) + 1e-12),
                         p_perm=(1 + (null >= obs_rho).sum()) / (1 + n_perm)))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ads", nargs="+", required=True, help="prepared_<trajectory>.h5ad")
    p.add_argument("--sceptic_obs", nargs="*", default=[], help="sceptic_obs_<trajectory>.csv")
    p.add_argument("--de_dir", default=None,
                   help="main DE dir; line_vs_primary up genes for the precursor overlap test")
    p.add_argument("--markers_tsv", default=None)
    p.add_argument("--state_min_z", type=float, default=0.5)
    p.add_argument("--flag_z", type=float, default=1.0)
    p.add_argument("--min_candidates", type=int, default=30)
    p.add_argument("--top_frac", type=float, default=0.05)
    p.add_argument("--n_perm", type=int, default=200,
                   help="random primary-cell sets for the rank-test null")
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    out = args.out_dir
    os.makedirs(out, exist_ok=True)
    modules = load_modules(args.markers_tsv)
    sceptic = pd.concat([pd.read_csv(f, index_col=0) for f in args.sceptic_obs]) \
        if args.sceptic_obs else pd.DataFrame()

    all_obs, prec_rows, prec_mod_rows, rank_rows = [], [], [], []
    for path in args.h5ads:
        adata = sc.read_h5ad(path)
        if "symbol" not in adata.var:
            adata.var["symbol"] = adata.var_names
        base = os.path.basename(path).removesuffix(".h5ad")
        name = base.removeprefix("prepared_") if base.startswith("prepared_") \
            else str(adata.obs["lineage"].iloc[0])
        species = str(adata.obs["species"].iloc[0])
        print(f"\n=== {name} ({adata.n_obs} cells) ===")
        used = score_modules(adata, modules)
        obs = call_states(adata.obs.copy(), used, args.state_min_z, args.flag_z)
        obs["trajectory"] = name
        ls = line_similarity(adata)
        if ls is not None:
            obs["line_similarity"] = ls
        if not sceptic.empty:
            obs = obs.join(sceptic[[c for c in sceptic.columns if c.startswith("sceptic_")]],
                           how="left")
        keep = ["trajectory", "species", "lineage", "condition", "sample_type", "source_file",
                "wolbachia_titer", "cell_state", "proliferating", "immune_active",
                "line_similarity"] + [c for c in obs if c.startswith(("score_", "z_", "sceptic_"))]
        obs[[c for c in keep if c in obs]].to_csv(os.path.join(out, f"states_{name}.csv.gz"))
        all_obs.append(obs[[c for c in keep if c in obs]])

        # marker dotplot by stage
        a2 = adata.copy()
        a2.obs["stage"] = pd.Categorical(a2.obs["sample_type"].astype(str),
                                         categories=[s for s in STAGES
                                                     if s in set(a2.obs["sample_type"].astype(str))])
        s2id = pd.Series(a2.var_names, index=a2.var["symbol"].astype(str).values)
        s2id = s2id[~s2id.index.duplicated()]
        gm = {m: [s2id[g] for g in used[m][:5]] for m in used}
        a2.var_names = [str(s) for s in a2.var["symbol"]]
        a2.var_names_make_unique()
        id2sym = dict(zip(adata.var_names, a2.var_names))
        gm = {m: [id2sym[i] for i in ids] for m, ids in gm.items()}
        sc.settings.figdir = out
        gm = {m: g for m, g in gm.items() if g}
        if gm:
            try:
                sc.pl.dotplot(a2, gm, groupby="stage", standard_scale="var", show=False,
                              save=f"_markers_{name}.pdf")
            except Exception as e:
                print(f"  WARNING: marker dotplot failed ({e})")

        # ── precursor analysis ────────────────────────────────────────────
        if "line_similarity" not in obs:
            continue
        prim = obs[obs["sample_type"].astype(str) == "primary_cells"]
        cand = prim["line_similarity"] > 0
        rule = "line_similarity > 0"
        if cand.sum() < args.min_candidates:
            thr = prim["line_similarity"].quantile(1 - args.top_frac)
            cand = prim["line_similarity"] >= thr
            rule = f"top {args.top_frac:.0%} line_similarity (>= {thr:.3f})"
        row = dict(trajectory=name, species=species, n_primary=len(prim),
                   n_candidates=int(cand.sum()), frac_candidates=float(cand.mean()),
                   frac_primary_closer_to_line=float((prim["line_similarity"] > 0).mean()),
                   candidate_rule=rule)
        for flag in ["proliferating", "immune_active"]:
            if flag in prim:
                row[f"{flag}_candidates"] = float(prim.loc[cand, flag].mean())
                row[f"{flag}_other_primary"] = float(prim.loc[~cand, flag].mean())
        if "sceptic_prob_cell_culture" in prim:
            row["sceptic_P_line_candidates"] = float(prim.loc[cand, "sceptic_prob_cell_culture"].mean())
            row["sceptic_P_line_other"] = float(prim.loc[~cand, "sceptic_prob_cell_culture"].mean())
        st = prim.loc[cand, "cell_state"].value_counts(normalize=True).round(3)
        row["candidate_states"] = "; ".join(f"{k}:{v}" for k, v in st.head(5).items())
        for m in used:
            c = f"score_{m}"
            if cand.sum() >= 3 and (~cand).sum() >= 3:
                u = mannwhitneyu(prim.loc[cand, c], prim.loc[~cand, c])
                prec_mod_rows.append(dict(trajectory=name, module=m,
                                          mean_candidates=prim.loc[cand, c].mean(),
                                          mean_other=prim.loc[~cand, c].mean(),
                                          diff=prim.loc[cand, c].mean() - prim.loc[~cand, c].mean(),
                                          p=u.pvalue))

        # exploratory marker genes: candidates vs other primary cells
        ap = adata[prim.index].copy()
        ap.obs["group"] = np.where(cand.values, "candidate", "other")
        if (ap.obs["group"] == "candidate").sum() >= 3:
            sc.tl.rank_genes_groups(ap, "group", groups=["candidate"], reference="other",
                                    method="wilcoxon", use_raw=False)
            rg = sc.get.rank_genes_groups_df(ap, "candidate")
            rg["symbol"] = rg["names"].map(dict(zip(adata.var_names, adata.var["symbol"].astype(str))))
            rg.to_csv(os.path.join(out, f"precursor_markers_{name}.csv"), index=False)
            de_p = os.path.join(args.de_dir or "", f"de_{species}_line_vs_primary.csv")
            if args.de_dir and os.path.exists(de_p):
                de = pd.read_csv(de_p, index_col=0)
                universe = set(de.index) & set(rg["names"])
                up_line = set(de.index[(de["padj"] < 0.05) & (de["log2FoldChange"] > 0)]) & universe
                up_cand = set(rg.loc[(rg["pvals_adj"] < 0.05) & (rg["logfoldchanges"] > 0.5),
                                     "names"]) & universe
                both = up_line & up_cand
                tab = [[len(both), len(up_cand - up_line)],
                       [len(up_line - up_cand), len(universe) - len(up_line | up_cand)]]
                orr, fp = fisher_exact(tab, alternative="greater")
                row.update(n_candidate_up_genes=len(up_cand), n_line_up_genes=len(up_line),
                           overlap=len(both), overlap_odds_ratio=orr, overlap_p=fp)
        prec_rows.append(row)

        # transcriptome-wide rank test (evaluation genes disjoint from HVGs)
        de_p = os.path.join(args.de_dir or "", f"de_{species}_line_vs_primary.csv")
        if args.de_dir and os.path.exists(de_p):
            de = pd.read_csv(de_p, index_col=0)
            k = int(cand.sum())
            masks = {"line_similarity_top": cand.values}
            for col, hi in [("score_proliferation", True), ("score_immune_AMP", False)]:
                if col in prim:
                    r = prim[col].rank(ascending=not hi, method="first")
                    masks[f"{'high' if hi else 'low'}_{col.replace('score_', '')}"] = (r <= k).values
            for r in rank_test(adata, prim.index, masks, de, args.n_perm):
                r.update(trajectory=name, species=species)
                rank_rows.append(r)

    obs = pd.concat(all_obs)
    obs["sample_type"] = pd.Categorical(obs["sample_type"].astype(str), categories=STAGES)

    # state composition
    comp = (obs.groupby(["trajectory", "sample_type", "cell_state"], observed=True).size()
            .unstack(fill_value=0))
    comp = comp.div(comp.sum(axis=1), axis=0)
    comp.to_csv(os.path.join(out, "state_composition.csv"))
    lineages = sorted(obs["trajectory"].unique())
    fig, axes = plt.subplots(1, len(lineages), figsize=(3 * len(lineages) + 3, 4.5),
                             sharey=True, squeeze=False)
    colors = dict(zip(comp.columns, sns.color_palette("tab20", comp.shape[1])))
    for ax, l in zip(axes[0], lineages):
        comp.loc[l].plot(kind="bar", stacked=True, ax=ax, legend=False, width=0.8,
                         color=[colors[c] for c in comp.columns])
        ax.set_title(l, fontsize=9); ax.set_xlabel(""); ax.tick_params(axis="x", rotation=30)
    axes[0, 0].set_ylabel("fraction of cells")
    axes[0, -1].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    fig.suptitle("Marker-based cell state by stage", fontsize=10)
    savefig(fig, os.path.join(out, "state_composition.pdf"))

    flags = [f for f in ["proliferating", "immune_active"] if f in obs]
    if flags:
        (obs.groupby(["trajectory", "sample_type"], observed=True)[flags].mean().round(3)
         .to_csv(os.path.join(out, "flags_by_stage.csv")))

    # module scores per sample + per-stage plot
    sc_cols = [c for c in obs if c.startswith("score_")]
    (obs.groupby(["trajectory", "sample_type", "condition", "source_file"], observed=True)[sc_cols]
     .median().round(3).to_csv(os.path.join(out, "module_scores_by_sample.csv")))
    long = obs.melt(id_vars=["trajectory", "sample_type"], value_vars=sc_cols,
                    var_name="module", value_name="score")
    long["module"] = long["module"].str.replace("score_", "", regex=False)
    g = sns.catplot(data=long, x="sample_type", y="score", hue="trajectory", col="module",
                    col_wrap=5, kind="box", showfliers=False, sharey=False, height=2.6,
                    aspect=1.1, order=STAGES)
    g.set_xticklabels(rotation=30)
    savefig(g.figure, os.path.join(out, "module_scores.pdf"))

    # line similarity by stage
    if "line_similarity" in obs:
        fig, ax = plt.subplots(figsize=(8, 4))
        sns.violinplot(data=obs, x="trajectory", y="line_similarity", hue="sample_type",
                       hue_order=STAGES, cut=0, inner=None, ax=ax, density_norm="width")
        ax.axhline(0, c="k", lw=0.6, ls="--")
        ax.set_ylabel("r(cell, line) - r(cell, primary)")
        ax.set_title("Similarity to the cell line vs. own primary culture")
        savefig(fig, os.path.join(out, "line_similarity_by_stage.pdf"))

    if prec_rows:
        pr = pd.DataFrame(prec_rows)
        pr.to_csv(os.path.join(out, "precursors_summary.csv"), index=False)
        print("\nPrecursor candidates in primary culture:\n" + pr.round(3).to_string(index=False))
    if rank_rows:
        rr = pd.DataFrame(rank_rows)[["trajectory", "species", "candidate_set", "n_cells",
                                      "n_eval_genes", "rho", "null_mean", "null_sd", "z", "p_perm"]]
        rr.to_csv(os.path.join(out, "precursor_rank_test.csv"), index=False)
        print("\nRank test vs. pseudobulk line-vs-primary log2FC (non-HVG genes):\n"
              + rr.round(3).to_string(index=False))
        fig, ax = plt.subplots(figsize=(7, 3.5))
        sns.barplot(data=rr, x="trajectory", y="z", hue="candidate_set", ax=ax)
        ax.axhline(0, c="k", lw=0.6)
        ax.set_ylabel("z vs. random primary-cell sets")
        ax.set_title("Primary-cell subsets: shift toward the cell-line expression change",
                     fontsize=9)
        ax.tick_params(axis="x", rotation=20)
        ax.legend(fontsize=7)
        savefig(fig, os.path.join(out, "precursor_rank_test.pdf"))
    if prec_mod_rows:
        pm = pd.DataFrame(prec_mod_rows)
        pm["padj"] = multipletests(pm["p"], method="fdr_bh")[1]
        pm.to_csv(os.path.join(out, "precursor_modules.csv"), index=False)
        piv = pm.pivot(index="module", columns="trajectory", values="diff")
        fig, ax = plt.subplots(figsize=(0.9 * piv.shape[1] + 3, 0.35 * len(piv) + 1.5))
        lim = np.nanmax(np.abs(piv.values))
        sns.heatmap(piv, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, annot=True, fmt=".2f",
                    annot_kws={"size": 7}, ax=ax, cbar_kws={"label": "score: candidates - other"})
        ax.set_title("Precursor candidates vs. other primary cells", fontsize=9)
        savefig(fig, os.path.join(out, "precursor_modules.pdf"))
    print("\nState composition:\n" + (comp * 100).round(1).to_string())
    print("\nDone ->", out)


if __name__ == "__main__":
    main()
