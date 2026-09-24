#!/usr/bin/env python3
"""
export_joint_tradeseq.py
========================
Step 4a (rule pseudotime_joint_export). Builds the joint tradeSeq input
for one pair of trajectories A and B (two lineages of one species, or one
per species): one counts matrix in Dmel gene space (1:1 reciprocal-best-hit
orthologs) holding the SAME subsampled cells used in each per-trajectory
tradeSeq fit. Each cell keeps its own trajectory's SCEPTIC pseudotime, and
the trajectory name (obs 'group') becomes the tradeSeq condition.

Genes kept: 1:1 orthologs present in both species' objects and detected in
>= --min_frac of cells in at least one trajectory.
"""

import os
import argparse

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

from pt_utils import load_orthologs, load_flybase_symbols


def subset_counts(h5ad, cells_csv, name, gene_map):
    a = sc.read_h5ad(h5ad)
    cells = pd.read_csv(cells_csv, index_col=0)
    a = a[cells.index]
    b = ad.AnnData(X=sp.csr_matrix(a.layers["counts"]), obs=cells.copy(),
                   var=pd.DataFrame(index=a.var_names))
    # native Dsim IDs are mapped; genes already in Dmel FBgn space pass through
    b.var_names = [gene_map.get(g, g) for g in b.var_names]
    b = b[:, b.var_names.isin(set(gene_map.values()))].copy()
    b.obs["group"] = name
    print(f"  {name}: {b.n_obs} cells x {b.n_vars} genes (Dmel ID space)")
    return b


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    for x in ["a", "b"]:
        p.add_argument(f"--{x}_h5ad", required=True, help="run_sceptic_stages.py h5ad")
        p.add_argument(f"--{x}_cells", required=True, help="tradeseq_inputs/cells.csv")
        p.add_argument(f"--{x}_name", required=True, help="trajectory name (tradeSeq condition)")
    p.add_argument("--ortholog_map", required=True)
    p.add_argument("--flybase_annotation", default=None)
    p.add_argument("--min_frac", type=float, default=0.05)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    d2m = load_orthologs(args.ortholog_map)
    ta = subset_counts(args.a_h5ad, args.a_cells, args.a_name, d2m)
    tb = subset_counts(args.b_h5ad, args.b_cells, args.b_name, d2m)

    shared = ta.var_names.intersection(tb.var_names)
    frac_a = np.asarray((ta[:, shared].X > 0).mean(axis=0)).ravel()
    frac_b = np.asarray((tb[:, shared].X > 0).mean(axis=0)).ravel()
    genes = shared[(frac_a >= args.min_frac) | (frac_b >= args.min_frac)]
    print(f"  shared 1:1 orthologs: {len(shared)}; detected in >= "
          f"{args.min_frac:.0%} of one trajectory: {len(genes)}")

    joint = ad.concat([ta[:, genes], tb[:, genes]], join="inner")
    mat = sp.csr_matrix(joint.X.T)
    mat.data = np.round(mat.data).astype(np.int64)
    scipy.io.mmwrite(os.path.join(args.out_dir, "counts_genesXcells.mtx"), mat)

    fb = load_flybase_symbols(args.flybase_annotation)
    pd.DataFrame({"gene": joint.var_names,
                  "label": [fb.get(g, g) for g in joint.var_names]}).to_csv(
        os.path.join(args.out_dir, "genes.tsv"), sep="\t", index=False)
    joint.obs[["group", "sceptic_pseudotime", "condition", "lineage",
               "sample_type"]].to_csv(os.path.join(args.out_dir, "cells.csv"))
    print(f"Wrote joint input: {joint.n_obs} cells x {joint.n_vars} genes -> {args.out_dir}")


if __name__ == "__main__":
    main()
