#!/usr/bin/env python3
"""
export_joint_tradeseq.py
========================
Step 4a (rule pseudotime_joint_export). Builds the cross-species tradeSeq
input for tradeseq_joint.R: one counts matrix in Dmel gene space (1:1
reciprocal-best-hit orthologs only) holding the SAME subsampled cells used
in each per-species tradeSeq fit. Each cell keeps its own species'
SCEPTIC pseudotime, and species becomes the tradeSeq condition.

Genes kept: 1:1 orthologs present in both species' objects and detected in
>= --min_frac of cells in at least one species.
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


def subset_counts(h5ad, cells_csv, species, gene_map=None):
    a = sc.read_h5ad(h5ad)
    cells = pd.read_csv(cells_csv, index_col=0)
    a = a[cells.index]
    b = ad.AnnData(X=sp.csr_matrix(a.layers["counts"]), obs=cells.copy(),
                   var=pd.DataFrame(index=a.var_names))
    if gene_map is not None:
        keep = b.var_names.isin(list(gene_map))
        b = b[:, keep].copy()
        b.var_names = [gene_map[g] for g in b.var_names]
    b.obs["species"] = species
    print(f"  {species}: {b.n_obs} cells x {b.n_vars} genes (Dmel ID space)")
    return b


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dmel_h5ad", required=True)
    p.add_argument("--dmel_cells", required=True)
    p.add_argument("--dsim_h5ad", required=True)
    p.add_argument("--dsim_cells", required=True)
    p.add_argument("--ortholog_map", required=True)
    p.add_argument("--flybase_annotation", default=None)
    p.add_argument("--min_frac", type=float, default=0.05)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    d2m = load_orthologs(args.ortholog_map)
    dmel_ids = set(d2m.values())
    mel = subset_counts(args.dmel_h5ad, args.dmel_cells, "Dmel")
    mel = mel[:, mel.var_names.isin(dmel_ids)].copy()
    sim = subset_counts(args.dsim_h5ad, args.dsim_cells, "Dsim", gene_map=d2m)

    shared = mel.var_names.intersection(sim.var_names)
    frac_mel = np.asarray((mel[:, shared].X > 0).mean(axis=0)).ravel()
    frac_sim = np.asarray((sim[:, shared].X > 0).mean(axis=0)).ravel()
    genes = shared[(frac_mel >= args.min_frac) | (frac_sim >= args.min_frac)]
    print(f"  shared 1:1 orthologs: {len(shared)}; detected in >= "
          f"{args.min_frac:.0%} of one species: {len(genes)}")

    joint = ad.concat([mel[:, genes], sim[:, genes]], join="inner")
    mat = sp.csr_matrix(joint.X.T)
    mat.data = np.round(mat.data).astype(np.int64)
    scipy.io.mmwrite(os.path.join(args.out_dir, "counts_genesXcells.mtx"), mat)

    fb = load_flybase_symbols(args.flybase_annotation)
    pd.DataFrame({"gene": joint.var_names,
                  "label": [fb.get(g, g) for g in joint.var_names]}).to_csv(
        os.path.join(args.out_dir, "genes.tsv"), sep="\t", index=False)
    joint.obs[["species", "sceptic_pseudotime", "condition", "lineage",
               "sample_type"]].to_csv(os.path.join(args.out_dir, "cells.csv"))
    print(f"Wrote joint input: {joint.n_obs} cells x {joint.n_vars} genes -> {args.out_dir}")


if __name__ == "__main__":
    main()
