#!/usr/bin/env python3
"""
merge_pseudotime.py
===================
Final step (rule pseudotime_merge). Writes
results/pseudotime/integrated_with_pseudotime.h5ad: a copy of
integrated.h5ad (unchanged X/.raw/obsm) with the per-species pseudotime
results added to .obs:

  species, lineage, stage_numeric, pt_in_trajectory,
  sceptic_pseudotime, sceptic_pred_stage, sceptic_prob_<stage>,
  dpt_pseudotime

Cells outside the trajectory (embryo cells whose tissue call is not
represented in culture) keep species/lineage/stage_numeric and
pt_in_trajectory = False; their pseudotime columns are NaN.
integrated.h5ad itself is left untouched.
"""

import argparse

import pandas as pd
import anndata as ad


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--integrated", required=True)
    p.add_argument("--cells_csv", nargs="+", required=True,
                   help="prepare_species.py --out_cells_csv, one per species")
    p.add_argument("--sceptic_obs_csv", nargs="+", required=True,
                   help="run_sceptic_stages.py --out_obs_csv, one per species")
    p.add_argument("--out_h5ad", required=True)
    args = p.parse_args()

    adata = ad.read_h5ad(args.integrated)
    cells = pd.concat([pd.read_csv(f, index_col=0) for f in args.cells_csv])
    pt = pd.concat([pd.read_csv(f, index_col=0) for f in args.sceptic_obs_csv])
    pt = pt.drop(columns=[c for c in ["species", "lineage", "stage_numeric"] if c in pt])

    add = cells[["species", "lineage", "stage_numeric", "pt_in_trajectory"]].join(pt, how="left")
    clash = [c for c in add.columns if c in adata.obs.columns]
    if clash:
        print(f"Overwriting existing obs columns: {clash}")
        adata.obs = adata.obs.drop(columns=clash)
    adata.obs = adata.obs.join(add, how="left")
    adata.obs["pt_in_trajectory"] = adata.obs["pt_in_trajectory"].fillna(False).astype(bool)
    for c in ["species", "lineage", "sceptic_pred_stage"]:
        adata.obs[c] = adata.obs[c].astype("category")

    print(adata.obs.groupby(["species", "sample_type"], observed=True)
          .agg(n=("pt_in_trajectory", "size"), in_traj=("pt_in_trajectory", "sum"),
               median_pt=("sceptic_pseudotime", "median")).to_string())
    adata.write(args.out_h5ad)
    print(f"Wrote {args.out_h5ad}")


if __name__ == "__main__":
    main()
