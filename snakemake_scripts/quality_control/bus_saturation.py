#!/usr/bin/env python
"""Sequencing saturation and a downsampling curve from a kallisto|bustools BUS file.

Reads `bustools text -p <output.unfiltered.bus>` on stdin
(columns: barcode, UMI, ec, count). Each BUS record is treated as one
molecule and `count` is the number of reads for that molecule.

Cells are barcodes with >= --min-umis molecules (a simple threshold, not
the knee-based cell calling in the filter script). The downsampling curve
uses the exact expectation for subsampling a fraction p of reads: a
molecule seen c times survives with probability 1 - (1 - p)^c.
"""
import argparse
import sys

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--min-umis", type=int, default=500,
                help="barcodes with >= this many molecules are called cells")
ap.add_argument("--summary", required=True)
ap.add_argument("--curve", required=True)
a = ap.parse_args()

# (barcode, reads-per-molecule) -> number of molecules, built in chunks
parts = []
for chunk in pd.read_csv(sys.stdin, sep="\t", header=None,
                         names=["bc", "umi", "ec", "count"],
                         usecols=["bc", "count"], chunksize=20_000_000):
    parts.append(chunk.groupby(["bc", "count"]).size())
hist = pd.concat(parts).groupby(level=[0, 1]).sum().rename("n").reset_index()

umis_per_bc = hist.groupby("bc")["n"].sum()
cells = umis_per_bc.index[umis_per_bc >= a.min_umis]
n_cells = len(cells)
if n_cells == 0:
    sys.exit(f"No barcodes with >= {a.min_umis} molecules in {a.sample}")

h = hist[hist["bc"].isin(cells)]
c, n = h["count"].to_numpy(), h["n"].to_numpy()
reads_all = int((hist["count"] * hist["n"]).sum())
reads, mols = int((c * n).sum()), int(n.sum())

pd.DataFrame([{
    "sample": a.sample,
    "min_umis": a.min_umis,
    "n_cells": n_cells,
    "reads_all_barcodes": reads_all,
    "frac_reads_in_cells": reads / reads_all,
    "reads_in_cells": reads,
    "umis_in_cells": mols,
    "mean_reads_per_cell": reads / n_cells,
    "median_umis_per_cell": umis_per_bc[cells].median(),
    "saturation": 1 - mols / reads,
}]).to_csv(a.summary, sep="\t", index=False)

rows = []
for p in np.geomspace(0.01, 1, 25):
    umis = (n * (1 - (1 - p) ** c)).sum()
    rows.append({"sample": a.sample, "fraction": p,
                 "reads_per_cell": p * reads / n_cells,
                 "umis_per_cell": umis / n_cells,
                 "saturation": 1 - umis / (p * reads)})
pd.DataFrame(rows).to_csv(a.curve, sep="\t", index=False)
