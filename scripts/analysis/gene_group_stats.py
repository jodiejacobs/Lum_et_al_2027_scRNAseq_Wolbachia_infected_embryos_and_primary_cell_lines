#!/usr/bin/env python
"""
Compute per-gene and per-cell count statistics for gene groups defined by
var_name prefixes in an AnnData object -- e.g. Wolbachia (GQ...) vs
Drosophila (FBgn...) gene IDs in a combined host+symbiont alignment.

Outputs three CSVs into --outdir:
    gene_stats.csv     per-gene total counts, n_cells_expressing, group
    group_summary.csv  per-group totals, % of library, and per-cell pct stats
    cell_stats.csv      per-cell total counts and per-group counts/pct

Usage:
    python gene_group_stats.py --input adata.h5ad --outdir stats/ \
        --prefix Wolbachia:GQ --prefix Drosophila:FB

    # if raw counts are in a layer instead of adata.raw:
    python gene_group_stats.py --input adata.h5ad --outdir stats/ \
        --source layer --layer-name counts
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse


def parse_prefixes(prefix_args):
    groups = {}
    for item in prefix_args:
        if ':' not in item:
            raise ValueError(f"--prefix must be LABEL:PREFIX, got '{item}'")
        label, prefix = item.split(':', 1)
        groups[label] = prefix
    return groups


def get_counts_matrix(adata, source, layer_name):
    """Return (X, var_names) for the requested counts source."""
    if source == 'raw':
        if adata.raw is None:
            sys.exit("Error: --source raw requested but adata.raw is None.")
        return adata.raw.X, adata.raw.var_names
    if source == 'layer':
        if layer_name not in adata.layers:
            sys.exit(
                f"Error: layer '{layer_name}' not found. "
                f"Available layers: {list(adata.layers.keys())}"
            )
        return adata.layers[layer_name], adata.var_names
    return adata.X, adata.var_names


def to_flat_array(mat):
    """Flatten a sparse or dense matrix-sum result to a 1D numpy array."""
    return np.asarray(mat).flatten()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', required=True, help='Path to .h5ad AnnData file')
    parser.add_argument('--outdir', required=True, help='Directory to write stats files into')
    parser.add_argument(
        '--prefix', action='append', default=[],
        help='LABEL:PREFIX pair for a gene group, e.g. Wolbachia:GQ. '
             'Repeat for multiple groups. Default: Wolbachia:GQ Drosophila:FB'
    )
    parser.add_argument(
        '--source', choices=['raw', 'layer', 'X'], default='raw',
        help="Where to pull counts from: adata.raw.X ('raw', default), "
             "adata.layers[--layer-name] ('layer'), or adata.X ('X')."
    )
    parser.add_argument(
        '--layer-name', default='counts',
        help="Layer name to use when --source layer (default: 'counts')"
    )
    args = parser.parse_args()

    prefix_args = args.prefix or ['Wolbachia:GQ', 'Drosophila:FB']
    groups = parse_prefixes(prefix_args)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)

    X, var_names = get_counts_matrix(adata, args.source, args.layer_name)
    var_names = pd.Index(var_names)

    # Sanity check: warn if this looks like scaled/normalized data (negative
    # values), which would make count-based stats meaningless.
    sample = X[:1000] if X.shape[0] > 1000 else X
    sample_arr = sample.toarray() if issparse(sample) else np.asarray(sample)
    if (sample_arr < 0).any():
        print(
            f"WARNING: negative values detected in the '{args.source}' matrix -- "
            "this looks like scaled/transformed data, not raw counts. "
            "Stats below may not be meaningful as counts. "
            "Consider --source raw or --source layer with --layer-name counts.",
            file=sys.stderr,
        )

    # Assign each gene to a group by prefix (first match wins; else 'other')
    gene_group = pd.Series('other', index=var_names)
    for label, prefix in groups.items():
        gene_group[var_names.str.startswith(prefix)] = label

    # --- per-gene stats ---
    total_counts = to_flat_array(X.sum(axis=0))
    n_cells_expr = to_flat_array((X > 0).sum(axis=0))
    pct_cells_expr = n_cells_expr / X.shape[0] * 100

    gene_stats = pd.DataFrame({
        'gene': var_names,
        'group': gene_group.values,
        'total_counts': total_counts,
        'n_cells_expressing': n_cells_expr,
        'pct_cells_expressing': pct_cells_expr,
    }).sort_values(['group', 'total_counts'], ascending=[True, False])

    gene_stats_path = outdir / 'gene_stats.csv'
    gene_stats.to_csv(gene_stats_path, index=False)
    print(f"Wrote per-gene stats: {gene_stats_path}")

    # --- per-cell stats ---
    cell_stats = pd.DataFrame(index=adata.obs_names)
    cell_stats['total_counts'] = to_flat_array(X.sum(axis=1))
    for label, prefix in groups.items():
        mask = var_names.str.startswith(prefix)
        group_counts = to_flat_array(X[:, mask].sum(axis=1))
        cell_stats[f'counts_{label}'] = group_counts
        cell_stats[f'pct_{label}'] = np.where(
            cell_stats['total_counts'] > 0,
            group_counts / cell_stats['total_counts'] * 100,
            np.nan,
        )

    cell_stats_path = outdir / 'cell_stats.csv'
    cell_stats.to_csv(cell_stats_path)
    print(f"Wrote per-cell stats: {cell_stats_path}")

    # --- group-level summary ---
    # Gene-level totals (sum of raw counts across all cells, per group)
    group_summary = gene_stats.groupby('group').agg(
        n_genes=('gene', 'count'),
        total_counts=('total_counts', 'sum'),
        mean_counts_per_gene=('total_counts', 'mean'),
    ).reset_index()
    library_total = total_counts.sum()
    group_summary['pct_of_library'] = group_summary['total_counts'] / library_total * 100

    # Per-cell pct stats (e.g. average of pct_GQ across cells), skipping the
    # NaNs assigned above for cells with zero total counts.
    pct_rows = []
    for label in groups:
        pct_col = cell_stats[f'pct_{label}']
        pct_rows.append({
            'group': label,
            'mean_pct_per_cell': pct_col.mean(),
            'median_pct_per_cell': pct_col.median(),
            'std_pct_per_cell': pct_col.std(),
        })
    pct_summary = pd.DataFrame(pct_rows)

    group_summary = group_summary.merge(pct_summary, on='group', how='left')

    summary_path = outdir / 'group_summary.csv'
    group_summary.to_csv(summary_path, index=False)
    print(f"Wrote group summary: {summary_path}")

    print("Done.")


if __name__ == '__main__':
    main()
