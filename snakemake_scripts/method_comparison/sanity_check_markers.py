"""
sanity_check_markers.py
========================
Quick biological sanity check for the atlas_<label> / embryo_<label> calls
written by annotate_with_flysta3d.py, annotate_with_flysta3d_ingest.py, or
map_cellline_to_embryo.py: does a known marker gene for each transferred
cell type actually look enriched in the cells that got that label?

High KNN/ingest confidence only tells you the projection was
self-consistent (query cells near each other got similar votes) -- it
doesn't tell you the label is biologically correct. This is the cheap
complementary check: give it 1+ canonical marker genes per expected cell
type and it reports whether each marker is actually higher, on average, in
cells carrying that label than in every other cell.

You'll need real (label_value, gene_symbol) pairs -- run
annotate_with_flysta3d.py / annotate_with_flysta3d_ingest.py without
--label_cols first (or check its log) to see your atlas's actual
annotation/tissue/germ_layer category names, then pick a marker or two you
already know for a few of them.

Run with:
    mamba activate scanpy
    python snakemake_scripts/method_comparison/sanity_check_markers.py \\
        --files results/embryo_annotated_ingest_test/*.h5ad \\
        --label_col atlas_annotation \\
        --markers neuron:elav muscle:Mef2 epithelium:crb mesoderm:twi \\
        --flybase_annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz \\
        --out_dir results/embryo_annotated_ingest_test/marker_check

Compare the two integration methods head-to-head by pointing this at both
output directories in turn (Harmony's results/embryo_annotated/*.h5ad vs.
ingest's results/embryo_annotated_ingest_test/*.h5ad, same --label_col,
same --markers) and diffing marker_summary.csv.
"""

import os
import glob
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from annotate_with_flysta3d import load_symbol_to_fbgn


def load_and_concat(paths, label_col):
    adatas = []
    for p in paths:
        a = sc.read_h5ad(p)
        if label_col not in a.obs.columns:
            print(f"  SKIP {p}: no '{label_col}' column")
            continue
        if a.raw is None:
            raise ValueError(
                f"{p} has no .raw -- expected filtered/annotated h5ad output "
                "with adata.raw set to pre-normalisation counts (same "
                "convention every other script in this pipeline relies on)."
            )
        raw = ad.AnnData(X=a.raw.X.copy(), obs=a.obs.copy(), var=a.raw.var.copy())
        raw.obs["source_file"] = os.path.splitext(os.path.basename(p))[0]
        adatas.append(raw)

    if not adatas:
        raise SystemExit(f"No input files had a '{label_col}' column -- check --label_col")

    combined = ad.concat(adatas, join="inner", index_unique="-")
    print(f"Loaded {combined.n_obs:,} cells x {combined.n_vars:,} genes from {len(adatas)} files")

    print("Normalising (1e4 per cell) + log1p ...")
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    return combined


def resolve_marker_genes(adata, markers, flybase_annotation):
    """markers: list of 'label:gene_symbol' strings -> {label: [var_names
    actually present in adata]}. Accepts either a gene symbol (resolved via
    --flybase_annotation) or a raw var_name (e.g. an FBgn ID) directly."""
    resolved = {}
    symbol_to_fbgn = None
    for entry in markers:
        if ":" not in entry:
            raise ValueError(f"--markers entry '{entry}' must be 'label:gene_symbol'")
        label, symbol = entry.split(":", 1)

        gene_id = None
        if symbol in adata.var_names:
            gene_id = symbol
        else:
            if symbol_to_fbgn is None and flybase_annotation:
                symbol_to_fbgn = load_symbol_to_fbgn(flybase_annotation)
            fbgn = symbol_to_fbgn.get(symbol) if symbol_to_fbgn else None
            if fbgn and fbgn in adata.var_names:
                gene_id = fbgn

        if gene_id is None:
            print(f"  WARNING: marker '{symbol}' for label '{label}' not found "
                  "in adata.var_names (checked as-is and via symbol->FBgn) -- "
                  "skipping")
            continue
        resolved.setdefault(label, []).append(gene_id)

    return resolved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--files", required=True, nargs="+",
                         help="Annotated h5ad file(s)/glob(s), e.g. "
                              "results/embryo_annotated_ingest_test/*.h5ad")
    parser.add_argument("--label_col", required=True,
                         help="e.g. atlas_annotation, atlas_tissue, "
                              "embryo_annotation")
    parser.add_argument("--markers", required=True, nargs="+",
                         help="label:gene_symbol pairs, e.g. neuron:elav "
                              "muscle:Mef2. Repeat a label for >1 marker.")
    parser.add_argument("--flybase_annotation", default=None,
                         help="Needed if adata.var_names are FBgn IDs but "
                              "your --markers are gene symbols.")
    parser.add_argument("--out_dir", required=True)
    args = parser.parse_args()

    paths = []
    for pattern in args.files:
        matches = glob.glob(pattern)
        paths.extend(matches if matches else [pattern])
    paths = sorted(set(paths))
    print(f"Input files ({len(paths)}):")
    for p in paths:
        print(f"  {p}")

    os.makedirs(args.out_dir, exist_ok=True)
    sc.settings.figdir = args.out_dir

    combined = load_and_concat(paths, args.label_col)
    combined.obs[args.label_col] = combined.obs[args.label_col].astype(str)

    resolved = resolve_marker_genes(combined, args.markers, args.flybase_annotation)
    if not resolved:
        raise SystemExit("None of the requested markers were found -- nothing to plot")
    print(f"\nResolved marker genes: {resolved}")

    all_genes = sorted({g for genes in resolved.values() for g in genes})

    print("\nRendering marker dotplot (rows = cell-type label, columns = "
          "every requested marker) ...")
    sc.pl.dotplot(
        combined, var_names=all_genes, groupby=args.label_col,
        standard_scale="var", save="_marker_dotplot.pdf",
    )

    # Per-marker summary: mean expression inside its own label vs. every
    # other cell, so you don't have to eyeball the dotplot to know if a
    # given marker/label pair is working.
    rows = []
    for label, genes in resolved.items():
        for gene in genes:
            col = combined[:, gene].X
            expr = np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
            in_label  = expr[combined.obs[args.label_col].values == label]
            out_label = expr[combined.obs[args.label_col].values != label]
            rows.append({
                "label": label,
                "gene": gene,
                "n_in_label": int(len(in_label)),
                "mean_in_label": float(in_label.mean()) if len(in_label) else np.nan,
                "mean_outside": float(out_label.mean()) if len(out_label) else np.nan,
            })

    summary = pd.DataFrame(rows)
    summary["enriched"] = summary["mean_in_label"] > summary["mean_outside"]
    summary_path = os.path.join(args.out_dir, "marker_summary.csv")
    summary.to_csv(summary_path, index=False)

    print(f"\n{summary.to_string(index=False)}")
    print(f"\nSaved: {summary_path}")

    n_bad = int((~summary["enriched"]).sum())
    if n_bad:
        print(f"\nWARNING: {n_bad}/{len(summary)} marker(s) are NOT higher "
              "inside their assigned label than outside -- check those "
              "label/marker pairs before trusting the projection for that "
              "cell type (could be the wrong marker, or a real projection "
              "problem).")
    else:
        print(f"\nAll {len(summary)} markers enriched in their assigned "
              "label -- good sign the projection is biologically sane.")


if __name__ == "__main__":
    main()
