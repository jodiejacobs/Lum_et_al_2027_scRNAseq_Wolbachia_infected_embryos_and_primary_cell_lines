"""
integrate_via_atlas_projection.py
==================================
Replaces the "rule integrate" joint Harmony/BBKNN re-clustering
(integrate_v2.py) with a single projection of EVERY cell -- embryo AND
primary cell line together -- onto the SAME frozen Flysta3D-v2 atlas
embedding annotate_with_flysta3d_ingest.py uses. One combined object,
atlas_<label> cell-type calls and atlas UMAP coordinates for every cell,
embryo and cell line included and directly comparable to each other for
the first time (no more two-step atlas -> embryo -> cell-line label relay
through map_cellline_to_embryo.py -- every cell here is projected straight
onto the atlas in one shot).

Do the query datasets need to be harmonised against EACH OTHER too?
---------------------------------------------------------------------
No, not for what this object is for. Reference projection means every
query cell's coordinates come from ONE fixed, external transform (the
atlas's own PCA/UMAP) -- no query file ever sees or influences another
query file's placement, so there is no "batch" axis between your own
samples for Harmony to remove here. That is the actual mechanism behind
why projection fixes the "clusters split by dataset" problem: there is no
joint fit left for a dataset identity to leak into. Every sample lands in
the same coordinate system because they're all measured against the same
fixed yardstick, not because they were corrected to agree with each other.

What this object is NOT for: the atlas was fit on ITS OWN biology (fly
developmental cell types from Flysta3D-v2), not on your experimental axis
(Wolbachia titer, infection condition, method/replicate). It has no reason
to be sensitive to the subtle within-cell-type expression shifts your
titer/condition analyses care about -- that signal was never part of what
the atlas's PCs were fit to capture, so it will not show up cleanly in
this embedding no matter how good the label transfer is. Use THIS object
(atlas_<label> + atlas UMAP) to answer "what cell type is this cell, and
how does cell-type composition differ across samples / cell lines vs. the
embryos they came from." Keep using integrate_v2.py's own Harmony
(method + replicate only) step, run on your cells alone, for "does
expression change with titer within a cell type" -- that embedding is
actually fit to be sensitive to your data's real structure. Two different
questions, two different embeddings; don't ask one to do both jobs.
Use analyze_titer_by_annotation.py --groupby atlas_annotation instead of
integrate_v2.py's own leiden-keyed analyze_titer_by_cluster on this object.

Schema compatibility with embryo_to_cellline_trajectory.py
--------------------------------------------------------------
This is meant to be a drop-in replacement for integrate_v2.py's output as
far as that script is concerned, so it also writes: is_embryo (bool, same
"embryo"-in-source_file convention the Snakefile's EMBRYO_SAMPLE_IDS split
already uses), cell_type_<label> (a copy of atlas_<label> -- integrate_v2.py
coalesces atlas_<label>/embryo_<label> into this name via its two-step
relay; every cell here already has atlas_<label> directly, so it's just
copied), a full-gene log1p-normalised .raw (needed for pseudobulk
correlation / marker-module scoring -- adata.X itself is restricted to the
atlas's HVG panel), and obsm['X_umap'] in addition to the explicitly named
obsm['X_umap_atlas']. Two columns integrate_v2.py's output has that this
one deliberately does NOT: 'leiden' (there's no clustering step here --
embryo_to_cellline_trajectory.py's Leiden-cluster-composition section skips
gracefully when it's absent, rather than being given a mislabeled copy of
atlas_<label>) and 'phase' (cell-cycle stage -- comes from rule
annotate_cell_cycle's cyclum-based results/annotated_h5ad/ output, a
different upstream path than the results/filtered_h5ad/ files this script
reads directly; that section skips gracefully too).

Run with:
    mamba activate scanpy
    python snakemake_scripts/method_comparison/integrate_via_atlas_projection.py \\
        --atlas resources/wcoembed_whole_embeding_downsampled_modified.h5ad \\
        --query results/filtered_h5ad/*.h5ad \\
        --out_path results/integrated/integrated.h5ad \\
        --fig_dir results/integrated/figures \\
        --label_cols annotation tissue germ_layer \\
        --flybase_annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz \\
        --ortholog_map reference/orthologs/dmel_dsim_orthologs_rbh.tsv

Reuses build_reference_embedding / project_query_onto_reference from
annotate_with_flysta3d_ingest.py (same frozen-atlas mechanism, same
tested code path) rather than duplicating the projection logic -- the
only thing genuinely new here is loading query files with ALL of their
obs columns kept (wolbachia_titer, condition, method, replicate, ...)
instead of just method/replicate, since this object is meant to be
analysis-ready, not just a label-transfer intermediate.
"""

import os
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse
import anndata as ad
import scanpy as sc

from annotate_with_flysta3d import (
    load_ortholog_map,
    load_atlas_reference,
    _looks_like_dsim,
    remap_dsim_to_dmel,
    _savefig,
)
from annotate_with_flysta3d_ingest import (
    build_reference_embedding,
    project_query_onto_reference,
)


def load_all_query_files(query_paths, dsim_to_dmel=None, dsim_ids=None, dmel_ids=None):
    """Like annotate_with_flysta3d.load_query_files, but keeps EVERY obs
    column from each original file (wolbachia_titer, condition, method,
    replicate, ...) instead of cherry-picking method/replicate -- this
    object is meant to be the final analysis-ready integrated object, not
    just a label-transfer intermediate."""
    adatas = []
    for path in query_paths:
        print(f"\n   Loading query: {path}")
        adata = sc.read_h5ad(path)
        print(f"   {adata.n_obs} cells x {adata.n_vars} genes")

        if adata.raw is None:
            raise ValueError(
                f"{path} has no .raw -- expected filtered h5ad output with "
                "adata.raw set to pre-normalisation counts."
            )

        raw_X = adata.raw.X
        raw_X = raw_X.tocsr() if scipy.sparse.issparse(raw_X) else scipy.sparse.csr_matrix(raw_X)
        raw_X.data = raw_X.data.astype(np.float32)

        basename = os.path.splitext(os.path.basename(path))[0]
        a = ad.AnnData(X=raw_X, obs=adata.obs.copy(), var=adata.raw.var.copy())
        a.obs["dataset"]     = "query"
        a.obs["source_file"] = basename

        if dsim_to_dmel:
            is_dsim, n_dsim, n_dmel = _looks_like_dsim(a.var_names, dsim_ids, dmel_ids)
            if is_dsim:
                print(f"   {basename}: detected as Dsim ({n_dsim} Dsim IDs vs "
                      f"{n_dmel} Dmel IDs) -- remapping var_names to Dmel "
                      "orthologs before atlas comparison")
                a = remap_dsim_to_dmel(a, dsim_to_dmel, label=basename)

        a.obs_names = [f"{basename}__{bc}" for bc in a.obs_names]
        adatas.append(a)

    print(f"\n-- Concatenating {len(adatas)} query files --")
    query = ad.concat(adatas, join="outer", index_unique=None)
    query.obs_names_make_unique()
    if scipy.sparse.issparse(query.X):
        query.X = query.X.tocsr()
    print(f"   Total query cells: {query.n_obs:,}")
    return query


def plot_diagnostics(query, label_cols, fig_dir):
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir
    query = query.copy()
    query.obsm["X_umap"] = query.obsm["X_umap_atlas"]

    # Embryo vs. cell line, using the same "embryo" substring-in-condition
    # convention the Snakefile's own EMBRYO_SAMPLE_IDS split already uses,
    # so this always agrees with how the rest of the pipeline classifies
    # samples.
    query.obs["_origin"] = np.where(
        query.obs["source_file"].str.contains("embryo", case=False), "embryo", "cell_line",
    )

    print(f"\n-- Diagnostic plots -- writing to {fig_dir}/ --")
    sc.pl.umap(query, color="_origin", save="_origin.pdf",
               title="Cell line vs. embryo, both projected onto the atlas")
    sc.pl.umap(query, color="source_file", save="_source_file.pdf",
               title="Every sample, atlas-projected")

    for col in label_cols:
        atlas_col = f"atlas_{col}"
        if atlas_col in query.obs.columns:
            sc.pl.umap(query, color=atlas_col, save=f"_{col}.pdf",
                       title=f"'{col}' transferred onto every cell")

    # Composition comparison: cell-type proportions per sample -- the
    # actual "where do cell lines land vs. the embryo they came from"
    # answer, as numbers rather than just a picture. Diff two source_file
    # rows in the CSV for any cell-line/parental-embryo pair you care about.
    for col in label_cols:
        atlas_col = f"atlas_{col}"
        if atlas_col not in query.obs.columns:
            continue
        counts = query.obs.groupby(["source_file", atlas_col], observed=True).size().unstack(fill_value=0)
        comp = counts.div(counts.sum(axis=1), axis=0)
        comp.to_csv(os.path.join(fig_dir, f"composition_by_sample_{col}.csv"))

        fig, ax = plt.subplots(figsize=(max(8, len(comp) * 0.6), 6))
        comp.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", legend=False)
        ax.set_ylabel(f"Fraction of cells ({atlas_col})")
        ax.set_title(f"Cell-type composition by sample -- {col}")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        plt.xticks(rotation=45, ha="right")
        _savefig(fig, os.path.join(fig_dir, f"composition_by_sample_{col}.pdf"))

    print(f"   Diagnostic plots complete -- see {fig_dir}/")
    print(f"   composition_by_sample_<col>.csv is the direct answer to \"where do "
          "cell lines land vs. the embryo they came from\" -- pick a cell-line "
          "source_file row and its parental-embryo source_file row and compare.")


def main():
    parser = argparse.ArgumentParser(
        description="Project EVERY query h5ad file (embryo AND primary cell "
                     "line) onto a FROZEN Flysta3D-v2 atlas PCA/UMAP in one "
                     "shot, replacing rule integrate's joint Harmony/BBKNN "
                     "re-clustering with pure reference projection."
    )
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--query", required=True, nargs="+",
                         help="ALL your filtered per-sample h5ad files -- "
                              "embryo AND primary cell line, e.g. "
                              "results/filtered_h5ad/*.h5ad")
    parser.add_argument("--out_path", required=True)
    parser.add_argument("--label_cols", nargs="+", default=None)
    parser.add_argument("--flybase_annotation", default=None)
    parser.add_argument("--ortholog_map", default=None)
    parser.add_argument("--k", type=int, default=30)
    parser.add_argument("--n_pcs", type=int, default=30)
    parser.add_argument("--n_top_genes", type=int, default=3000)
    parser.add_argument("--subsample_ref", type=int, default=None)
    parser.add_argument("--fig_dir", default=None)
    args = parser.parse_args()

    query_paths = []
    for pattern in args.query:
        matches = glob.glob(pattern)
        query_paths.extend(matches if matches else [pattern])
    query_paths = sorted(set(query_paths))
    print(f"Query files ({len(query_paths)}):")
    for p in query_paths:
        print(f"  {p}")

    ref = load_atlas_reference(
        args.atlas, args.label_cols,
        flybase_annotation=args.flybase_annotation,
        subsample_ref=args.subsample_ref,
    )

    dsim_to_dmel, dsim_ids, dmel_ids = (
        load_ortholog_map(args.ortholog_map) if args.ortholog_map else ({}, set(), set())
    )
    query = load_all_query_files(
        query_paths, dsim_to_dmel=dsim_to_dmel, dsim_ids=dsim_ids, dmel_ids=dmel_ids,
    )

    # Full-gene, log1p-normalised copy BEFORE the projection step restricts
    # query down to the atlas's HVG panel -- stashed as .raw afterwards, the
    # same convention integrate_v2.py's own preprocess() uses (adata.raw set
    # right after normalize_total+log1p, before HVG subsetting). Needed for
    # anything downstream that wants full-gene expression (pseudobulk
    # correlation, marker-module scoring in embryo_to_cellline_trajectory.py)
    # rather than just the ~n_top_genes atlas HVGs project_query_onto_reference
    # restricts adata.X to.
    print("\n-- Building full-gene log1p-normalised .raw layer --")
    raw_full = query.copy()
    sc.pp.normalize_total(raw_full, target_sum=1e4)
    sc.pp.log1p(raw_full)

    ref = build_reference_embedding(
        ref, n_pcs=args.n_pcs, n_top_genes=args.n_top_genes, n_neighbors=args.k,
    )

    projected = project_query_onto_reference(ref, query, args.label_cols, k=args.k)

    # Reattach by obs_names (not position) -- robust regardless of any
    # internal reordering project_query_onto_reference's gene-padding/concat
    # steps may have done.
    projected.raw = raw_full[projected.obs_names].copy()

    # Name these explicitly rather than leaving them as the scanpy-default
    # 'X_umap'/'X_pca' keys, so it's unambiguous once this sits next to any
    # other embedding (e.g. integrate_v2.py's own Harmony one) on the same
    # object -- but ALSO keep the conventional 'X_umap' key so generic
    # scanpy-based tooling (sc.pl.umap with no explicit basis, etc.) and
    # embryo_to_cellline_trajectory.py's UMAP-overview section work
    # unmodified against this object.
    projected.obsm["X_umap_atlas"] = projected.obsm["X_umap"]
    projected.obsm["X_pca_atlas"]  = projected.obsm["X_pca"]

    # Backward-compat columns for scripts written against integrate_v2.py's
    # output schema (e.g. embryo_to_cellline_trajectory.py):
    #   - is_embryo: boolean, same "embryo" substring-in-source_file
    #     convention the Snakefile's own EMBRYO_SAMPLE_IDS split already
    #     uses, so this always agrees with how the rest of the pipeline
    #     classifies samples.
    #   - cell_type_<label>: integrate_v2.py coalesces atlas_<label>
    #     (embryo cells) and embryo_<label> (cell line cells, transferred
    #     via map_cellline_to_embryo.py) into one cell_type_<label> column
    #     per cell. Every cell here already has atlas_<label> directly (all
    #     cells were projected straight onto the atlas, no relay needed),
    #     so cell_type_<label> is simply a copy -- this just saves
    #     downstream scripts' detect_label_cols()-style helpers from falling
    #     back with a warning.
    # Intentionally NOT adding a fake 'leiden' column (there was no Leiden
    # clustering step here) -- embryo_to_cellline_trajectory.py's Leiden-
    # cluster-composition section already skips gracefully when 'leiden' is
    # absent, and mislabeling atlas_<label> as 'leiden' would be misleading
    # to anyone inspecting this object later. Similarly, cell-cycle 'phase'
    # is NOT present -- that comes from rule annotate_cell_cycle's
    # cyclum-based output (results/annotated_h5ad/), a different upstream
    # path than the results/filtered_h5ad/ files this script reads directly;
    # the trajectory script's cell-cycle section also skips gracefully.
    projected.obs["is_embryo"] = projected.obs["source_file"].astype(str).str.contains(
        "embryo", case=False,
    )
    if args.label_cols:
        for col in args.label_cols:
            atlas_col = f"atlas_{col}"
            if atlas_col in projected.obs.columns:
                projected.obs[f"cell_type_{col}"] = projected.obs[atlas_col]

    out_dir = os.path.dirname(args.out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    projected.write(args.out_path)
    print(f"\nWrote combined atlas-projected object -> {args.out_path}  "
          f"({projected.n_obs:,} cells x {projected.n_vars:,} genes)")

    if args.fig_dir:
        plot_diagnostics(projected, args.label_cols, args.fig_dir)

    print("\n" + "=" * 60)
    print("COMPLETE (atlas-projected integration)")
    print("=" * 60)
    print(f"-> {args.out_path}")
    print("Every cell -- embryo and cell line -- now carries atlas_<label> / "
          "cell_type_<label> / atlas_<label>_confidence, is_embryo, a "
          "full-gene .raw layer, and obsm['X_umap']/['X_umap_atlas'] in the "
          "same frozen coordinate system -- schema-compatible with "
          "embryo_to_cellline_trajectory.py as-is. Two things it won't have: "
          "'leiden' (no clustering step here -- that analysis section skips "
          "gracefully) and 'phase' (cell-cycle -- comes from a different "
          "upstream rule; that section skips gracefully too). Use "
          "analyze_titer_by_annotation.py --groupby atlas_annotation for "
          "titer analysis instead of integrate_v2.py's leiden-keyed "
          "analyze_titer_by_cluster.")


if __name__ == "__main__":
    main()
