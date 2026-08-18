"""
annotate_with_flysta3d_ingest.py
=================================
Same job as annotate_with_flysta3d.py -- transfer Flysta3D-v2 atlas
cell-type labels onto your filtered per-sample h5ad files -- but via
scanpy.tl.ingest projection instead of joint Harmony re-integration.

Why this exists
----------------
annotate_with_flysta3d.py (and map_cellline_to_embryo.py, integrate_by_ref.py,
integrate_v2.py) all use the same recipe: concatenate reference + query,
jointly re-select HVGs, refit PCA from scratch, refit Harmony, then
re-cluster. That's *joint re-integration* -- both the atlas and your cells
move, and however much of the atlas-vs-query technical gap Harmony fails to
remove survives into the embedding as its own axis, which is what shows up
as clusters splitting by dataset instead of by cell type. That gap is
larger than a typical batch effect here: Flysta3D-v2 is a Stereo-seq
(BGI DNB) *spatial* atlas, not droplet scRNA-seq, so on top of the usual
method/replicate batch structure you're also correcting across assay
chemistries (Stereo-seq vs. 10x vs. PIPseq vs. kallisto|bustools).

This script instead does *reference projection*: the atlas's PCA + neighbor
graph + UMAP are fit ONCE, on the atlas alone, and then frozen. Every query
cell is placed into that already-existing space via scanpy.tl.ingest --
which applies the atlas's stored PCA loadings (a fixed linear transform,
not a refit) and then walks each query cell onto the atlas's existing UMAP
manifold using the atlas's own fitted neighbor graph. Atlas cells and the
atlas UMAP itself never move. Because nothing about the query is allowed to
reshape the embedding, "dataset" cannot become its own axis of separation
the way it can under Harmony -- a query cell only ends up wherever it's
transcriptionally most similar to atlas cells, honestly. This also makes
adding a new sample later an O(1) operation: project it onto the same
frozen reference, no need to recompute anything for samples you already
projected (unlike the Harmony arm, where adding one new file changes the
joint fit for everyone).

Trade-off: ingest does NOT do any explicit batch correction -- there's no
Harmony step massaging the query to look more like the reference. If the
technical gap between an assay and the atlas is large enough (e.g. very
different capture sensitivity or library size distribution), query cells
can still land unevenly. Check figures/umap_dataset.pdf here: if projected
cells scatter through the same atlas regions as real atlas cells of the
same type, this is working. If they still clump off to one side, that's
the sign to move up to scVI/scANVI + scArches (explicit batch modeling,
frozen shared latent space) rather than tuning this further.

Everything else -- gene ID harmonisation (Dsim->Dmel ortholog remap, atlas
symbol->FBgn remap, positional-index recovery, ATAC-only zero-count cell
dropping), query file loading/concatenation, and writing labelled copies of
your ORIGINAL per-sample h5ad files back out -- is unchanged from
annotate_with_flysta3d.py and is imported from it directly rather than
duplicated, so both scripts stay in sync on that logic.

Run with (same flags as annotate_with_flysta3d.py, minus --harmony_vars
which is accepted-but-ignored so the Snakemake rule doesn't need editing --
see below):
    mamba activate scanpy
    python snakemake_scripts/method_comparison/annotate_with_flysta3d_ingest.py \\
        --atlas resources/wcoembed_whole_embeding_downsampled_modified.h5ad \\
        --query results/filtered_h5ad/*embryos*.h5ad \\
        --out_dir results/embryo_annotated \\
        --label_cols annotation tissue germ_layer \\
        --flybase_annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz \\
        --ortholog_map reference/orthologs/dmel_dsim_orthologs_rbh.tsv \\
        --fig_dir results/embryo_annotated/figures_ingest

To swap rule annotate_with_atlas over to this script with no Snakefile
edits, in config.yaml:
    annotate_atlas_script: "snakemake_scripts/method_comparison/annotate_with_flysta3d_ingest.py"

Same --label_cols discovery behaviour as annotate_with_flysta3d.py: run
without --label_cols first to print the atlas's candidate obs columns and
exit before doing any expensive work.
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

# Reuse the atlas-loading / gene-ID-harmonisation / query-loading / output
# -writing machinery from the Harmony version instead of duplicating it --
# none of that logic changes between the two integration methods, only how
# reference and query get combined into one embedding does.
from annotate_with_flysta3d import (
    load_ortholog_map,
    load_atlas_reference,
    load_query_files,
    write_annotated_outputs,
    _savefig,
)


# -----------------------------------------------------------------------------
# Step 3 (replaces joint_preprocess_and_harmony) -- fit the reference space
# ONCE, on the atlas alone. This embedding is frozen from here on: nothing
# downstream is allowed to move these coordinates.
# -----------------------------------------------------------------------------

def build_reference_embedding(ref, n_pcs=30, n_top_genes=3000, n_neighbors=30):
    print("\n-- Building frozen reference embedding (atlas only, no query) --")
    ref = ref.copy()

    if scipy.sparse.issparse(ref.X):
        ref.X = ref.X.tocsr()
    else:
        ref.X = scipy.sparse.csr_matrix(ref.X)

    # HVGs on raw counts (seurat_v3), same flavor the Harmony version uses,
    # but computed on the reference ALONE -- this fixed gene panel is what
    # every query file gets projected through below.
    print(f"   Highly variable genes (seurat_v3, atlas-only, top {n_top_genes}) ...")
    sc.pp.highly_variable_genes(ref, flavor="seurat_v3", n_top_genes=n_top_genes)

    print("   Normalising (1e4 per cell) + log1p ...")
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)

    ref = ref[:, ref.var["highly_variable"]].copy()
    print(f"   Reference HVG panel: {ref.n_vars:,} genes")

    # Deliberately NOT scaling (no sc.pp.scale) here. ingest re-applies the
    # reference's own PCA loadings directly to the query's log-normalised
    # values; scaling the reference with its own per-gene mean/std and then
    # NOT being able to re-apply those exact same per-gene stats to the
    # query (ingest doesn't do this for you) is a known footgun -- this
    # mirrors scanpy's own "Integrating data using ingest" tutorial recipe.
    print(f"   PCA ({n_pcs} components) ...")
    sc.pp.pca(ref, n_comps=n_pcs, svd_solver="arpack")

    print(f"   Neighbors (k={n_neighbors}) + UMAP ...")
    sc.pp.neighbors(ref, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(ref)

    print(f"   Reference embedding fit: {ref.n_obs:,} cells x {ref.n_vars:,} HVGs, "
          f"{n_pcs} PCs -- this UMAP is now frozen.")
    return ref


# -----------------------------------------------------------------------------
# Step 4 (replaces knn_label_transfer) -- project query cells into the
# frozen reference space via scanpy.tl.ingest, then compute a KNN
# confidence score in that same space (ingest itself only writes a hard
# label, no confidence, per label_col).
# -----------------------------------------------------------------------------

def project_query_onto_reference(ref, query, label_cols, k=30):
    """Project ALL query cells into ref's frozen PCA/UMAP space.

    Nothing about `ref` is modified by this call -- ref.obsm['X_pca'] /
    ref.obsm['X_umap'] and every reference cell's coordinates are untouched.
    Returns a NEW AnnData (the query cells) carrying obsm['X_pca'] /
    obsm['X_umap'] in the reference's own coordinate system, plus
    atlas_<col> / atlas_<col>_confidence obs columns for every label_col.
    """
    ref_genes = ref.var_names
    n_present = int(query.var_names.isin(ref_genes).sum())
    print(f"\n-- Projecting {query.n_obs:,} query cells onto the frozen "
          f"reference embedding --")
    print(f"   Query genes overlapping reference HVG panel: "
          f"{n_present}/{len(ref_genes)}")
    if n_present < 0.5 * len(ref_genes):
        print("   WARNING: fewer than half the reference HVGs are present in "
              "the query -- projection quality may be poor. Check gene ID "
              "harmonisation (Dsim ortholog remap / FBgn-vs-symbol) above.")

    # ingest applies the reference's PCA loadings matrix directly, so query
    # var_names must exactly match ref_genes, same order. Reindex, zero-
    # filling any reference HVG this particular query set doesn't have
    # (e.g. genes that only survived filtering in some samples).
    common = ref_genes.intersection(query.var_names)
    q = query[:, common].copy()
    missing = ref_genes.difference(common)
    if len(missing):
        print(f"   Zero-filling {len(missing)} reference HVGs absent from "
              "this query set")
        pad = ad.AnnData(
            X=scipy.sparse.csr_matrix((q.n_obs, len(missing)), dtype=np.float32),
            obs=pd.DataFrame(index=q.obs_names.copy()),
            var=pd.DataFrame(index=missing),
        )
        q = ad.concat([q, pad], axis=1, join="outer")
    q = q[:, ref_genes].copy()
    q.obs = query.obs.copy()

    if scipy.sparse.issparse(q.X):
        q.X = q.X.tocsr()
    else:
        q.X = scipy.sparse.csr_matrix(q.X)

    print("   Normalising (1e4 per cell) + log1p (same recipe as reference) ...")
    sc.pp.normalize_total(q, target_sum=1e4)
    sc.pp.log1p(q)

    print(f"   Running sc.tl.ingest (label_cols={label_cols}) ...")
    sc.tl.ingest(q, ref, obs=label_cols, embedding_method=["umap", "pca"])

    # ingest overwrites obs[col] in place with the transferred value; rename
    # to atlas_<col> to match what map_cellline_to_embryo.py / integrate_v2.py
    # expect from this rule.
    q.obs = q.obs.rename(columns={col: f"atlas_{col}" for col in label_cols})

    _add_knn_confidence(q, ref, label_cols, k=k)

    return q


def _add_knn_confidence(query, ref, label_cols, k=30):
    """ingest only writes a hard label transfer -- compute a majority-vote
    confidence score (winning fraction among k nearest reference neighbours
    in the same frozen PCA space ingest used) so downstream figures /
    filtering that already expect atlas_<label>_confidence keep working."""
    from sklearn.neighbors import NearestNeighbors

    print(f"\n-- KNN confidence scoring (k={k}, in reference PCA space) --")
    nbrs = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=-1)
    nbrs.fit(ref.obsm["X_pca"])
    _, indices = nbrs.kneighbors(query.obsm["X_pca"])

    for col in label_cols:
        labels_arr = ref.obs[col].astype(str).values
        confidence = []
        for row_idx in indices:
            neigh = labels_arr[row_idx]
            counts = pd.Series(neigh).value_counts()
            confidence.append(counts.iloc[0] / k)
        query.obs[f"atlas_{col}_confidence"] = confidence

        print(f"\n   Transferred 'atlas_{col}' distribution:")
        print(query.obs[f"atlas_{col}"].value_counts().head(20).to_string())
        mean_conf = float(np.mean(confidence))
        low_conf = float(np.mean(np.array(confidence) < 0.5))
        print(f"   Mean confidence: {mean_conf:.3f}")
        if low_conf > 0.2:
            print(f"   WARNING: {low_conf*100:.1f}% of cells have confidence "
                  f"< 0.5 for '{col}' -- treat this column cautiously, or "
                  "increase --k, or check gene-ID harmonisation upstream")


# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

def plot_diagnostics(ref, query, label_cols, fig_dir):
    """QC plots for the projection. Unlike the Harmony version's diagnostics,
    the UMAP axes here ARE the frozen reference atlas's own UMAP -- these
    plots show literally where query cells land on the atlas's own map, not
    a freshly recomputed joint embedding."""
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    print(f"\n-- Diagnostic plots -- writing to {fig_dir}/ --")

    combined = ad.concat(
        [ref, query], join="outer", index_unique=None,
        label="_source", keys=["reference", "query"],
    )
    combined.obs_names_make_unique()
    combined.obsm["X_umap"] = np.concatenate(
        [ref.obsm["X_umap"], query.obsm["X_umap"]], axis=0
    )

    sc.pl.umap(combined, color="_source", save="_dataset.pdf",
               title="Atlas (frozen ref) vs. projected query cells")

    ref_names   = ref.obs_names
    query_names = query.obs_names

    for col in label_cols:
        atlas_col = f"atlas_{col}"
        conf_col  = f"{atlas_col}_confidence"

        display_col = f"_display_{col}"
        combined.obs[display_col] = pd.Series(index=combined.obs_names, dtype=object)
        combined.obs.loc[ref_names, display_col]   = ref.obs[col].astype(str).values
        combined.obs.loc[query_names, display_col] = query.obs[atlas_col].astype(str).values
        combined.obs[display_col] = combined.obs[display_col].fillna("NA")

        sc.pl.umap(combined, color=display_col, save=f"_{col}.pdf",
                   title=f"'{col}': atlas ground truth + ingest-projected query")

        conf_vals = query.obs[conf_col].astype(float)

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(conf_vals.dropna(), bins=30, color="#2196F3", edgecolor="black", alpha=0.8)
        ax.axvline(conf_vals.mean(), color="red", linestyle="--",
                   label=f"mean={conf_vals.mean():.2f}")
        ax.set_xlabel(f"KNN confidence ({col})")
        ax.set_ylabel("Query cells")
        ax.set_title(f"Ingest label-transfer confidence -- {col}")
        ax.legend()
        _savefig(fig, os.path.join(fig_dir, f"confidence_hist_{col}.pdf"))

        conf_df = pd.DataFrame({
            "confidence":  conf_vals.values,
            "source_file": query.obs["source_file"].values,
        })
        samples = sorted(conf_df["source_file"].unique())
        fig, ax = plt.subplots(figsize=(max(8, len(samples) * 1.2), 5))
        ax.boxplot([conf_df.loc[conf_df["source_file"] == s, "confidence"].dropna().values
                    for s in samples],
                   labels=samples, showfliers=False)
        ax.set_ylabel(f"KNN confidence ({col})")
        ax.set_title(f"Ingest confidence by sample -- {col}")
        ax.set_ylim(0, 1.05)
        plt.xticks(rotation=45, ha="right")
        _savefig(fig, os.path.join(fig_dir, f"confidence_by_sample_{col}.pdf"))

        conf_df.groupby("source_file")["confidence"].agg(
            ["mean", "median", "std", "count"]
        ).to_csv(os.path.join(fig_dir, f"confidence_summary_{col}.csv"))

    print(f"   Diagnostic plots complete -- see {fig_dir}/")
    print(f"   Compare {fig_dir}/umap_dataset.pdf against the Harmony version's "
          "equivalent plot: cells should scatter through the same atlas "
          "regions as real atlas cells of the same type, not clump off to "
          "one side by dataset.")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Project query h5ad files onto a FROZEN Flysta3D-v2 "
                     "atlas PCA/UMAP via scanpy.tl.ingest -- no joint "
                     "Harmony refit, reference cells and the reference UMAP "
                     "never move. Same CLI and outputs as "
                     "annotate_with_flysta3d.py (atlas_<label> / "
                     "atlas_<label>_confidence obs columns) so it drops "
                     "into rule annotate_with_atlas via the "
                     "annotate_atlas_script config key with no Snakefile "
                     "edits -- run both and compare figures/umap_dataset.pdf to "
                     "see which actually separates by cell type instead of "
                     "by dataset."
    )
    parser.add_argument("--atlas", required=True,
                         help="Path to the downloaded Flysta3D-v2 h5ad, e.g. "
                              "resources/wcoembed_whole_embeding_downsampled_modified.h5ad")
    parser.add_argument("--query", required=True, nargs="+",
                         help="Filtered per-sample h5ad file(s), e.g. "
                              "results/filtered_h5ad/*.h5ad")
    parser.add_argument("--out_dir", required=True,
                         help="Directory to write annotated copies of each "
                              "query file (same filenames, new atlas_* obs "
                              "columns added)")
    parser.add_argument("--label_cols", nargs="+", default=None,
                         help="obs column(s) in the atlas to transfer (e.g. "
                              "annotation tissue germ_layer). Omit to print "
                              "candidate columns from the atlas and exit.")
    parser.add_argument("--flybase_annotation", default=None,
                         help="reference/fbgn_annotation_ID_fb_2025_04.tsv.gz "
                              "-- used to remap atlas gene symbols to FBgn IDs "
                              "if the atlas isn't already FBgn-indexed.")
    parser.add_argument("--ortholog_map", type=str, default=None,
                         help="TSV with Dsim/Dmel FlyBase ID columns (reciprocal "
                              "best hit orthologs). If given, query files "
                              "detected as Dsim have their var_names remapped "
                              "to the orthologous Dmel FlyBase ID before "
                              "projection.")
    parser.add_argument("--k", type=int, default=30,
                         help="Reference neighbours used both for the "
                              "neighbor graph the UMAP is built from and for "
                              "the post-hoc KNN confidence score.")
    parser.add_argument("--n_pcs", type=int, default=30)
    parser.add_argument("--n_top_genes", type=int, default=3000,
                         help="HVGs selected on the REFERENCE ALONE (unlike "
                              "the Harmony version, which selects them "
                              "jointly with the query) -- this fixed gene "
                              "panel is what the query is projected through.")
    parser.add_argument("--harmony_vars", nargs="+", default=None,
                         help="Accepted for CLI compatibility with "
                              "annotate_with_flysta3d.py / rule "
                              "annotate_with_atlas's shell command (which "
                              "always passes --harmony_vars) -- IGNORED, "
                              "this script does no Harmony correction at all.")
    parser.add_argument("--subsample_ref", type=int, default=None,
                         help="Subsample the atlas to this many cells before "
                              "building the reference embedding -- use for a "
                              "fast first pass.")
    parser.add_argument("--fig_dir", type=str, default=None,
                         help="If given, write QC plots here. Omit to skip "
                              "plotting.")

    args = parser.parse_args()

    if args.harmony_vars:
        print(f"NOTE: --harmony_vars {args.harmony_vars} was passed but is "
              "ignored -- this script (ingest projection) does no Harmony "
              "correction.")

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
    query, originals = load_query_files(
        query_paths, dsim_to_dmel=dsim_to_dmel, dsim_ids=dsim_ids, dmel_ids=dmel_ids,
    )

    ref = build_reference_embedding(
        ref, n_pcs=args.n_pcs, n_top_genes=args.n_top_genes, n_neighbors=args.k,
    )

    projected = project_query_onto_reference(
        ref, query, args.label_cols, k=args.k,
    )

    if args.fig_dir:
        plot_diagnostics(ref, projected, args.label_cols, args.fig_dir)

    atlas_cols = [c for c in projected.obs.columns if c.startswith("atlas_")]
    transferred_df = projected.obs[atlas_cols].copy()
    transferred_df.index = projected.obs_names

    write_annotated_outputs(originals, transferred_df, args.out_dir)

    print("\n" + "=" * 60)
    print("COMPLETE (ingest projection)")
    print("=" * 60)
    print(f"Annotated files -> {args.out_dir}/")
    print("Every query cell's atlas_<label> / atlas_<label>_confidence was "
          "assigned by projecting it onto the FROZEN atlas UMAP. Point rule "
          "integrate / integrate_uninfected at this directory the same way "
          "you already do for the Harmony version.")


if __name__ == "__main__":
    main()
