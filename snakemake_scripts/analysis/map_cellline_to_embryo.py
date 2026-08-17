"""
map_cellline_to_embryo.py
==========================
Maps the primary cell line samples (everything in samples.csv that ISN'T an
embryo sample) onto the Flysta3D-v2-annotated embryo cells (rule
annotate_with_atlas / annotate_with_flysta3d.py output), transferring each
embryo cell's atlas-derived cell-type label(s) onto its nearest cell line
cells via Harmony + KNN.

Why map cell lines onto embryos rather than re-running the atlas transfer
directly on the cell lines
-------------------------------------------------------------------------
The Flysta3D-v2 atlas is a whole-embryo spatiotemporal atlas -- every
reference cell comes from a staged embryo, and its "cell_type"/"tissue"/
"germ_layer" labels describe embryonic cell identities. Cultured primary
cell lines are not embryos: transferring atlas labels straight onto them
would force each cultured cell into the closest embryonic identity in the
atlas's own batch structure, which conflates two different batch-correction
problems in one Harmony run (atlas vs. embryo-sample vs. cell-line-sample).
Instead we transfer through your OWN already-annotated embryo cells (one
correction: your embryo samples vs. your cell line samples, both captured
the same way you always capture data), so cell type calls line up with
what's actually already known about your own embryo samples' composition,
not the atlas's raw embedding.

Design (mirrors annotate_with_flysta3d.py's own atlas-vs-query recipe)
-------------------------------------------------------------------------
1. Load all "reference" files (per-sample output of rule annotate_with_atlas
   / annotate_with_flysta3d.py -- your OWN embryo h5ads, each already
   carrying atlas_<label> / atlas_<label>_confidence obs columns) in one
   pass and concatenate them, tagged dataset="embryo_reference".
2. Load all "query" files (filtered per-sample h5ads for the remaining,
   non-embryo / cell-line samples) in one pass and concatenate them, tagged
   dataset="cellline_query". D. simulans query samples (detected by gene ID
   overlap against --ortholog_map) have their var_names remapped onto the
   orthologous D. melanogaster FlyBase ID first, so every sample -- embryo
   reference or cell line query, Dmel or Dsim -- ends up in the same Dmel
   FBgn gene-ID space.
3. Concatenate reference + query on shared genes, normalise + log1p,
   restrict to jointly-computed HVGs, scale, PCA, Harmony (batch =
   embryo_reference vs. cellline_query, plus method/replicate if present).
4. KNN in Harmony PCA space: each cell line query cell gets the
   majority-vote label (+ confidence = winning fraction of k neighbours)
   from its k nearest EMBRYO REFERENCE cells, for every atlas_<label>
   column found on the reference files.
5. Split the labelled query cells back out by source_file and write each
   ORIGINAL cell line h5ad back out with new embryo_<label> /
   embryo_<label>_confidence obs columns added -- nothing else changes.

Run with:
    mamba activate scanpy
    python snakemake_scripts/analysis/map_cellline_to_embryo.py \\
        --reference results/embryo_annotated/*.h5ad \\
        --query results/filtered_h5ad/JW18wMel*.h5ad \\
                results/filtered_h5ad/ubkhc-[0-9]*.h5ad \\
                results/filtered_h5ad/Dsim*.h5ad \\
        --out_dir results/celllines_mapped_to_embryo \\
        --label_cols cell_type tissue germ_layer \\
        --ortholog_map reference/orthologs/dmel_dsim_orthologs_rbh.tsv

--label_cols should match whatever --label_cols you actually used for rule
annotate_with_atlas (the atlas_<label> columns it wrote). Omit --label_cols
to auto-detect every atlas_<label> column present on the reference files.
"""

import os
import glob
import argparse

import numpy as np
import pandas as pd
import scipy.sparse
import anndata as ad
import scanpy as sc
import harmonypy as hm


# -----------------------------------------------------------------------------
# Dsim -> Dmel ortholog remapping (same recipe as integrate_v2.py /
# annotate_with_flysta3d.py). Cell line samples include D. simulans lines
# (Dsim-Merrill23, Dsim6B, Dsim6B-wMel -- quantified against the Dsim genome,
# so var_names are Dsim NCBI/Gnomon IDs), which must land in the same Dmel
# FBgn gene-ID space as the (already Dmel-indexed, post rule
# annotate_with_atlas) embryo reference files before Harmony/KNN.
# -----------------------------------------------------------------------------

def load_ortholog_map(path):
    """Load a Dsim -> Dmel FlyBase ID reciprocal-best-hit ortholog table.

    Expects tab-separated columns Dsim, Dmel (extra columns like pident/
    evalue/bitscore are ignored). Rows where a Dsim or Dmel ID appears more
    than once are dropped so the mapping stays strictly 1:1.

    Returns (dsim_to_dmel dict, dsim_ids set, dmel_ids set).
    """
    df = pd.read_csv(path, sep="\t")
    n_raw = len(df)
    df = df.drop_duplicates(subset="Dsim", keep=False)
    df = df.drop_duplicates(subset="Dmel", keep=False)
    n_kept = len(df)
    if n_kept < n_raw:
        print(f"  Ortholog map: dropped {n_raw - n_kept}/{n_raw} non-1:1 rows from {path}")
    dsim_to_dmel = dict(zip(df["Dsim"], df["Dmel"]))
    print(f"  Loaded {len(dsim_to_dmel)} 1:1 Dsim->Dmel orthologs from {path}")
    return dsim_to_dmel, set(df["Dsim"]), set(df["Dmel"])


def _looks_like_dsim(var_names, dsim_ids, dmel_ids):
    """Decide whether a sample's gene IDs belong to the Dsim or Dmel FlyBase
    ID namespace, by counting how many var_names land in each side of the
    ortholog table."""
    vs = set(var_names)
    n_dsim = len(vs & dsim_ids)
    n_dmel = len(vs & dmel_ids)
    return n_dsim > n_dmel, n_dsim, n_dmel


def remap_dsim_to_dmel(adata, dsim_to_dmel, label=""):
    """Rename a Dsim sample's var_names (Dsim FlyBase/NCBI IDs) to the
    orthologous Dmel FlyBase ID, dropping genes with no 1:1 ortholog."""
    mapped = adata.var_names.map(dsim_to_dmel)
    keep = mapped.notna().values
    n_total, n_kept = adata.n_vars, int(keep.sum())
    print(f"  [{label}] Dsim->Dmel remap: {n_kept}/{n_total} genes have a "
          f"1:1 Dmel ortholog (kept); {n_total - n_kept} dropped (no ortholog)")
    adata = adata[:, keep].copy()
    adata.var_names = mapped[keep].astype(str).values
    adata.var_names_make_unique()
    return adata


# -----------------------------------------------------------------------------
# Step 1 -- Load the annotated embryo reference files
# -----------------------------------------------------------------------------

def load_embryo_reference(ref_paths, label_cols):
    """Load all annotated embryo h5ads (rule annotate_with_atlas output),
    extract raw counts, and concatenate. Auto-detects atlas_<label> obs
    columns if label_cols isn't given."""
    adatas = []
    detected_cols = None

    for path in ref_paths:
        print(f"\n   Loading embryo reference: {path}")
        adata = sc.read_h5ad(path)
        print(f"   {adata.n_obs} cells x {adata.n_vars} genes")

        atlas_cols = [c for c in adata.obs.columns
                      if c.startswith("atlas_") and not c.endswith("_confidence")]
        if not atlas_cols:
            raise ValueError(
                f"{path} has no atlas_<label> obs columns -- expected rule "
                "annotate_with_atlas / annotate_with_flysta3d.py output "
                "(embryo samples labelled with the Flysta3D-v2 atlas)."
            )
        if detected_cols is None:
            detected_cols = set(atlas_cols)
        else:
            detected_cols &= set(atlas_cols)

        if adata.raw is None:
            raise ValueError(
                f"{path} has no .raw -- expected filter_h5ad output with "
                "adata.raw set to pre-normalisation counts."
            )
        raw_X = adata.raw.X
        if not scipy.sparse.issparse(raw_X):
            raw_X = scipy.sparse.csr_matrix(raw_X)
        raw_X = raw_X.tocsr()
        raw_X.data = raw_X.data.astype(np.float32)

        basename = os.path.splitext(os.path.basename(path))[0]
        keep_label_cols = label_cols or atlas_cols
        obs_cols = [c for c in keep_label_cols if c in adata.obs.columns]
        a = ad.AnnData(
            X=raw_X,
            obs=adata.obs[obs_cols].copy(),
            var=adata.raw.var.copy(),
        )
        a.obs["dataset"]     = "embryo_reference"
        a.obs["source_file"] = basename
        if "method" in adata.obs.columns:
            a.obs["method"] = adata.obs["method"].values
        if "replicate" in adata.obs.columns:
            a.obs["replicate"] = adata.obs["replicate"].values
        a.obs_names = [f"{basename}__{bc}" for bc in a.obs_names]

        adatas.append(a)

    if label_cols:
        missing = [c for c in label_cols if detected_cols is not None and c not in detected_cols]
        if missing:
            raise ValueError(
                f"--label_cols {missing} not found as atlas_<label> columns "
                f"on every reference file. Columns common to all reference "
                f"files: {sorted(detected_cols) if detected_cols else []}"
            )
        final_cols = label_cols
    else:
        final_cols = sorted(detected_cols) if detected_cols else []
        print(f"\n   No --label_cols given -- auto-detected atlas_<label> "
              f"columns common to all reference files: {final_cols}")
    if not final_cols:
        raise ValueError("No atlas_<label> columns found on any reference file.")

    print(f"\n-- Concatenating {len(adatas)} embryo reference files --")
    ref = ad.concat(adatas, join="outer", index_unique=None)
    ref.obs_names_make_unique()
    if scipy.sparse.issparse(ref.X):
        ref.X = ref.X.tocsr()
    print(f"   Total embryo reference cells: {ref.n_obs:,}")

    return ref, final_cols


# -----------------------------------------------------------------------------
# Step 2 -- Load cell line query files, remap Dsim, concatenate
# -----------------------------------------------------------------------------

def load_cellline_query(query_paths, dsim_to_dmel=None, dsim_ids=None, dmel_ids=None):
    adatas = []
    originals = {}  # basename -> original AnnData (kept in memory for re-merge)

    for path in query_paths:
        print(f"\n   Loading cell line query: {path}")
        adata = sc.read_h5ad(path)
        print(f"   {adata.n_obs} cells x {adata.n_vars} genes")

        if adata.raw is None:
            raise ValueError(
                f"{path} has no .raw -- expected filter_h5ad output with "
                "adata.raw set to pre-normalisation counts."
            )

        raw_X = adata.raw.X
        if not scipy.sparse.issparse(raw_X):
            raw_X = scipy.sparse.csr_matrix(raw_X)
        raw_X = raw_X.tocsr()
        raw_X.data = raw_X.data.astype(np.float32)

        basename = os.path.splitext(os.path.basename(path))[0]
        a = ad.AnnData(
            X=raw_X,
            obs=pd.DataFrame(index=adata.obs_names.copy()),
            var=adata.raw.var.copy(),
        )
        a.obs["dataset"]     = "cellline_query"
        a.obs["source_file"] = basename
        if "method" in adata.obs.columns:
            a.obs["method"] = adata.obs["method"].values
        if "replicate" in adata.obs.columns:
            a.obs["replicate"] = adata.obs["replicate"].values

        if dsim_to_dmel:
            is_dsim, n_dsim, n_dmel = _looks_like_dsim(a.var_names, dsim_ids, dmel_ids)
            if is_dsim:
                print(f"   {basename}: detected as Dsim ({n_dsim} Dsim IDs vs "
                      f"{n_dmel} Dmel IDs among var_names) -- remapping "
                      "var_names to Dmel orthologs before embryo comparison")
                a = remap_dsim_to_dmel(a, dsim_to_dmel, label=basename)

        a.obs_names = [f"{basename}__{bc}" for bc in a.obs_names]

        adatas.append(a)
        originals[basename] = adata

    print(f"\n-- Concatenating {len(adatas)} cell line query files --")
    query = ad.concat(adatas, join="outer", index_unique=None)
    query.obs_names_make_unique()
    if scipy.sparse.issparse(query.X):
        query.X = query.X.tocsr()
    print(f"   Total cell line query cells: {query.n_obs:,}")

    return query, originals


# -----------------------------------------------------------------------------
# Step 3 -- Joint preprocessing + Harmony
# -----------------------------------------------------------------------------

def joint_preprocess_and_harmony(ref, query, harmony_vars,
                                  n_pcs=30, n_top_genes=3000):
    print("\n-- Joint preprocessing --")

    combined = ad.concat(
        [ref, query],
        join="outer",
        index_unique=None,
        label="_source",
        keys=["reference", "query"],
    )
    combined.obs_names_make_unique()

    n_shared = len(set(ref.var_names) & set(query.var_names))
    print(f"   Shared genes (embryo ref ^ cell line query): {n_shared:,}")
    if n_shared < 500:
        raise ValueError(
            f"Only {n_shared} genes shared between embryo reference and cell "
            "line query after gene ID harmonisation -- check that both are "
            "Dmel-FBgn-indexed (Dsim samples remapped via --ortholog_map) "
            "before proceeding."
        )

    if scipy.sparse.issparse(combined.X):
        combined.X = combined.X.tocsr()
    else:
        combined.X = scipy.sparse.csr_matrix(combined.X)

    ref_mask = combined.obs["dataset"].values == "embryo_reference"
    print(f"   Combined: {combined.n_obs:,} cells x {combined.n_vars:,} genes")
    print(f"   Embryo reference: {ref_mask.sum():,}  Cell line query: {(~ref_mask).sum():,}")

    print("   Highly variable genes (seurat_v3 on raw counts, "
          "batch_key='dataset') ...")
    sc.pp.highly_variable_genes(
        combined, flavor="seurat_v3", n_top_genes=n_top_genes,
        batch_key="dataset", subset=False,
    )
    print(f"   HVGs selected: {combined.var['highly_variable'].sum()}")

    print("   Normalising (1e4 per cell) + log1p ...")
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)

    combined = combined[:, combined.var["highly_variable"]].copy()

    print("   Scaling (max_value=10) ...")
    sc.pp.scale(combined, max_value=10)

    print(f"   PCA ({n_pcs} components) ...")
    sc.tl.pca(combined, n_comps=n_pcs, svd_solver="arpack")

    missing_vars = [v for v in harmony_vars if v not in combined.obs.columns]
    if missing_vars:
        print(f"   NOTE: harmony_vars {missing_vars} not present in every "
              "cell -- filling with 'NA'")
    for v in harmony_vars:
        if v not in combined.obs.columns:
            combined.obs[v] = "NA"
        combined.obs[v] = combined.obs[v].astype(str).fillna("NA")

    print(f"   Harmony correction on: {harmony_vars} ...")
    ho = hm.run_harmony(
        combined.obsm["X_pca"], combined.obs, harmony_vars,
        max_iter_harmony=30, random_state=42,
    )
    combined.obsm["X_pca_harmony"] = ho.Z_corr.T

    print("   Harmony complete")
    return combined, ref_mask


# -----------------------------------------------------------------------------
# Step 4 -- KNN label transfer (multiple label columns from the same neighbours)
# -----------------------------------------------------------------------------

def knn_label_transfer(combined, ref_mask, ref_obs, label_cols, k=30):
    from sklearn.neighbors import NearestNeighbors

    print(f"\n-- KNN label transfer (k={k}) --")

    ref_pca   = combined.obsm["X_pca_harmony"][ref_mask]
    query_pca = combined.obsm["X_pca_harmony"][~ref_mask]
    query_names = combined.obs_names[~ref_mask]

    ref_names_in_combined = combined.obs_names[ref_mask]
    ref_obs_aligned = ref_obs.loc[ref_names_in_combined]

    nbrs = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=-1)
    nbrs.fit(ref_pca)
    distances, indices = nbrs.kneighbors(query_pca)

    result = pd.DataFrame(index=query_names)
    for col in label_cols:
        labels_arr = ref_obs_aligned[col].astype(str).values
        transferred, confidence = [], []
        for row_idx in indices:
            neigh = labels_arr[row_idx]
            counts = pd.Series(neigh).value_counts()
            transferred.append(counts.index[0])
            confidence.append(counts.iloc[0] / k)
        out_col = col.replace("atlas_", "embryo_", 1) if col.startswith("atlas_") else f"embryo_{col}"
        result[out_col]                = transferred
        result[f"{out_col}_confidence"] = confidence

        print(f"\n   Transferred '{out_col}' distribution:")
        print(pd.Series(transferred).value_counts().head(20).to_string())
        print(f"   Mean confidence: {np.mean(confidence):.3f}")
        low_conf = np.mean(np.array(confidence) < 0.5)
        if low_conf > 0.2:
            print(f"   WARNING: {low_conf*100:.1f}% of cells have confidence "
                  f"< 0.5 for '{out_col}' -- treat this column cautiously, or "
                  "increase --k, or check batch correction quality")

    return result


# -----------------------------------------------------------------------------
# Step 5 -- Merge labels back onto each ORIGINAL cell line query file
# -----------------------------------------------------------------------------

def write_mapped_outputs(originals, transferred_df, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    split = transferred_df.index.to_series().str.split("__", n=1, expand=True)
    transferred_df = transferred_df.copy()
    transferred_df["_basename"] = split[0].values
    transferred_df["_barcode"]  = split[1].values

    label_cols = [c for c in transferred_df.columns
                  if c not in ("_basename", "_barcode")]

    for basename, adata in originals.items():
        sub = transferred_df[transferred_df["_basename"] == basename]
        sub = sub.set_index("_barcode")
        sub = sub.reindex(adata.obs_names)

        n_missing = sub[label_cols[0]].isna().sum()
        if n_missing:
            print(f"   WARNING: {basename}: {n_missing}/{adata.n_obs} cells "
                  "have no transferred label -- left as NaN")

        for col in label_cols:
            adata.obs[col] = sub[col].values

        out_path = os.path.join(out_dir, f"{basename}.h5ad")
        adata.write(out_path)
        print(f"   -> {out_path}  ({adata.n_obs} cells, "
              f"{sub[label_cols[0]].notna().sum()} labelled)")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Map primary cell line samples onto already atlas-"
                     "annotated embryo cells via Harmony + KNN, transferring "
                     "embryo cell-type labels onto the cell lines."
    )
    parser.add_argument("--reference", required=True, nargs="+",
                         help="Annotated embryo h5ad file(s) -- rule "
                              "annotate_with_atlas / annotate_with_flysta3d.py "
                              "output, e.g. results/embryo_annotated/*.h5ad")
    parser.add_argument("--query", required=True, nargs="+",
                         help="Filtered per-sample h5ad file(s) for the "
                              "remaining (non-embryo / cell line) samples, "
                              "e.g. results/filtered_h5ad/{cell line samples}.h5ad")
    parser.add_argument("--out_dir", required=True,
                         help="Directory to write mapped copies of each cell "
                              "line query file (same filenames, new "
                              "embryo_* obs columns added)")
    parser.add_argument("--label_cols", nargs="+", default=None,
                         help="atlas_<label> column(s) (as written by rule "
                              "annotate_with_atlas) to transfer, given without "
                              "the 'atlas_' prefix is NOT required -- pass the "
                              "exact atlas_<label> column name(s), e.g. "
                              "atlas_cell_type atlas_tissue. Omit to "
                              "auto-detect every atlas_<label> column common "
                              "to all reference files.")
    parser.add_argument("--ortholog_map", type=str, default=None,
                         help="TSV with Dsim/Dmel FlyBase ID columns "
                              "(reciprocal best hit orthologs). If given, "
                              "cell line query files detected as Dsim (Dsim-"
                              "Merrill23, Dsim6B, Dsim6B-wMel, ...) have their "
                              "var_names remapped to the orthologous Dmel "
                              "FlyBase ID before comparison against the "
                              "(Dmel-indexed) embryo reference.")
    parser.add_argument("--k", type=int, default=30,
                         help="Number of nearest embryo-reference neighbours "
                              "for the majority-vote label transfer")
    parser.add_argument("--n_pcs", type=int, default=30)
    parser.add_argument("--n_top_genes", type=int, default=3000,
                         help="Jointly-computed HVGs (seurat_v3) used for "
                              "PCA/Harmony")
    parser.add_argument("--harmony_vars", nargs="+",
                         default=["dataset", "method"],
                         help="obs columns Harmony corrects for. 'dataset' "
                              "(embryo_reference vs. cellline_query) should "
                              "always be included.")

    args = parser.parse_args()

    ref_paths = []
    for pattern in args.reference:
        matches = glob.glob(pattern)
        ref_paths.extend(matches if matches else [pattern])
    ref_paths = sorted(set(ref_paths))
    print(f"Embryo reference files ({len(ref_paths)}):")
    for p in ref_paths:
        print(f"  {p}")

    query_paths = []
    for pattern in args.query:
        matches = glob.glob(pattern)
        query_paths.extend(matches if matches else [pattern])
    query_paths = sorted(set(query_paths))
    print(f"\nCell line query files ({len(query_paths)}):")
    for p in query_paths:
        print(f"  {p}")

    ref, label_cols = load_embryo_reference(ref_paths, args.label_cols)
    ref_obs = ref.obs.copy()

    dsim_to_dmel, dsim_ids, dmel_ids = (
        load_ortholog_map(args.ortholog_map) if args.ortholog_map else ({}, set(), set())
    )
    query, originals = load_cellline_query(
        query_paths, dsim_to_dmel=dsim_to_dmel, dsim_ids=dsim_ids, dmel_ids=dmel_ids,
    )

    combined, ref_mask = joint_preprocess_and_harmony(
        ref, query, harmony_vars=args.harmony_vars,
        n_pcs=args.n_pcs, n_top_genes=args.n_top_genes,
    )

    transferred_df = knn_label_transfer(
        combined, ref_mask, ref_obs, label_cols, k=args.k,
    )

    write_mapped_outputs(originals, transferred_df, args.out_dir)

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Mapped cell line files -> {args.out_dir}/")
    print("Point rule integrate at these files (alongside the annotated "
          "embryo files) so every sample carries a cell-type label -- "
          "atlas_<label> for embryos, embryo_<label> for cell lines.")


if __name__ == "__main__":
    main()
