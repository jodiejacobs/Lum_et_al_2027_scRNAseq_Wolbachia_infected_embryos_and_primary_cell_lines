"""
annotate_with_flysta3d.py
==========================
Transfers Drosophila embryo cell-type labels from the Flysta3D-v2 atlas
(Wang et al. 2025, Cell 188:4734-4753 -- https://db.cngb.org/stomics/flysta3d-v2/)
onto your own filtered per-sample h5ad files, BEFORE they go into your
existing multi-condition integrate.py run.

Why annotate per-sample files rather than the already-integrated object
-------------------------------------------------------------------------
Cell identity is a property of the cell, not of your experiment. Doing the
transfer here means:
  - atlas_cell_type rides through integrate.py's own concat/Harmony step as
    a fixed obs column (same mechanism wolbachia_titer already uses),
    instead of needing a second, more complicated Harmony correction that
    would have to reconcile three batch structures at once (atlas source,
    your capture method, your infection condition/timepoint).
  - integrate.py's Harmony step keeps correcting for method + replicate
    only, exactly as documented in its own docstring -- cell type doesn't
    leak into that correction.
  - You get one annotation per cell, computed once, that infected and
    uninfected cells in every downstream analysis can all reference.

Design (mirrors integrate_by_ref.py's reference+query Harmony/KNN recipe)
-------------------------------------------------------------------------
1. Load the Flysta3D-v2 atlas once. Detect whether its var_names are
   FlyBase FBgn IDs or gene symbols (Flysta3D is generally symbol-indexed);
   if symbols, remap to FBgn using your existing fbgn_annotation_ID table
   so genes line up with your own FBgn-indexed kallisto|bustools output.
2. Load ALL your query h5ad files' raw counts in one pass and concatenate
   them (tagged by source_file) -- one joint Harmony fit sees the full
   diversity of your batches at once, rather than fitting separately per
   file (which would give inconsistent KNN neighborhoods file-to-file).
3. Concatenate reference + query on shared genes, normalize + log1p,
   restrict to jointly-computed HVGs, scale, PCA, Harmony (batch = atlas
   vs. query, plus method/replicate if present).
4. KNN in Harmony PCA space: each query cell gets the majority-vote label
   (+ confidence = winning fraction of k neighbours) from its k nearest
   ATLAS cells, for every --label_cols column you ask for.
5. Split the labelled query cells back out by source_file and write each
   ORIGINAL h5ad back out with the new atlas_<label> / atlas_<label>_conf
   obs columns added -- nothing else about your file changes.

Memory note: the atlas file is large (the "downsampled" whole-embryo
co-embed is ~12 GB). This script stays sparse everywhere except the final
HVG-subset scale/PCA step (same as your own preprocess() in integrate.py).
Use --subsample_ref for a fast first pass before committing a big SLURM job.

Run with:
    mamba activate scanpy
    python snakemake_scripts/method_comparison/annotate_with_flysta3d.py \\
        --atlas resources/wcoembed_whole_embeding_downsampled_modified.h5ad \\
        --query results/filtered_h5ad/*embryos*.h5ad \\
        --out_dir results/embryo_annotated \\
        --label_cols cell_type tissue germ_layer \\
        --flybase_annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz \\
        --ortholog_map reference/orthologs/dmel_dsim_orthologs_rbh.tsv

Only pass your EMBRYO query files here (this atlas is embryo-specific --
transferring embryo cell-type labels onto cultured primary cell lines
directly doesn't make biological sense). Primary cell line samples get
their cell-type label via map_cellline_to_embryo.py instead, which maps
them onto these already-annotated embryo cells.

First run without --label_cols to see what's actually in the atlas obs --
the script will print every candidate column and exit before doing any
expensive work.
"""

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scipy.sparse
import anndata as ad
import scanpy as sc
import harmonypy as hm


FBGN_RE = re.compile(r"^FBgn\d{7,8}$")


# -----------------------------------------------------------------------------
# Dsim -> Dmel ortholog remapping (same recipe as integrate_v2.py). The
# Flysta3D-v2 atlas is Dmel-indexed (FBgn, after harmonise_atlas_gene_ids
# below), but several of our own query samples are D. simulans (quantified
# against the Dsim genome, so their var_names are Dsim NCBI/Gnomon IDs, e.g.
# "LOC120284240" -- a different namespace entirely). Left unmapped, every
# Dsim query cell would look like it has ~0 genes in common with the atlas
# once restricted to shared var_names, and the label transfer below would be
# meaningless for those cells. Remap each Dsim query file's var_names onto
# the orthologous Dmel FlyBase ID (reciprocal-best-hit table) before
# concatenation, exactly like integrate_v2.py does for the final integration
# step -- this keeps every query file, Dmel or Dsim, embryo or cell line, in
# the same Dmel FBgn gene-ID space the atlas already uses.
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
    # Index.map() returns an Index; Index.notna() already returns a plain
    # ndarray (unlike Series.notna()), so no trailing .values here.
    keep = mapped.notna()
    n_total, n_kept = adata.n_vars, int(keep.sum())
    print(f"  [{label}] Dsim->Dmel remap: {n_kept}/{n_total} genes have a "
          f"1:1 Dmel ortholog (kept); {n_total - n_kept} dropped (no ortholog)")
    adata = adata[:, keep].copy()
    adata.var_names = mapped[keep].astype(str).values
    adata.var_names_make_unique()
    return adata


# -----------------------------------------------------------------------------
# Gene ID harmonisation: Flysta3D is generally symbol-indexed; your own
# kallisto|bustools output is FBgn-indexed (see swap_gene_id_to_symbol.py --
# that script is only run for the *optional* symbol-keyed index, the default
# t2g.txt/index.idx used by map_pipseq/map_10x keys gene_id = FBgn).
# -----------------------------------------------------------------------------

def load_symbol_to_fbgn(flybase_path):
    """Build a gene_symbol -> primary FBgn# map from a FlyBase
    fbgn_annotation_ID_fb_*.tsv.gz file (same file already used by
    nmf_programs.py's load_gene_mapping, just inverted: column 0 is
    gene_symbol, column 2 is primary_FBgn#).

    Ambiguous symbols (mapping to >1 distinct FBgn) are dropped rather than
    picking one arbitrarily -- mirrors the drop_duplicates(keep=False)
    pattern already used for the Dsim/Dmel ortholog table in integrate.py.
    """
    import gzip
    from io import StringIO

    with gzip.open(flybase_path, "rt") as f:
        lines = [line for line in f if not line.startswith("#")]
    df = pd.read_csv(StringIO("".join(lines)), sep="\t", header=None)
    df = df[[0, 2]].rename(columns={0: "symbol", 2: "fbgn"}).dropna()

    n_raw = len(df)
    df = df.drop_duplicates(subset="symbol", keep=False)
    n_kept = len(df)
    if n_kept < n_raw:
        print(f"  FlyBase symbol->FBgn map: dropped {n_raw - n_kept}/{n_raw} "
              f"ambiguous (non-unique) symbol rows from {flybase_path}")

    mapping = dict(zip(df["symbol"], df["fbgn"]))
    print(f"  Loaded {len(mapping)} unique symbol->FBgn mappings")
    return mapping


def harmonise_atlas_gene_ids(atlas, flybase_annotation=None):
    """Make atlas.var_names FBgn IDs so they line up with your query files.

    Priority:
      1. var_names already look like FBgn IDs -> leave as-is.
      2. atlas.var has a column that looks like it holds FBgn/FlyBase IDs
         (name contains 'fbgn'/'flybase'/'gene_id', values match FBgn regex)
         -> use that column directly.
      3. Fall back to symbol->FBgn remap via --flybase_annotation.
    """
    frac_fbgn = atlas.var_names.str.match(FBGN_RE).mean()
    print(f"  Atlas var_names matching FBgn pattern: {frac_fbgn*100:.1f}%")
    print(f"  Atlas var_names sample (first 20): {atlas.var_names[:20].tolist()}")
    print(f"  Atlas var.columns: {atlas.var.columns.tolist()}")
    if frac_fbgn > 0.5:
        print("  Atlas already FBgn-indexed -- no remap needed")
        return atlas

    for col in atlas.var.columns:
        if re.search(r"fbgn|flybase|gene_?id", col, re.IGNORECASE):
            vals = atlas.var[col].astype(str)
            frac = vals.str.match(FBGN_RE).mean()
            if frac > 0.5:
                print(f"  Using atlas.var['{col}'] as FBgn ID "
                      f"({frac*100:.1f}% match FBgn pattern)")
                atlas = atlas.copy()
                atlas.var_names = vals.values
                atlas.var_names_make_unique()
                return atlas

    if flybase_annotation is None:
        raise ValueError(
            "Atlas var_names don't look like FBgn IDs and no --flybase_annotation "
            "was given to remap them. Either pass --flybase_annotation "
            "reference/fbgn_annotation_ID_fb_2025_04.tsv.gz, or inspect "
            "atlas.var.columns yourself and tell the script which column holds "
            "FlyBase IDs."
        )

    print("  Atlas appears symbol-indexed -- remapping to FBgn via "
          f"{flybase_annotation}")
    symbol_to_fbgn = load_symbol_to_fbgn(flybase_annotation)
    mapped = atlas.var_names.map(symbol_to_fbgn)
    # Index.map() returns an Index; Index.notna() already returns a plain
    # ndarray (unlike Series.notna()), so no trailing .values here.
    keep = mapped.notna()
    n_total, n_kept = atlas.n_vars, int(keep.sum())
    print(f"  Symbol->FBgn remap: {n_kept}/{n_total} atlas genes matched "
          f"({n_total - n_kept} dropped -- no unique FBgn for that symbol)")
    if n_kept < 1000:
        sample_symbols = list(symbol_to_fbgn.keys())[:20]
        print(f"  Atlas var_names sample (unmatched): {atlas.var_names[:20].tolist()}")
        print(f"  FlyBase table symbol sample (what we're matching against): "
              f"{sample_symbols}")
        raise ValueError(
            f"Only {n_kept} atlas genes remapped to FBgn -- something is "
            "wrong with the symbol matching (check for a species-prefix "
            "like 'Dmel\\\\' or case mismatches in atlas.var_names, or "
            "whether atlas.var_names are actually symbols at all -- compare "
            "the two samples printed above)."
        )
    atlas = atlas[:, keep].copy()
    atlas.var_names = mapped[keep].astype(str).values
    atlas.var_names_make_unique()
    return atlas


def _recover_var_index(var, label="var"):
    """Some h5ad exports (e.g. R/Seurat -> h5ad round-trips, or a "_modified"
    re-export like this atlas) end up with a purely positional integer var
    index (0, 1, 2, ...) while the real gene identifiers survive only as an
    ordinary column -- commonly literally named "index", a naming collision
    with anndata's own special "_index" marker used on disk. Left alone,
    var_names is just row numbers and every downstream gene-ID match (FBgn
    pattern, symbol->FBgn remap) silently returns ~0 hits instead of
    erroring loudly.

    Detects a numeric-looking index and, if found, promotes whichever
    column looks like real (non-numeric) gene IDs back to the index.
    """
    idx_is_numeric = pd.Series(var.index.astype(str)).str.match(r"^\d+$").mean() > 0.9
    if not idx_is_numeric:
        return var

    for col in var.columns:
        vals = var[col].astype(str)
        if vals.str.match(r"^\d+$").mean() < 0.5:
            print(f"   {label}: index is positional (0,1,2,...) -- "
                  f"promoting column '{col}' (e.g. {vals.iloc[0]!r}) to be "
                  "the real var_names")
            # drop=True: a leftover column literally named "index" (as seen
            # in this atlas export) would otherwise collide with anndata's
            # own on-disk "_index" convention when the object is written
            # back out later.
            var = var.set_index(col, drop=True)
            var.index.name = None
            return var

    print(f"   WARNING: {label} index looks positional (0,1,2,...) and no "
          f"column looks like real gene IDs either (columns: "
          f"{list(var.columns)}) -- leaving var_names as row numbers, gene "
          "ID matching below will fail")
    return var


# -----------------------------------------------------------------------------
# Step 1 -- Load atlas reference
# -----------------------------------------------------------------------------

def load_atlas_reference(atlas_path, label_cols, flybase_annotation=None,
                          subsample_ref=None, random_state=42):
    print(f"\n-- Loading Flysta3D-v2 atlas: {atlas_path} --")
    atlas = sc.read_h5ad(atlas_path)
    print(f"   {atlas.n_obs:,} cells x {atlas.n_vars:,} genes")

    if not label_cols:
        candidates = [c for c in atlas.obs.columns
                      if atlas.obs[c].dtype == object
                      or str(atlas.obs[c].dtype).startswith("category")]
        print("\n   No --label_cols given. Candidate annotation columns in "
              "the atlas obs (dtype object/category):")
        for c in candidates:
            n_uniq = atlas.obs[c].nunique()
            print(f"     {c!r}  ({n_uniq} unique values)")
        print("\n   All obs columns:")
        print(f"     {list(atlas.obs.columns)}")
        raise SystemExit(
            "\nPick one or more of the columns above and re-run with "
            "--label_cols <col1> [<col2> ...]."
        )

    missing = [c for c in label_cols if c not in atlas.obs.columns]
    if missing:
        raise ValueError(
            f"--label_cols {missing} not found in atlas obs. "
            f"Available: {list(atlas.obs.columns)}"
        )

    for c in label_cols:
        print(f"\n   Atlas '{c}' distribution (top 15):")
        print(atlas.obs[c].astype(str).value_counts().head(15).to_string())

    if subsample_ref and atlas.n_obs > subsample_ref:
        print(f"\n   Subsampling atlas to {subsample_ref:,} cells "
              f"(--subsample_ref) for a faster first pass")
        sc.pp.subsample(atlas, n_obs=subsample_ref, random_state=random_state)

    # Prefer raw counts for renormalisation consistency with the query data.
    if atlas.raw is not None:
        X = atlas.raw.X
        var = atlas.raw.var.copy()
        print("   Using atlas.raw.X as counts source")
    elif "counts" in atlas.layers:
        X = atlas.layers["counts"]
        var = atlas.var.copy()
        print("   Using atlas.layers['counts'] as counts source")
    else:
        X = atlas.X
        var = atlas.var.copy()
        print("   WARNING: no .raw or 'counts' layer found in the atlas -- "
              "using .X as-is. If this atlas object is already "
              "normalised/log1p'd, renormalising it here is not strictly "
              "correct (values won't be integer counts), but the relative "
              "structure used for PCA/Harmony/KNN is only mildly affected. "
              "Check atlas.X.max() / whether values look like log1p counts "
              "if this matters for your analysis.")

    # This atlas export (and apparently others like it -- see
    # _recover_var_index docstring) can carry a purely positional var index
    # with the real gene IDs stranded in an ordinary column instead; fix
    # that up before var_names is used for anything downstream.
    var = _recover_var_index(var, label="atlas var")

    if scipy.sparse.issparse(X):
        X = X.tocsr()
    else:
        X = scipy.sparse.csr_matrix(X)
    X.data = X.data.astype(np.float32)

    # The Flysta3D-v2 export co-embeds scRNA-seq and scATAC-seq cells in the
    # same object (see the `assay` obs column). The ATAC-only barcodes carry
    # nCount_RNA == 0 / nFeature_RNA == 0 -- i.e. an all-zero row in the
    # counts matrix we just extracted. Left in, an all-zero row survives
    # normalize_total (0/0 -> 0, harmless) but after sc.pp.scale() every gene
    # in that row becomes a strong, *uniform* z-score of -mean/std rather
    # than "no information" -- which pulls all such cells together into
    # tight, artificial clusters sitting off to the side of the real
    # biological manifold. That matches the "small mixed clusters
    # surrounding the actual atlas" you're seeing, and they're composed of
    # atlas cells because they *are* atlas cells -- just ones with no RNA
    # signal to cluster on. Drop them here before anything downstream sees
    # them.
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    keep_cell = total_counts > 0
    n_dropped = int((~keep_cell).sum())
    if n_dropped:
        print(f"   Dropping {n_dropped:,} / {X.shape[0]:,} atlas cells with "
              f"zero total counts (likely co-embedded ATAC-only barcodes, "
              f"not RNA transcriptomes)")
        if "assay" in atlas.obs.columns:
            dropped_assay_counts = (
                atlas.obs.loc[~keep_cell, "assay"].value_counts()
            )
            kept_assay_counts = (
                atlas.obs.loc[keep_cell, "assay"].value_counts()
            )
            print("   Dropped cells by `assay`:")
            for val, n in dropped_assay_counts.items():
                print(f"     {val}: {n:,}")
            print("   Kept cells by `assay`:")
            for val, n in kept_assay_counts.items():
                print(f"     {val}: {n:,}")
        X = X[keep_cell]
        atlas = atlas[keep_cell].copy()
    else:
        print("   No zero-count atlas cells found (nothing dropped)")

    ref = ad.AnnData(X=X, obs=atlas.obs[label_cols].copy(), var=var)
    ref.obs_names = atlas.obs_names
    ref.obs["dataset"] = "flysta3d_atlas"
    ref = harmonise_atlas_gene_ids(ref, flybase_annotation)

    print(f"\n   Reference ready: {ref.n_obs:,} cells x {ref.n_vars:,} genes")
    return ref


# -----------------------------------------------------------------------------
# Step 2 -- Load all query files' raw counts, tag by source_file
# -----------------------------------------------------------------------------

def load_query_files(query_paths, dsim_to_dmel=None, dsim_ids=None, dmel_ids=None):
    adatas = []
    originals = {}  # basename -> original AnnData (kept in memory for re-merge)

    for path in query_paths:
        print(f"\n   Loading query: {path}")
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
        a.obs["dataset"]     = "query"
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
                      "var_names to Dmel orthologs before atlas comparison")
                a = remap_dsim_to_dmel(a, dsim_to_dmel, label=basename)

        a.obs_names = [f"{basename}__{bc}" for bc in a.obs_names]

        adatas.append(a)
        originals[basename] = adata

    print(f"\n-- Concatenating {len(adatas)} query files --")
    query = ad.concat(adatas, join="outer", index_unique=None)
    query.obs_names_make_unique()
    if scipy.sparse.issparse(query.X):
        query.X = query.X.tocsr()
    print(f"   Total query cells: {query.n_obs:,}")

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
    print(f"   Shared genes (ref ^ query): {n_shared:,}")
    if n_shared < 500:
        raise ValueError(
            f"Only {n_shared} genes shared between atlas and query after gene "
            "ID harmonisation -- check that both are FBgn-indexed for the "
            "same species/genome build before proceeding."
        )

    if scipy.sparse.issparse(combined.X):
        combined.X = combined.X.tocsr()
    else:
        combined.X = scipy.sparse.csr_matrix(combined.X)

    ref_mask = combined.obs["dataset"].values == "flysta3d_atlas"
    print(f"   Combined: {combined.n_obs:,} cells x {combined.n_vars:,} genes")
    print(f"   Reference: {ref_mask.sum():,}  Query: {(~ref_mask).sum():,}")

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
              "cell (e.g. atlas has no 'method') -- filling with 'NA'")
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

    # ref_obs is indexed by the *original* atlas obs_names; combined's
    # reference rows are in the same order atlas cells were concatenated in,
    # so re-derive per-row labels via the reference obs_names actually kept
    # in `combined`.
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
        result[f"atlas_{col}"]           = transferred
        result[f"atlas_{col}_confidence"] = confidence

        print(f"\n   Transferred '{col}' distribution:")
        print(pd.Series(transferred).value_counts().head(20).to_string())
        print(f"   Mean confidence: {np.mean(confidence):.3f}")
        low_conf = np.mean(np.array(confidence) < 0.5)
        if low_conf > 0.2:
            print(f"   WARNING: {low_conf*100:.1f}% of cells have confidence "
                  f"< 0.5 for '{col}' -- treat this column cautiously, or "
                  "increase --k, or check batch correction quality")

    return result


# -----------------------------------------------------------------------------
# Diagnostics -- QC plots showing the label transfer is doing something sane
# -----------------------------------------------------------------------------

def _savefig(fig, path):
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"   Saved: {path}")


def plot_diagnostics(combined, ref_mask, ref_obs, transferred_df, label_cols,
                      fig_dir, k):
    """QC plots for the atlas -> query label transfer:
      - UMAP colored by dataset (atlas vs. query) -- the sanity check is
        that query cells scatter through the same regions as atlas cells
        rather than clumping off in their own separate island, which would
        mean Harmony didn't actually integrate the two batches.
      - UMAP colored by each label_col (atlas ground truth on reference
        cells, KNN-transferred value on query cells, same color scale) --
        query cells should land in/near the same-colored atlas regions if
        the transfer is picking up real biology.
      - Confidence score histogram + per-source-file confidence boxplot for
        each label_col -- flags samples or whole label columns where the
        KNN vote was weak (transfer is a guess, not a confident call).
    """
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    print(f"\n-- Diagnostic plots (k={k}) -- writing to {fig_dir}/ --")
    print("   Computing neighbors + UMAP on combined Harmony embedding ...")
    sc.pp.neighbors(combined, use_rep="X_pca_harmony", n_neighbors=30)
    sc.tl.umap(combined)

    combined.obs["_dataset_display"] = np.where(ref_mask, "flysta3d_atlas", "query")
    sc.pl.umap(combined, color="_dataset_display", save="_dataset.pdf",
               title="Atlas vs. query cells (post-Harmony)")

    ref_names_in_combined   = combined.obs_names[ref_mask]
    query_names_in_combined = combined.obs_names[~ref_mask]
    ref_obs_aligned = ref_obs.loc[ref_names_in_combined]

    for col in label_cols:
        atlas_col = f"atlas_{col}"
        conf_col  = f"{atlas_col}_confidence"

        display_col = f"_display_{col}"
        combined.obs[display_col] = pd.Series(index=combined.obs_names, dtype=object)
        combined.obs.loc[ref_names_in_combined, display_col] = \
            ref_obs_aligned[col].astype(str).values
        combined.obs.loc[query_names_in_combined, display_col] = \
            transferred_df.loc[query_names_in_combined, atlas_col].astype(str).values
        combined.obs[display_col] = combined.obs[display_col].fillna("NA")

        sc.pl.umap(combined, color=display_col, save=f"_{col}.pdf",
                   title=f"'{col}': atlas ground truth + KNN-transferred query")

        conf_vals = transferred_df.loc[query_names_in_combined, conf_col].astype(float)

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(conf_vals.dropna(), bins=30, color="#2196F3", edgecolor="black", alpha=0.8)
        ax.axvline(conf_vals.mean(), color="red", linestyle="--",
                   label=f"mean={conf_vals.mean():.2f}")
        ax.set_xlabel(f"KNN confidence ({col})")
        ax.set_ylabel("Query cells")
        ax.set_title(f"Label transfer confidence -- {col} (k={k})")
        ax.legend()
        _savefig(fig, os.path.join(fig_dir, f"confidence_hist_{col}.pdf"))

        conf_df = pd.DataFrame({
            "confidence":  conf_vals.values,
            "source_file": combined.obs.loc[query_names_in_combined, "source_file"].values,
        })
        samples = sorted(conf_df["source_file"].unique())
        fig, ax = plt.subplots(figsize=(max(8, len(samples) * 1.2), 5))
        ax.boxplot([conf_df.loc[conf_df["source_file"] == s, "confidence"].dropna().values
                    for s in samples],
                   labels=samples, showfliers=False)
        ax.set_ylabel(f"KNN confidence ({col})")
        ax.set_title(f"Label transfer confidence by sample -- {col}")
        ax.set_ylim(0, 1.05)
        plt.xticks(rotation=45, ha="right")
        _savefig(fig, os.path.join(fig_dir, f"confidence_by_sample_{col}.pdf"))

        conf_df.groupby("source_file")["confidence"].agg(
            ["mean", "median", "std", "count"]
        ).to_csv(os.path.join(fig_dir, f"confidence_summary_{col}.csv"))

    print(f"   Diagnostic plots complete -- see {fig_dir}/")


# -----------------------------------------------------------------------------
# Step 5 -- Merge labels back onto each ORIGINAL query file and write out
# -----------------------------------------------------------------------------

def write_annotated_outputs(originals, transferred_df, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # transferred_df index is "{basename}__{barcode}" (unique_only) or
    # "{basename}__{barcode}-N" if obs_names_make_unique appended a suffix.
    # Recover (basename, barcode) by splitting on the first "__".
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
                  "have no transferred label (dropped somewhere upstream, "
                  "e.g. all-zero after gene intersection) -- left as NaN")

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
        description="Transfer Flysta3D-v2 atlas cell-type labels onto your "
                     "filtered per-sample h5ad files via Harmony + KNN, before "
                     "your own integrate.py run."
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
                              "cell_type tissue germ_layer). Omit to print "
                              "candidate columns from the atlas and exit.")
    parser.add_argument("--flybase_annotation", default=None,
                         help="reference/fbgn_annotation_ID_fb_2025_04.tsv.gz "
                              "-- used to remap atlas gene symbols to FBgn IDs "
                              "if the atlas isn't already FBgn-indexed.")
    parser.add_argument("--ortholog_map", type=str, default=None,
                         help="TSV with Dsim/Dmel FlyBase ID columns (reciprocal "
                              "best hit orthologs; same file passed to "
                              "integrate_v2.py's --ortholog_map). If given, "
                              "query files detected as Dsim (by gene ID overlap "
                              "with this table) have their var_names remapped "
                              "to the orthologous Dmel FlyBase ID before "
                              "comparison against the Dmel-indexed atlas. If "
                              "omitted, no remapping is done and Dsim query "
                              "files will share almost no genes with the atlas.")
    parser.add_argument("--k", type=int, default=30,
                         help="Number of nearest atlas neighbours for the "
                              "majority-vote label transfer (default 30; the "
                              "atlas has many more diverse cell types than "
                              "your own uninfected reference, so a larger k "
                              "than integrate_by_ref.py's default of 15 is "
                              "used here)")
    parser.add_argument("--n_pcs", type=int, default=30)
    parser.add_argument("--n_top_genes", type=int, default=3000,
                         help="Jointly-computed HVGs (seurat_v3) used for "
                              "PCA/Harmony")
    parser.add_argument("--harmony_vars", nargs="+",
                         default=["dataset", "method"],
                         help="obs columns Harmony corrects for. 'dataset' "
                              "(atlas vs. query) should always be included.")
    parser.add_argument("--subsample_ref", type=int, default=None,
                         help="Subsample the atlas to this many cells before "
                              "integration -- use for a fast first pass.")
    parser.add_argument("--fig_dir", type=str, default=None,
                         help="If given, write QC plots here: UMAP colored "
                              "by atlas-vs-query, UMAP colored by each "
                              "transferred label, and confidence-score "
                              "histograms/boxplots per label per sample. "
                              "Adds a neighbors+UMAP computation on top of "
                              "the Harmony step. Omit to skip plotting.")

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
    ref_obs = ref.obs.copy()

    dsim_to_dmel, dsim_ids, dmel_ids = (
        load_ortholog_map(args.ortholog_map) if args.ortholog_map else ({}, set(), set())
    )
    query, originals = load_query_files(
        query_paths, dsim_to_dmel=dsim_to_dmel, dsim_ids=dsim_ids, dmel_ids=dmel_ids,
    )

    combined, ref_mask = joint_preprocess_and_harmony(
        ref, query, harmony_vars=args.harmony_vars,
        n_pcs=args.n_pcs, n_top_genes=args.n_top_genes,
    )

    transferred_df = knn_label_transfer(
        combined, ref_mask, ref_obs, args.label_cols, k=args.k,
    )

    if args.fig_dir:
        plot_diagnostics(
            combined, ref_mask, ref_obs, transferred_df, args.label_cols,
            args.fig_dir, k=args.k,
        )

    write_annotated_outputs(originals, transferred_df, args.out_dir)

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Annotated files -> {args.out_dir}/")
    print("Point rule integrate / integrate_uninfected at this directory "
          "instead of results/filtered_h5ad/ to carry atlas_* columns "
          "through your existing Harmony/BBKNN integration.")


if __name__ == "__main__":
    main()
