"""Shared helpers for the snakemake_scripts/pseudotime/ scripts."""

import os

import pandas as pd
import matplotlib.pyplot as plt


def savefig(fig, path):
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


def load_flybase_symbols(flybase_annotation):
    """FBgn -> gene symbol from FlyBase fbgn_annotation_ID*.tsv.gz."""
    if not flybase_annotation or not os.path.exists(flybase_annotation):
        return {}
    df = pd.read_csv(flybase_annotation, sep="\t", comment="#", header=None,
                     usecols=[0, 2], dtype=str).dropna()
    return dict(zip(df[2], df[0]))


def load_orthologs(path):
    """Strict 1:1 Dsim -> Dmel map from the RBH table (columns Dsim, Dmel)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df.drop_duplicates("Dsim", keep=False).drop_duplicates("Dmel", keep=False)
    print(f"  {len(df)} 1:1 Dsim->Dmel orthologs from {path}")
    return dict(zip(df["Dsim"], df["Dmel"]))


def load_gene_labels(flybase_annotation, ortholog_map, species):
    """gene id -> readable label (Dmel symbol; Dsim IDs via their 1:1 Dmel ortholog)."""
    fb = load_flybase_symbols(flybase_annotation)
    if species == "Dmel" or not ortholog_map:
        return fb
    # works whether Dsim genes are native IDs or already remapped to FBgn
    return {**fb, **{d: fb.get(m, m) for d, m in load_orthologs(ortholog_map).items()}}


def remap_to_dmel(adata, to_dmel, label=""):
    """Rename var_names to their 1:1 Dmel FBgn ortholog and drop genes without
    one. Same logic as remap_dsim_to_dmel() in
    method_comparison/annotate_with_flysta3d.py (used by rule integrate), so
    Dsim genes are annotated the same way here as in integrated.h5ad."""
    mapped = adata.var_names.map(to_dmel)
    keep = mapped.notna()
    print(f"  [{label}] ortholog remap: {int(keep.sum())}/{adata.n_vars} genes have a "
          f"1:1 Dmel ortholog (kept)")
    adata = adata[:, keep].copy()
    adata.var_names = mapped[keep].astype(str).values
    return adata
