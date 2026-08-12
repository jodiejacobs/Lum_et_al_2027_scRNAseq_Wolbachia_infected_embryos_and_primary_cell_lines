# resources/

External reference data too large to keep in git. Download to your cluster's
scratch/data disk and either place it here (as a symlink is fine) or point
`flysta3d_atlas` in `config/config.yaml` at wherever it actually lives.

## Flysta3D-v2 atlas

Wang et al. 2025, *Cell* 188(17):4734-4753 -- *A Drosophila single-cell 3D
spatiotemporal multi-omics atlas unveils panoramic key regulators of
cell-type differentiation*. doi:10.1016/j.cell.2025.05.047

Used by `snakemake_scripts/method_comparison/annotate_with_flysta3d.py` /
`rule annotate_with_atlas` to transfer embryo cell-type labels onto our
filtered per-sample h5ad files before integration.

Download the whole-embryo co-embedded object (the file the pipeline expects
by default):

```
wget https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000398/wcoembed_whole_embeding_downsampled_modified.h5ad \
    -O resources/wcoembed_whole_embeding_downsampled_modified.h5ad
```

~12 GB. Database: https://db.cngb.org/stomics/flysta3d-v2/
CNSA project accession: CNP0005060

Other files on the same download page (per-timepoint scStereo-seq/scRNA-seq
objects, germ-layer-split co-embeds, ATAC objects) are listed at
https://db.cngb.org/stomics/flysta3d-v2/download if a different subset is
more useful later.
