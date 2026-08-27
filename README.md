# Methods for the single-cell analysis of differential expression in *Wolbachia*-infected *Drosophila* embryos and primary cell lines

In prep.

Snakemake pipeline for scRNA-seq (PIPseq and 10x Genomics) processing, QC, host+*Wolbachia* joint quantification, cell-type annotation, and cross-condition analysis of *Wolbachia*-infected and uninfected *Drosophila* embryos and cultured primary cell lines. Supports three host species (*D. melanogaster*, *D. simulans*, *D. willistoni*) and four *Wolbachia* strains (wMel, wRi_Riv84, wRi_M23, wWil).

## Requirements

- `snakemake` >= 9.0, run from a `mamba`/`conda` environment (`mamba activate snakemake`)
- A Slurm cluster (the Snakefile's per-rule `resources:` blocks assume `--executor slurm`)
- Four additional conda/mamba environments, built from `config/envs/*.yml` and pointed to by path in `config/config.yaml`:
  - `scanpy_env` (`config/envs/scanpy_env.yml`) — scanpy/anndata QC, filtering, integration, and downstream analysis scripts
  - `kallisto_env` (`config/envs/kallisto.yaml`) — `kb-python` / kallisto\|bustools pseudoalignment
  - `cyclum_env` (`config/envs/cyclum_env.yml`) — Cyclum cell-cycle scoring (legacy step, see below)
  - `sra_tools_env` (`config/envs/sra-tools.yaml`) — `bwa`/`samtools`, used for the 16S read-alignment branch

## Setup

1. Build the four environments above and set their paths under `scanpy_env` / `cyclum_env` / `kallisto_env` / `sra_tools_env` in `config/config.yaml`.
2. Fill in `config/samples.csv` — one row per FASTQ pair, no header, columns `condition, genome, seq_platform, replicate, R1, R2`. `genome` must match a key defined in `config/config.yaml` (see `genome_components`, below). `condition` names containing `embryo` are routed through the embryo arm of the pipeline; everything else is treated as a cultured primary cell line.
3. Pre-build the combined host+*Wolbachia*+16S kallisto\|bustools references for every `genome` key used in `samples.csv` — this happens **offline**, not as part of the Snakemake DAG. See `snakemake_scripts/alignment/build_dmel_dsim_dwil_transcriptomes.sh` and `build_combined_host_wolbachia_references.sh`, and point `config.yaml`'s `<genome_key>:` entries at the resulting directories.
4. Pre-build the Dsim→Dmel ortholog map (`ortholog_map` in `config.yaml`) via the reciprocal-best-hit scripts in `snakemake_scripts/reference/` (see `README_dsim_dmel_orthologs.md` there).
5. Download the Flysta3D-v2 whole-embryo atlas (~12 GB) per `resources/README.md` and point `flysta3d_atlas` at it (defaults to `resources/wcoembed_whole_embeding_downsampled_modified.h5ad`).

## Running

```bash
mamba activate snakemake
snakemake --executor slurm --default-resources slurm_partition=medium slurm_time="2:00:00" runtime=120 mem_mb=8000 -j 16 -n   # dry run
```

Drop `-n` to actually launch jobs. `rule all` is the default target; `rule clean_only` is a lighter alternative target (per-sample annotated/QC'd h5ad + 16S read counts only, skips integration/atlas/trajectory).

## Pipeline structure

The current `rule all` DAG (see `pipeline_rulegraph.svg` / regenerate with `snakemake --rulegraph | dot -Tsvg > pipeline_rulegraph.svg`) has two independent arms that both feed the final target:

**h5ad arm** — per sample, then joint:

1. **`map_pipseq`** / **`map_10x`** — `kb count` (kallisto\|bustools) against that sample's combined host+*Wolbachia*+16S reference → raw per-sample `results/h5ad_results/{sample_id}.h5ad`. Platform is picked automatically per sample from the `_pipseq`/`_10x` suffix Snakemake derives from `samples.csv`'s `seq_platform` column.
2. **`filter_h5ad`** — QC filtering and per-cell *Wolbachia* titer calculation (using the host/symbiont rRNA gene lists resolved per sample from `genome_components`) → `results/filtered_h5ad/{sample_id}.h5ad`.
3. **`integrate`** — projects every filtered sample (embryo **and** cell line, together) onto the frozen Flysta3D-v2 atlas embedding in one shot (`snakemake_scripts/method_comparison/integrate_via_atlas_projection.py`) → `results/integrated/integrated.h5ad`. This is the current default integration path; it replaced an earlier Harmony/BBKNN re-clustering step (`rule integrate_harmony`, removed — see git history).
4. From `results/integrated/integrated.h5ad`:
   - **`titer_by_annotation_atlas`** — *Wolbachia* titer/infection-rate stats and plots grouped by the transferred atlas cell-type annotation (not a Leiden cluster, since none is computed on this embedding) → `results/integrated/figures_atlas/`.
   - **`embryo_to_cellline_trajectory`** — exploratory comparison of the cultured primary cell lines against the embryonic tissue they were derived from: composition, diversity, transfer-confidence, pseudobulk correlation, marker-module scoring, cell-cycle shift, *Wolbachia*-effect, species, and cluster-composition analyses → `results/trajectory_analysis/`.

**16S read-alignment arm** — independent of the h5ad arm above, run directly from raw FASTQs:

5. **`bwa_index_symbiont_genome`** — one shared BWA index per *Wolbachia* strain genome.
6. **`count_16s_reads`** — aligns R2 to the sample's own strain genome and reports reads on that strain's 16S rRNA locus vs. total/mapped reads → `results/rRNA_analysis/read_counts/{sample_id}/{gene}_read_counts.txt`.

### Standalone rules not on the `rule all` critical path

- **`annotate_with_atlas`** and **`map_celllines_to_embryo`** — the original two-step atlas label-transfer approach (Harmony+KNN of embryo samples onto the atlas, then cell lines onto the annotated embryo cells). `integrate` no longer depends on either — every sample now projects onto the atlas directly — but both rules still run standalone if you target their outputs (`results/embryo_annotated/`, `results/celllines_mapped_to_embryo/`) for diagnostics or comparison.
- **`annotate_cell_cycle`** (Cyclum) and **`combine_files_by_condition_platform`** — per-sample cell-cycle scoring and a per-condition Harmony-based replicate merge. Neither feeds into `filter_h5ad`/`integrate`; they're only reachable via `rule clean_only`'s `results/annotated_h5ad/{sample_id}.h5ad` target.

## Repository layout

```
Snakefile                  # pipeline definition (see rule docstrings for per-rule rationale)
config/
  config.yaml              # sample sheet path, env paths, reference paths, per-rule Slurm resources
  samples.csv              # condition, genome, seq_platform, replicate, R1, R2
  envs/*.yml               # conda/mamba env specs (scanpy, kallisto_bustools, cyclum, sra-tools)
snakemake_scripts/
  alignment/                 # offline reference-build scripts (host+symbiont kallisto|bustools indices)
  filtering/                 # QC filtering + titer calculation
  analysis/                  # cell cycle, condition-combine, cell-line<->embryo mapping, trajectory, NMF programs
  method_comparison/          # atlas label transfer, atlas-projection integration, titer-by-annotation, cluster/pathway, pseudotime
  reference/                  # ortholog-map (reciprocal-best-hit) build scripts, rRNA gene-list finders
  plotting/, rRNA_analysis/   # QC and 16S/coverage plotting helpers
resources/                  # large external references (Flysta3D-v2 atlas) — see resources/README.md
results/                    # pipeline outputs (see below)
```

## Output layout

- `results/h5ad_results/` — raw per-sample kallisto\|bustools output
- `results/filtered_h5ad/` — QC-filtered, titer-annotated per-sample h5ad
- `results/embryo_annotated/`, `results/celllines_mapped_to_embryo/` — standalone atlas label-transfer outputs (see above)
- `results/integrated/integrated.h5ad` — atlas-projected integration of all samples; `figures/` and `figures_atlas/` hold its QC and titer-by-annotation plots
- `results/trajectory_analysis/` — embryo→cell-line comparison plots and their matching CSVs
- `results/rRNA_analysis/read_counts/` — per-sample 16S vs. total/mapped read counts
- `results/Figures/` — manuscript figure sources
