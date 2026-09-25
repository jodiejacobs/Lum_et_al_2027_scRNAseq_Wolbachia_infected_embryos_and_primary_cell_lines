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
2. Fill in `config/samples.csv` — one row per FASTQ pair, no header, columns `condition, genome, seq_platform, replicate, R1, R2, sample_type`. `genome` must match a key defined in `config/config.yaml` (see `genome_components`, below). `sample_type` must be one of `embryo`, `primary_cells` (freshly dissociated, briefly cultured -- an intermediate between whole embryos and established lines), or `cell_culture` (established, continuously-cultured lines e.g. JW18wMel, ubkhc, Dsim6B). Rows with `sample_type` `embryo` are routed through the embryo arm of the pipeline; `primary_cells` and `cell_culture` both go through the cell-line arm (same embryo-bridge label transfer for both) but stay distinguishable downstream via `obs['sample_type']` on the integrated object -- see `config/condition_sample_type.tsv`, regenerated from this column on every Snakefile parse.
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

![Pipeline rule graph](pipeline_rulegraph.svg)

The current `rule all` DAG (see `pipeline_rulegraph.svg` / regenerate with `snakemake --rulegraph | dot -Tsvg > pipeline_rulegraph.svg`) has two independent arms that both feed the final target:

**h5ad arm** — per sample, then joint:

1. **`map_pipseq`** / **`map_10x`** — `kb count` (kallisto\|bustools) against that sample's combined host+*Wolbachia*+16S reference → raw per-sample `results/h5ad_results/{sample_id}.h5ad`. Platform is picked automatically per sample from the `_pipseq`/`_10x` suffix Snakemake derives from `samples.csv`'s `seq_platform` column.
2. **`filter_h5ad`** — QC filtering and per-cell *Wolbachia* titer calculation (using the host/symbiont rRNA gene lists resolved per sample from `genome_components`) → `results/filtered_h5ad/{sample_id}.h5ad`.
3. **`integrate`** — projects every filtered sample (embryo **and** cell line, together) onto the frozen Flysta3D-v2 atlas embedding in one shot (`snakemake_scripts/method_comparison/integrate_via_atlas_projection.py`) → `results/integrated/integrated.h5ad`. This is the current default integration path; it replaced an earlier Harmony/BBKNN re-clustering step (`rule integrate_harmony`, removed — see git history).
4. From `results/integrated/integrated.h5ad`:
   - **`titer_by_annotation_atlas`** — *Wolbachia* titer/infection-rate stats and plots grouped by the transferred atlas cell-type annotation (not a Leiden cluster, since none is computed on this embedding) → `results/integrated/figures_atlas/`.
   - **`embryo_to_cellline_trajectory`** — exploratory comparison of the cultured primary cell lines against the embryonic tissue they were derived from: composition, diversity, transfer-confidence, pseudobulk correlation, marker-module scoring, cell-cycle shift, *Wolbachia*-effect, species, and cluster-composition analyses → `results/trajectory_analysis/`.

**Pseudotime arm** (`rule pseudotime_*`, `snakemake_scripts/pseudotime/`): embryo → primary cells → immortalized cell line, one trajectory per matched lineage. Adapted from the SCEPTIC + tradeSeq workflow in [Jacobs et al. 2026](https://github.com/jodiejacobs/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq). Lineages are defined by `pseudotime_lineages` in `config.yaml`:

| Trajectory | Embryo | Primary cells | Cell line(s) |
|---|---|---|---|
| `Dmel_Ubkhc` | ubkhc-wMel-embryos | wMel_Ubkhc_primary_cells | ubkhc |
| `Dmel_JupiterGFP` | JupiterGFP-wMel-embryos | wMel_JupiterGFP_primary_cells | JW18wMel |
| `Dsim_w` | Dsim-w-wRi-embryos | wRi_Dsimw-_primary_cells | Dsim6B, Dsim6B-wMel |
| `Dsim_Merrill23` | Dsim-Merrill23-wRi-embryos | wRi_Merrill23_primary_cells | Dsim-Merrill23 |

- **`pseudotime_prepare`** (`prepare_species.py`): pulls one lineage's counts from `results/filtered_h5ad/` and its atlas labels from `integrated.h5ad`. Dsim is remapped to 1:1 Dmel orthologs first (`pseudotime_ortholog_species`), so every trajectory uses Dmel FlyBase annotation, as in `integrated.h5ad`. It then fits a trajectory-specific PCA; the atlas embedding isn't used because it isn't fit to capture culture adaptation. Embryo cells are kept only if their atlas cell type makes up at least 1% of that lineage's cultured cells. It also runs UMAP, leiden, a diffusion map, DPT rooted in embryo cells, and PAGA as an unsupervised check.
- **`pseudotime_sceptic`** (`run_sceptic_stages.py`): supervised SCEPTIC (xgboost) on the three stage labels (embryo 0, primary_cells 1, cell_culture 2), so every trajectory runs from 0 to 2. Readouts are the confusion matrix, SCEPTIC vs DPT, *Wolbachia* titer vs pseudotime, and atlas identity/confidence along pseudotime. The step also exports a stratified subsample (1,000 cells per sample) for tradeSeq.
- **`pseudotime_tradeseq`** (`tradeseq_species.R`): associationTest, startVsEndTest, and stage-transition tests (embryo→primary, primary→cell line), per trajectory.
- **`pseudotime_nmf`** (`nmf_along_pseudotime.py`): per-trajectory NMF programs and their usage along pseudotime.
- **`pseudotime_joint_export`** / **`pseudotime_joint_tradeseq`** (`export_joint_tradeseq.py`, `tradeseq_joint.R`): for each pair of trajectories, one tradeSeq fit on 1:1 orthologs with trajectory as the condition. Outputs are conditionTest results plus a shape-only similarity (Pearson r between the two smoothers).
- **`pseudotime_compare`** (`compare_species.py`): per pair, shared vs trajectory-specific dynamic genes, logFC concordance per transition, GSEA NES, joint gene classes, and NMF program matching (Hungarian matching on Jaccard of top genes). Pairs default to all within-species and all cross-species lineage pairs (6); override with `pseudotime_comparisons`.
> **Default scope.** The kNN stage-mixing check (`knn_stage_mixing_<lineage>.csv`) showed embryo, primary and cell-line cells share <0.5% of neighbours: there are no intermediate cells to order, so a continuous trajectory isn't supported. By default (`pseudotime_continuous: false`) the pipeline runs prepare, SCEPTIC (kept as a per-cell stage-similarity score), the pseudobulk DE (main analysis), the atlas composition by stage, and the merge. tradeSeq, joint tradeSeq, NMF along pseudotime, the pairwise comparisons and the DE-vs-tradeSeq concordance run only with `pseudotime_continuous: true`.

- **`pseudotime_composition`** (`composition_by_stage.py`): per lineage, Flysta3D atlas cell-type composition at each stage (whole embryo, primary cells, cell line; confidence ≥ 0.5) and log2 enrichment vs. that lineage's embryo, i.e. which embryonic cell types carry through into culture → `results/pseudotime/composition/`.
- **`pseudotime_de`** (`de_stages.py`): discrete counterpart to the pseudotime. One pseudobulk profile per sample (same trajectory cells), PyDESeq2 with lineages as replicates: per species (`~lineage + stage`), pooled across species on shared genes, and a species × stage interaction (Dsim change minus Dmel change) for embryo→primary, primary→line and embryo→line. Also per-lineage sign consistency of pooled DE genes and preranked GSEA on every contrast → `results/pseudotime/de/`. Needs `pydeseq2==0.5.2` in `scanpy_env` (0.5.3+ forces numpy>=2, which breaks numba).
- **`pseudotime_de_downsampled`**: the same DE with every cell downsampled to `de_downsample_counts` UMIs (default 2,000) before pseudobulk, because primary cells were sequenced ~3–4× deeper than the lines and embryos. `robustness_vs_reference.csv` gives, per contrast, how many of the main run's DE genes survive depth matching (split by up/down) and fold-change correlation → `results/pseudotime/de_downsampled/`.
- **`pseudotime_de_celltype`**: the same models within one atlas cell type at a time (auto-selected types present at all three stages in ≥3 lineages; `celltype_selection.csv`), separating "cell types are lost in culture" from "cells change expression". One subfolder per cell type, each with `robustness_vs_reference.csv` → `results/pseudotime/de_celltype/`.
- **`pseudotime_cell_states`** (`cell_states.py`): marker-module scores per cell (plasmatocyte, crystal cell, lamellocyte, fat body, neural, neuroblast, muscle, epidermis, germline; plus proliferation, AMP/immune, injury JAK/JNK, apoptosis, hemocyte core) and a marker-based `cell_state`, because the embryo-atlas labels mislabel cultured cells (lines are called "CNS" but lack neural markers). Also finds candidate cell-line precursors in each primary culture (cells closer to the line than to their own culture's average) and tests whether their markers overlap the pseudobulk line-vs-primary up genes. A rank-based test (`precursor_rank_test.csv`) correlates each candidate set's mean expression difference vs. the rest of the primary culture with the line-vs-primary log2FC across non-HVG genes (Spearman), against a null of `precursor_n_perm` random equal-size primary-cell sets; high-proliferation and low-AMP sets of the same size serve as comparisons → `results/pseudotime/cell_states/`.
- **`pseudotime_de_state`**: the cell-type-matched DE, using marker-based `cell_state` instead of atlas labels → `results/pseudotime/de_state/`. Not in `rule all`: its automatic state choice picks states that are not meaningful here; run it by name if needed.
- **`pseudotime_infection`** (`infection_by_stage.py`): Wolbachia gene fraction, Wolbachia and other-bacteria 16S reads per sample and per cell by stage, and whether the primary-culture AMP program tracks Wolbachia or bacterial reads → `results/pseudotime/infection/`. `sample_summary.csv` compares observed infection (Wolbachia gene fraction ≥ 1e-3) with `infection_status` in the config and flags mismatches; set it for conditions whose FASTQ names are misleading (the ubkhc line is uninfected). The 16S outputs are not used for interpretation.
- **`pseudotime_de_concordance`** (`de_vs_tradeseq.py`): per trajectory, agreement between the pseudobulk DE contrasts and the tradeSeq transition tests (fold-change rank correlation, overlap, sign agreement) → `results/pseudotime/de/concordance/`.
- **`pseudotime_merge`** (`merge_pseudotime.py`): `results/pseudotime/integrated_with_pseudotime.h5ad`, which is `integrated.h5ad` plus the `species`, `lineage`, `stage_numeric`, `pt_in_trajectory`, `sceptic_pseudotime`, `sceptic_pred_stage`, `sceptic_prob_*`, and `dpt_pseudotime` columns in `.obs`. Each cell's pseudotime comes from its own lineage's trajectory.

Requires `sceptic` + `xgboost` (Python) and `tradeSeq`, `SingleCellExperiment`, `BiocParallel`, `patchwork` (and optionally `ggrepel`) (R) in `scanpy_env`, the same setup as Jacobs et al. 2026. Stage is confounded with sequencing run, so treat the stage axis as stage + batch.

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
  pseudotime/                # embryo -> primary -> cell line pseudotime (SCEPTIC, DPT/PAGA, tradeSeq, NMF, species comparison)
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
- `results/pseudotime/<trajectory>/` — per-lineage prepared/SCEPTIC h5ad, `figures/`, `tradeseq/`, `nmf/`; `results/pseudotime/pairs/<A>__vs__<B>/` — joint tradeSeq (`tradeseq/`) and comparison tables/plots (`compare/`) per pair; `results/pseudotime/integrated_with_pseudotime.h5ad` — merged object
- `results/rRNA_analysis/read_counts/` — per-sample 16S vs. total/mapped read counts
- `results/Figures/` — manuscript figure sources
