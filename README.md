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
- **`pseudotime_cell_states`** (`cell_states.py`): marker-module scores per cell (plasmatocyte, crystal cell, lamellocyte, fat body, neural, neuroblast, muscle, epidermis, germline; plus proliferation, AMP/immune, injury JAK/JNK, apoptosis, hemocyte core) and a marker-based `cell_state`, because the embryo-atlas labels mislabel cultured cells (lines are called "CNS" but lack neural markers). Also finds candidate cell-line precursors in each primary culture (cells closer to the line than to their own culture's average) and tests whether their markers overlap the pseudobulk line-vs-primary up genes. A rank-based test (`precursor_rank_test.csv`) correlates each candidate set's mean expression difference vs. the rest of the primary culture with the line-vs-primary log2FC across non-HVG genes (Spearman), against a null of `precursor_n_perm` random equal-size primary-cell sets; high-proliferation and low-AMP sets of the same size serve as comparisons. The test runs on all genes and again without proliferation and ribosome-biogenesis genes (`gene_set = no_prolif_ribo`: the proliferation module plus leading-edge genes of line-up GSEA terms for cell cycle, DNA replication, ribosome biogenesis, translation and chromatin; list in `precursor_rank_test_excluded_genes.csv`), and a third time also without mitochondrial/respiration genes (`no_prolif_ribo_mito`). Each test reports a random-set null (`z`, `p_perm`) and a null matched on per-cell UMI depth deciles (`z_depth`, `p_perm_depth`); `depth_ratio` is the set's median UMIs over the rest → `results/pseudotime/cell_states/`.
- **`pseudotime_de_state`**: the cell-type-matched DE, using marker-based `cell_state` instead of atlas labels → `results/pseudotime/de_state/`. Not in `rule all`: its automatic state choice picks states that are not meaningful here; run it by name if needed.
- **`pseudotime_infection_de`** (`infection_line_de.py`): pseudobulk DE between an infected line and its uninfected parent (`infection_line_pairs`, default Dsim6B-wMel vs Dsim6B; 2 vs 2 samples), GSEA, a panel of immune/pathway genes, and whether the infection log2FC tracks the line-vs-primary log2FC → `results/pseudotime/infection_de/`. Caveat: Dsim6B-wMel expresses *white* and has lost the plasmatocyte/basement-membrane genes that Dsim6B expresses, so this pair differs in more than infection.
- **`atlas_stats`** (`snakemake_scripts/analysis/atlas_stats.py`): tests for the atlas-projection claims with conditions as replicates: Shannon entropy embryo vs cell line (exact Mann-Whitney), within- vs between-group pseudobulk correlation (exhaustive label permutation), rank of each cell line's parental embryo (`atlas_parents`), tissue specificity, and Wolbachia titer among infected cells (`atlas_infected_lines`) → `results/atlas_stats/`.
- **`pseudotime_primary_proliferative`** (`primary_proliferative.py`): clusters each primary cell line alone and takes its proliferative cells from the `cell_states` `proliferating` flag, whose proliferation z-score spans all stages of the lineage, so primary cells are called proliferating on the same scale as embryos and cell lines (`primary_prolif_definition: shared`). The older within-primary definition (`within`: Leiden clusters or cells with z ≥ 1 inside the primary cell line) always returns the top tail of the culture and is kept only as a sensitivity check. Lineages with fewer than `primary_prolif_min_cells` proliferative cells are listed in `summary.csv` and skipped. Compares it with the rest of the primary cell line: UMIs, genes detected, % mito, doublet score, cell-cycle phase, Wilcoxon markers (cell-cycle genes flagged), marker-state mix, pseudobulk similarity to the cell line (permutation null), Wolbachia titer (Mann-Whitney and OLS with log UMI depth), and embryonic origin from atlas labels and from kNN transfer to the lineage's own embryo cells using only stage-stable genes. GSEA (`gsea/`) on the Wilcoxon z-scores per lineage and combined across lineages (Stouffer), each with and without cell-cycle genes. Also scores each cell for the shared cell-line signature (up minus down genes changing in both species), OXPHOS, injury/JNK and AMP programs, and writes `species_comparison*.csv` (2 lineages per species; descriptive) → `results/pseudotime/primary_proliferative/`.
- **`pseudotime_wolbachia_load`** (`wolbachia_load_primary.py`): tests whether Wolbachia load in primary culture (per-cell Wolbachia gene UMI fraction) goes with a more cell-line-like, proliferative, or apoptotic state. Between lineages it is descriptive (n = 4; species and strain confounded); within each primary cell line it reports partial Spearman correlations controlling for host UMIs and cell state, BH-adjusted and combined across lineages (Stouffer), plus load quintiles → `results/pseudotime/wolbachia_load/`.
- **`pseudotime_wolbachia_states`** (`wolbachia_target_states.py`): tests Wolbachia load against three immortalization-relevant states defined beforehand from host genes only: dividing (S or G2/M score vs expression-matched random gene sets, one scale across stages), cell-line-like (out-of-fold logit of a cell line vs primary classifier on depth-downsampled counts, excluding cell-cycle, stress, AMP/injury and Wolbachia-responsive genes), and death resistance (Diap1/Diap2/Buffy/Bruce minus rpr/hid/grim/skl). Writes a validation table before joining Wolbachia data, then fixed tests: partial Spearman with log host UMIs, group comparisons, Stouffer across primary cell lines → `results/pseudotime/wolbachia_states/`.
- **`pseudotime_culture_stress`** (`culture_stress.py`): tests whether Wolbachia could stabilize early culture. Pre-specified programs (injury/JAK-JNK, antimicrobial peptides, NRF2 redox, iron/heme, starvation, OXPHOS) are compared embryo vs primary per lineage (Cohen's d with bootstrap CI; pseudobulk log2FC) and related to primary Wolbachia load across lineages (exact Spearman, n = 4, descriptive), within cultures (partial Spearman with log host UMIs), and between the infected Dsim6B-wMel and cured Dsim6B lines → `results/pseudotime/culture_stress/`.
- **`pseudotime_rtk_switch`** (`rtk_switch.py`): the Egfr → Pvr receptor switch. Per cell: receptor class (Egfr only / Pvr only / both / neither) and switch index (z Pvr − z Egfr). Tests (1) the stage shift per lineage (Cohen's d, Fisher, pseudobulk log2 Pvr/Egfr); (2) ligand–receptor co-detection, depth-adjusted (Pvf1–3 with Pvr; spi/Krn/vn/grk and rho/ru with Egfr) and sender × receiver state scores with permutation p; (3) hemocyte vs mesoderm identity of Pvr+ cells; (4) Pvr/Egfr in the primary proliferative population; (5) MAPK-target score by receptor class; (6) Dsim vs Dmel (Egfr retained into primary culture; interaction DE); (7) Wolbachia UMIs vs the switch within samples (partial Spearman with log host UMIs, load quintiles, receptor class), across lineages (n = 4, descriptive) and Dsim6B-wMel vs cured Dsim6B → `results/pseudotime/rtk_switch/`.
- **`pseudotime_species_candidates`** (`species_candidates.py`): candidate genes for why D. simulans immortalizes more readily. Ranks culture-emergent species differences (embryos similar; every Dsim primary cell line differs from every Dmel one; flagged if maintained in all cell lines, including uninfected Dsim lines), annotates shared cell-line signature membership and proliferative-population markers, and flags genes whose gene model or ortholog mapping needs checking → `results/pseudotime/species_candidates/`. The primary-vs-embryo interaction is not used because the Dsim and Dmel embryo samples appear stage-mismatched.
- **`pseudotime_silencing`** (`silencing_h3k27me3.py`, runs only when `h3k27me3_beds` is set): fraction of genes with ≥50% of the gene body in H3K27me3 domains, for genes switched off in all lines, consistently down, unchanged and consistently up (line vs primary); Fisher test vs unchanged and a logistic OR adjusted for gene length and primary expression. BEDs must be dm6; the log prints how many Hox genes are marked as a positive control. Also a trend test (`h3k27me3_trend.csv`) and GO enrichment of marked vs unmarked switched-off genes (`go_off_*.csv`) → `results/pseudotime/silencing/`.
- **`pseudotime_cell_states` also writes `composition_tests.csv`** (stage effect on each marker state/flag: paired t-test on per-lineage logits, BH) and BH-adjusted rank-test p-values (`padj`, `padj_depth`).
- **`pseudotime_cell_states` also writes `pvr_mapk_by_sample.csv/.pdf`**: per sample, Spearman rho (raw and controlling for UMI depth) of Pvr/Pvf2 expression with a MAPK-target score (sty, pnt, kek1, Mkp3, aos, CG6006) and the proliferation score; MAPK score in Pvr+ vs Pvr- cells; Pvr/MAPK in the most proliferative primary cells.
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
