# Snakefile for processing scRNA-seq data with Snakemake and kallisto bustools
# mamba activate snakemake #Needs snakemake>=9.0
# snakemake --executor slurm --default-resources slurm_partition=medium slurm_time="2:00:00" runtime=120 mem_mb=8000 -j 16 -n 

import pandas as pd
import os
import re

# Configuration
configfile: "config/config.yaml"
SCANPY_ENV = config["scanpy_env"]
CYCLUM_ENV = config["cyclum_env"]
KALLISTO_ENV = config["kallisto_env"]
SRA_TOOLS_ENV = config["sra_tools_env"]


# Load samples information
# Columns: 0=condition, 1=genome, 2=seq_platform, 3=replicate, 4=R1, 5=R2,
# 6=sample_type (embryo / primary_cells / cell_culture -- explicit column
# added when primary_cells was introduced as a third category between
# whole embryos and established, continuously-cultured cell lines; see
# EMBRYO_SAMPLE_IDS / PRIMARY_CELLS_SAMPLE_IDS / CELL_CULTURE_SAMPLE_IDS below)
samples_df = pd.read_csv(config["samples_file"], header=None, sep=',')
samples_df = samples_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

# Change index format to use hyphen between condition and replicate
samples_df.index = [(row[0] + "-" + str(row[3]) + "_" + row[2]) for _, row in samples_df.iterrows()]

# Fail fast with a clear message if two rows collapse to the same sample_id
# (this happens when condition + replicate + seq_platform is duplicated in
# samples.csv, e.g. two rows with the same condition/replicate/platform but
# different genome). A duplicated index turns samples_df.loc[sample_id] into
# a DataFrame instead of a Series everywhere below, which surfaces as
# confusing errors like "unhashable type: 'Series'".
_dupe_ids = sorted(set(samples_df.index[samples_df.index.duplicated(keep=False)]))
if _dupe_ids:
    raise ValueError(
        f"Duplicate sample_id(s) in {config['samples_file']}: {_dupe_ids}. "
        "Each condition/replicate/seq_platform combination must be unique — "
        "check for duplicate rows or a reused replicate number for that condition/platform."
    )

# Get actual sample IDs that exist
SAMPLE_IDS = samples_df.index.tolist()

# Validate the explicit sample_type column (column 6) -- added when
# primary_cells was introduced as a third category between whole embryos
# and established cell_culture lines. This used to be inferred by checking
# whether "embryo" appeared in the condition string (column 0), which had
# no way to represent a third category; samples.csv now says explicitly
# what each row is, so a typo in the condition name (e.g. a sample missing
# "embryo"/"primary_cells" in its name) can't silently misclassify it.
_ALLOWED_SAMPLE_TYPES = {"embryo", "primary_cells", "cell_culture"}
_bad_sample_types = sorted(set(samples_df[6].astype(str).str.strip()) - _ALLOWED_SAMPLE_TYPES)
if _bad_sample_types:
    raise ValueError(
        f"Unrecognized sample_type value(s) in {config['samples_file']}: {_bad_sample_types}. "
        f"Column 7 (0-indexed 6) of every row must be one of {sorted(_ALLOWED_SAMPLE_TYPES)}."
    )

# Split samples by sample_type (column 6). EMBRYO_SAMPLE_IDS routes embryo
# samples through the Flysta3D-v2 atlas cell-type transfer (rule
# annotate_with_atlas); PRIMARY_CELLS_SAMPLE_IDS and CELL_CULTURE_SAMPLE_IDS
# are the finer 3-way split (used for per-sample-type metadata/analysis --
# see condition_sample_type.tsv below); CELLLINE_SAMPLE_IDS keeps its
# original meaning (everything non-embryo) since routing through the
# embryo-reference transfer (rule map_celllines_to_embryo) is the same for
# primary_cells and cell_culture samples -- both need the embryo bridge,
# not the atlas directly -- before both arms feed into rule integrate.
EMBRYO_SAMPLE_IDS        = [s for s in SAMPLE_IDS if samples_df.loc[s][6] == "embryo"]
PRIMARY_CELLS_SAMPLE_IDS = [s for s in SAMPLE_IDS if samples_df.loc[s][6] == "primary_cells"]
CELL_CULTURE_SAMPLE_IDS  = [s for s in SAMPLE_IDS if samples_df.loc[s][6] == "cell_culture"]
CELLLINE_SAMPLE_IDS      = PRIMARY_CELLS_SAMPLE_IDS + CELL_CULTURE_SAMPLE_IDS

# condition -> sample_type lookup, regenerated from samples.csv on every
# Snakefile parse (dry runs included) so it can never drift out of sync.
# Single source of truth for the standalone scripts (integrate_via_atlas_
# projection.py, embryo_to_cellline_trajectory.py, analyze_titer_by_
# annotation.py) that operate on globbed h5ad files after the fact and
# have no access to samples.csv itself -- see --condition_sample_type in
# rule integrate below.
CONDITION_SAMPLE_TYPE_PATH = "config/condition_sample_type.tsv"
_condition_sample_type = samples_df.drop_duplicates(subset=[0]).set_index(0)[6]
_condition_sample_type.rename_axis("condition").rename("sample_type").to_csv(
    CONDITION_SAMPLE_TYPE_PATH, sep="\t"
)

CONDITIONS = samples_df[0].unique().tolist()
REPLICATES = samples_df[3].unique().tolist()
SEQUENCING_PLATFORMS = samples_df[2].unique().tolist()

# Create a dictionary to track which condition-seq_platform combinations exist
# and what replicates are available for each
CONDITION_PLATFORM_REPS = {}
for sample_id in SAMPLE_IDS:
    condition = samples_df.loc[sample_id][0]
    seq_platform = samples_df.loc[sample_id][2]
    replicate = samples_df.loc[sample_id][3]

    key = (condition, seq_platform)
    if key not in CONDITION_PLATFORM_REPS:
        CONDITION_PLATFORM_REPS[key] = []
    CONDITION_PLATFORM_REPS[key].append(replicate)

# Get list of condition-seq_platform combinations that actually exist
CONDITION_PLATFORM_COMBOS = list(CONDITION_PLATFORM_REPS.keys())

# Find conditions that have BOTH 10x and pipseq data for method comparison
CONDITIONS_WITH_BOTH_METHODS = []
for condition in CONDITIONS:
    has_10x = any(p == '10x' for c, p in CONDITION_PLATFORM_COMBOS if c == condition)
    has_pipseq = any(p == 'pipseq' for c, p in CONDITION_PLATFORM_COMBOS if c == condition)
    if has_10x and has_pipseq:
        CONDITIONS_WITH_BOTH_METHODS.append(condition)

print(f"Conditions with both methods for comparison: {CONDITIONS_WITH_BOTH_METHODS}")

# Sample lookup by sample name; samples_df.loc['sample_name'][3] gives the R1 path
print(f"Conditions: {CONDITIONS}")
print(f"Replicates: {REPLICATES}")
print(f"Sequencing Platforms: {SEQUENCING_PLATFORMS}")
print(f"Sample IDs: {SAMPLE_IDS}")
print(f"Condition-Platform combinations: {CONDITION_PLATFORM_COMBOS}")
print(f"Embryo sample IDs: {EMBRYO_SAMPLE_IDS}")
print(f"Primary cells sample IDs: {PRIMARY_CELLS_SAMPLE_IDS}")
print(f"Cell culture sample IDs: {CELL_CULTURE_SAMPLE_IDS}")
print(f"Cell line sample IDs (primary_cells + cell_culture): {CELLLINE_SAMPLE_IDS}")

# Helper function to get fastq files for a sample
def get_fastq_files(sample_id):
    """Get R1 and R2 fastq files for a sample."""
    sample_info = samples_df.loc[sample_id]
    r1_path = sample_info[4]  # R1 path is in column 4 (0-indexed)
    r2_path = sample_info[5]  # R2 path is in column 5 (0-indexed)
    return r1_path, r2_path

# Helper function to get the reference genome name for a sample (looked up in config.yaml)
def get_genome(sample_id):
    """Get the genome name for a sample; used as a key into config.yaml for the kallisto index path."""
    sample_info = samples_df.loc[sample_id]
    return sample_info[1]  # Genome is in column 1 (0-indexed)

# Helper function to get replicates for a condition-seq_platform combo
def get_replicates_for_combo(condition, seq_platform):
    """Get the list of replicates that exist for a given condition-seq_platform combo."""
    key = (condition, seq_platform)
    return CONDITION_PLATFORM_REPS.get(key, [])

# Helper functions to resolve the rRNA gene lists used for the Wolbachia titer
# calculation, for any host species x Wolbachia strain combination -- looks
# up genome_components[genome] (host, symbiont) in config.yaml, then that
# host's / that strain's gene list path. Missing entries return None; the
# filter script treats a missing/absent file as "skip titer for this sample"
# rather than crashing, so new genome keys can be added incrementally.
def get_host_rrna_genes(sample_id):
    genome = get_genome(sample_id)
    components = config.get("genome_components", {}).get(genome)
    if not components:
        return ""
    return config.get("host_rrna_genes", {}).get(components["host"], "")

def get_symbiont_rrna_genes(sample_id):
    genome = get_genome(sample_id)
    components = config.get("genome_components", {}).get(genome)
    if not components:
        return ""
    return config.get("symbiont_rrna_genes", {}).get(components["symbiont"], "")

def get_symbiont_strain(sample_id):
    """This sample's Wolbachia strain key (e.g. 'wMel'), via genome_components."""
    genome = get_genome(sample_id)
    components = config.get("genome_components", {}).get(genome)
    return components["symbiont"] if components else None

def get_symbiont_fasta(sample_id):
    """This sample's own Wolbachia strain genome fasta -- used as the BWA
    reference for count_16s_reads so non-wMel strains align (and are
    therefore counted) against their own genome instead of a wMel-only
    reference. Falls back to config['ref_fasta'] if the strain can't be
    resolved."""
    strain = get_symbiont_strain(sample_id)
    return config.get("wolbachia_genome", {}).get(strain, {}).get("fasta", config["ref_fasta"])

# Directory-name wildcard <-> Wolbachia strain key, precomputed once, used by
# rule bwa_index_symbiont_genome (Snakemake doesn't allow output: to be a
# function, only input:, so the per-strain genome directory name -- already
# unique per strain -- is used as the rule's wildcard instead of the strain
# key itself).
def _wolbachia_dir_name(strain_key):
    return os.path.basename(os.path.dirname(config["wolbachia_genome"][strain_key]["fasta"]))

WOLBACHIA_DIR_BY_STRAIN = {s: _wolbachia_dir_name(s) for s in config.get("wolbachia_genome", {})}
STRAIN_BY_WOLBACHIA_DIR = {v: k for k, v in WOLBACHIA_DIR_BY_STRAIN.items()}
WOLBACHIA_REF_ROOT = (
    os.path.dirname(os.path.dirname(config["wolbachia_genome"]["wMel"]["fasta"]))
    if config.get("wolbachia_genome") else ""
)

def get_bwa_index_flag(sample_id):
    """Path to this sample's Wolbachia strain's shared BWA-index-done flag
    (rule bwa_index_symbiont_genome) -- lets count_16s_reads depend on the
    index instead of every per-sample job racing to build it themselves
    when several samples share the same strain."""
    strain = get_symbiont_strain(sample_id)
    strain_dir = WOLBACHIA_DIR_BY_STRAIN.get(strain)
    if not strain_dir:
        return []
    return WOLBACHIA_REF_ROOT + "/" + strain_dir + "/.bwa_index.done"

_SYMBIONT_16S_REGION_CACHE = {}

def get_symbiont_16s_region(sample_id):
    """Resolve 'locus_tag::seqid:start-end' for this sample's Wolbachia
    strain's 16S rRNA gene: looks up the single locus tag configured per
    strain in config['symbiont_16s_gene'], then finds its coordinates by
    scanning that strain's own GTF (wolbachia_genome[strain].gtf). Falls
    back to config['rRNA_16S_region'] (or its wMel-only default) if the
    strain/locus_tag/GTF aren't resolvable or the locus_tag isn't found, so
    an unconfigured strain degrades gracefully instead of crashing the DAG.
    Result is cached per (strain, locus_tag) so the GTF is only scanned
    once even though every sample of that strain calls this.
    """
    strain = get_symbiont_strain(sample_id)
    fallback = config.get("rRNA_16S_region", "GQX67_05945::NZ_CP046925.1:1167785-1169290")
    locus_tag = config.get("symbiont_16s_gene", {}).get(strain)
    gtf = config.get("wolbachia_genome", {}).get(strain, {}).get("gtf")
    if not (strain and locus_tag and gtf):
        return fallback

    cache_key = (strain, locus_tag)
    if cache_key in _SYMBIONT_16S_REGION_CACHE:
        return _SYMBIONT_16S_REGION_CACHE[cache_key]

    region = fallback
    needle = locus_tag
    try:
        with open(gtf) as fh:
            for line in fh:
                if line.startswith("#") or needle not in line:
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                attrs = fields[8]
                # NCBI cmsearch/Infernal (Rfam) rRNA features don't carry a
                # locus_tag attribute at all -- their gene_id/transcript_id
                # are the locus tag prefixed with "gene-"/"rna-" (e.g.
                # gene_id "gene-WRI_RS06005" instead of locus_tag
                # "WRI_RS06005"), so an exact `locus_tag "X"` / `gene_id "X"`
                # match silently misses them. Match with an optional
                # gene-/rna- prefix so both annotation styles resolve.
                if (f'locus_tag "{locus_tag}"' in attrs
                        or re.search(
                            rf'(?:gene_id|transcript_id) "(?:gene-|rna-)?{re.escape(locus_tag)}"',
                            attrs)):
                    region = f"{locus_tag}::{fields[0]}:{fields[3]}-{fields[4]}"
                    break
        if region == fallback:
            print(f"WARNING: locus_tag {locus_tag!r} not found in {gtf} -- "
                  f"using fallback 16S region for strain {strain}")
    except FileNotFoundError:
        print(f"WARNING: GTF not found for strain {strain} ({gtf}) -- "
              f"using fallback 16S region")

    _SYMBIONT_16S_REGION_CACHE[cache_key] = region
    return region

# Main rule that defines the final output
rule all:
    input:
        # Integration
        "results/integrated/integrated.h5ad",
        # Wolbachia titer by transferred atlas cell-type annotation
        "results/integrated/figures_atlas/.titer_by_annotation.done",
        # Embryo -> cell line trajectory/identity analysis
        "results/trajectory_analysis/.done",
        expand("results/rRNA_analysis/read_counts/{sample_id}/{gene}_read_counts.txt",
               sample_id=SAMPLE_IDS,
               gene=config.get("target_genes", ["GQX67_05945"]))

rule clean_only:
    input:
        expand("results/annotated_h5ad/{sample_id}.h5ad", sample_id=SAMPLE_IDS),
        expand("results/rRNA_analysis/read_counts/{sample_id}/{gene}_read_counts.txt",
                sample_id=SAMPLE_IDS,
                gene=config.get("target_genes", ["GQX67_05945"]))

# ─────────────────────────────────────────────────────────────────────────────
# Reference resolution: all combined host/Wolbachia/16S kallisto|bustools
# references are pre-built (see snakemake_scripts/alignment/
# build_dmel_dsim_dwil_transcriptomes.sh and
# build_combined_host_wolbachia_references.sh for the offline build commands).
# map_pipseq/map_10x pull index.idx/t2g.txt straight from the genome
# directories listed in config.yaml (config[genome_key]) as plain, already-
# existing input files -- there is no in-DAG build rule, so Snakemake never
# tries to (re)build a reference; it just errors if a config path is missing.
# ─────────────────────────────────────────────────────────────────────────────

# Establish rule precedencex
ruleorder: map_pipseq > combine_files_by_condition_platform
ruleorder: map_10x > combine_files_by_condition_platform

# Process PIPseq samples with kallisto bustools
rule map_pipseq:
    input:
        read1 = lambda wildcards: get_fastq_files(wildcards.sample_id)[0],
        read2 = lambda wildcards: get_fastq_files(wildcards.sample_id)[1],
        kallisto_index = lambda wildcards: os.path.join(config[get_genome(wildcards.sample_id)], "index.idx"),
        transcripts_to_genes = lambda wildcards: os.path.join(config[get_genome(wildcards.sample_id)], "t2g.txt")
    output:
        h5ad = "results/h5ad_results/{sample_id}.h5ad",
        bus = "results/pipseq/{sample_id}/output.unfiltered.bus",
        ec = "results/pipseq/{sample_id}/matrix.ec",
        transcripts = "results/pipseq/{sample_id}/transcripts.txt"
    params:
        sample_id = "{sample_id}",
        outdir = "results/pipseq/{sample_id}",
        genome = lambda wildcards: get_genome(wildcards.sample_id)
    wildcard_constraints:
        sample_id = ".*_pipseq"  # Only match samples ending with _pipseq
    log:
        "logs/pipseq/{sample_id}.log"
    threads:
        config["pipseeker_threads"]
    resources:
        slurm_partition = config["pipseeker_partition"],
        mem_mb = config["pipseeker_mem"],
        slurm_time = config["pipseeker_time"]
    shell:
        """
        exec > {log} 2>&1
        echo "Starting PIPseq processing for {params.sample_id}"
        echo "Input files: {input.read1}, {input.read2}"
        echo "Output directory: {params.outdir}"

        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {KALLISTO_ENV}

        kb count \
            -i {input.kallisto_index} \
            --keep-tmp \
            -g {input.transcripts_to_genes} \
            -x 0,0,16:0,16,28:1,0,0 \
            -o {params.outdir} \
            -t {threads} \
            --h5ad \
            {input.read1} {input.read2}

        echo "Moving h5ad file to final location"
        # Move the h5ad file to the expected location
        mv {params.outdir}/counts_unfiltered/adata.h5ad {output.h5ad}
        echo "PIPseq processing complete for {params.sample_id}"
        """

# Process 10X samples with kallisto bustools
rule map_10x:
    input:
        read1 = lambda wildcards: get_fastq_files(wildcards.sample_id)[0],
        read2 = lambda wildcards: get_fastq_files(wildcards.sample_id)[1],
        kallisto_index = lambda wildcards: os.path.join(config[get_genome(wildcards.sample_id)], "index.idx"),
        transcripts_to_genes = lambda wildcards: os.path.join(config[get_genome(wildcards.sample_id)], "t2g.txt")
    output:
        h5ad = "results/h5ad_results/{sample_id}.h5ad",
        bus = "results/10x/{sample_id}/output.unfiltered.bus",
        ec = "results/10x/{sample_id}/matrix.ec",
        transcripts = "results/10x/{sample_id}/transcripts.txt"
    params:
        sample_id = "{sample_id}",
        outdir = "results/10x/{sample_id}",
        genome = lambda wildcards: get_genome(wildcards.sample_id)
    wildcard_constraints:
        sample_id = ".*_10x"  # Only match samples ending with _10x
    log:
        "logs/10x/{sample_id}.log"
    threads:
        config["cellranger_threads"]
    resources:
        slurm_partition = config["cellranger_partition"],
        mem_mb = config["cellranger_mem"],
        slurm_time = config["cellranger_time"]
    shell:
        """
        exec > {log} 2>&1
        echo "Starting 10x processing for {params.sample_id}"
        echo "Input files: {input.read1}, {input.read2}"
        echo "Output directory: {params.outdir}"

        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {KALLISTO_ENV}

        kb count \
            --kallisto /private/home/jomojaco/kallisto/build/src/kallisto \
            -i {input.kallisto_index} \
            -g {input.transcripts_to_genes} \
            --keep-tmp \
            -x 10xv3 \
            -o {params.outdir} \
            -t {threads} \
            --h5ad \
            {input.read1} {input.read2}

        echo "Moving h5ad file to final location"
        # Move the h5ad file to the expected locatio
        mv {params.outdir}/counts_unfiltered/adata.h5ad {output.h5ad}
        echo "10x processing complete for {params.sample_id}"
        """

# Filter h5ad output and output qc:
rule filter_h5ad:
    input:
        h5ad = "results/h5ad_results/{sample_id}.h5ad",
        host_rrna_genes = lambda wildcards: get_host_rrna_genes(wildcards.sample_id),
        symbiont_rrna_genes = lambda wildcards: get_symbiont_rrna_genes(wildcards.sample_id)
    output:
        filtered_h5ad = "results/filtered_h5ad/{sample_id}.h5ad"
    params:
        script = config["filter_script"]
    log:
        "logs/filter/{sample_id}.log"
    threads:
        config["filter_threads"]
    resources:
        slurm_partition = config["filter_partition"],
        mem_mb = config["filter_mem"],
        slurm_time = config["filter_time"]
    shell:
        """
        exec > {log} 2>&1
        echo "Starting filtering for {wildcards.sample_id}"
        echo "Input file: {input.h5ad}"

        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --input {input.h5ad} \
            --output {output.filtered_h5ad} \
            --host_rrna_genes "{input.host_rrna_genes}" \
            --symbiont_rrna_genes "{input.symbiont_rrna_genes}"

        echo "Compressing original h5ad file"
        gzip {input.h5ad}
        echo "Filtering complete for {wildcards.sample_id}"
        """

# Cell cycle Annotation
rule annotate_cell_cycle: # This needs the cyclum conda environment
    input:
        h5ad = "results/filtered_h5ad/{sample_id}.h5ad" # Output of filtered script 
    output:
        annotated_h5ad = "results/annotated_h5ad/{sample_id}.h5ad"
    params:
        script = config["cell_cycle_script"]
    log:
        "logs/annotate/{sample_id}.log"
    threads:
        config["cell_cycle_threads"]
    resources:
        slurm_partition = config["cell_cycle_partition"],
        mem_mb = config["cell_cycle_mem"],
        slurm_time = config["cell_cycle_time"]
    shell:
        """
        exec > {log} 2>&1
        echo "Starting cell cycle annotation for {wildcards.sample_id}"
        echo "Input file: {input.h5ad}"
        
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {CYCLUM_ENV}

        python {params.script} \
            --input {input.h5ad} \
            --output {output.annotated_h5ad}
        
        echo "Compressing filtered h5ad file"
        gzip {input.h5ad}
        echo "Cell cycle annotation complete for {wildcards.sample_id}"
        """

rule combine_files_by_condition_platform:
    input:
        input_files = lambda wildcards: expand("results/annotated_h5ad/{condition}-{replicate}_{seq_platform}.h5ad", 
                                              condition=wildcards.condition, 
                                              replicate=get_replicates_for_combo(wildcards.condition, wildcards.seq_platform),
                                              seq_platform=wildcards.seq_platform)
    output:
        combined = "results/combined/{condition}_{seq_platform}.h5ad"
    params:
        combine_script = config["combine_script"],
        sample = '{condition}_{seq_platform}',  # Remove {replicate} since you're combining replicates
        fig_dir = 'results/combined/{condition}_{seq_platform}'  # Remove {replicate}
    log:
        "logs/combine/{condition}_{seq_platform}.log"
    threads: 
        config["combine_threads"]
    resources:
        slurm_partition = config["combine_partition"],
        mem_mb = config["combine_mem"],
        slurm_time = config["combine_time"]
    shell:
        """
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        python {params.combine_script} --files {input.input_files} --out_path {output.combined} --fig_dir {params.fig_dir} --sample {params.sample} --batch_key batch --min_cells 3 --min_genes 200
        """

##################################################################
# rRNA analysis rules
##################################################################

# One-time BWA index per Wolbachia strain genome, shared by every sample
# mapped to that strain. Built as its own rule (rather than indexing
# ad hoc inside count_16s_reads' shell block) so Snakemake's own DAG
# dedups it -- without this, multiple samples of the same strain running
# concurrently under SLURM would all try to `bwa index` the same fasta at
# once and race/corrupt each other's index files.
rule bwa_index_symbiont_genome:
    input:
        fasta = lambda wildcards: config["wolbachia_genome"][STRAIN_BY_WOLBACHIA_DIR[wildcards.strain_dir]]["fasta"]
    output:
        flag = touch(WOLBACHIA_REF_ROOT + "/{strain_dir}/.bwa_index.done")
    wildcard_constraints:
        strain_dir = "|".join(STRAIN_BY_WOLBACHIA_DIR.keys())
    log:
        "logs/reference/bwa_index_{strain_dir}.log"
    resources:
        slurm_partition = "medium", mem_mb = 8000, slurm_time = "1:00:00"
    shell:
        """
        exec > {log} 2>&1
        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {SRA_TOOLS_ENV}
        bwa index {input.fasta}
        """

##################################################################
# Flysta3D-v2 atlas cell-type label transfer (embryo samples only)
##################################################################
# Runs BEFORE integrate: transfers cell-type labels from the Wang et al.
# 2025 Cell / Flysta3D-v2 Drosophila embryo atlas
# (https://db.cngb.org/stomics/flysta3d-v2/) onto every EMBRYO sample's
# filtered h5ad in one joint Harmony+KNN run (see docstring in
# snakemake_scripts/method_comparison/annotate_with_flysta3d.py for why
# this happens per-sample, before your own condition integration, rather
# than after it). Restricted to EMBRYO_SAMPLE_IDS -- the atlas is a
# whole-embryo atlas, so cell-type identification against it only makes
# sense for embryo samples; primary cell line samples get their cell type
# via rule map_celllines_to_embryo below instead, which maps them onto
# these already-annotated embryo cells rather than the atlas directly.
# D. simulans embryo samples (var_names in the Dsim NCBI/Gnomon namespace)
# are remapped to Dmel FlyBase orthologs internally via --ortholog_map, so
# they compare correctly against the Dmel-indexed atlas.
#
# First run with atlas_label_cols left empty in config.yaml to print the
# atlas's candidate obs columns and stop; fill in atlas_label_cols with the
# real column name(s) (e.g. ["cell_type"]) once you know them, then rerun.
# rule integrate (via rule map_celllines_to_embryo) automatically picks up
# whatever atlas_<label> columns this rule writes -- no further wiring
# needed once atlas_label_cols is set.
rule annotate_with_atlas:
    input:
        files              = expand("results/filtered_h5ad/{sample_id}.h5ad", sample_id=EMBRYO_SAMPLE_IDS),
        atlas              = config.get("flysta3d_atlas", "resources/wcoembed_whole_embeding_downsampled_modified.h5ad"),
        flybase_annotation = config["flybase_annotation"],
        orthologs          = config["ortholog_map"],
    output:
        annotated = expand("results/embryo_annotated/{sample_id}.h5ad", sample_id=EMBRYO_SAMPLE_IDS),
    params:
        script       = config.get("annotate_atlas_script",
                           "snakemake_scripts/method_comparison/annotate_with_flysta3d.py"),
        out_dir      = "results/embryo_annotated",
        fig_dir      = "results/embryo_annotated/figures",
        k            = config.get("atlas_k", 30),
        n_pcs        = config.get("atlas_n_pcs", 30),
        harmony_vars = config.get("atlas_harmony_vars", ["dataset", "method"]),
        # Precomputed as plain strings so the shell block below stays a
        # straight template -- easier to read/debug than inline conditionals.
        label_cols_flag = (
            "--label_cols " + " ".join(config["atlas_label_cols"])
            if config.get("atlas_label_cols") else ""
        ),
        subsample_flag = (
            f"--subsample_ref {config['atlas_subsample_ref']}"
            if config.get("atlas_subsample_ref") else ""
        ),
    log:
        "logs/annotate_with_atlas/annotate_with_atlas.log"
    threads:
        config.get("annotate_atlas_threads", 16)
    resources:
        slurm_partition = config.get("annotate_atlas_partition", "medium"),
        mem_mb          = config.get("annotate_atlas_mem", 500000),
        slurm_time      = config.get("annotate_atlas_time", "12:00:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Starting Flysta3D-v2 atlas label transfer (embryo samples only)"

        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --atlas {input.atlas} \
            --query {input.files} \
            --out_dir {params.out_dir} \
            --fig_dir {params.fig_dir} \
            --flybase_annotation {input.flybase_annotation} \
            --ortholog_map {input.orthologs} \
            --k {params.k} \
            --n_pcs {params.n_pcs} \
            --harmony_vars {params.harmony_vars} \
            {params.label_cols_flag} \
            {params.subsample_flag}

        echo "Atlas label transfer complete"
        """

##################################################################
# Map primary cell line samples onto the annotated embryo cells
##################################################################
# Runs AFTER annotate_with_atlas, BEFORE integrate: transfers the
# atlas-derived cell-type label(s) from the now-annotated embryo cells onto
# every remaining (non-embryo) sample -- the cultured primary cell lines --
# via the same Harmony+KNN recipe (see docstring in
# snakemake_scripts/analysis/map_cellline_to_embryo.py for why this maps
# through your own embryo cells rather than the atlas directly). D.
# simulans cell line samples are remapped to Dmel FlyBase orthologs
# internally via --ortholog_map, same as the embryo arm above.
rule map_celllines_to_embryo:
    input:
        reference = expand("results/embryo_annotated/{sample_id}.h5ad", sample_id=EMBRYO_SAMPLE_IDS),
        query     = expand("results/filtered_h5ad/{sample_id}.h5ad", sample_id=CELLLINE_SAMPLE_IDS),
        orthologs = config["ortholog_map"],
    output:
        mapped = expand("results/celllines_mapped_to_embryo/{sample_id}.h5ad", sample_id=CELLLINE_SAMPLE_IDS),
    params:
        script       = config.get("map_cellline_script",
                           "snakemake_scripts/analysis/map_cellline_to_embryo.py"),
        out_dir      = "results/celllines_mapped_to_embryo",
        fig_dir      = "results/celllines_mapped_to_embryo/figures",
        k            = config.get("cellline_embryo_k", 30),
        n_pcs        = config.get("cellline_embryo_n_pcs", 30),
        harmony_vars = config.get("cellline_embryo_harmony_vars", ["dataset", "method"]),
        # map_cellline_to_embryo.py's --label_cols expects the actual
        # atlas_<label> column names rule annotate_with_atlas wrote (not
        # the bare atlas_label_cols entries used for that rule's own
        # --label_cols flag), so prefix each one with "atlas_" here.
        label_cols_flag = (
            "--label_cols " + " ".join(f"atlas_{c}" for c in config["atlas_label_cols"])
            if config.get("atlas_label_cols") else ""
        ),
    log:
        "logs/map_celllines_to_embryo/map_celllines_to_embryo.log"
    threads:
        config.get("cellline_embryo_threads", 16)
    resources:
        slurm_partition = config.get("cellline_embryo_partition", "medium"),
        mem_mb          = config.get("cellline_embryo_mem", 250000),
        slurm_time      = config.get("cellline_embryo_time", "8:00:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Mapping primary cell line samples onto the annotated embryo reference"

        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --reference {input.reference} \
            --query {input.query} \
            --out_dir {params.out_dir} \
            --fig_dir {params.fig_dir} \
            --ortholog_map {input.orthologs} \
            --k {params.k} \
            --n_pcs {params.n_pcs} \
            --harmony_vars {params.harmony_vars} \
            {params.label_cols_flag}

        echo "Cell line -> embryo mapping complete"
        """

##################################################################
# Integration via frozen-atlas projection (now the default path)
##################################################################
# Produces the canonical results/integrated/integrated.h5ad -- replaces the
# old Harmony/BBKNN re-clustering (formerly rule integrate_harmony, removed
# from this Snakefile -- see git history if you need it back) with
# integrate_via_atlas_projection.py: every cell -- embryo AND primary cell
# line -- is projected onto the SAME frozen Flysta3D-v2 atlas embedding in
# one shot, instead of Harmony re-clustering
# results/embryo_annotated + results/celllines_mapped_to_embryo. See that
# script's docstring for why cross-sample harmonisation isn't needed for
# this (no joint fit for a dataset axis to leak into) and what this object
# is/isn't good for (cell-type identity and composition, not
# titer-sensitive fine structure -- see analyze_titer_by_annotation.py and
# rule titer_by_annotation_atlas below for that). Takes ALL filtered
# samples directly -- no longer depends on rule annotate_with_atlas or
# rule map_celllines_to_embryo (their outputs are unused by this rule; both
# rules still exist and still run standalone if requested, but are no
# longer on the critical path to results/integrated/integrated.h5ad).
#
# The output object is schema-compatible with
# rule embryo_to_cellline_trajectory's expectations (is_embryo,
# cell_type_<label>, full-gene .raw, obsm['X_umap']) -- see
# integrate_via_atlas_projection.py's docstring for the two columns it
# deliberately does NOT provide (leiden, phase) and why those analysis
# sections skip gracefully rather than error.
rule integrate:
    input:
        files              = expand("results/filtered_h5ad/{sample_id}.h5ad", sample_id=SAMPLE_IDS),
        atlas              = config.get("flysta3d_atlas", "resources/wcoembed_whole_embeding_downsampled_modified.h5ad"),
        flybase_annotation = config["flybase_annotation"],
        orthologs          = config["ortholog_map"],
    output:
        integrated = "results/integrated/integrated.h5ad"
    params:
        script   = config.get("integrate_atlas_script",
                       "snakemake_scripts/method_comparison/integrate_via_atlas_projection.py"),
        fig_dir  = "results/integrated/figures",
        k        = config.get("atlas_k", 30),
        n_pcs    = config.get("atlas_n_pcs", 30),
        label_cols_flag = (
            "--label_cols " + " ".join(config["atlas_label_cols"])
            if config.get("atlas_label_cols") else ""
        ),
        subsample_flag = (
            f"--subsample_ref {config['atlas_subsample_ref']}"
            if config.get("atlas_subsample_ref") else ""
        ),
        condition_sample_type = CONDITION_SAMPLE_TYPE_PATH,
    log:
        "logs/integrate/integrate.log"
    threads:
        config.get("integrate_atlas_threads", 16)
    resources:
        slurm_partition = config.get("integrate_atlas_partition", "medium"),
        mem_mb          = config.get("integrate_atlas_mem", 500000),
        slurm_time      = config.get("integrate_atlas_time", "12:00:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Starting atlas-projected integration (all samples)"

        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --atlas {input.atlas} \
            --query {input.files} \
            --out_path {output.integrated} \
            --fig_dir {params.fig_dir} \
            --flybase_annotation {input.flybase_annotation} \
            --ortholog_map {input.orthologs} \
            --k {params.k} \
            --n_pcs {params.n_pcs} \
            --condition_sample_type {params.condition_sample_type} \
            {params.label_cols_flag} \
            {params.subsample_flag}

        echo "Atlas-projected integration complete"
        """

##################################################################
# Wolbachia titer by transferred cell-type annotation (atlas path)
##################################################################
# Runs analyze_titer_by_annotation.py on rule integrate's output (atlas-
# projected integration) -- titer/infection-rate stats and plots grouped by
# atlas_<label> instead of a Leiden cluster (see that script's docstring for
# why: Leiden on this embedding would cluster the ATLAS's own developmental
# biology, not your titer axis). --groupby defaults to the first entry of
# atlas_label_cols (atlas_<that column>); override with titer_groupby in
# config.yaml if you want a different transferred column (e.g. atlas_tissue).
rule titer_by_annotation_atlas:
    input:
        h5ad = rules.integrate.output.integrated
    output:
        flag = touch("results/integrated/figures_atlas/.titer_by_annotation.done")
    params:
        script      = config.get("titer_by_annotation_script",
                          "snakemake_scripts/method_comparison/analyze_titer_by_annotation.py"),
        fig_dir     = "results/integrated/figures_atlas",
        groupby     = config.get("titer_groupby",
                          f"atlas_{config['atlas_label_cols'][0]}" if config.get("atlas_label_cols") else "atlas_annotation"),
        condition_col = config.get("titer_condition_col", "condition"),
        origin_col  = config.get("titer_origin_col", "sample_type"),
        sample      = "wolbachia_infection",
    log:
        "logs/integrate/titer_by_annotation.log"
    threads:
        config.get("titer_by_annotation_threads", 4)
    resources:
        slurm_partition = config.get("titer_by_annotation_partition", "medium"),
        mem_mb          = config.get("titer_by_annotation_mem", 32000),
        slurm_time      = config.get("titer_by_annotation_time", "2:00:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Starting titer-by-annotation analysis"

        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --adata {input.h5ad} \
            --groupby {params.groupby} \
            --condition_col {params.condition_col} \
            --origin_col {params.origin_col} \
            --fig_dir {params.fig_dir} \
            --sample {params.sample}

        echo "Titer-by-annotation analysis complete"
        """

##################################################################
# Embryo -> cell line trajectory/identity analysis
##################################################################
# Exploratory pass over the integrated object: how do the cultured primary
# cell lines relate to the embryonic tissues they were derived from? Runs
# composition, diversity, confidence, pseudobulk correlation, marker-module
# scoring, cell cycle shift, Wolbachia effects, species, and cluster
# composition analyses -- see the module docstring in
# snakemake_scripts/analysis/embryo_to_cellline_trajectory.py for the full
# rationale behind each one. Every plot has a matching CSV of the
# underlying numbers.
rule embryo_to_cellline_trajectory:
    input:
        h5ad = rules.integrate.output.integrated
    output:
        flag = touch("results/trajectory_analysis/.done")
    params:
        script          = config.get("trajectory_script",
                               "snakemake_scripts/analysis/embryo_to_cellline_trajectory.py"),
        fig_dir         = "results/trajectory_analysis",
        conf_threshold  = config.get("trajectory_conf_threshold", 0.5),
        min_cells       = config.get("trajectory_min_cells", 20),
        top_n_markers   = config.get("trajectory_top_n_markers", 50),
        skip_markers    = "--skip_markers" if config.get("trajectory_skip_markers", False) else "",
        skip_pseudobulk = "--skip_pseudobulk" if config.get("trajectory_skip_pseudobulk", False) else "",
        skip_umap       = "--skip_umap" if config.get("trajectory_skip_umap", False) else "",
    log:
        "logs/trajectory_analysis/trajectory_analysis.log"
    threads:
        config.get("trajectory_threads", 8)
    resources:
        slurm_partition = config.get("trajectory_partition", "medium"),
        mem_mb          = config.get("trajectory_mem", 64000),
        slurm_time      = config.get("trajectory_time", "4:00:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Starting embryo -> cell line trajectory analysis"

        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --input {input.h5ad} \
            --fig_dir {params.fig_dir} \
            --conf_threshold {params.conf_threshold} \
            --min_cells {params.min_cells} \
            --top_n_markers {params.top_n_markers} \
            {params.skip_markers} \
            {params.skip_pseudobulk} \
            {params.skip_umap}

        echo "Trajectory analysis complete"
        """
        
# Count reads aligning to Wolbachia 16S rRNA (GQX67_05945) vs total reads per sample
rule count_16s_reads:
    input:
        r2             = lambda wildcards: get_fastq_files(wildcards.sample_id)[1],
        ref            = lambda wildcards: get_symbiont_fasta(wildcards.sample_id),
        bwa_index_flag = lambda wildcards: get_bwa_index_flag(wildcards.sample_id)
    output:
        counts = "results/rRNA_analysis/read_counts/{sample_id}/{gene}_read_counts.txt"
    wildcard_constraints:
        # Only allow sample_ids that currently exist in the dataframe
        sample_id = "|".join(SAMPLE_IDS)
    params:
        # Resolved per-sample from that sample's own Wolbachia strain (see
        # get_symbiont_16s_region) -- e.g. wMel samples get wMel's
        # GQX67_RS05935 region, wRi_M23 samples get M23_00679's, etc.
        # {wildcards.gene} is just the generic path label from
        # config['target_genes']; the region (and the gene label actually
        # written into the output row below) is the real per-strain locus.
        region = lambda wildcards: get_symbiont_16s_region(wildcards.sample_id)
    log:
        "logs/count_16s/{sample_id}_{gene}.log"
    threads: 16
    resources:
        slurm_partition = "medium",
        mem_mb          = 32000,
        slurm_time      = "2:00:00"
    shell:
        """
        exec > {log} 2>&1
        echo "Counting 16S vs total reads for {wildcards.sample_id} (region: {params.region})"

        source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
        conda activate {SRA_TOOLS_ENV}

        mkdir -p results/rRNA_analysis/read_counts/{wildcards.sample_id}

        SORTED=results/rRNA_analysis/read_counts/{wildcards.sample_id}/all_aligned.bam

        # Align R2 to this sample's own Wolbachia strain genome (not a
        # wMel-only reference -- see get_symbiont_fasta), sort, index
        bwa mem -t {threads} {input.ref} {input.r2} | \
            samtools view -Sb | \
            samtools sort -@ {threads} -o $SORTED
        samtools index $SORTED

        # {params.region} is "<locus_tag>::<seqid>:<start>-<end>" for this
        # sample's own strain -- use the locus_tag piece as the reported
        # gene label since it's the gene actually analysed for this sample,
        # not the generic {wildcards.gene} path label
        GENE_LABEL="{params.region}"
        GENE_LABEL="${{GENE_LABEL%%::*}}"

        # Reads mapped to the 16S region (primary alignments only, excludes
        # unmapped/secondary/supplementary via -F 0x904)
        SIXTEEN_S=$(samtools view -c -F 0x904 $SORTED "{params.region}")

        # Total reads in the input (all records in the BAM)
        TOTAL=$(samtools view -c $SORTED)

        # Mapped reads only (optional denominator if you prefer mapped over total)
        MAPPED=$(samtools view -c -F 0x4 $SORTED)

        # Write tab-delimited summary
        printf "sample\\tgene\\tregion\\tsixteenS_reads\\ttotal_reads\\tmapped_reads\\n" > {output.counts}
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \
            "{wildcards.sample_id}" "$GENE_LABEL" "{params.region}" \
            "$SIXTEEN_S" "$TOTAL" "$MAPPED" >> {output.counts}

        echo "Done. 16S=$SIXTEEN_S  total=$TOTAL  mapped=$MAPPED"
        """
