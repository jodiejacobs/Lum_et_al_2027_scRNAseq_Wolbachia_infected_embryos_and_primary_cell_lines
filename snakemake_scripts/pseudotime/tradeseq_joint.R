#!/usr/bin/env Rscript
# tradeseq_joint.R
# ================
# Step 4b (rule pseudotime_joint_tradeseq). One tradeSeq model for a pair of
# trajectories A and B (two lineages of one species, or one per species) on
# 1:1 ortholog counts (Dmel gene space), with the trajectory (cells.csv
# 'group') as the condition. Each cell keeps its own trajectory's SCEPTIC
# pseudotime. All trajectories share the 0-2 stage scale, so no joint
# embedding is needed.
#
# conditionTest asks whether a gene's smoother differs between A and B.
# That covers expression-level offsets as well as shape. Level offsets can
# come from annotation, mapping, or batch differences, so this script also
# reports a shape-only measure: the Pearson r between the two smoothed
# curves on the shared pseudotime grid.
#
# Usage:
#   Rscript tradeseq_joint.R --indir results/pseudotime/pairs/<A>__vs__<B>/tradeseq_inputs \
#       --outdir results/pseudotime/pairs/<A>__vs__<B>/tradeseq --nknots 6 --nworkers 16
#
# Outputs: tradeseq_joint_sce.rds, tradeseq_condition_test.csv,
#   tradeseq_joint_association.csv, tradeseq_joint_smooth_tidy.csv.gz,
#   tradeseq_shape_similarity.csv, shape_similarity_hist.pdf,
#   condition_test_top9_smoothers.pdf

suppressPackageStartupMessages({
    library(tradeSeq)
    library(SingleCellExperiment)
    library(BiocParallel)
    library(Matrix)
    library(ggplot2)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) == 0) {
        if (is.null(default)) stop(sprintf("Required argument %s not provided", flag))
        return(default)
    }
    args[idx + 1]
}
in_dir    <- parse_arg("--indir")
out_dir   <- parse_arg("--outdir")
n_knots   <- as.integer(parse_arg("--nknots", "6"))
n_workers <- as.integer(parse_arg("--nworkers", "8"))
n_points  <- as.integer(parse_arg("--npoints", "50"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("[1/4] Loading joint inputs\n")
counts <- as(readMM(file.path(in_dir, "counts_genesXcells.mtx")), "CsparseMatrix")
genes  <- read.delim(file.path(in_dir, "genes.tsv"), stringsAsFactors = FALSE)
cells  <- read.csv(file.path(in_dir, "cells.csv"), row.names = 1, stringsAsFactors = FALSE)
rownames(counts) <- genes$gene
colnames(counts) <- rownames(cells)
label_of <- setNames(genes$label, genes$gene)
# trajectory A is whichever group appears first in cells.csv
grp_levels <- unique(cells$group)
if (length(grp_levels) != 2) stop("cells.csv must contain exactly 2 groups")
species  <- factor(cells$group, levels = grp_levels)
name_a <- grp_levels[1]; name_b <- grp_levels[2]
cat(sprintf("  %d genes x %d cells\n", nrow(counts), ncol(counts)))
print(table(species, cells$sample_type))

pt_mat <- matrix(cells$sceptic_pseudotime, ncol = 1,
                 dimnames = list(rownames(cells), "pseudotime"))
wt_mat <- matrix(1, nrow = nrow(cells), ncol = 1, dimnames = list(rownames(cells), "w1"))

cat(sprintf("\n[2/4] fitGAM with conditions = trajectory (%d knots, %d workers)\n",
            n_knots, n_workers))
set.seed(42)
sce <- fitGAM(counts = counts, pseudotime = pt_mat, cellWeights = wt_mat,
              conditions = species, nknots = n_knots, parallel = TRUE,
              BPPARAM = MulticoreParam(workers = n_workers), verbose = FALSE)
saveRDS(sce, file.path(out_dir, "tradeseq_joint_sce.rds"))

finish <- function(df, name) {
    df <- as.data.frame(df)
    df$gene  <- rownames(df)
    df$label <- label_of[df$gene]
    df$padj  <- p.adjust(df$pvalue, method = "BH")
    df$sig   <- !is.na(df$padj) & df$padj < 0.05
    df <- df[order(df$waldStat, decreasing = TRUE), ]
    write.csv(df, file.path(out_dir, name), row.names = FALSE)
    cat(sprintf("  %s: %d / %d genes padj < 0.05\n", name, sum(df$sig), nrow(df)))
    df
}

cat("\n[3/4] Tests\n")
cond <- finish(conditionTest(sce, global = TRUE, pairwise = FALSE),
               "tradeseq_condition_test.csv")
assoc <- tryCatch(
    finish(associationTest(sce, global = TRUE, lineages = TRUE),
           "tradeseq_joint_association.csv"),
    error = function(e) { cat("  associationTest failed:", conditionMessage(e), "\n"); NULL })

cat("\n[4/4] Shape similarity of the two smoothers\n")
sm <- predictSmooth(sce, gene = rownames(counts), nPoints = n_points, tidy = TRUE)
gz <- gzfile(file.path(out_dir, "tradeseq_joint_smooth_tidy.csv.gz"), "w")
write.csv(sm, gz, row.names = FALSE); close(gz)

cat("  predictSmooth condition values:", paste(unique(sm$condition), collapse = ", "), "\n")
sm$ly <- log1p(sm$yhat)
wide  <- split(sm, sm$gene)
shape <- do.call(rbind, lapply(wide, function(d) {
    is_mel <- as.character(d$condition) == name_a
    is_sim <- as.character(d$condition) == name_b
    a <- d$ly[is_mel][order(d$time[is_mel])]
    b <- d$ly[is_sim][order(d$time[is_sim])]
    r <- if (sd(a) > 0 && sd(b) > 0) cor(a, b) else NA_real_
    data.frame(gene = d$gene[1], shape_r = r,
               mean_log_ratio_a_vs_b = mean(a - b),
               range_a = diff(range(a)), range_b = diff(range(b)))
}))
shape$label <- label_of[shape$gene]
ct <- cond[, c("gene", "waldStat", "padj", "sig")]
names(ct) <- c("gene", "condition_waldStat", "condition_padj", "condition_sig")
shape <- merge(shape, ct, by = "gene")
shape <- shape[order(shape$shape_r), ]
write.csv(shape, file.path(out_dir, "tradeseq_shape_similarity.csv"), row.names = FALSE)

p <- ggplot(shape, aes(shape_r, fill = condition_sig)) +
    geom_histogram(bins = 50, alpha = 0.8, position = "identity") +
    scale_fill_manual(values = c("FALSE" = "#95a5a6", "TRUE" = "#c0392b"),
                      labels = c("conditionTest ns", "conditionTest padj < 0.05")) +
    labs(x = sprintf("Pearson r, %s vs %s smoothed curve (log1p)", name_a, name_b), y = "Genes", fill = NULL,
         title = "Shape similarity of pseudotime dynamics") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out_dir, "shape_similarity_hist.pdf"), p, width = 7, height = 5)

top9 <- head(cond$gene[cond$sig], 9)
if (length(top9) > 0) {
    plots <- lapply(top9, function(g)
        plotSmoothers(sce, counts = counts, gene = g, alpha = 0.3) +
            labs(title = label_of[[g]], x = "SCEPTIC pseudotime") +
            theme_bw(base_size = 9))
    ggsave(file.path(out_dir, "condition_test_top9_smoothers.pdf"),
           wrap_plots(plots, ncol = 3), width = min(3, length(plots)) * 5,
           height = ceiling(length(plots) / 3) * 5)
}
cat("\n=== joint tradeSeq complete ===\n")
