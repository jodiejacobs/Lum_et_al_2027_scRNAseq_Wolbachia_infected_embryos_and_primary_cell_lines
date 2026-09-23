#!/usr/bin/env Rscript
# tradeseq_species.R
# ==================
# Step 3 (rule pseudotime_tradeseq). tradeSeq GAMs along SCEPTIC pseudotime
# (0 = embryo, 1 = primary cells, 2 = cell line) for one species.
# Adapted from Jacobs et al. 2026 run_tradeseq.R. Changes: sparse MatrixMarket
# input, readable gene labels, and stage-transition tests.
#
# Tests
#   associationTest                       : any change along pseudotime
#   startVsEndTest (global)               : embryo end vs cell-line end
#   startVsEndTest(pseudotimeValues=...)  : embryo -> primary  (min(pt) vs 1)
#                                           primary -> cell line (1 vs max(pt))
#
# Usage:
#   Rscript tradeseq_species.R --indir results/pseudotime/Dmel/tradeseq_inputs \
#       --outdir results/pseudotime/Dmel/tradeseq --nknots 6 --nworkers 16
#
# Outputs (in --outdir):
#   tradeseq_sce.rds, tradeseq_association.csv, tradeseq_startvsend.csv,
#   tradeseq_embryo_to_primary.csv, tradeseq_primary_to_cellline.csv,
#   tradeseq_smooth_predictions.csv, tradeseq_heatmap_genes_ordered.csv,
#   tradeseq_association_volcano.pdf, tradeseq_smooth_heatmap.pdf,
#   tradeseq_smooth_curves_top9.pdf, tradeseq_<transition>_barplot.pdf

suppressPackageStartupMessages({
    library(tradeSeq)
    library(SingleCellExperiment)
    library(BiocParallel)
    library(Matrix)
    library(ggplot2)
    library(patchwork)
})
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

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
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Load
# ─────────────────────────────────────────────────────────────────────────────
cat("[1/5] Loading inputs from", in_dir, "\n")
counts <- as(readMM(file.path(in_dir, "counts_genesXcells.mtx")), "CsparseMatrix")
genes  <- read.delim(file.path(in_dir, "genes.tsv"), stringsAsFactors = FALSE)
cells  <- read.csv(file.path(in_dir, "cells.csv"), row.names = 1, stringsAsFactors = FALSE)
rownames(counts) <- genes$gene
colnames(counts) <- rownames(cells)
label_of <- setNames(ifelse(is.na(genes$label) | genes$label == "", genes$gene, genes$label),
                     genes$gene)
cat(sprintf("  %d genes x %d cells\n", nrow(counts), ncol(counts)))
print(table(cells$sample_type))

pt  <- cells$sceptic_pseudotime
names(pt) <- rownames(cells)
pt_mat <- matrix(pt, ncol = 1, dimnames = list(names(pt), "pseudotime"))
wt_mat <- matrix(1, nrow = length(pt), ncol = 1, dimnames = list(names(pt), "w1"))

# ─────────────────────────────────────────────────────────────────────────────
# Fit
# ─────────────────────────────────────────────────────────────────────────────
cat(sprintf("\n[2/5] fitGAM (%d genes, %d knots, %d workers)\n",
            nrow(counts), n_knots, n_workers))
set.seed(42)
sce <- fitGAM(counts = counts, pseudotime = pt_mat, cellWeights = wt_mat,
              nknots = n_knots, parallel = TRUE,
              BPPARAM = MulticoreParam(workers = n_workers), verbose = FALSE)
saveRDS(sce, file.path(out_dir, "tradeseq_sce.rds"))

# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────
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

cat("\n[3/5] Tests\n")
assoc <- finish(associationTest(sce, global = TRUE, lineages = FALSE),
                "tradeseq_association.csv")
sve   <- finish(startVsEndTest(sce, global = TRUE, lineages = FALSE),
                "tradeseq_startvsend.csv")

t_min <- min(pt); t_max <- max(pt)
transitions <- list(
    embryo_to_primary   = c(t_min, 1),
    primary_to_cellline = c(1, t_max)
)
trans_res <- list()
for (nm in names(transitions)) {
    tv <- transitions[[nm]]
    if (tv[1] >= tv[2]) {
        cat(sprintf("  Skipping %s: pseudotime range does not span %s\n", nm, paste(tv, collapse = "-")))
        next
    }
    trans_res[[nm]] <- finish(
        startVsEndTest(sce, global = TRUE, lineages = FALSE, pseudotimeValues = tv),
        sprintf("tradeseq_%s.csv", nm))
}

# ─────────────────────────────────────────────────────────────────────────────
# Smooth predictions + plots (same figures as 2026 run_tradeseq.R)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[4/5] Smooth predictions\n")
top_genes <- head(assoc$gene[assoc$sig], 300)
if (length(top_genes) < 5) top_genes <- head(assoc$gene, 50)
yhat <- predictSmooth(sce, gene = top_genes, nPoints = 100, tidy = FALSE)
write.csv(yhat, file.path(out_dir, "tradeseq_smooth_predictions.csv"))

cat("\n[5/5] Plots\n")
save_plot <- function(p, name, w, h) {
    ggsave(file.path(out_dir, name), p, width = w, height = h, limitsize = FALSE)
    cat("  Saved:", name, "\n")
}

# Volcano
assoc$neg_log10p <- -log10(assoc$padj + 1e-300)
p <- ggplot(assoc, aes(waldStat, neg_log10p, colour = sig)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_colour_manual(values = c("FALSE" = "#bdc3c7", "TRUE" = "#c0392b"),
                        labels = c("ns", "padj < 0.05")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    labs(x = "Wald statistic", y = expression(-log[10](p[adj])), colour = NULL,
         title = "associationTest: genes dynamic along embryo -> cell line pseudotime") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
if (has_ggrepel)
    p <- p + ggrepel::geom_text_repel(data = head(assoc[assoc$sig, ], 20),
                                      aes(label = label), size = 2.5, colour = "black")
save_plot(p, "tradeseq_association_volcano.pdf", 8, 6)

# Heatmap ordered by peak pseudotime
hm_genes <- intersect(head(assoc$gene[assoc$sig], 100), rownames(yhat))
if (length(hm_genes) >= 5) {
    mat_z <- t(scale(t(as.matrix(yhat[hm_genes, , drop = FALSE]))))
    mat_z[is.nan(mat_z)] <- 0
    peak  <- apply(mat_z, 1, which.max)
    ord   <- order(peak)
    mat_z <- mat_z[ord, , drop = FALSE]; peak <- peak[ord]
    pt_grid <- seq(t_min, t_max, length.out = ncol(mat_z))
    df_long <- data.frame(
        gene = factor(rep(label_of[rownames(mat_z)], times = ncol(mat_z)),
                      levels = rev(unique(label_of[rownames(mat_z)]))),
        pt   = rep(pt_grid, each = nrow(mat_z)),
        expr = as.vector(mat_z))
    p <- ggplot(df_long, aes(pt, gene, fill = expr)) + geom_tile() +
        scale_fill_viridis_c(name = "Scaled\nexpr") +
        geom_vline(xintercept = 1, linetype = "dashed", colour = "white") +
        scale_x_continuous(expand = c(0, 0), name = "SCEPTIC pseudotime (0 embryo, 1 primary, 2 cell line)") +
        theme_minimal(base_size = 8) +
        theme(axis.text.y = element_text(size = ifelse(length(hm_genes) > 60, 4, 6)),
              panel.grid = element_blank()) +
        labs(title = sprintf("Top %d dynamic genes ordered by peak pseudotime", length(hm_genes)), y = NULL)
    save_plot(p, "tradeseq_smooth_heatmap.pdf", 10, max(5, length(hm_genes) * 0.12 + 1))

    peak_pt <- pt_grid[peak]
    out <- assoc[match(rownames(mat_z), assoc$gene), c("gene", "label", "waldStat", "padj")]
    out$peak_pseudotime <- round(peak_pt, 3)
    out$peak_stage <- cut(peak_pt, breaks = c(-Inf, 0.5, 1.5, Inf),
                          labels = c("embryo", "primary_cells", "cell_culture"))
    write.csv(out, file.path(out_dir, "tradeseq_heatmap_genes_ordered.csv"), row.names = FALSE)
}

# Top 9 smoothers
top9 <- head(assoc$gene[assoc$sig], 9)
if (length(top9) == 0) top9 <- head(assoc$gene, 9)
plots <- lapply(top9, function(g)
    plotSmoothers(sce, counts = counts, gene = g) +
        labs(title = label_of[[g]], x = "SCEPTIC pseudotime", y = "log(expression + 1)") +
        theme_bw(base_size = 9) + theme(legend.position = "none"))
save_plot(wrap_plots(plots, ncol = 3), "tradeseq_smooth_curves_top9.pdf",
          min(3, length(plots)) * 5, ceiling(length(plots) / 3) * 5)

# Transition barplots (top 30 by Wald, coloured by direction)
for (nm in names(trans_res)) {
    d <- head(trans_res[[nm]][trans_res[[nm]]$sig, ], 30)
    if (nrow(d) == 0) next
    lfc <- grep("^logFC", colnames(d), value = TRUE)[1]
    d$direction <- ifelse(d[[lfc]] > 0, "Up", "Down")
    d$label <- factor(make.unique(d$label), levels = rev(make.unique(d$label)))
    p <- ggplot(d, aes(waldStat, label, fill = direction)) + geom_col(alpha = 0.8) +
        scale_fill_manual(values = c(Up = "#c0392b", Down = "#2980b9")) +
        labs(x = "Wald statistic", y = NULL, fill = NULL, title = gsub("_", " ", nm)) +
        theme_bw(base_size = 10) + theme(legend.position = "top")
    save_plot(p, sprintf("tradeseq_%s_barplot.pdf", nm), 7, max(4, nrow(d) * 0.25 + 2))
}

cat("\n=== tradeSeq complete ===\n")
