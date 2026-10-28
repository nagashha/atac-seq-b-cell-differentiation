# ---------------------------------------------------------------------------
# Fit the model, stabilise the variance, and run exploratory QC.
#
# In:  se_filtered.rds
# Out: ddsA.rds, vsdA.rds, qc/*.png
# ---------------------------------------------------------------------------
suppressMessages({
  library(DESeq2); library(SummarizedExperiment)
  library(vsn); library(hexbin); library(pheatmap)
  library(RColorBrewer); library(ggplot2); library(matrixStats)
})

se <- readRDS("se_filtered.rds")
dir.create("qc", showWarnings = FALSE)

# Design A, the main analysis: ~ cell_type on all seven samples.
# A paired design (~ Donor + cell_type) is not used as the main model because
# only 2 of the 5 donors contributed both cell types, which leaves too few
# residual degrees of freedom. 04_validation_sensitivity.R fits the paired
# model on those two donors as a sensitivity check.
ddsA <- DESeq(DESeqDataSet(se, design = ~ cell_type))
cat("size factors:\n"); print(round(sizeFactors(ddsA), 3))

# blind = FALSE lets the transformation use the experimental design.
vsdA <- vst(ddsA, blind = FALSE)

# --- variance stabilisation ------------------------------------------------
# The red running-mean line should flatten after the transformation. Before
# it, variance rises with the mean, which would let the few highest-count
# regions dominate every distance and every principal component.
png("qc/meanSd_before_vst.png", 900, 700, res = 120)
print(meanSdPlot(log2(counts(ddsA, normalized = TRUE) + 1), plot = FALSE)$gg +
      ggtitle("meanSdPlot: before VST"))
dev.off()

png("qc/meanSd_after_vst.png", 900, 700, res = 120)
print(meanSdPlot(assay(vsdA), plot = FALSE)$gg + ggtitle("meanSdPlot: after VST"))
dev.off()

# --- principal components --------------------------------------------------
# DESeq2's plotPCA silently uses only the ntop = 500 most variable regions, so
# the percentages it prints describe variance within that subset and not
# within the whole matrix. Both are computed here, and both are reported,
# because the two numbers answer different questions.
mat <- assay(vsdA)

pc_all  <- prcomp(t(mat))
var_all <- 100 * pc_all$sdev^2 / sum(pc_all$sdev^2)

top500  <- head(order(rowVars(mat), decreasing = TRUE), 500)
pc_500  <- prcomp(t(mat[top500, ]))
var_500 <- 100 * pc_500$sdev^2 / sum(pc_500$sdev^2)

cat(sprintf("PC1: %.1f%% of variance across all %d regions\n", var_all[1], nrow(mat)))
cat(sprintf("PC1: %.1f%% of variance across the 500 most variable regions\n", var_500[1]))

pcaData <- data.frame(PC1 = pc_500$x[, 1], PC2 = pc_500$x[, 2],
                      cell_type = vsdA$cell_type, Donor = vsdA$Donor)
png("qc/pca.png", 1000, 800, res = 120)
print(ggplot(pcaData, aes(PC1, PC2, colour = cell_type, shape = Donor)) +
        geom_point(size = 5) +
        xlab(sprintf("PC1: %.0f%% variance", var_500[1])) +
        ylab(sprintf("PC2: %.0f%% variance", var_500[2])) +
        labs(caption = sprintf(
          "500 most variable regions (DESeq2 plotPCA default); across all %d regions PC1 = %.0f%%",
          nrow(mat), var_all[1])) +
        theme_bw(base_size = 14))
dev.off()

# --- sample-to-sample distances -------------------------------------------
# An unsupervised view: samples should group by cell type, not by donor.
sampleDists <- dist(t(mat))
dmat <- as.matrix(sampleDists)
labs <- paste(vsdA$cell_type, vsdA$Donor, sep = " | ")
rownames(dmat) <- labs; colnames(dmat) <- labs
png("qc/sample_distance_heatmap.png", 900, 800, res = 120)
pheatmap(dmat,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()

saveRDS(ddsA, "ddsA.rds")
saveRDS(vsdA, "vsdA.rds")
