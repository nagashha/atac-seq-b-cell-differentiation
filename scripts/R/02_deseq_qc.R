# Fit the model, stabilise variance, and run exploratory QC.
suppressMessages({
  library(DESeq2); library(SummarizedExperiment)
  library(vsn); library(hexbin); library(pheatmap)
  library(RColorBrewer); library(ggplot2)
})

se <- readRDS("se_filtered.rds")
dir.create("qc", showWarnings = FALSE)

# Design A: ~ cell_type on all 7 samples (main analysis).
# A paired design (~ Donor + cell_type) is not used here because only 2 of the
# 5 donors contributed both cell types; see 04_validation_sensitivity.R.
ddsA <- DESeqDataSet(se, design = ~ cell_type)
ddsA <- DESeq(ddsA)
print(round(sizeFactors(ddsA), 3))

# blind = FALSE lets the transformation account for the experimental design
vsdA <- vst(ddsA, blind = FALSE)

# meanSdPlot: the red running-mean line should flatten after transformation
png("qc/meanSd_before_vst.png", 900, 700, res = 120)
print(meanSdPlot(log2(counts(ddsA, normalized = TRUE) + 1), plot = FALSE)$gg +
      ggtitle("meanSdPlot: before VST"))
dev.off()

png("qc/meanSd_after_vst.png", 900, 700, res = 120)
print(meanSdPlot(assay(vsdA), plot = FALSE)$gg + ggtitle("meanSdPlot: after VST"))
dev.off()

# PCA. NOTE: plotPCA uses ntop = 500 by default, so the reported percentages
# describe variance among the 500 most variable peaks, not the whole matrix.
pcaData <- plotPCA(vsdA, intgroup = c("cell_type", "Donor"), returnData = TRUE)
pv <- round(100 * attr(pcaData, "percentVar"))
png("qc/pca.png", 1000, 800, res = 120)
print(ggplot(pcaData, aes(PC1, PC2, color = cell_type, shape = Donor)) +
        geom_point(size = 5) +
        xlab(paste0("PC1: ", pv[1], "% variance")) +
        ylab(paste0("PC2: ", pv[2], "% variance")) +
        theme_bw(base_size = 14))
dev.off()
cat("PC1 =", pv[1], "%  PC2 =", pv[2], "% (of the top 500 variable peaks)\n")

# Sample-to-sample distances
sampleDists <- dist(t(assay(vsdA)))
mat  <- as.matrix(sampleDists)
labs <- paste(vsdA$cell_type, vsdA$Donor, sep = " | ")
rownames(mat) <- labs; colnames(mat) <- labs
png("qc/sample_distance_heatmap.png", 900, 800, res = 120)
pheatmap(mat, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()

saveRDS(ddsA, "ddsA.rds"); saveRDS(vsdA, "vsdA.rds")
