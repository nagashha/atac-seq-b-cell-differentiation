# Differential accessibility testing with LFC shrinkage.
suppressMessages({ library(DESeq2); library(apeglm); library(EnhancedVolcano) })

ddsA <- readRDS("ddsA.rds")

# Plasmablast vs Bulk_B (reference). Positive log2FC = more open in plasmablasts.
# alpha = 0.05 matches the DAR criterion, so independent filtering is optimised
# for the threshold actually used.
res <- results(ddsA, contrast = c("cell_type", "Plasmablast", "Bulk_B"), alpha = 0.05)
summary(res)
print(resultsNames(ddsA))

# Low-count peaks show inflated fold changes from sampling noise alone.
# apeglm shrinks those toward zero and leaves well-supported peaks intact.
resS <- lfcShrink(ddsA, coef = "cell_type_Plasmablast_vs_Bulk_B", type = "apeglm")

# DAR criteria: FDR < 0.05 AND |log2FC| >= 1
df <- as.data.frame(resS); df$PeakID <- rownames(df)
df <- df[!is.na(df$padj), ]
dars <- df[df$padj < 0.05 & abs(df$log2FoldChange) >= 1, ]
up   <- dars[dars$log2FoldChange >=  1, ]
down <- dars[dars$log2FoldChange <= -1, ]

cat("Peaks tested:", nrow(df), "\n")
cat("Total DARs:",   nrow(dars), "\n")
cat("  more open in plasmablasts:", nrow(up), "\n")
cat("  more open in bulk B      :", nrow(down), "\n")

png("qc/volcano.png", 1100, 950, res = 120)
print(EnhancedVolcano(df, lab = df$PeakID, x = "log2FoldChange", y = "padj",
      pCutoff = 0.05, FCcutoff = 1,
      xlab = bquote(~Log[2]~"fold change"),
      ylab = bquote(~-Log[10]~italic(FDR)),
      pointSize = 1.5, labSize = 3.0, colAlpha = 0.6))
dev.off()

# NOTE: plotMA colours by `alpha` only; it does not additionally apply the
# |log2FC| >= 1 threshold, so blue points are not identical to the DAR set.
png("qc/ma_plot.png", 950, 750, res = 120)
plotMA(resS, alpha = 0.05, ylim = c(-6, 6))
dev.off()

df <- df[order(df$padj), ]
write.table(df,   "DA_results_all.tsv",    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dars[order(dars$padj), ], "DARs_significant.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(resS, "resShrink.rds")
