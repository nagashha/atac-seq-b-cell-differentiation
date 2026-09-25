# Biological positive control and paired-design sensitivity analysis.
suppressMessages({ library(DESeq2); library(apeglm) })

se   <- readRDS("se_filtered.rds")
resS <- readRDS("resShrink.rds")
df   <- as.data.frame(resS); df$PeakID <- rownames(df)

# --- 1. Positive control: PRDM1 (BLIMP1), the plasma-cell master regulator.
# "The locus" here is the consensus interval peak_123058-peak_123067.
prdm1 <- paste0("peak_", 123058:123067)
v <- df[rownames(df) %in% prdm1, c("baseMean", "log2FoldChange", "padj")]
v$DAR <- ifelse(!is.na(v$padj) & v$padj < 0.05 & abs(v$log2FoldChange) >= 1, "DAR", "-")
print(round(v[, 1:3], 3)); print(table(v$DAR))

# --- 2. Sensitivity analysis: true paired design on the matched donors only.
matched <- colData(se)$Donor %in% c("D1001", "D1002")
seB <- se[, matched]
colData(seB)$Donor <- droplevels(colData(seB)$Donor)
print(table(colData(seB)$Donor, colData(seB)$cell_type))

ddsB <- DESeqDataSet(seB, design = ~ Donor + cell_type)
ddsB <- DESeq(ddsB)
resB <- lfcShrink(ddsB, coef = "cell_type_Plasmablast_vs_Bulk_B", type = "apeglm")

dfB <- as.data.frame(resB); dfB$PeakID <- rownames(dfB)
dfB <- dfB[!is.na(dfB$padj), ]
darsB <- dfB[dfB$padj < 0.05 & abs(dfB$log2FoldChange) >= 1, ]

# --- 3. Compare. If the unblocked design were detecting donor artefacts,
# controlling for donor should dissolve or reverse a substantial fraction.
darsA <- df[!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) >= 1, ]
ov <- intersect(rownames(darsA), rownames(darsB))
cat("Design A DARs:", nrow(darsA), "  Design B DARs:", nrow(darsB), "\n")
cat("Shared:", length(ov),
    sprintf(" (%.1f%% of B)\n", 100 * length(ov) / nrow(darsB)))

common <- intersect(rownames(darsA), rownames(dfB))
agree  <- sign(darsA[common, "log2FoldChange"]) == sign(dfB[common, "log2FoldChange"])
cat(sprintf("log2FC sign agreement: %.1f%% (n=%d)\n", 100 * mean(agree), length(common)))

write.table(darsB[order(darsB$padj), ], "DARs_designB_paired.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
