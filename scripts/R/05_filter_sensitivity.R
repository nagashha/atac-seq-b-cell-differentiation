# ---------------------------------------------------------------------------
# Filter-threshold sensitivity analysis: N = 3 versus N = 4.
#
# The project slides specify a low-count filter of >= 10 counts in at least
# N = 3 samples (the size of the smaller group); the assignment PDF specifies
# N = 4. The main analysis uses N = 3. Rather than argue for that choice, this
# script runs both thresholds to the end and reports what actually differs.
#
# Everything downstream of the filter is identical between the two arms, so
# any difference in the result is attributable to the filter alone.
#
# In:  counts_clean.txt, clinical_table.txt
# Out: DARs_N4_threshold.tsv, filter_comparison.rds
# ---------------------------------------------------------------------------
suppressMessages({
  library(DESeq2); library(SummarizedExperiment); library(matrixStats)
})

SHRINK <- Sys.getenv("SHRINK_TYPE", "apeglm")

counts  <- as.matrix(read.delim("counts_clean.txt", row.names = 1,
                                check.names = FALSE))
coldata <- read.delim("clinical_table.txt", row.names = 1)
stopifnot(all(colnames(counts) == rownames(coldata)))

coldata$cell_type <- factor(coldata$cell_type, levels = c("Bulk_B", "Plasmablast"))
coldata$Donor     <- factor(coldata$Donor)
se_all <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

cat("consensus regions:", nrow(se_all), " samples:", ncol(se_all), "\n\n")

# --- one arm: filter at N, fit ~ cell_type, shrink, call DARs ---------------
run_arm <- function(se, N, label) {
  keep <- rowSums(assay(se, "counts") >= 10) >= N
  se_f <- se[keep, ]
  cat("=== ", label, " (N = ", N, ") ===\n", sep = "")
  cat(sprintf("regions retained : %d (%.1f%% of %d removed)\n",
              nrow(se_f), 100 * (nrow(se) - nrow(se_f)) / nrow(se), nrow(se)))

  dds <- DESeq(DESeqDataSet(se_f, design = ~ cell_type), quiet = TRUE)

  # alpha = 0.05 so that independent filtering is optimised for the same
  # threshold used to call DARs, rather than DESeq2's default 0.1.
  results(dds, contrast = c("cell_type", "Plasmablast", "Bulk_B"), alpha = 0.05)
  resS <- lfcShrink(dds, coef = "cell_type_Plasmablast_vs_Bulk_B",
                    type = SHRINK, quiet = TRUE)

  df <- as.data.frame(resS); df$PeakID <- rownames(df)
  df <- df[!is.na(df$padj), ]
  dars <- df[df$padj < 0.05 & abs(df$log2FoldChange) >= 1, ]

  cat("regions tested   :", nrow(df), "\n")
  cat("DARs             :", nrow(dars),
      " (up", sum(dars$log2FoldChange >= 1),
      "/ down", sum(dars$log2FoldChange <= -1), ")\n")

  # PCA on the full retained matrix as well as DESeq2's top-500 default.
  vsd <- vst(dds, blind = FALSE); m <- assay(vsd)
  v_all <- prcomp(t(m))$sdev^2; v_all <- 100 * v_all / sum(v_all)
  top   <- head(order(rowVars(m), decreasing = TRUE), 500)
  v_500 <- prcomp(t(m[top, ]))$sdev^2; v_500 <- 100 * v_500 / sum(v_500)
  cat(sprintf("PC1              : %.1f%% (all %d) | %.1f%% (top 500)\n\n",
              v_all[1], nrow(m), v_500[1]))

  list(df = df, dars = dars, n_retained = nrow(se_f),
       pc1_all = v_all[1], pc1_500 = v_500[1])
}

a3 <- run_arm(se_all, 3, "main analysis")
a4 <- run_arm(se_all, 4, "assignment-PDF threshold")

# --- concordance -----------------------------------------------------------
cat("=== N = 3 versus N = 4 ===\n")
d3 <- rownames(a3$dars); d4 <- rownames(a4$dars)
ov <- intersect(d3, d4)
cat(sprintf("DARs: N=3 %d, N=4 %d, shared %d (%.1f%% of N=3, %.1f%% of N=4)\n",
            length(d3), length(d4), length(ov),
            100 * length(ov) / length(d3), 100 * length(ov) / length(d4)))

common <- intersect(rownames(a3$df), rownames(a4$df))
l3 <- a3$df[common, "log2FoldChange"]; l4 <- a4$df[common, "log2FoldChange"]
cat(sprintf("regions tested in both : %d\n", length(common)))
cat(sprintf("log2FC Pearson r %.4f | Spearman rho %.4f | sign agreement %.2f%%\n",
            cor(l3, l4), cor(l3, l4, method = "spearman"),
            100 * mean(sign(l3) == sign(l4))))
cat(sprintf("|difference| median %.4f, max %.3f\n",
            median(abs(l3 - l4)), max(abs(l3 - l4))))
cat(sprintf("DARs unique to N=3 : %d | unique to N=4 : %d\n",
            length(setdiff(d3, d4)), length(setdiff(d4, d3))))

write.table(a4$dars[order(a4$dars$padj), ], "DARs_N4_threshold.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(a3 = a3, a4 = a4), "filter_comparison.rds")
