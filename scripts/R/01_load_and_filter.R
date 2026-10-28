# ---------------------------------------------------------------------------
# Assemble a SummarizedExperiment from the count matrix and apply the
# low-count filter.
#
# Out: se_filtered.rds
# ---------------------------------------------------------------------------
suppressMessages({
  library(SummarizedExperiment)
  library(DESeq2)
})

counts  <- as.matrix(read.delim("counts_clean.txt", row.names = 1,
                                check.names = FALSE))
coldata <- read.delim("clinical_table.txt", row.names = 1)

# The single most important check in the whole analysis. If the count columns
# and the metadata rows are in different orders, DESeq2 compares the wrong
# groups, returns a plausible-looking result and never raises an error.
stopifnot(all(colnames(counts) == rownames(coldata)))
cat(">>> column order match VERIFIED\n")

# Bulk_B as the reference level, so a positive log2 fold change means
# "more accessible in plasmablasts".
coldata$cell_type <- factor(coldata$cell_type, levels = c("Bulk_B", "Plasmablast"))
coldata$Donor     <- factor(coldata$Donor)

se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

# ---------------------------------------------------------------------------
# Low-count filter: at least 10 counts in at least N samples.
#
# N defaults to the size of the smaller experimental group, so that a region
# accessible throughout one entire cell type survives even when it is absent
# from the other. Note what the rule does NOT do: it does not require the N
# supporting samples to come from the same group, so a region supported by any
# mixture of N samples is also retained.
#
# The project slides specify N = 3 and the assignment PDF specifies N = 4.
# N = 3 is used here; 05_filter_sensitivity.R runs both to the end and reports
# how little turns on the choice.
# ---------------------------------------------------------------------------
N <- as.integer(Sys.getenv("MIN_SAMPLES", min(table(coldata$cell_type))))
cat("low-count filter: >= 10 counts in >= ", N, " samples\n", sep = "")

keep <- rowSums(assay(se, "counts") >= 10) >= N
cat(sprintf("regions before %d, after %d (%.1f%% removed)\n",
            nrow(se), sum(keep), 100 * (nrow(se) - sum(keep)) / nrow(se)))

se <- se[keep, ]
saveRDS(se, "se_filtered.rds")
