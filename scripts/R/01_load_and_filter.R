# Assemble a SummarizedExperiment and apply the low-count filter.
suppressMessages({ library(SummarizedExperiment); library(DESeq2) })

counts  <- as.matrix(read.delim("counts_clean.txt", row.names = 1, check.names = FALSE))
coldata <- read.delim("clinical_table.txt", row.names = 1)

# Critical: count columns and metadata rows must be in identical order.
# A mismatch here would silently compare the wrong groups and never error.
stopifnot(all(colnames(counts) == rownames(coldata)))
cat(">>> Column order match VERIFIED\n")

# Bulk_B as reference => positive log2FC means "more open in plasmablasts"
coldata$cell_type <- factor(coldata$cell_type, levels = c("Bulk_B", "Plasmablast"))
coldata$Donor     <- factor(coldata$Donor)

se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

# Low-count filter: >=10 counts in at least N samples.
# N = size of the smallest experimental group, so a peak present throughout
# one entire cell type survives even if absent from the other. Note this does
# not require the N samples to come from the same group.
N <- min(table(coldata$cell_type))
cat("Smallest group size -> N =", N, "\n")

keep <- rowSums(assay(se, "counts") >= 10) >= N
cat("Before filter:", nrow(se), " after:", sum(keep),
    sprintf(" (%.1f%% removed)\n", 100 * (nrow(se) - sum(keep)) / nrow(se)))

se <- se[keep, ]
saveRDS(se, "se_filtered.rds")
