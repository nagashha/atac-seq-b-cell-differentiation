# Build the sample metadata table from the official SRA run table.
#
# The CSV contains quoted fields with embedded commas, so read.csv is used
# rather than awk -F, which would split those fields incorrectly.

sra <- read.csv("SraRunTable.csv", stringsAsFactors = FALSE)
cat("SraRunTable rows:", nrow(sra), "\n")

# Take the sample order from the count matrix header. This order is
# authoritative: DESeq2 requires colData rows to match count columns exactly.
counts_header <- strsplit(readLines("counts_clean.txt", n = 1), "\t")[[1]]
samples <- counts_header[-1]
cat("Count matrix samples:", paste(samples, collapse = ", "), "\n\n")

clin <- sra[match(samples, sra$Run),
            c("Run", "Donor", "cell_type", "lineage", "treatment")]

# Prefix donor IDs with a letter: R factor levels should not start with a digit
clin$Donor     <- paste0("D", clin$Donor)
clin$cell_type <- ifelse(clin$cell_type == "Plasmablasts", "Plasmablast", clin$cell_type)
colnames(clin)[1] <- "sample"

print(clin)

stopifnot(all(clin$sample == samples))
stopifnot(!any(is.na(clin$cell_type)))
cat("\n>>> Order and integrity checks PASSED\n")
cat("Treatment (should be single):", unique(clin$treatment), "\n")
cat("Lineage (should be single):",   unique(clin$lineage), "\n")
cat("Group sizes:\n");       print(table(clin$cell_type))
cat("Donor x cell type:\n"); print(table(clin$Donor, clin$cell_type))

write.table(clin, "clinical_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)
