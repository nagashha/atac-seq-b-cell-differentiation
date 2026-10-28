# Filter-threshold sensitivity: N = 3 versus N = 4

Run after submission, to replace an argument with a measurement.

## Why this was run

The low-count filter keeps a consensus region if it has at least 10 counts in at least *N* samples. The project slides specify N = 3, the size of the smaller experimental group; the assignment PDF specifies N = 4. The submitted analysis used N = 3, and the report defends that choice by reasoning: N = 3 is the size of the plasmablast group, so a region accessible throughout one entire cell type survives even when absent from the other.

The reasoning is sound but it is still only reasoning. This check runs both thresholds to the end. Everything downstream of the filter is identical between the two arms, so any difference is attributable to the filter alone.

## What was run

`scripts/R/05_filter_sensitivity.R`, from `counts_clean.txt` and `clinical_table.txt`, both arms with `~ cell_type` on all seven samples, `results(alpha = 0.05)` and the same LFC shrinkage, differential regions called at FDR < 0.05 and |log2FC| ≥ 1.

One caveat on reproduction. This re-run used DESeq2 1.42.0 with DESeq2's built-in `normal` shrinkage, because `apeglm` could not be installed in the environment available at the time. The submitted analysis used DESeq2 1.50.2 with `apeglm`. The N = 3 arm here therefore returns 14,992 differential regions where the submitted analysis reported 14,662 — a 2.2% difference attributable to the shrinkage estimator, not to the data. Both arms below use the same estimator, so the comparison between them is unaffected.

The re-run does reproduce the filtering step exactly: 101,593 regions retained at N = 3, of which 101,592 are tested. That one-region gap is a single region with `padj = NA` from Cook's-distance outlier filtering, which is what the report assumed without having verified it.

## Results

| | N = 3 | N = 4 |
|---|---|---|
| Regions retained of 155,948 | 101,593 (34.9% removed) | 83,434 (46.5% removed) |
| Regions tested | 101,592 | 80,198 |
| Differential regions | 14,992 | 12,932 |
| more accessible in plasmablasts | 8,154 | 6,443 |
| more accessible in bulk B | 6,838 | 6,489 |
| PC1, 500 most variable regions | 86.6% | 88.3% |
| PC1, all retained regions | 44.7% | 46.6% |

Concordance across the 80,197 regions tested under both thresholds:

| | |
|---|---|
| Shared differential regions | 12,373 — 95.7% of the N = 4 set, 82.5% of the N = 3 set |
| log2FC correlation | Pearson *r* = 0.9999, Spearman ρ = 0.9999 |
| log2FC difference | median 0.025, maximum 0.595 |
| Direction agreement, all tested regions | 98.36% |
| Direction agreement, regions differential in either arm | **100.00%** (n = 12,992) |

## Reading it

**The stricter threshold does not reveal anything the looser one missed.** 95.7% of the N = 4 differential regions are also N = 3 differential regions. Only 559 are unique to N = 4, and those arise because removing 21,395 regions from the matrix changes the multiple-testing burden and the independent-filtering threshold, not because new signal appeared.

**The 2,619 regions unique to N = 3 are exactly what the threshold is there to decide about.** Their median baseMean is 12.4, against 31.8 for the differential set as a whole. These are low-count regions — accessible in the three plasmablasts and near-zero in bulk B, which is the pattern N = 3 is designed to retain and N = 4 discards. Whether they are signal or noise is the substance of the disagreement between the slides and the PDF, and the data here cannot settle it; what the data do show is that including them changes no effect-size estimate anywhere else.

**Effect sizes are essentially identical.** Pearson *r* = 0.9999 over 80,197 regions, with a median absolute difference of 0.025 in log2 fold change. The 1,316 regions where the two arms disagree on direction all have |log2FC| ≤ 0.042, which is noise around zero; not one of them is differential in either arm. Restricted to regions that are differential in either arm, direction agreement is 100%.

**PC1 is 87% of the variance among the 500 most variable regions, not of the whole matrix.** Across all 101,593 retained regions PC1 accounts for 44.7%. Both numbers are honest and they answer different questions; the report quotes the first, which is DESeq2's `plotPCA` default, and should say which it is. The ordering is unchanged under either threshold.

## Conclusion

The filter threshold is not a load-bearing choice. Every biological conclusion in the report — the direction of the differential set, the IRF/ETS motif split, the target genes — rests on regions that are differential under both thresholds. Had the assignment's N = 4 been used as the main analysis, the report would have said the same things with slightly smaller numbers.
