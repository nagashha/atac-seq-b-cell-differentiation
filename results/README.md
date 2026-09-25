# Results

Small output tables. Large intermediates (BAM, count matrix, full DA table)
are not committed — see the top-level README.

| File | Contents |
|---|---|
| `clinical_table.txt` | Sample metadata: donor, cell type, lineage, treatment |
| `DARs_significant.tsv` | 14,662 DARs, main analysis (`~ cell_type`, 7 samples) |
| `DARs_designB_paired.tsv` | 3,696 DARs, paired sensitivity analysis (4 matched samples) |
| `genes_up_promoter.txt` | 731 genes with a promoter-proximal DAR opening in plasmablasts |
| `genes_down_promoter.txt` | 733 genes with a promoter-proximal DAR staying open in bulk B |
| `target_genes.txt` | The six target genes selected for downstream analysis |
| `target_gene_DARs.tsv` | Every DAR assigned to one of the six targets, with distance and effect size |
| `homer_up_top15.tsv` | Top known motifs in regions opening in plasmablasts |
| `homer_down_top15.tsv` | Top known motifs in regions staying open in bulk B |

DAR columns: `baseMean`, `log2FoldChange` (apeglm-shrunk), `lfcSE`, `pvalue`,
`padj`, `PeakID`. Positive log2FC means more accessible in plasmablasts.
