# A peak-derived background for GO enrichment

Run after submission, to turn a documented limitation into a number.

## Why this was run

GO enrichment in the submitted analysis was run through the g:Profiler web interface against its default background, the whole annotated genome. The report lists this as a limitation, correctly: the query lists can only contain genes with a consensus region near their TSS, so the default background includes genes the experiment had no way of reporting. That inflates the apparent enrichment of anything the assay is biased toward.

The report also speculated that the effect "is most likely to inflate broad immune categories rather than to manufacture specific ones." That was not knowable without running it, and the claim should not have been made.

## What the right background is

Every protein-coding gene whose TSS lies within 5 kb of any consensus region — the same window used to call a peak-to-gene assignment credible for the query lists. In other words, the set of genes the experiment could in principle have reported.

`scripts/03_annotation_and_motifs.sh` now derives this with `bedtools window`, and the list is in [`../results/gprofiler_background.txt`](../results/gprofiler_background.txt).

## Result

| | |
|---|---|
| Protein-coding genes in the Ensembl 110 annotation | 19,419 |
| With a TSS within 5 kb of a consensus region | **16,776 (86.4%)** |
| Excluded from the proper background | 2,643 (13.6%) |

## Reading it

The universe shrinks by 13.6%, not by half. This bounds the problem: because 86.4% of protein-coding genes have an accessible region near their promoter in at least one of these seven samples, the default genome background and the peak-derived background differ far less than they would in a sparser assay. A 13.6% reduction in the universe moves enrichment *p*-values, but it cannot plausibly reverse an enrichment reported at the magnitudes seen here.

That this is a small correction is itself a finding about ATAC-seq in primary immune cells: the accessible fraction of promoters is large. It does not make the original choice right, and the direction of the shift still cannot be stated without re-running the enrichment.

## To re-run the enrichment

The g:Profiler REST API was not reachable from the environment where this was prepared, so the enrichment itself has not been recomputed. Doing it takes about a minute in the web interface:

1. Open <https://biit.cs.ut.ee/gprofiler/gost>
2. Paste `results/genes_up_promoter.txt` into the query box; organism *Homo sapiens*
3. Open **Advanced options** → **Statistical domain scope** → **Custom** and paste the contents of `results/gprofiler_background.txt`
4. Run, then repeat for `results/genes_down_promoter.txt`

What to compare: whether the top terms keep their identity and their ordering, and by how much the *p*-values move. If the top terms are stable, the submitted enrichment stands and the limitation is closed with a measurement rather than a caveat.
