# Post 4: what changes if rRNA loci are ignored?

**Question.** `featureCounts -t gene` does not produce rows for the rRNA loci (they are typed
`ncRNA_gene`, see the repo's "Known limits"). If they are simply left out of the count matrix —
not filtered afterwards, just never counted — what happens to normalisation and to DE calls?

**Short answer.** This is a sensitivity analysis using a *synthetic* rRNA total (our own
current-annotation rRNA recount is known to badly undercount, see `NOTES` in the wider
investigation — not usable as a real per-sample rRNA count). Treating the paper's published rRNA
fraction as a synthetic addition to the counted coding-gene library shows that omitting it
rescales every coding gene's within-library *share* by a sample-specific constant (panel A).
DESeq2's median-of-ratios normalisation is largely unchanged by this, because the added row is one
among thousands. Total-count normalisation is more sensitive, since its size factor depends
directly on the matrix column sum — omitting rRNA can increase *or* decrease the number of
`padj < 0.05` calls relative to the randomized group labels (panel B). Some splits produce
hundreds of such calls under either treatment, showing rRNA omission is not the only source of
structure among these replicates — one of them (replicate 6) carries 23.7% rRNA by itself.

![rRNA ignored](rRNA_ignored.png)

## What the analysis does (`rRNA_ignored.R`)
Run from the repository root; uses the 14 ExpA/ExpB coding-gene counts (`work/STAR_both/fC/`).

**The synthetic rRNA row.** `rrna <- colSums(coding) * p/(1-p)`, where `p` is the paper's
published rRNA fraction (Table S2C) for that sample. Solving `R = C·p/(1-p)` gives
`R/(C+R) = p` by construction — this is not a recovered rRNA read count, it is *the amount of
rRNA that would have to exist for the counted coding total to represent fraction `p` of
coding+rRNA*. It also treats `colSums(coding)` as if it were the entire non-rRNA library, which
ignores unassigned, multimapping, intergenic and other-noncoding-biotype reads — a simplifying
assumption, not a full reconstruction of the true library composition.

**Panel A.** One library (replicate 11, 31% rRNA per the paper). log2 share of the library per
coding gene, computed two ways: with the synthetic rRNA row in the total, and with rRNA ignored
(coding genes only). Raw coding-gene counts never change — only their computed share of the total
does, and it shifts by a constant `log2(1/(1-p))` ≈ 0.54 log2 for this sample. That constant is
sample-specific: it scales with each library's own rRNA fraction, not a fixed number across the
dataset.

**Panel B.** Same batch (ExpA), same genotype — there is no *designed* group difference between
the two labels, so a `padj < 0.05` call is false with respect to the randomized labels (this does
not mean the replicates have no real biological or technical heterogeneity among themselves —
replicate 6, 23.7% rRNA, is included in every split below, by construction). 8 sampled 3-vs-3
splits are drawn from the 7 ExpA replicates: replicate 6 plus 2 others in group 2, 3 of the
remaining 4 in group 1 (1 of the 7 unused per split; `set.seed(1)` makes the selection
reproducible).

Each split is run 4 ways, crossing rRNA-in-the-size-factors vs. ignored with DESeq2's own
median-of-ratios vs. simple total-count scaling — but DE is always tested on the *same*
coding-only matrix; only the size factors differ between the "in the matrix" and "ignored" runs.
This isolates the normalisation effect: any difference in DE counts comes from the size factors
themselves, not from DESeq2's dispersion-trend fitting or independent filtering seeing a different
set of rows (which adding a single extra rRNA row could otherwise influence). Panel B's `dB` table
is written to `rRNA_ignored_panelB_results.csv`.

## Why it happens
- Ignoring rRNA removes real reads from the total without removing them from any individual
  gene's count — every other gene's share of what's left is mechanically larger. On the log scale
  that is a constant shift within one sample (post 3 already showed correlation-based QC cannot
  see a constant shift like this) — but the size of that constant differs sample to sample,
  scaling with each sample's own rRNA fraction.
- DESeq2's median-of-ratios size factor is a per-gene ratio to a reference, so it is comparatively
  insulated from a single dominant biotype being in or out of the normalisation matrix.
- Total-count scaling uses the sum of the whole library directly, so removing a 20-30%-rRNA
  contribution changes the size factor itself — panel B's total-count facet shows more movement
  between "in the matrix" and "ignored" than the median-of-ratios facet does, in both directions.
- Splits that put replicate 6 in an unlucky group are already bad (hundreds of `padj < 0.05`
  calls) regardless of whether rRNA is in the size factors or not — panel B's highest-count split
  shows that. Ignoring rRNA is not the sole source of that structure; it shifts the count, it does
  not create the underlying problem.

## What this means for QC
- "The GTF/annotation doesn't have rRNA gene models" is not a reason to skip counting rRNA reads
  separately — the omission is not inert, it is a real (if here synthetically modelled)
  sample-specific compositional rescaling of every other gene's apparent expression.
- DESeq2's default normalisation is more robust to this specific omission than total-count
  scaling, but robustness to one failure mode is not the same as immunity to a contaminated
  replicate in the design (see the high-count splits in both facets of panel B).
- The rRNA fraction per sample (post 2's direct measurement) is still the check that catches
  this; a good genome-wide correlation or a reasonable-looking DESeq2 run will not flag it.

## Limits
- **This is a synthetic sensitivity analysis, not a reconstruction of the real rRNA count
  matrix.** The rRNA row is back-calculated from the paper's published percentage and the counted
  coding total, not recovered from the reads. An empirical version would need real per-locus rRNA
  counts in the same units as the coding matrix, and the paper's own Table S2C denominator
  (reads vs. read pairs vs. aligned vs. input) confirmed against how featureCounts was run here —
  neither is currently available (our own current-annotation recount undercounts badly enough to
  invert the sample ranking, see the wider investigation's notes).
- `colSums(coding)` is treated as the entire non-rRNA library, which excludes unassigned reads,
  multimappers, intergenic reads and other noncoding biotypes not counted by `-t gene` — a
  simplifying assumption, not a complete accounting of the library.
- Panel A's 0.54 log2 shift is specific to one library (31% rRNA); it scales with each library's
  own rRNA fraction, not a fixed constant across samples.
- Panel B draws 8 of the 15 possible 3-vs-3-from-7 splits, not an exhaustive enumeration — a
  teaching-scale illustration, not a full power analysis. Replicate 6 is in every split by
  construction, so panel B specifically probes "what happens when a known-contaminated replicate
  is in the design", not a randomly representative sample of ExpA as a whole.
- The strongest result here is not that ignoring rRNA necessarily changes DE calls in one
  direction; it is that total-count normalisation is composition-sensitive to a single dominant
  biotype, while median-of-ratios is comparatively robust to this specific one-row perturbation.
