# Post 7: significance vs. effect size

**Question.** Post 5 counted a gene as DE when `padj < 0.05`. Many analyses call DE on effect size
instead, e.g. `|log2FC| > 1` (2-fold) or `|log2FC| > 2` (4-fold), or combine both. Using DESeq2's
median-of-ratios, the normalization that was stable in post 5, how many genes does each rule call
in comparisons with no designed difference?

**Short answer.**

- **A fold-change cutoff alone is mostly noise.** In the 8 same-batch splits, 2,098-2,448 genes
  pass `|log2FC| > 1` and 589-826 pass `|log2FC| > 2`, including splits with 0 genes at
  `padj < 0.05`. Almost all of them have a handful of reads (median baseMean 0.7-0.9), where
  0 vs. 2 reads already is a "4-fold change" (panel A).
- **Shrunken fold changes (apeglm) remove that noise**: 1-102 genes at `> 1`, 0-10 at `> 2`.
- **Combining significance with effect size** (panel B) keeps only the large, well-supported
  changes: `padj < 0.05 & |log2FC| > 1` gives 0-140 genes with the raw log2FC and 0-93 with the
  shrunken one; `& |log2FC| > 2` gives 0-11 (raw) and 0-8 (shrunken).
- **Split 6 (post 6) stays the outlier under every rule that accounts for uncertainty**, but most
  of its 689 genes change by less than 2-fold: its signal is many modest shifts, not a few large
  ones.

![effect size](effect_size.png)

## What the analysis does (`effect_size.R`)
Run from the repository root; needs DESeq2 and apeglm; takes about 15 minutes (apeglm).
Same coding-gene counts, synthetic rRNA feature, 8 splits (`set.seed(1)`) and 4 size-factor
treatments as post 5 (median-of-ratios or total-count scaling, rRNA feature included or ignored
for size factors; DE always tested on the coding-only matrix). Each fit is scored under ten rules:

| Rule | What it uses |
|---|---|
| `padj < 0.05` | post 5's rule (Wald test, H0: log2FC = 0) |
| `\|log2FC\| > 1`, `\|log2FC\| > 2` (raw) | maximum-likelihood log2FC, no uncertainty |
| `\|log2FC\| > 1`, `\|log2FC\| > 2` (shrunken) | apeglm-shrunken log2FC (`lfcShrink`), pulls noisy estimates towards 0 |
| `padj < 0.05 & \|log2FC\| > 1`, `& \|log2FC\| > 2` (raw) | significance plus an effect-size filter on the raw log2FC |
| `padj < 0.05 & \|log2FC\| > 1`, `& \|log2FC\| > 2` (shrunken) | significance plus an effect-size filter on the shrunken log2FC |
| test `\|log2FC\| > 1`, `padj < 0.05` | `results(lfcThreshold = 1)`, H0: \|log2FC\| <= 1 |

**The figure** shows median-of-ratios only. Including or ignoring the rRNA feature gives the same
median-of-ratios size factors (post 5), so each split is drawn once. Panel A: significance vs.
effect-size cutoffs alone (raw and shrunken). Panel B: `padj < 0.05` alone and combined with
`|log2FC| > 1` and `> 2`, raw and shrunken. One line per split; split 6 highlighted.

Outputs: `effect_size_results.csv` (genes called per split, normalization, rRNA treatment and
rule), `raw_log2FC_expression.csv` (how well expressed the genes passing a raw cutoff are),
`log2FC_shift.csv` (see "Side note" below).

## Results (median-of-ratios)
Genes called, min-max over the 8 splits:

| Rule | genes | split 6 |
|---|---|---|
| `padj < 0.05` | 0-689 | 689 |
| `\|log2FC\| > 1` (raw) | 2,098-2,448 | 2,448 |
| `\|log2FC\| > 2` (raw) | 589-826 | 622 |
| `\|log2FC\| > 1` (shrunken) | 1-102 | 102 |
| `\|log2FC\| > 2` (shrunken) | 0-10 | 10 |
| `padj < 0.05 & \|log2FC\| > 1` (raw) | 0-140 | 140 |
| `padj < 0.05 & \|log2FC\| > 1` (shrunken) | 0-93 | 93 |
| `padj < 0.05 & \|log2FC\| > 2` (raw) | 0-11 | 11 |
| `padj < 0.05 & \|log2FC\| > 2` (shrunken) | 0-8 | 8 |
| test `\|log2FC\| > 1` | 0-3 | 3 |

Genes passing a raw cutoff: median baseMean 0.7-0.9 reads (all genes: about 410); 88-97% of the
`|log2FC| > 1` genes and 96-100% of the `|log2FC| > 2` genes have baseMean below 10.

The median-of-ratios counts are identical with and without the synthetic rRNA feature under every
rule, except one gene once (`|log2FC| > 1`, 2,214 vs. 2,215): adding one row changes the size
factors by about 1e-6, enough to move a gene that sits exactly at the cutoff.

## Why it happens
- Raw log2FC estimates from genes with almost no reads are extreme by construction. They dominate
  any raw fold-change cutoff, which is why thousands of genes pass in comparisons with no designed
  difference.
- Shrinkage and significance both account for that uncertainty, so they agree: the same few
  splits, led by split 6, carry the real signal.
- Effect-size filters on top of `padj` answer a different question ("changed, and by how much?").
  For split 6 they show that its 689 genes are mostly modest changes (median |log2FC| 0.71, post 6).
- On top of `padj`, raw and shrunken filters mostly agree (140 vs. 93 genes for split 6 at `> 1`):
  significance has already removed the noisy low-count genes; shrinkage then trims the estimates of
  the remaining ones slightly, so genes just above the cutoff fall below it. Post 9 shows how.

## What this means for QC
- "Fold change > 2" without a measure of uncertainty is not a DE call. In 3-vs-3 designs it is
  dominated by genes with a few reads.
- Use shrunken log2FC for ranking and effect-size filtering, combine a fold-change filter with
  `padj`, or use `lfcThreshold` when the question is "changed by more than X-fold".
- Report the rule together with the DE count: the same data give 0-3, 0-93, 0-140, 0-689 or
  2,000+ "DE genes" depending on it.

**Next posts:** what median-of-ratios actually computes (post 8), how apeglm shrinks a fold change
(post 9), and DESeq2 (+ apeglm) vs. rlog + t-test on the same splits (post 10).

## Side note: total-count scaling shifts effect sizes
With total-count scaling (not in the figure), ignoring the synthetic rRNA feature shifts every
gene's log2FC by about the same amount: median over genes 0.09-0.12 log2 per split (IQR over genes
about 0.02), close to the 0.11-0.14 expected from the size factors alone (`log2FC_shift.csv`).
Fixed fold-change cutoffs are directly exposed to such a shift: the raw `|log2FC|` counts change in
all 8 splits (e.g. 2,345 -> 2,093 genes), the shrunken `> 1` count in 7, `padj < 0.05` and both
`padj & |log2FC| > 1` rules in 4 (post 5), both `padj & |log2FC| > 2` rules in 2, the `lfcThreshold`
test in 1. With median-of-ratios the shift is 0. This is why the
figure uses median-of-ratios.

## Limits
- Same synthetic rRNA feature and the same 8 sampled splits as post 5 (not all 60 rep-6 splits);
  see posts 4-5 for the two-component model and its assumptions.
- One shrinkage method (apeglm; ashr and normal not compared), one significance level. apeglm
  reports optimizer warnings ("line search routine failed") for a few dozen genes across the 32
  fits; they do not change the counts above.
- Gene counts were made unstranded (`featureCounts` default `-s 0`) although the libraries are
  stranded (paper section 2.1); post 11 tests `-s 0/1/2` and reruns this post on stranded counts.
