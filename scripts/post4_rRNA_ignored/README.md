# Post 4: what changes if rRNA loci are ignored?

**Question.** `featureCounts -t gene` does not produce rows for the rRNA loci (they are typed
`ncRNA_gene`, see the repo's "Known limits"). If they are simply left out of the count matrix —
not filtered afterwards, just never counted — what happens to normalisation and to DE calls?

**Short answer.** Leaving rRNA out does not just drop a few rows. It silently shrinks the library
total, so every remaining gene's apparent share of the library shifts up by the same constant
amount (panel A). That constant shift feeds into size-factor estimation and, in a same-batch,
same-genotype comparison where every DE call is by definition a false positive, it inflates false
DE counts — more so under naive total-count scaling than under DESeq2's median-of-ratios
(panel B).

![rRNA ignored](rRNA_ignored.png)

## What the analysis does (`rRNA_ignored.R`)
Run from the repository root; uses the 14 ExpA/ExpB coding-gene counts (`work/STAR_both/fC/`)
and adds rRNA back as one extra row, at each library's published rRNA fraction (Table S2C).

**Panel A.** One library (replicate 11, 31% rRNA per the paper). Density of each coding gene's
log10 share of the library, computed two ways: with rRNA counted in the total, and with rRNA
ignored (coding genes only). The whole distribution shifts right by log2(1/(1-0.31)) ≈ 0.54 log2
when rRNA is ignored — a constant, not a per-gene effect.

**Panel B.** Same batch (ExpA), same genotype — no true biological difference is possible, so any
`padj < 0.05` call is a false positive. Replicate 6 (23.7% rRNA) is always placed in group 2;
8 random 3-vs-3 splits of the other ExpA replicates are drawn. Each split is run 4 ways: rRNA in
the matrix vs. ignored, crossed with DESeq2's own median-of-ratios size factors vs. simple
total-count scaling. One line per split connects its "in the matrix" and "ignored" result.

## Why it happens
- Ignoring rRNA removes real reads from the total without removing them from any individual
  gene's count — every other gene's share of what's left is mechanically larger. On the log scale
  that is a constant shift (post 3 already showed correlation-based QC cannot see a constant
  shift like this).
- DESeq2's median-of-ratios size factor is a per-gene ratio to a reference, so it is comparatively
  insulated from a single dominant biotype being in or out of the matrix.
- Total-count scaling uses the sum of the whole library directly, so removing a 20-30%-rRNA
  contribution changes the size factor itself, not just one row — the false-DE inflation is
  visible mainly on this side of panel B.
- Splits that put replicate 6 in an unlucky group are already bad (hundreds of false DE genes)
  before the rRNA-ignored question is even asked — panel B's flattest, highest split shows that.
  Ignoring rRNA on top of that is not what creates the false positives, but it is not neutral
  either, especially for total-count scaling.

## What this means for QC
- "The GTF/annotation doesn't have rRNA gene models" is not a reason to skip counting rRNA reads
  separately — the omission is not inert, it is a small constant bias baked into every gene's
  apparent expression.
- DESeq2's default normalisation is more robust to this specific omission than total-count
  scaling, but robustness to one failure mode is not the same as immunity to a contaminated
  replicate in the design (see the high-false-DE splits in both facets of panel B).
- The rRNA fraction per sample (post 2's direct measurement) is still the check that catches
  this; a good genome-wide correlation or a reasonable-looking DESeq2 run will not flag it.

## Limits
- Panel A's shift is computed for one library (31% rRNA); the exact shift scales with each
  library's own rRNA fraction, not a fixed constant across samples.
- Panel B uses only ExpA (14 replicates minus rep 6, 8 of the possible 3-vs-3 splits including
  rep 6), not an exhaustive enumeration — a teaching-scale illustration of the effect, not a full
  power analysis.
- rRNA is added back as a single synthetic row (its total published fraction), not as real
  per-locus counts — real rRNA reads would spread across the 4 annotated loci differently, though
  the library-total effect in panel A does not depend on that detail.
