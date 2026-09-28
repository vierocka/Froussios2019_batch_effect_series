# Post 4: what does an omitted feature do to normalization?

**Question.** `featureCounts -t gene` does not produce rows for the rRNA loci (they are typed
`ncRNA_gene`, see the repo's "Known limits"). If 31% of a library is rRNA but rRNA never enters
the count matrix, what happens to the coding genes that are counted?

**Short answer.** The coding-gene counts themselves do not change. What changes is the denominator
used to calculate each gene's *share* of the library — and every coding gene's calculated share
rises by the same amount, because the total it's divided by got smaller.

This is a sensitivity analysis using the paper's published rRNA fraction, not a reconstruction of
the real rRNA count matrix — our own current-annotation recount is known to badly undercount rRNA
(see the wider investigation's notes), so it isn't usable as a real per-sample count here.

Under a simplified two-component model — counted coding genes plus rRNA — the published 31.21%
rRNA fraction (replicate 11) corresponds to a 0.54 log2 rescaling of coding-gene shares when rRNA
is omitted from the denominator:

**A gene representing 1% of the complete library becomes 1.45% of the coding-only total. It did
not gain a single read; its calculated share rose by 45%.**

![rRNA ignored](rRNA_ignored.png)

## What the analysis does (`rRNA_ignored.R`)
Run from the repository root; uses the 14 ExpA/ExpB coding-gene counts (`work/STAR_both/fC/`) and
one library, replicate 11 (31% rRNA per the paper's Table S2C).

**The synthetic rRNA total.** `rrna <- coding_total * p/(1-p)`, where `p` is the paper's published
rRNA fraction for this sample. Solving the two-component model `R = C·p/(1-p)` gives
`R/(C+R) = p` by construction — this is not a recovered rRNA read count, it is the amount of rRNA
that would have to exist, under that simplified model, for the counted coding total to represent
fraction `p` of coding+rRNA. It also treats the coding total as the entire non-rRNA library, which
ignores unassigned, multimapping, intergenic and other-noncoding reads — a simplifying assumption,
not a full reconstruction of library composition.

**The figure.** log2 share of the library per coding gene, computed two ways: with the synthetic
rRNA total in the denominator, and with rRNA ignored (coding genes only). Raw coding-gene counts
never change — only their computed share of the total does, by a constant `log2(1/(1-p))` for this
sample. That constant is sample-specific: it scales with each library's own published rRNA
fraction, not a fixed number across the dataset.

## Why it happens
- Ignoring rRNA removes real reads from the total without removing them from any individual gene's
  count — every other gene's share of what's left is mechanically larger. On the log scale that is
  a constant shift within one sample (post 3 already showed correlation-based QC cannot see a
  constant shift like this) — but the size of that constant differs sample to sample, scaling with
  each sample's own published rRNA fraction.
- **Total-count normalization is directly sensitive to this denominator** — its size factor is the
  column sum itself, which a 31%-rRNA omission changes substantially.
- **Median-of-ratios asks what happens to the typical gene**, not to the column sum, so it should
  be much less sensitive to one omitted dominant feature — a per-gene ratio to a reference isn't
  moved much by a single row being present or absent among thousands.

**Next post: does that difference in sensitivity actually change DE calls?** (post 5)

## What this means for QC
- "The GTF/annotation doesn't have rRNA gene models" is not a reason to skip counting rRNA reads
  separately — a dominant, uncounted feature is a real (if here synthetically modelled)
  sample-specific compositional rescaling of every other gene's apparent expression.
- The rRNA fraction per sample (the paper's own Table S2C, used directly in post 2) is still the
  measurement that reveals this; a good genome-wide correlation will not flag it (post 3).

## Limits
- **This is a synthetic sensitivity analysis, not a reconstruction of the real rRNA count
  matrix.** The rRNA total is back-calculated from the paper's published percentage and the
  counted coding total under a simplified two-component model, not recovered from the reads. An
  empirical version would need real per-locus rRNA counts in the same units as the coding matrix,
  and the paper's own Table S2C denominator (reads vs. read pairs vs. aligned vs. input) confirmed
  against how featureCounts was run here — neither is currently available (our own
  current-annotation recount undercounts badly enough to invert the sample ranking, see the wider
  investigation's notes).
- The coding total is treated as the entire non-rRNA library, which excludes unassigned reads,
  multimappers, intergenic reads and other noncoding biotypes not counted by `-t gene` — a
  simplifying assumption, not a complete accounting of the library.
- The 0.54 log2 shift is specific to one library (31% rRNA); it scales with each library's own
  published rRNA fraction, not a fixed constant across samples.
