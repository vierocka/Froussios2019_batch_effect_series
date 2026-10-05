# Roadmap

Map of the repository. New files are added post by post; this file is updated with them.

```
Froussios2019_batch_effect_series/
├── README.md                          purpose, run order, known limits
├── ROADMAP.md                         this file
├── config.sh                          shared paths and thread count (override with environment variables)
├── .gitignore                         keeps raw data, BAMs, work/ and AED files out of git
├── work/                              not in git: generated FASTQ, reference, BAMs and counts (see README "Generated data")
│
├── publication/
│   ├── Froussios2019_Bioinformatics_btz089.pdf   the paper (open access, CC BY 4.0)
│   └── Froussios2019_S2_Table.csv                paper's Supplementary Table S2 (alignment and rRNA summaries)
│
├── data/                              only the small tables the posts use
│   ├── run_list.tsv                   17 ENA runs: accession, batch, ERCC mix, FTP paths, MD5
│   ├── replicate_map.tsv              run -> paper's replicate number / experiment (differs from ENA sample titles)
│   └── paper_TableS2C_rRNA.csv        rRNA reads and % per replicate, as published (Table S2C)
│
└── scripts/
    ├── sra_download/
    │   └── download_sra.sh            download the 17 runs with sra-tools; skips finished runs
    ├── reference/
    │   └── prep_reference.sh          TAIR10 genome + annotation (Ensembl Plants 63) and STAR index
    │
    ├── post2_rRNA_correlation/        rRNA content vs. genome-wide correlation
    │   ├── trim_star_featurecounts.sh   Trimmomatic -> STAR -> featureCounts gene counts
    │   ├── count_rRNA.sh                align with the paper's STAR settings, count reads on rRNA genes
    │   ├── rRNA_correlation_blindspot.R rlog correlation with/without rRNA genes; makes the figure
    │   ├── rRNA_correlation_blindspot.png   figure of the post
    │   └── README.md                    explanation, run order and limits (incl. the rRNA recount caveat)
    │
    ├── post3_compositional_shift/     teaching simulation
    │   ├── simulation.R               two correlated log2 vectors; 10 extreme genes changed, 1000 loops
    │   ├── compositional_shift.png    figure of the post
    │   └── README.md                  explanation of the simulation and what it means for QC
    │
    ├── post4_rRNA_ignored/            what does an omitted feature do to normalization?
    │   ├── rRNA_ignored.R             one library (rep 11), two-component share rescaling
    │   ├── rRNA_ignored.png           figure of the post
    │   └── README.md                  explanation of the synthetic two-component model and its limits
    │
    ├── post5_normalization_false_DE/  can normalization choice change false-DE counts?
    │   ├── false_DE_by_normalization.R          8 same-batch 3 vs 3 splits (rep 6), DE with size factors from a synthetic rRNA-added vs. coding-only matrix
    │   ├── false_DE_by_normalization_results.csv   underlying numbers
    │   ├── false_DE_by_normalization.png        figure of the post
    │   └── README.md                            results, interpretation and limits
    │
    ├── post6_split_structure/         why does one same-batch 3-vs-3 split give 689 DE genes?
    │   ├── split_structure.R          split 6 DE, all 70 ExpA splits, PCA, leave-one-out/swaps; makes the figure
    │   ├── split6_DE_genes.csv        the 689 DE genes (padj < 0.05)
    │   ├── all_70_splits.csv          DE count and rRNA separation for every 3-vs-3 split of ExpA
    │   ├── split_structure.png        figure of the post
    │   └── README.md                  results, interpretation and limits
    │
    └── post7_effect_size/             significance vs. effect size as DE rules
        ├── effect_size.R              post 5's 8 splits x 4 size-factor treatments, ten DE call rules
        ├── effect_size_results.csv    genes called per split, treatment and rule
        ├── log2FC_shift.csv           log2FC shift per split when the rRNA feature is ignored
        ├── raw_log2FC_expression.csv  expression of the genes passing a raw log2FC cutoff
        ├── effect_size.png            figure of the post
        └── README.md                  results, interpretation and limits
```

## Posts
| Post | Topic | Where |
|---|---|---|
| 1 | QC observations from the paper text and supplement | `publication/` |
| 2 | rRNA fraction vs. genome-wide correlation | `scripts/post2_rRNA_correlation/` |
| 3 | 10 extreme genes barely move a correlation | `scripts/post3_compositional_shift/` |
| 4 | what does an omitted feature do to normalization? | `scripts/post4_rRNA_ignored/` |
| 5 | can normalization choice change false-DE counts? | `scripts/post5_normalization_false_DE/` |
| 6 | why does one same-batch 3-vs-3 split give 689 DE genes? | `scripts/post6_split_structure/` |
| 7 | significance vs. effect size: padj < 0.05, \|log2FC\| > 1 / > 2 and combinations | `scripts/post7_effect_size/` |
| 8 (next) | what does DESeq2's median-of-ratios actually do? | `scripts/post8_median_of_ratios/` (planned) |
| 9 | how does apeglm shrink a fold change? | `scripts/post9_apeglm_shrinkage/` (planned) |
| 10 | DESeq2 (+ apeglm) vs. rlog + t-test on the same splits | `scripts/post10_deseq2_vs_rlog_ttest/` (planned) |
| 11 | how can 31% rRNA become 0.006% in a different counting approach? (+ strandedness) | `scripts/post11_rRNA_recount/` (planned) |

## Order of use
`sra_download` -> `reference` -> `post2_rRNA_correlation` (trim/STAR/count, then rRNA count, then figure). `post3_compositional_shift` runs on its own. `post4_rRNA_ignored` and `post5_normalization_false_DE` both only need `post2`'s coding-gene counts (`work/STAR_both/fC/`) and `data/paper_TableS2C_rRNA.csv`; post 5 additionally needs DESeq2. `post6_split_structure` and `post7_effect_size` need the same counts (post 6 also ggrepel, a few minutes with `THREADS`; post 7 also apeglm, about 15 minutes).

The teaching arc:
- **Posts 2-3, detection:** rRNA content varies about 25-fold among "identical" replicates, and
  genome-wide correlation cannot see it (post 3 shows why).
- **Posts 4-5, consequences:** omitting a dominant feature rescales every gene's share (post 4);
  whether that changes DE calls depends on the normalization (post 5).
- **Posts 6-7, reading DE results:** a split with hundreds of DE genes reflects structure among the
  replicates, not rRNA (post 6); the DE rule (significance, raw or shrunken effect size, or both)
  changes the count by orders of magnitude (post 7).
- **Posts 8-10, under the hood** of the tools used in posts 5-7: what median-of-ratios computes and
  why one extra row barely moves it (post 8); how apeglm shrinks fold changes (post 9); what changes
  when the same splits are tested with rlog + t-test instead of DESeq2 (post 10).
- **Post 11, back upstream:** why a published 31.21% rRNA fraction becomes 0.006% in a recount,
  and whether counting unstranded (`-s 0`) changed any result in posts 2-10.

## Ready to publish

**Post 6: "Same genotype, same batch, no designed difference: why does one 3-vs-3 split give
689 DE genes?"** Done: `scripts/post6_split_structure/` (script, figure, README with all numbers).
Short version: the 7 ExpA replicates carry a gradient that follows replicate numbering (PC1 64%),
not rRNA; replacing rep 6 with rep 7 gives 1,389 DE genes, swapping rep 3 and rep 7 gives 5.
Open follow-ups:
- GO enrichment of the 689 genes (top genes ATHB-2, HAT1, ATHB4 suggest light / time of day).
- Sample metadata: ArrayExpress E-MTAB-5446 SDRF for harvest order, plate, extraction and library
  dates; if absent, ask the authors.
- Does ExpB show the same kind of gradient (PCA of replicates 8-14)?
- Side finding for a sentence somewhere: chloroplast (Pt) genes take 58-61% of the counted
  coding-gene reads in every ExpA sample (Mt 0.3%); a second dominant feature, stable across
  samples, relevant to total-count scaling (posts 4-5).

**Post 7: "Significance vs. effect size: how many genes does each DE rule call when there is no
designed difference?"** Done: `scripts/post7_effect_size/` (script, figure, README with all
numbers). DESeq2 median-of-ratios (the stable normalization from post 5), same 8 splits as post 5,
ten DE call rules. Figure: A = `padj` vs. raw and shrunken `|log2FC| > 1 / > 2`; B = `padj < 0.05`
alone and combined with raw and shrunken `|log2FC| > 1 / > 2`.
Short version: raw fold-change cutoffs alone flag 2,098-2,448 (`> 1`) and 589-826 (`> 2`) genes per
split, almost all with a handful of reads; shrinkage (apeglm), `padj & |log2FC|` (raw or shrunken)
and `lfcThreshold` give 0-140 and keep split 6 as the outlier, whose 689 genes are mostly < 2-fold changes. Side note
in the README: under total-count scaling, ignoring the rRNA feature shifts every log2FC by about
0.1 log2, so raw cutoff counts change in all 8 splits.
Open follow-ups:
- All 60 rep-6 splits instead of 8. (apeglm vs. normal vs. ashr shrinkage: post 9.)
- Same rules on post 6's all-70-split table: which rule best separates split-6-like structure
  from noise?

## Planned posts, in order

**Post 8: "What does DESeq2's median-of-ratios actually do?"** Planned folder:
`scripts/post8_median_of_ratios/` (own script + README + figure). Explains the normalization that
was stable in post 5 and used in posts 6-7.
- Step by step on the ExpA counts: per-gene geometric mean across samples (the pseudo-reference;
  genes with a zero in any sample drop out: how many?), each sample's ratios to it, the median of
  those ratios = the size factor.
- Compare with total-count size factors (column sums), sample by sample; replicate 6 (23.7% rRNA)
  is where they should differ most once the synthetic rRNA row is added.
- Why one dominant row barely moves it: one more ratio among about 20,000; post 7 measured the
  change at about 1e-6. Contrast: total-count scaling moves by the row's whole share.
- Chloroplast genes take 58-61% of the counted reads (post 6 side finding): what that does to
  total-count vs. median-of-ratios size factors.
- Assumption and when it breaks: most genes are not DE, and there is no global shift in total
  mRNA. Spike-ins (ERCC, present in these libraries but not in our reference) are the usual check.
- Figure idea: A = per-sample distribution of gene ratios to the pseudo-reference, median marked;
  B = size factors from median-of-ratios vs. total counts, with and without the rRNA row.

**Post 9: "How does apeglm shrink a fold change?"** Planned folder:
`scripts/post9_apeglm_shrinkage/`. Explains why raw `|log2FC| > 1` calls 2,000+ genes per
no-difference split in post 7 but shrunken `|log2FC| > 1` only 1-102.
- MA plots (log2FC vs. mean expression), raw vs. shrunken, for one null split and split 6.
- Shrinkage vs. information: low-count / high-variance genes are pulled to 0, well-measured genes
  keep their estimate. Pick a few example genes (e.g. 0 vs. 2 reads: raw 4-fold, shrunken ~0; a
  well-expressed split-6 gene: almost unchanged) and show raw and shrunken values with their SE.
- How it works, in plain terms: the estimate is the posterior mode under a heavy-tailed prior on
  log2FC whose scale is estimated from all genes; it changes the effect size, not the Wald p-value
  (`padj` stays the same; apeglm's s-values are an alternative).
- Compare apeglm with `type = "normal"` and `"ashr"` (install ashr).

**Post 10: "DESeq2 (+ apeglm) vs. rlog + t-test: same splits, different answers?"** Planned folder:
`scripts/post10_deseq2_vs_rlog_ttest/`.
- Same splits as post 5 (8) and post 6 (all 70). Per split: DESeq2 `padj < 0.05`; DESeq2 `padj` &
  shrunken `|log2FC| > 1`; rlog (and vst) + per-gene Welch t-test + BH; limma on rlog as the
  moderated middle ground.
- Questions: how many genes each calls in no-designed-difference splits; overlap of the gene lists
  in split 6; which genes only the t-test calls (low counts? genes with one extreme replicate?);
  `rlog(blind = TRUE)` vs. `blind = FALSE`.
- Point: with n = 3 per group a gene-by-gene t-test estimates each variance from 4 degrees of
  freedom; DESeq2 shares information across genes (dispersion trend). rlog/vst are meant for
  visualization and clustering, not for DE testing (DESeq2 vignette); show what goes wrong when
  they are used for it.

**Post 11: "How can a sample reported to contain 31% rRNA become 0.006% in a different counting
approach?"** Planned folder: `scripts/post11_rRNA_recount/` (own script + README + figure).
Placed last because it needs the BAMs on the HPC (realignment / recounting) and reruns the earlier
posts on stranded counts.
The direct rRNA recount (`count_rRNA.sh`: current TAIR10.63 annotation, STAR
`--outFilterMultimapNmax 2`, featureCounts without `-M`) is not just a methods footnote; it is
arguably a stronger result than the synthetic perturbations of posts 4-5.

- Paper (Table S2C): replicate 11 has 31.21% rRNA, ranks **highest** of the 14 ExpA/ExpB samples.
- Our current-annotation/default-counting recount: 0.006%, ranks **lowest**.
- Not a quantitative underestimate: the sample ordering is not preserved (rep 11 flips from most-
  to least-contaminated).
- Shows that "I counted the annotated rRNA genes" is not a valid rRNA QC measurement when the
  annotation contains too few rRNA loci (rDNA is highly repetitive; most true copies are
  collapsed or unassembled in the reference), the aligner drops reads mapping to many loci
  (`--outFilterMultimapNmax`), and the counter drops multimappers by default (no `-M`).
- Source material exists in the wider investigation's notes (`NOTES.md`,
  `froussios2019_paperrepro_04_rRNA_multimapper_ramses.sh`, the lenient-multimap rerun that tests
  how much of the gap that closes); adapt it to this repo's post style, do not write from scratch.
- Minimum figure: per-sample rRNA % from the paper vs. our recount(s), on a log axis, with ranks,
  so the flip of rep 11 is visible. Ideally a ladder of counting choices (see below) showing which
  step recovers how much of the 31%.
- Open question to settle first: what is Table S2C's denominator (reads, read pairs, aligned
  reads, input)? A % is only comparable if numerator and denominator are in the same units.

**Post 11 to-do: strandedness (`featureCounts -s`) as one of the counting choices.**
The paper used "Illumina TruSeq Stranded Total RNA with Ribo-Zero Plant" (section 2.1), i.e. a
reverse-stranded library, but `trim_star_featurecounts.sh` counts genes with the featureCounts
default `-s 0` (unstranded), while `count_rRNA.sh` uses `-s 2`. Posts 2 and 4-7 all use the
unstranded gene counts. Plan:
1. **Confirm the strandedness from the data, not only from the kit name.** Cheapest: run
   featureCounts with `-s 1` and `-s 2` on the same BAMs; for a TruSeq Stranded library the `-s 2`
   assigned fraction should be far higher than `-s 1`. Optionally RSeQC `infer_experiment.py` on 2-3
   BAMs (needs a BED12 made from the GFF3). Report the fractions per sample.
2. **Gene counts with `-s 0`, `-s 1`, `-s 2`** from the existing BAMs (no realignment; BAMs are on
   the HPC, `work/STAR_both/*.bam`). Write each to its own folder (`work/STAR_both/fC_s0|s1|s2/`)
   and add a `STRAND` setting to `config.sh` (default 0, so the published posts still reproduce).
   Per sample: assigned %, and how many genes change by more than 10% / 2-fold between `-s 0` and
   `-s 2` (expected: genes overlapping antisense transcripts or neighbouring genes on the other
   strand lose reads they should never have had).
3. **rRNA counts with `-s 0`, `-s 1`, `-s 2`** (`count_rRNA.sh`) as one rung of the post 11 ladder:
   paper 31.21% -> current recount (0.006%) -> `-s` variants -> featureCounts `-M` (multimappers)
   -> higher `--outFilterMultimapNmax` -> annotation-independent count (SortMeRNA/BBDuk vs. SILVA).
   For each rung: rep 11's %, and the Spearman correlation of the 14 samples' ranking with Table
   S2C. Expected: `-s` alone explains little of the gap (the main losses are missing rDNA loci and
   discarded multimappers), but a wrong `-s` would make it worse; check, don't assume.
4. **Downstream check: rerun posts 2 and 4-10 on the `-s 2` gene counts** and report whether any
   conclusion changes (post 2 correlation level, post 4's 0.54 log2 shift, post 5's 4-of-8 splits,
   post 6's 689 DE genes and PC1 gradient, post 7's ranges, posts 8-10). If nothing material
   changes, that is one sentence in post 11; if something does, it gets its own correction note in the affected
   post's README.

## Investigation backlog

Every result in posts 2-11 depends on choices made upstream. Each axis below is a candidate
sensitivity analysis (or post) asking: does the conclusion survive a different reasonable choice?

**Normalisations**
- Beyond median-of-ratios and total-count scaling (post 5): edgeR TMM, upper-quartile, CPM/TPM,
  and ERCC spike-in based scaling (the runs carry ERCC mix 1/2, see `data/run_list.tsv`).
- Removal of unwanted variation (RUVg with ERCCs as controls, RUVs with replicates), and
  `~ batch` (ExpA/ExpB) in the design vs. not.
- Post 5 done exhaustively: all 60 rep-6 splits instead of 8, ExpB as well, and with rep 11 in.
- Real per-sample rRNA counts (once post 11 has a usable recount) in place of the synthetic feature.

**Filtering**
- Low-count gene filters: `rowSums > 0` (current) vs. `edgeR::filterByExpr` vs. >= 10 counts in
  >= n samples; DESeq2 independent filtering on vs. off.
- Sample filtering: with vs. without rep 11 (excluded by the paper) and rep 6 (kept); what
  rRNA-fraction threshold would a reasonable exclusion rule use, and is it applied consistently?
- Read-level filtering: Trimmomatic settings, STAR `--outFilterMultimapNmax` (2 vs. 10 vs.
  higher), mismatch limits, MAPQ filters downstream.

**Annotations**
- Ensembl Plants 63 (4 rRNA loci typed `ncRNA_gene`) vs. Araport11 vs. the original TAIR10 GFF:
  how many rRNA loci (45S, 5S, organellar) does each contain, and under which feature type /
  biotype? The `grep ribos | grep ncRNA_gene` in `count_rRNA.sh` may miss rRNA typed otherwise.
- Organellar rRNA (chloroplast 16S/23S/4.5S/5S, mitochondrial): is it annotated and counted,
  and should it count as "rRNA" in the same sense as Table S2C?
- Do gene-level conclusions (posts 2, 5) change with the annotation version?

**Pipelines**
- Reproduce the paper's own pipeline as described in its Methods vs. ours (Trimmomatic -> STAR ->
  featureCounts) vs. a standard one (e.g. nf-core/rnaseq).
- Alignment-free quantification (Salmon, kallisto) with and without decoys.
- Different aligner (HISAT2) with the same counting step.

**Counting methods**
- featureCounts `-M` (multimappers), `-M --fraction`, `-O` (overlapping features), `-t exon` vs.
  `-t gene`, `-p` with/without `-B -P`.
- Strandedness (`-s 0/1/2`): planned as part of post 11, see "Post 11 to-do" above.
- HTSeq-count (union / intersection modes) vs. featureCounts.
- Annotation-independent rRNA measurement: SortMeRNA or BBDuk against SILVA/Rfam rRNA,
  FastQ Screen; this does not depend on how many rDNA copies the reference has.

**References**
- TAIR10 (rDNA arrays largely collapsed or missing) vs. newer near-complete Col-0 assemblies
  (e.g. Col-CEN, Naish et al. 2021); check how much of the NOR 45S rDNA arrays each contains.
- TAIR10 plus a full rDNA unit added as an extra contig, so rRNA reads have somewhere to map.
- Toplevel FASTA includes chloroplast and mitochondrion; confirm, and test with/without them.
- Ensembl Plants release (63 vs. older) and the ERCC sequences: are spike-ins in the reference?
