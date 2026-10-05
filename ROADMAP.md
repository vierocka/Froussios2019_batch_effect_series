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
    └── post7_effect_size/             does normalization also shift effect sizes (log2FC)?
        ├── effect_size.R              post 5's 8 splits x 4 size-factor treatments, five DE call rules
        ├── effect_size_results.csv    genes called per split, treatment and rule
        ├── log2FC_shift.csv           log2FC shift per split when the rRNA feature is ignored
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
| 7 | does normalization only move genes across padj < 0.05, or also shift effect sizes (log2FC)? | `scripts/post7_effect_size/` |
| 8 (next) | how can 31% rRNA become 0.006% in a different counting approach? | `scripts/post8_rRNA_recount/` (planned) |

## Order of use
`sra_download` -> `reference` -> `post2_rRNA_correlation` (trim/STAR/count, then rRNA count, then figure). `post3_compositional_shift` runs on its own. `post4_rRNA_ignored` and `post5_normalization_false_DE` both only need `post2`'s coding-gene counts (`work/STAR_both/fC/`) and `data/paper_TableS2C_rRNA.csv`; post 5 additionally needs DESeq2. `post6_split_structure` and `post7_effect_size` need the same counts (post 6 also ggrepel, a few minutes with `THREADS`; post 7 also apeglm, about 15 minutes).

The teaching arc across posts 3-8: post 3 shows correlation can miss a dominant minority of
features; post 4 shows omitting those features changes relative composition and normalization;
post 5 shows normalization choice can alter downstream DE calls (padj < 0.05); post 6 looks at
the one split that gives hundreds of DE genes under every normalization, i.e. structure among
"identical" replicates; post 7 asks whether normalization also shifts effect sizes (log2FC);
post 8 asks why a published 31.21% rRNA signal becomes 0.006% in a modern recount.

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

**Post 7: "Does normalization only move genes across the padj < 0.05 line, or does it also shift
the estimated effect sizes?"** Done: `scripts/post7_effect_size/` (script, figure, README with all
numbers). Same 8 splits and 4 size-factor treatments as post 5, seven DE call rules.
Short version: with total-count scaling, ignoring the rRNA feature shifts every gene's log2FC by
about 0.1 log2 (median-of-ratios: 0), so raw `|log2FC|` cutoff counts change in all 8 splits. Raw
cutoffs alone are mostly noise from genes with a few reads (2,093-2,494 genes pass `|log2FC| > 1`
in splits with no designed difference); shrunken log2FC (apeglm), `padj & |log2FC|` and
`lfcThreshold` give 0-140 genes and keep split 6 as the outlier.
Open follow-ups:
- All 60 rep-6 splits instead of 8; ashr and normal shrinkage next to apeglm (ashr not installed).
- Same rules on post 6's all-70-split table: which rule best separates split-6-like structure
  from noise?

## Next post

**Post 8: "How can a sample reported to contain 31% rRNA become 0.006% in a different counting
approach?"** Planned folder: `scripts/post8_rRNA_recount/` (own script + README + figure).
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

**Post 8 to-do: strandedness (`featureCounts -s`) as one of the counting choices.**
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
3. **rRNA counts with `-s 0`, `-s 1`, `-s 2`** (`count_rRNA.sh`) as one rung of the post 8 ladder:
   paper 31.21% -> current recount (0.006%) -> `-s` variants -> featureCounts `-M` (multimappers)
   -> higher `--outFilterMultimapNmax` -> annotation-independent count (SortMeRNA/BBDuk vs. SILVA).
   For each rung: rep 11's %, and the Spearman correlation of the 14 samples' ranking with Table
   S2C. Expected: `-s` alone explains little of the gap (the main losses are missing rDNA loci and
   discarded multimappers), but a wrong `-s` would make it worse; check, don't assume.
4. **Downstream check: rerun posts 2 and 4-7 on the `-s 2` gene counts** and report whether any
   conclusion changes (post 2 correlation level, post 4's 0.54 log2 shift, post 5's 4-of-8 splits,
   post 6's 689 DE genes and PC1 gradient, post 7's ranges). If nothing material changes, that is
   one sentence in post 8; if something does, it gets its own correction note in the affected
   post's README.

## Investigation backlog

Every result in posts 2-8 depends on choices made upstream. Each axis below is a candidate
sensitivity analysis (or post) asking: does the conclusion survive a different reasonable choice?

**Normalisations**
- Beyond median-of-ratios and total-count scaling (post 5): edgeR TMM, upper-quartile, CPM/TPM,
  and ERCC spike-in based scaling (the runs carry ERCC mix 1/2, see `data/run_list.tsv`).
- Removal of unwanted variation (RUVg with ERCCs as controls, RUVs with replicates), and
  `~ batch` (ExpA/ExpB) in the design vs. not.
- Post 5 done exhaustively: all 60 rep-6 splits instead of 8, ExpB as well, and with rep 11 in.
- Real per-sample rRNA counts (once post 8 has a usable recount) in place of the synthetic feature.

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
- Strandedness (`-s 0/1/2`): planned as part of post 8, see "Post 8 to-do" above.
- HTSeq-count (union / intersection modes) vs. featureCounts.
- Annotation-independent rRNA measurement: SortMeRNA or BBDuk against SILVA/Rfam rRNA,
  FastQ Screen; this does not depend on how many rDNA copies the reference has.

**References**
- TAIR10 (rDNA arrays largely collapsed or missing) vs. newer near-complete Col-0 assemblies
  (e.g. Col-CEN, Naish et al. 2021); check how much of the NOR 45S rDNA arrays each contains.
- TAIR10 plus a full rDNA unit added as an extra contig, so rRNA reads have somewhere to map.
- Toplevel FASTA includes chloroplast and mitochondrion; confirm, and test with/without them.
- Ensembl Plants release (63 vs. older) and the ERCC sequences: are spike-ins in the reference?
