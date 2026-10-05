# Froussios et al. 2019 – batch-effect series

Purpose. 
This repository is not a critique of the original study. It uses that study's openly available, well-documented data (17 biological replicates of one genotype) as a teaching case for batch effects in RNA-seq: what they are, why caution is needed, how to detect them, how to handle them, and what they change downstream. The original paper had a different aim (benchmarking differential-expression tools).

Reanalysis of the Arabidopsis Col-0 RNA-seq benchmarking dataset of Froussios et al. 2019
(*Bioinformatics* 35(18):3372, doi:10.1093/bioinformatics/btz089; ArrayExpress E-MTAB-5446, ENA ERP021226).
Paper (open access, CC BY 4.0): https://doi.org/10.1093/bioinformatics/btz089 ; also on ResearchGate: https://www.researchgate.net/publication/335922274
Each LinkedIn post of the series has its own folder; files are added post by post.

```
publication/                            paper PDF and its Supplementary Table S2
data/                                   run list, corrected replicate numbering, paper's rRNA fractions (Table S2C)
config.sh                               shared paths/threads (override with environment variables)
scripts/sra_download/                   download the 17 runs with sra-tools
scripts/reference/                      TAIR10 (Ensembl Plants 63) genome, annotation, STAR index
scripts/post2_rRNA_correlation/         trimming, STAR, featureCounts, rRNA counts; rRNA % vs. correlation
scripts/post3_compositional_shift/      simulation: 10 extreme genes vs. a correlation
scripts/post4_rRNA_ignored/             what does an omitted feature do to library shares?
scripts/post5_normalization_false_DE/   can normalization choice change false-DE counts?
scripts/post6_split_structure/          why does one same-batch 3-vs-3 split give 689 DE genes?
scripts/post7_effect_size/              does normalization also shift effect sizes (log2FC)?
```

| Post | Topic | Code |
|---|---|---|
| 1 | QC observations from the paper text and supplement | none |
| 2 | rRNA fraction vs. genome-wide correlation | `scripts/post2_rRNA_correlation/` |
| 3 | 10 extreme genes barely move a correlation (simulation) | `scripts/post3_compositional_shift/` |
| 4 | what does an omitted feature do to normalization? | `scripts/post4_rRNA_ignored/` |
| 5 | can normalization choice change false-DE counts? | `scripts/post5_normalization_false_DE/` |
| 6 | why does one same-batch 3-vs-3 split give 689 DE genes? | `scripts/post6_split_structure/` |
| 7 | does normalization also shift effect sizes (log2FC)? | `scripts/post7_effect_size/` |

Each post folder has a README with the question, method, results and limits. Planned posts
(next: how a reported 31% rRNA becomes 0.006% in a different counting approach) and open
questions are in `ROADMAP.md`.

## Run order
```bash
bash scripts/sra_download/download_sra.sh
bash scripts/reference/prep_reference.sh
ADAPTERS=/path/to/TruSeq3-PE.fa bash scripts/post2_rRNA_correlation/trim_star_featurecounts.sh
bash scripts/post2_rRNA_correlation/count_rRNA.sh
Rscript scripts/post2_rRNA_correlation/rRNA_correlation_blindspot.R
Rscript scripts/post3_compositional_shift/simulation.R
Rscript scripts/post4_rRNA_ignored/rRNA_ignored.R
Rscript scripts/post5_normalization_false_DE/false_DE_by_normalization.R
Rscript scripts/post6_split_structure/split_structure.R
Rscript scripts/post7_effect_size/effect_size.R
```
Requires sra-tools, Trimmomatic, STAR, Subread (featureCounts) on `PATH`; R with DESeq2, ggplot2, patchwork, MASS, ggrepel (post 6), apeglm (post 7).
Long steps skip outputs that already exist and can be resubmitted on an HPC scheduler.
The R scripts of posts 2-5 run in a minute or two; post 6 takes a few minutes (set `THREADS` for its 70 DESeq2
fits), post 7 about 15 minutes (apeglm shrinkage).

## Generated data
Reads, reference, BAMs and count tables are not in the repository (`work/` is git-ignored); the
scripts regenerate them from public sources into `work/` (or `$WORK_DIR`):

| Needed by | Output | Produced by |
|---|---|---|
| all | `work/fastq/` (17 runs, ENA ERP021226) | `scripts/sra_download/download_sra.sh` (accessions and MD5s in `data/run_list.tsv`) |
| all | `work/ref/` (TAIR10, Ensembl Plants 63, STAR index) | `scripts/reference/prep_reference.sh` |
| posts 2, 4-7 | `work/STAR_both/fC/*.count` (gene counts) | `scripts/post2_rRNA_correlation/trim_star_featurecounts.sh` |
| post 2 | `work/rRNA_qc/counts/*_rRNA.count` (4 rRNA loci) | `scripts/post2_rRNA_correlation/count_rRNA.sh` |

Post 3 needs no data. The figures and result tables (`*.csv`) committed in the post folders were
made from these outputs.

## Replicate numbering
ENA `sample_title` numbering differs from the paper's replicate numbering. Use `data/replicate_map.tsv`
(paper's replicate 11 = ERR1811891; replicate 6 = ERR1811898).

## Known limits
- `featureCounts -t gene` does not count the rRNA loci (typed `ncRNA_gene`); they are counted separately in `count_rRNA.sh`.
- The rRNA recount uses the current annotation (4 rRNA loci) and `--outFilterMultimapNmax 2`, so it undercounts severely relative to the paper's Table S2C (replicate 11: about 0.006% vs. 31.21%). Post 2's "incl. 4 rRNA" correlation line therefore carries almost no rRNA signal (see `scripts/post2_rRNA_correlation/README.md`); posts 4, 5 and 7 use the published fractions instead.
- The libraries are stranded (TruSeq Stranded Total RNA, paper section 2.1), but the gene counts
  were made unstranded (`featureCounts` default `-s 0`); `count_rRNA.sh` uses `-s 2`. Post 8
  will test `-s 0/1/2` and rerun posts 2 and 4-7 on stranded counts (plan in `ROADMAP.md`).
- Posts 4, 5 and 7 add a **synthetic** rRNA feature calibrated to the published percentages; they
  are sensitivity analyses, not reconstructions of the real rRNA counts.
