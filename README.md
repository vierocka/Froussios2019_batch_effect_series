# Froussios et al. 2019 – batch-effect series

Reanalysis of the Arabidopsis Col-0 RNA-seq benchmarking dataset of Froussios et al. 2019
(*Bioinformatics* 35(18):3372, doi:10.1093/bioinformatics/btz089; ArrayExpress E-MTAB-5446, ENA ERP021226).
Each LinkedIn post of the series has its own folder; files are added post by post.

```
publication/                  paper PDF and its Supplementary Table S2
data/                         run list, corrected replicate numbering, paper's rRNA fractions (Table S2C)
config.sh                     shared paths/threads (override with environment variables)
scripts/sra_download/         download the 17 runs with sra-tools
scripts/reference/            TAIR10 (Ensembl Plants 63) genome, annotation, STAR index
scripts/post2_rRNA_correlation/   trimming, STAR, featureCounts, rRNA counts, figure
scripts/post3_compositional_shift/   simulation for teaching
```

| Post | Topic | Code |
|---|---|---|
| 1 | QC observations from the paper text and supplement | none |
| 2 | rRNA fraction vs. genome-wide correlation | `scripts/post2_rRNA_correlation/` |
| 3 | 10 extreme genes barely move a correlation (simulation) | `scripts/post3_compositional_shift/` |

## Run order
```bash
bash scripts/sra_download/download_sra.sh
bash scripts/reference/prep_reference.sh
ADAPTERS=/path/to/TruSeq3-PE.fa bash scripts/post2_rRNA_correlation/trim_star_featurecounts.sh
bash scripts/post2_rRNA_correlation/count_rRNA.sh
Rscript scripts/post2_rRNA_correlation/rRNA_correlation_blindspot.R
Rscript scripts/post3_compositional_shift/simulation.R
```
Requires sra-tools, Trimmomatic, STAR, Subread (featureCounts) on `PATH`; R with DESeq2, ggplot2, patchwork, MASS.
Long steps skip outputs that already exist and can be resubmitted on an HPC scheduler.

## Replicate numbering
ENA `sample_title` numbering differs from the paper's replicate numbering. Use `data/replicate_map.tsv`
(paper's replicate 11 = ERR1811891; replicate 6 = ERR1811898).

## Known limits
- `featureCounts -t gene` does not count the rRNA loci (typed `ncRNA_gene`); they are counted separately in `count_rRNA.sh`.
- The rRNA recount uses the current annotation (4 rRNA loci) and `--outFilterMultimapNmax 2`, so it undercounts relative to the paper's Table S2C.
