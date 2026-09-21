# Roadmap

Map of the repository. New files are added post by post; this file is updated with them.

```
Froussios2019_batch_effect_series/
├── README.md                          purpose, run order, known limits
├── ROADMAP.md                         this file
├── config.sh                          shared paths and thread count (override with environment variables)
├── .gitignore                         keeps raw data, BAMs, work/ and AED files out of git
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
    │   └── rRNA_correlation_blindspot.png   figure of the post
    │
    └── post3_compositional_shift/     teaching simulation
        ├── simulation.R               two correlated log2 vectors; 10 extreme genes changed, 1000 loops
        ├── compositional_shift.png    figure of the post
        └── README.md                  explanation of the simulation and what it means for QC
```

## Posts
| Post | Topic | Where |
|---|---|---|
| 1 | QC observations from the paper text and supplement | `publication/` |
| 2 | rRNA fraction vs. genome-wide correlation | `scripts/post2_rRNA_correlation/` |
| 3 | 10 extreme genes barely move a correlation | `scripts/post3_compositional_shift/` |

## Order of use
`sra_download` -> `reference` -> `post2_rRNA_correlation` (trim/STAR/count, then rRNA count, then figure). `post3_compositional_shift` runs on its own.
