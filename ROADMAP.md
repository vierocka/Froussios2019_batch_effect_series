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
    │
    └── post4_rRNA_ignored/            sensitivity analysis using the paper's published rRNA fractions
        ├── rRNA_ignored.R             panel A: two-component share rescaling; panel B: 8 same-batch 3 vs 3 splits, DE with size factors from a synthetic rRNA-added vs. coding-only matrix
        ├── rRNA_ignored_panelB_results.csv   panel B's underlying numbers
        └── rRNA_ignored.png          figure of the post
```

## Posts
| Post | Topic | Where |
|---|---|---|
| 1 | QC observations from the paper text and supplement | `publication/` |
| 2 | rRNA fraction vs. genome-wide correlation | `scripts/post2_rRNA_correlation/` |
| 3 | 10 extreme genes barely move a correlation | `scripts/post3_compositional_shift/` |
| 4 | sensitivity analysis using the paper's published rRNA fractions | `scripts/post4_rRNA_ignored/` |

## Order of use
`sra_download` -> `reference` -> `post2_rRNA_correlation` (trim/STAR/count, then rRNA count, then figure). `post3_compositional_shift` runs on its own.

## Ideas for future posts

**Post 5 (candidate): "How can 31% rRNA become 0.006%? Annotation and multimapping can make
contamination disappear."** The failed direct rRNA recount (current TAIR10.63 annotation, default
STAR/featureCounts settings) isn't just a methods footnote — it's arguably a stronger, more
dramatic result than post 4's synthetic-perturbation analysis, and could stand as its own post.

- Paper (Table S2C): replicate 11 has 31.21% rRNA, ranks **highest** of the 14 ExpA/ExpB samples.
- Our own current-annotation/default-counting recount: 0.006%, ranks **lowest**.
- Not a quantitative underestimate — a complete failure to preserve sample ordering (rep 11 flips
  from most- to least-contaminated).
- Demonstrates "I counted the annotated rRNA genes" is not a valid rRNA QC measurement when: the
  annotation contains too few rRNA loci (rDNA is highly repetitive, most true copies are
  collapsed/unassembled in the reference); the aligner excludes reads mapping to many loci
  (`--outFilterMultimapNmax`); and the counting step excludes multimappers by default (no `-M`).
- Source material already exists in the wider investigation's notes (`NOTES.md`,
  `froussios2019_paperrepro_04_rRNA_multimapper_ramses.sh` for the lenient-multimap rerun that
  tests how much of the gap that closes) — would need adapting into this repo's teaching-post
  style (own script + README + figure), not written from scratch.
