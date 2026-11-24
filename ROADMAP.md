# Roadmap

| Item | Topic | Content | Post |
|---|---|---|---|
| [`00_preprocessing/`](00_preprocessing/) | preprocessing | reads to gene and rRNA counts; strandedness check | – |
| [`01_paper_and_dataset/`](01_paper_and_dataset/) | paper and dataset | QC observations from the paper; the dataset and its design | 1, 6.5 |
| [`02_rRNA_correlation/`](02_rRNA_correlation/) | rRNA vs. correlation | genome-wide correlation does not reveal rRNA content | 2 |
| [`03_compositional_shift/`](03_compositional_shift/) | compositional shift | ten extreme genes barely move a correlation (simulation) | 3 |
| [`04_rRNA_ignored/`](04_rRNA_ignored/) | rRNA ignored | what an omitted feature does to library shares | 4 |
| [`05_normalization_false_DE/`](05_normalization_false_DE/) | normalization and false DE | normalization choice changes false-DE counts; series summary | 5, 5.5 |
| [`06_split_structure/`](06_split_structure/) | split structure | one same-batch 3-vs-3 split gives 689 DE genes | 6A-6D |
| [`07_effect_size/`](07_effect_size/) | effect size | significance vs. effect size as DE rules | 7 |
| [`08_hidden_degradation/`](08_hidden_degradation/) | hidden degradation | TIN, rRNA, singletons and false DE within ExpA | – |
| [`09_batch_equals_condition/`](09_batch_equals_condition/) | batch = condition | confounded 4-vs-4 design; Shiny apps in R and Python | – |
| [`10_mislabelling/`](10_mislabelling/) | mislabelling | one swapped sample per group hides all DE | – |
| [`workflow/`](workflow/) | count-data workflow | Snakemake pipeline from download to counts and TIN; subway map | – |
| [`lib/`](lib/) | shared code | GO enrichment in R and Python | – |
| [`data/`](data/) | metadata | run list, replicate numbering, published rRNA fractions | – |
| [`publication/`](publication/) | paper | Froussios et al. 2019 PDF and Supplementary Table S2 | – |
