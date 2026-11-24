# 00 Preprocessing: reads to counts

| | |
|---|---|
| **Aim** | Turn the 17 public runs into gene counts and rRNA counts; check which strand setting fits the libraries. |
| **Input** | 17 paired-end runs, ENA ERP021226 ([`data/run_list.tsv`](../data/run_list.tsv)); TAIR10 genome and annotation, Ensembl Plants 63. |
| **Scripts** | All HPC steps (00a-00e) also as one Snakemake workflow: [`workflow/Snakefile`](../workflow/Snakefile) ([subway map](../workflow/subway_map.png)). · [`00a_download_sra.sh`](scripts/00a_download_sra.sh) · [`00b_prep_reference.sh`](scripts/00b_prep_reference.sh) · [`00c_trim_star_featurecounts.sh`](scripts/00c_trim_star_featurecounts.sh) · [`00d_count_rRNA.sh`](scripts/00d_count_rRNA.sh) · [`00e_count_strandedness.sh`](scripts/00e_count_strandedness.sh) · [`00f_strandedness_comparison.R`](scripts/00f_strandedness_comparison.R) / [`.py`](scripts/00f_strandedness_comparison.py) |
| **Assumptions** | TruSeq3-PE adapters. Gene counts: `featureCounts -t gene -p`, unstranded (`-s 0`), each mate counted (no `--countReadPairs`); used by all topics. rRNA loci are typed `ncRNA_gene`, so they are counted separately (`00d`, paper's STAR settings, `-s 2`, 4 annotated loci, no multimappers). |
| **Statistical approach** | None (read processing). `00f`: assigned fraction per strand setting, per-gene log2 ratio `-s 2` / `-s 0`, topic 02's mean correlation on both. |
| **Output** | `pipeline_data/` (not in git): `fastq/`, `ref/`, `STAR_both/` (BAMs, gene counts in `fC/`), `STAR_Fonly/`, `STAR_Ronly/` (reads whose mate was removed by trimming), `rRNA_qc/counts/` (rRNA counts). `00e`/`00f`: `pipeline_data/STAR_both/fC_s0, fC_s1, fC_s2, fC_s2_pairs/`, `table00a_strandedness_assignment.csv`, `table00b_strandedness_gene_ratio.csv`, `table00c_strandedness_topic02_correlation.csv`, `figure00a_strandedness.png`. |
| **Conclusion** | The direct rRNA count undercounts severely (replicate 11: about 0.006% vs. 31.21% published), so topics 04-07 use the published fractions. Strandedness (`-s 2`) and fragment counting are open: `00e`/`00f` need the BAMs and have not been run. |
