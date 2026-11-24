# 01 The paper and the dataset

| | |
|---|---|
| **Aim** | Describe the dataset and the QC observations that can be made from the paper and its supplement alone. |
| **Input** | Froussios et al. 2019, *Bioinformatics* 35(18):3372 ([`publication/`](../publication/)), Supplementary Table S2. |
| **Scripts** | none |
| **Assumptions** | 17 biological replicates of *A. thaliana* Col-0 in three experiments (ExpA 1-7, ExpB 8-14, ExpC 15-17), one growth condition, so any group difference between them is undesigned. Replicate numbers follow the paper ([`data/replicate_map.tsv`](../data/replicate_map.tsv)), not ENA sample titles. |
| **Statistical approach** | None (reading the paper). |
| **Output** | Post 1: [`post01.txt`](results/post01.txt). Post 6.5 (dataset reminder): [`post06.5_dataset.txt`](results/post06.5_dataset.txt). |
| **Conclusion** | The deposited material has no ERCC analysis; replicate 11 (31.21% rRNA) was excluded for poor correlation while replicate 6 (23.7% rRNA) was kept; no read trimming is described. None of this shows the paper's conclusions are wrong; it shows that "high-quality dataset" is not one property. |
