"""Python equivalent of 02a_rRNA_correlation.R.

rRNA fraction per sample (paper Table S2C) vs. mean correlation with the other samples. Run from
the repository root. Needs pipeline_data/STAR_both/fC/ (00c_trim_star_featurecounts.sh) and
pipeline_data/rRNA_qc/counts/ (00d_count_rRNA.sh). pydeseq2 has no rlog, so the correlation is computed on
the blind vst instead. Output: 02_rRNA_correlation/results/python/.
"""
import os
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet

warnings.filterwarnings("ignore")
WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
OUT = Path("02_rRNA_correlation/results/python")
OUT.mkdir(parents=True, exist_ok=True)


def read_count(f):
    d = pd.read_csv(f, sep="\t", comment="#")
    return pd.Series(d.iloc[:, -1].values, index=d["Geneid"])


# 14 ExpA + ExpB samples; replicate numbers follow the paper, not ENA sample titles
smp = pd.read_csv("data/replicate_map.tsv", sep="\t")
smp = smp[smp.paper_experiment.isin(["ExpA", "ExpB"])]
pct = pd.read_csv("data/paper_TableS2C_rRNA.csv").set_index("replicate")["pct_rRNA"]
smp["pct_rRNA"] = smp.paper_replicate.map(pct)
smp = smp.sort_values("pct_rRNA").reset_index(drop=True)
smp["flag"] = np.select([smp.paper_replicate == 11, smp.paper_replicate == 6], ["rep 11\n(excluded)", "rep 6\n(kept)"], "")

# gene counts (no rRNA rows) plus the 4 rRNA genes from the rRNA count files
rrna_genes = ["AT2G01010", "AT2G01020", "AT3G41768", "AT3G41979"]
gene = pd.DataFrame({r: read_count(WORK / "STAR_both/fC" / f"{r}.count") for r in smp.run})
rrna = pd.DataFrame({r: read_count(WORK / "rRNA_qc/counts" / f"{r}_rRNA.count")[rrna_genes] for r in smp.run})
assert not set(rrna_genes) & set(gene.index)
counts = pd.concat([gene, rrna]).astype(int)
counts = counts[counts.sum(axis=1) > 0]

meta = pd.DataFrame({"x": ["a"] * counts.shape[1]}, index=counts.columns)
dds = DeseqDataSet(counts=counts.T, metadata=meta, design="~1", quiet=True)
dds.vst(use_design=False)
vs = pd.DataFrame(dds.layers["vst_counts"].T, index=counts.index, columns=counts.columns)


def mean_cor(genes):
    cm = np.corrcoef(vs.loc[genes].T.values)
    np.fill_diagonal(cm, np.nan)
    return pd.Series(np.nanmean(cm, axis=1), index=vs.columns)


rng = np.random.default_rng(42)
lab = [f"All {len(vs):,} genes (incl. 4 rRNA)", f"{len(vs) - 4:,} genes (4 rRNA removed)",
       "Random 5,000 genes (avg. of 50 draws)"]
curves = [mean_cor(vs.index), mean_cor(vs.index.difference(rrna_genes)),
          pd.concat([mean_cor(rng.choice(vs.index, 5000, replace=False)) for _ in range(50)], axis=1).mean(axis=1)]

pal = {"ExpA": "#2C7FB8", "ExpB": "#D95F0E"}
x = np.arange(len(smp))
fig, (a, b) = plt.subplots(2, 1, figsize=(8.5, 9.5))
a.bar(x, smp.pct_rRNA, width=0.7, color=smp.paper_experiment.map(pal))
for xi, (v, f) in enumerate(zip(smp.pct_rRNA, smp.flag)):
    a.text(xi, v, f, ha="center", va="bottom", fontsize=8, fontweight="bold")
a.set_xticks(x, smp.run, rotation=45, ha="right", fontsize=8)
a.set_ylim(0, smp.pct_rRNA.max() * 1.25)
a.set_ylabel("% of reads mapping to rRNA\n(paper's Table S2C)")
a.set_title(f"A. rRNA content varies ~{round(smp.pct_rRNA.max() / smp.pct_rRNA.min())}-fold across 14 "
            "'same-condition' replicates", loc="left", fontweight="bold", fontsize=11)
a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in pal.values()], labels=list(pal), title="Experiment",
         loc="upper left", frameon=False)
for curve, l, col in zip(curves, lab, ["#333333", "#E7298A", "#66A61E"]):
    b.plot(x, curve[smp.run].values, marker="o", color=col, label=l)
b.set_xticks(x, smp.run, rotation=45, ha="right", fontsize=8)
b.set_ylim(min(c.min() for c in curves) - 0.002, 1)
b.set(xlabel="sample (ordered by rRNA%, same order as panel A)", ylabel="mean vst correlation\nwith the other 13 samples")
b.set_title("B. ...but genome-wide correlation barely moves either way", loc="left", fontweight="bold", fontsize=11)
b.legend(fontsize=9, frameon=False, ncol=2, loc="lower left")
fig.suptitle("Froussios et al. 2019 (Bioinformatics), Arabidopsis Col-0, ExpA+ExpB (n=14)", fontsize=10,
             style="italic", x=0.02, ha="left")
fig.tight_layout()
fig.savefig(OUT / "figure02a_rRNA_vs_correlation.png", dpi=300, facecolor="white")
print(pd.DataFrame({"run": smp.run, "pct_rRNA": smp.pct_rRNA, "all": curves[0][smp.run].values,
                    "no_rRNA": curves[1][smp.run].values, "random5000": curves[2][smp.run].values}).round(4).to_string(index=False))
