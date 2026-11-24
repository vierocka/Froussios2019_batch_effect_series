"""Python equivalent of 00f_strandedness_comparison.R.

Gene counts with featureCounts -s 0 (fC/), -s 1 and -s 2, and -s 2 with --countReadPairs. Run from
the repository root after 00e_count_strandedness.sh (pipeline_data/STAR_both/fC_*/). pydeseq2 has no rlog,
so step 4 uses the blind vst. Output: 00_preprocessing/results/python/.
1. assigned fraction per sample and variant (which strand setting fits the library?)
2. check: fC_s0 must equal fC/
3. gene level: how many genes change between -s 0 and -s 2, and which ones
4. topic 02's result (mean correlation with the other 13 samples) on -s 0 vs. -s 2 counts
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
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
OUT = Path("00_preprocessing/results/python")
OUT.mkdir(parents=True, exist_ok=True)
cnt = lambda v: WORK / "STAR_both" / v
variants = [v for v in ("fC_s0", "fC_s1", "fC_s2", "fC_s2_pairs") if cnt(v).is_dir()]
assert {"fC_s0", "fC_s2"} <= set(variants), "run 00e_count_strandedness.sh first"

smp = pd.read_csv("data/replicate_map.tsv", sep="\t")
pct = pd.read_csv("data/paper_TableS2C_rRNA.csv").set_index("replicate")["pct_rRNA"]
smp["pct_rRNA"] = smp.paper_replicate.map(pct)


def mat(v, runs):
    out = {}
    for r in runs:
        d = pd.read_csv(cnt(v) / f"{r}.count", sep="\t", comment="#")
        out[r] = pd.Series(d.iloc[:, -1].values, index=d["Geneid"])
    return pd.DataFrame(out)


# 1. featureCounts summaries
rows = []
for v in variants:
    for r, rep, exp in zip(smp.run, smp.paper_replicate, smp.paper_experiment):
        s = pd.read_csv(cnt(v) / f"{r}.count.summary", sep="\t")
        n = pd.Series(s.iloc[:, 1].values, index=s.iloc[:, 0])
        rows.append(dict(variant=v, run=r, replicate=rep, experiment=exp, total=n.sum(), assigned=n["Assigned"],
                         no_features=n["Unassigned_NoFeatures"], ambiguity=n["Unassigned_Ambiguity"]))
summ = pd.DataFrame(rows)
summ["pct_assigned"] = 100 * summ.assigned / summ.total
summ["pct_ambiguity"] = 100 * summ.ambiguity / summ.total
summ.to_csv(OUT / "table00a_strandedness_assignment.csv", index=False)
print("assigned % per variant (median, min, max over the 17 samples):")
print(summ.groupby("variant").pct_assigned.agg(["median", "min", "max"]).round(1).to_string())
print("ambiguity % per variant (median):")
print(summ.groupby("variant").pct_ambiguity.median().round(2).to_string())
a = summ.pivot(index="run", columns="variant", values="assigned")
if "fC_s1" in variants:
    q = a.fC_s2 / a.fC_s1
    print(f"assigned -s 2 / -s 1: median {q.median():.1f} (range {q.min():.1f}-{q.max():.1f}); "
          "> 1 means reverse-stranded (read 1 antisense)")
if "fC_s2_pairs" in variants:
    q = a.fC_s2_pairs / a.fC_s2
    print(f"assigned fragments (--countReadPairs) / assigned reads, -s 2: median {q.median():.3f} "
          f"(range {q.min():.3f}-{q.max():.3f})")

# 2. fC_s0 must reproduce fC/, the counts used by topics 02-07
pub, s0, s2 = mat("fC", smp.run), mat("fC_s0", smp.run), mat("fC_s2", smp.run)
same = pub.equals(s0)
print(f"\nfC_s0 identical to fC/: {same}")
if not same:
    warnings.warn("fC_s0 differs from fC/: check featureCounts version and BAMs before comparing")

# 3. gene level, 14 ExpA/ExpB samples (as in topic 02)
ab = list(smp.run[smp.paper_experiment != "ExpC"])
keep = (s0[ab].sum(axis=1) > 0) | (s2[ab].sum(axis=1) > 0)
a0, a2 = s0.loc[keep, ab], s2.loc[keep, ab]
lr = np.log2((a2 + 1) / (a0 + 1))
e10 = a0 >= 10
per_sample = pd.DataFrame({"replicate": smp.set_index("run").paper_replicate[ab], "total_s0": a0.sum(),
                           "total_s2": a2.sum(), "genes_expressed_s0": e10.sum(),
                           "genes_change_10pct": (e10 & (lr.abs() > np.log2(1.1))).sum(),
                           "genes_change_2fold": (e10 & (lr.abs() > 1)).sum(), "genes_gain_2fold": (e10 & (lr > 1)).sum()})
print("\nper sample, genes with >= 10 reads at -s 0 that change between -s 0 and -s 2:")
print(per_sample.to_string())
m0, m2 = a0.mean(axis=1), a2.mean(axis=1)
g = np.log2((m2 + 1) / (m0 + 1))
gene_ratio = pd.DataFrame({"gene": a0.index, "mean_s0": m0, "mean_s2": m2, "log2_s2_vs_s0": g}).sort_values("log2_s2_vs_s0")
gene_ratio.to_csv(OUT / "table00b_strandedness_gene_ratio.csv", index=False)
e = m0 >= 10
print(f"\ngenes (mean >= 10 reads at -s 0): {e.sum()}; lose > 10%: {(e & (g < -np.log2(1.1))).sum()}; "
      f"lose > 2-fold: {(e & (g < -1)).sum()}; gain > 10%: {(e & (g > np.log2(1.1))).sum()}; gain > 2-fold: {(e & (g > 1)).sum()}")
print("genes losing most (antisense reads counted on them when unstranded):")
print(gene_ratio[gene_ratio.mean_s0 >= 10].head(10).to_string(index=False))
print("genes gaining most (unstranded counting discarded their reads as ambiguous with an antisense neighbour):")
print(gene_ratio[(gene_ratio.mean_s2 >= 10) & (gene_ratio.log2_s2_vs_s0 > 0)]
      .sort_values("log2_s2_vs_s0", ascending=False).head(10).to_string(index=False))

# 4. topic 02's panel B (genes only) on -s 0 vs. -s 2
o = smp[smp.paper_experiment != "ExpC"].sort_values("pct_rRNA")


def mean_cor(m):
    m = m.loc[m[o.run].sum(axis=1) > 0, o.run].astype(int)
    dds = DeseqDataSet(counts=m.T, metadata=pd.DataFrame({"x": ["a"] * m.shape[1]}, index=m.columns),
                       design="~1", quiet=True)
    dds.vst(use_design=False)
    cm = np.corrcoef(dds.layers["vst_counts"])
    np.fill_diagonal(cm, np.nan)
    return np.nanmean(cm, axis=1)


cors = pd.DataFrame({"run": o.run.values, "replicate": o.paper_replicate.values,
                     "experiment": o.paper_experiment.values, "pct_rRNA": o.pct_rRNA.values,
                     "s0": mean_cor(s0), "s2": mean_cor(s2)})
cors.to_csv(OUT / "table00c_strandedness_topic02_correlation.csv", index=False)
print(f"\ntopic 02 mean vst correlation: -s 0 {cors.s0.min():.4f}-{cors.s0.max():.4f}, -s 2 {cors.s2.min():.4f}-"
      f"{cors.s2.max():.4f}; Spearman (rank of samples) {spearmanr(cors.s0, cors.s2)[0]:.2f}")
print(f"Spearman with rRNA %: -s 0 {spearmanr(cors.s0, cors.pct_rRNA)[0]:.2f}, -s 2 {spearmanr(cors.s2, cors.pct_rRNA)[0]:.2f}")

# figure
lab = {"fC_s0": "-s 0 (unstranded, fC/)", "fC_s1": "-s 1", "fC_s2": "-s 2"}
pal = {"fC_s0": "#333333", "fC_s1": "#66A61E", "fC_s2": "#D95F0E"}
fig = plt.figure(figsize=(13, 10))
ax = fig.add_axes([0.07, 0.58, 0.4, 0.33])
for v in [v for v in lab if v in variants]:
    d = summ[summ.variant == v].sort_values("replicate")
    ax.plot(d.replicate, d.pct_assigned, marker="o", color=pal[v], label=lab[v])
ax.set_xticks(range(1, 18))
ax.set(xlabel="replicate (paper numbering)", ylabel="% of reads assigned to a gene")
ax.set_title("A. Which strand setting fits the library?", loc="left", fontweight="bold")
ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.15))
ax = fig.add_axes([0.56, 0.58, 0.4, 0.33])
d = gene_ratio[(gene_ratio.mean_s0 >= 1) | (gene_ratio.mean_s2 >= 1)]
ax.scatter(np.log2(d.mean_s0 + 1), d.log2_s2_vs_s0, s=1, alpha=0.3, color="0.2")
for yv, ls in ((-1, ":"), (0, "-"), (1, ":")):
    ax.axhline(yv, ls=ls, color="0.5")
ax.set(xlabel="log2 mean count at -s 0 (+1), 14 samples", ylabel="log2 (-s 2 / -s 0)")
ax.set_title(f"B. Per gene: stranded vs. unstranded counts\n{e.sum()} genes with >= 10 reads; "
             f"{(e & (g.abs() > 1)).sum()} change more than 2-fold", loc="left", fontsize=10)
ax = fig.add_axes([0.07, 0.08, 0.89, 0.36])
x = np.arange(len(cors))
ax.plot(x, cors.s0, marker="o", color="#333333", label="-s 0 (fC/)")
ax.plot(x, cors.s2, marker="o", color="#D95F0E", label="-s 2")
ax.set_xticks(x, cors.replicate)
ax.set(xlabel="replicate, ordered by rRNA % (as in topic 02)", ylabel="mean vst correlation\nwith the other 13 samples")
ax.set_title("C. Topic 02 panel B on unstranded vs. stranded counts", loc="left", fontweight="bold")
ax.legend(frameon=False)
fig.suptitle("Froussios et al. 2019, TruSeq Stranded Total RNA: featureCounts -s 0 vs. -s 1 vs. -s 2",
             fontweight="bold", x=0.02, ha="left")
fig.savefig(OUT / "figure00a_strandedness.png", dpi=300, facecolor="white")
print("written: figure00a, table00a-c in", OUT)
