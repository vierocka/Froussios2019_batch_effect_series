"""Python equivalent of 08c_DE_with_without_TIN.R (pydeseq2 instead of DESeq2).

False DE in all ExpA splits (7 replicates: 35 3-vs-4 splits), without (~ group) and with TIN as a covariate (~ TIN + group).
Run from the repository root after 08b_degradation_metrics.py. Output: 08_hidden_degradation/results/python/.
The 5 splits shown in the figures are the ones drawn by the R script (read from its table08d), so both
versions show the same splits. EXCLUDE_REPS (e.g. "6") leaves replicates out (run 08b with the same setting
first); with 6 replicates the splits are 3 vs. 3 and the output goes to results/without_rep<N>/python/.
"""
import itertools
import os
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
CPUS = int(os.environ.get("THREADS", "4"))
EXCL = [int(x) for x in os.environ.get("EXCLUDE_REPS", "").split(",") if x]
RDIR = Path("08_hidden_degradation/results") / ("without_rep" + "_".join(map(str, EXCL)) if EXCL else "")
OUT = RDIR / "python"
ORANGE, BLUE = "#D95F0E", "#2C7FB8"
FIVE = pd.read_csv(RDIR / "table08d_DE_5_random_splits.csv").group1.tolist()   # splits drawn by the R script

met = pd.read_csv(OUT / "table08a_degradation_metrics.csv")


def read_fc(run):
    d = pd.read_csv(WORK / "STAR_both/fC" / f"{run}.count", sep="\t", comment="#")
    return d.set_index("Geneid").iloc[:, -1]


counts = pd.DataFrame({f"r{rep}": read_fc(run) for run, rep in zip(met.run, met.replicate)})
counts = counts[counts.sum(axis=1) > 0]
tin = pd.Series(met.TIN_median.values, index=[f"r{r}" for r in met.replicate])


def n_de(g1, g2, with_tin):
    ids = list(g1) + list(g2)
    meta = pd.DataFrame({"group": ["g1"] * len(g1) + ["g2"] * len(g2),
                         "TIN": (tin[ids] - tin[ids].mean()).values}, index=ids)
    dds = DeseqDataSet(counts=counts[ids].T, metadata=meta, design="~TIN + group" if with_tin else "~group",
                       quiet=True, n_cpus=CPUS)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "g2", "g1"], quiet=True, n_cpus=CPUS)
    st.summary()
    return int((st.results_df.padj < 0.05).sum())


reps = list(counts.columns)
K = len(reps) // 2
g1s = list(itertools.combinations(reps, K))
if 2 * K == len(reps):
    g1s = [g for g in g1s if reps[0] in g]          # drop mirror images
DESIGN = f"{K}-vs-{len(reps) - K}"
rows = []
for k, g1 in enumerate(g1s, 1):
    g2 = [r for r in reps if r not in g1]
    rows.append(dict(split=k, group1=",".join(r[1:] for r in g1), group2=",".join(r[1:] for r in g2),
                     rep6_in="not included" if "r6" not in reps else "group of 3" if "r6" in g1 else "group of 4",
                     TIN_group1=tin[list(g1)].mean(), TIN_group2=tin[g2].mean(),
                     abs_dTIN=abs(tin[list(g1)].mean() - tin[g2].mean()),
                     DE_without_TIN=n_de(g1, g2, False), DE_with_TIN=n_de(g1, g2, True)))
    print(f"  {k}/{len(g1s)} {rows[-1]['group1']} vs {rows[-1]['group2']}: {rows[-1]['DE_without_TIN']} / {rows[-1]['DE_with_TIN']}",
          flush=True)
res = pd.DataFrame(rows)
res.to_csv(OUT / "table08c_DE_all_splits.csv", index=False)

five = res[res.group1.isin(FIVE)].copy()
five["label"] = five.group1 + " vs. " + five.group2
five.to_csv(OUT / "table08d_DE_5_random_splits.csv", index=False)
print(five[["label", "rep6_in", "abs_dTIN", "DE_without_TIN", "DE_with_TIN"]].to_string(index=False))
print("\nall splits, median (mean) DE genes:")
print(res.groupby("rep6_in")[["DE_without_TIN", "DE_with_TIN"]].agg(lambda x: f"{np.median(x):g} ({x.mean():.0f})"))
print(f"Spearman rho, DE vs. |dTIN|: without TIN {spearmanr(res.DE_without_TIN, res.abs_dTIN).statistic:.2f}, "
      f"with TIN {spearmanr(res.DE_with_TIN, res.abs_dTIN).statistic:.2f}")


def plot_splits(col, title, fname):
    fig, axes = plt.subplots(1, 5, figsize=(13, 5), sharey=True)
    rng = np.random.default_rng(1)
    for ax, (_, row) in zip(axes, five.iterrows()):
        for gi, grp in enumerate([row.group1, row.group2]):
            for r in grp.split(","):
                x = gi + rng.uniform(-0.18, 0.18)
                ax.scatter(x, tin[f"r{r}"], s=260, color=ORANGE if r == "6" else BLUE)
                ax.text(x, tin[f"r{r}"], r, color="white", ha="center", va="center", fontsize=9, fontweight="bold")
        ax.set_xticks([0, 1], ["group 1", "group 2"]); ax.set_xlim(-0.6, 1.6)
        ax.set_title(f"{row.label}\n{row[col]} DE genes", fontsize=10, fontweight="bold")
    axes[0].set_ylabel("sample TIN (median)")
    fig.suptitle(title, fontweight="bold", x=0.01, ha="left")
    left = f", replicate {', '.join(map(str, EXCL))} left out" if EXCL else ""
    fig.text(0.01, 0.9, f"5 random {DESIGN} splits of {len(reps)} ExpA replicates (same genotype and batch){left}; "
             "DE = padj < 0.05", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.9)); fig.savefig(OUT / fname, dpi=300); plt.close(fig)


plot_splits("DE_without_TIN", "DE without TIN in the model (~ group)", "figure08d_splits_without_TIN.png")
plot_splits("DE_with_TIN", "DE with TIN as a covariate (~ TIN + group)", "figure08e_splits_with_TIN.png")

fig, axes = plt.subplots(1, 2, figsize=(9, 5), sharey=True)
for ax, (col, model) in zip(axes, [("DE_without_TIN", "~ group"), ("DE_with_TIN", "~ TIN + group")]):
    for where, c in [("group of 3", ORANGE), ("group of 4", BLUE), ("not included", BLUE)]:
        s = res[res.rep6_in == where]
        if len(s):
            ax.scatter(s.abs_dTIN, s[col] + 1, color=c, alpha=0.8, label=f"replicate 6 in the {where}")
    ax.set_yscale("log"); ax.set_title(model); ax.set_xlabel("|difference in mean TIN between the groups|")
axes[0].set_ylabel("DE genes + 1 (log scale)")
if not EXCL:
    axes[0].legend(fontsize=8)
fig.suptitle(f"All {len(res)} {DESIGN} splits of ExpA" + (f" without replicate {', '.join(map(str, EXCL))}" if EXCL else ""), fontweight="bold", x=0.01, ha="left")
fig.tight_layout(); fig.savefig(OUT / "figure08f_DE_vs_TIN_imbalance.png", dpi=300); plt.close(fig)
