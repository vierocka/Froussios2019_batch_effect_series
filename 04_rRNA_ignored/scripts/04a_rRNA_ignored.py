"""Python equivalent of 04a_rRNA_ignored.R.

What does an omitted feature do to normalization? Run from the repository root.
Output: 04_rRNA_ignored/results/python/.
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
OUT = Path("04_rRNA_ignored/results/python")
OUT.mkdir(parents=True, exist_ok=True)


def read_count(run):
    d = pd.read_csv(WORK / "STAR_both/fC" / f"{run}.count", sep="\t", comment="#")
    return pd.Series(d.iloc[:, -1].values, index=d["Geneid"])


# coding-gene counts (featureCounts -t gene has no rRNA rows) for the 14 ExpA/ExpB samples
smp = pd.read_csv("data/replicate_map.tsv", sep="\t")
smp = smp[smp.paper_experiment != "ExpC"]
pct = pd.read_csv("data/paper_TableS2C_rRNA.csv").set_index("replicate")["pct_rRNA"]
smp["p"] = smp.paper_replicate.map(pct) / 100
coding = pd.DataFrame({r: read_count(r) for r in smp.run})
coding = coding[coding.sum(axis=1) > 0]

# SYNTHETIC rRNA total, not a recovered read count: rrna = coding_total * p / (1 - p) solves
# R / (C + R) = p by construction (two-component model: counted coding genes plus rRNA). The
# paper's published fraction (Table S2C) is used because the direct rRNA recount
# (00d_count_rRNA.sh) undercounts rRNA severely.
s = smp.run[smp.paper_replicate == 11].iloc[0]
p11 = smp.p[smp.paper_replicate == 11].iloc[0]
rrna_s = round(coding[s].sum() * p11 / (1 - p11))

# One library (replicate 11): log2 share per coding gene with the synthetic rRNA total in the
# denominator and with rRNA ignored. Raw counts never change; their share shifts by log2(1/(1-p)).
x = coding[s][coding[s] > 0].astype(float)
with_rrna = np.log2(x / (x.sum() + rrna_s))
ignored = np.log2(x / x.sum())
grid = np.linspace(min(with_rrna.min(), ignored.min()) - 1, max(with_rrna.max(), ignored.max()) + 1, 512)

fig, ax = plt.subplots(figsize=(7.5, 6.5))
ax.plot(grid, gaussian_kde(with_rrna)(grid), color="#D95F0E", lw=2,
        label=f"rRNA in library total (rRNA = {100 * p11:.0f}%)")
ax.plot(grid, gaussian_kde(ignored)(grid), color="#2C7FB8", lw=2, label="rRNA ignored (coding genes only)")
ax.set(xlabel="log2 share of library, per coding gene", ylabel="density")
fig.suptitle(f"One library (rep 11, {100 * p11:.0f}% rRNA, synthetic):\n"
             f"all coding-gene library shares shift by {np.log2(1 / (1 - p11)):.2f} log2",
             fontweight="bold", x=0.02, ha="left")
ax.set_title(f"A gene worth 1% of the complete library becomes {1 / (1 - p11):.2f}% of the coding-only total\n"
             f"(its share rises {100 * (1 / (1 - p11) - 1):.0f}%; it did not gain a single read)",
             fontsize=9, loc="left")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), frameon=False)
fig.tight_layout()
fig.savefig(OUT / "figure04a_rRNA_ignored.png", dpi=300, facecolor="white")
print(f"rep 11: p = {p11:.4f}, shift {np.log2(1 / (1 - p11)):.3f} log2; written {OUT}")
