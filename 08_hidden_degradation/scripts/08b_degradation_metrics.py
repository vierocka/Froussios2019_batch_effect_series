"""Python equivalent of 08b_degradation_metrics.R.

Degradation metrics of the 7 ExpA replicates: TIN per sample, TIN vs. rRNA %, TIN vs. singleton reads.
Run from the repository root. Output: 08_hidden_degradation/results/python/.
EXCLUDE_REPS (e.g. "6") leaves replicates out; the output then goes to results/without_rep<N>/python/.
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
EXCL = [int(x) for x in os.environ.get("EXCLUDE_REPS", "").split(",") if x]
OUT = Path("08_hidden_degradation/results") / ("without_rep" + "_".join(map(str, EXCL)) if EXCL else "") / "python"
OUT.mkdir(parents=True, exist_ok=True)
ORANGE, BLUE = "#D95F0E", "#2C7FB8"

rmap = pd.read_csv("data/replicate_map.tsv", sep="\t")
rmap = rmap[(rmap.paper_experiment == "ExpA") & ~rmap.paper_replicate.isin(EXCL)].copy()
pct = pd.read_csv("data/paper_TableS2C_rRNA.csv").set_index("replicate")["pct_rRNA"]
rmap["pct_rRNA"] = rmap.paper_replicate.map(pct)


def tin_summary(run):
    s = pd.read_csv(WORK / "TIN" / f"{run}_Aligned.sortedByCoord.out.summary.txt", sep="\t")
    return pd.Series({"TIN_mean": s.iloc[0, 1], "TIN_median": s.iloc[0, 2], "TIN_sd": s.iloc[0, 3]})


def star_input(folder, run):
    for line in open(WORK / folder / f"{run}_Log.final.out"):
        if "Number of input reads" in line:
            return float(line.split("|")[1])


d = rmap[["run", "paper_replicate", "pct_rRNA"]].rename(columns={"paper_replicate": "replicate"}).set_index("run")
d = d.join(pd.DataFrame({r: tin_summary(r) for r in d.index}).T)
d["pairs"] = [star_input("STAR_both", r) for r in d.index]
d["fwd_only"] = [star_input("STAR_Fonly", r) for r in d.index]
d["rev_only"] = [star_input("STAR_Ronly", r) for r in d.index]
tot = d.pairs + d.fwd_only + d.rev_only
d["pct_singletons"] = 100 * (d.fwd_only + d.rev_only) / tot
d["pct_fwd_only"] = 100 * d.fwd_only / tot
d = d.sort_values("replicate").reset_index()
d.to_csv(OUT / "table08a_degradation_metrics.csv", index=False)
print(d.to_string(index=False))


def cors(x, y):
    s = spearmanr(x, y); p = pearsonr(x, y)
    return dict(spearman_rho=s.statistic, spearman_p=s.pvalue, pearson_r=p.statistic, pearson_p=p.pvalue)


no6 = d.replicate != 6
ct = pd.DataFrame([dict(comparison="TIN median vs. rRNA %", **cors(d.TIN_median, d.pct_rRNA))]
    + ([dict(comparison="TIN median vs. rRNA % without replicate 6", **cors(d.TIN_median[no6], d.pct_rRNA[no6]))]
       if (~no6).any() else [])
    + [dict(comparison="TIN median vs. singleton %", **cors(d.TIN_median, d.pct_singletons)),
    dict(comparison="rRNA % vs. singleton %", **cors(d.pct_rRNA, d.pct_singletons))])
ct.to_csv(OUT / "table08b_correlations.csv", index=False)
print(ct.round(3).to_string(index=False))

# figure a: per-transcript TIN per sample
fig, ax = plt.subplots(figsize=(7, 5.5))
vals = [pd.read_csv(WORK / "TIN" / f"{r}_Aligned.sortedByCoord.out.tin.xls", sep="\t")["TIN"] for r in d.run]
vals = [v[v > 0] for v in vals]
bp = ax.boxplot(vals, tick_labels=d.replicate, patch_artist=True, widths=0.6, flierprops=dict(markersize=1, alpha=0.2))
for patch, rep in zip(bp["boxes"], d.replicate):
    patch.set_facecolor(ORANGE if rep == 6 else BLUE)
ax.scatter(range(1, len(d) + 1), d.TIN_median, marker="D", s=40, color="white", edgecolor="black", zorder=3)
ax.set(xlabel="ExpA replicate", ylabel="transcript integrity number (TIN)")
fig.suptitle("TIN per sample", fontweight="bold", x=0.02, ha="left")
ax.set_title("boxes: all transcripts with TIN > 0; diamond: median (tin.py summary)", fontsize=10, loc="left")
fig.tight_layout(); fig.savefig(OUT / "figure08a_TIN_per_sample.png", dpi=300); plt.close(fig)


def scatter(x, xlab, title, row, fname):
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    b = np.polyfit(d[x], d.TIN_median, 1)
    xs = np.linspace(d[x].min(), d[x].max(), 50)
    ax.plot(xs, np.polyval(b, xs), "--", color="grey", lw=0.8)
    ax.scatter(d[x], d.TIN_median, s=180, c=[ORANGE if r == 6 else BLUE for r in d.replicate])
    for xi, yi, r in zip(d[x], d.TIN_median, d.replicate):
        ax.text(xi, yi, str(r), color="white", ha="center", va="center", fontsize=8, fontweight="bold")
    ax.set(xlabel=xlab, ylabel="TIN (median over transcripts)")
    fig.suptitle(title, fontweight="bold", x=0.02, ha="left")
    ax.set_title(f"Spearman rho = {row.spearman_rho:.2f} (p = {row.spearman_p:.2g}); "
                 f"Pearson r = {row.pearson_r:.2f} (p = {row.pearson_p:.2g}); n = {len(d)}", fontsize=9, loc="left")
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=300); plt.close(fig)


scatter("pct_rRNA", "rRNA % (paper's Table S2C)", "TIN vs. rRNA", ct.iloc[0], "figure08b_TIN_vs_rRNA.png")
scatter("pct_singletons", "singleton reads, % of mapped fragments (STAR logs)", "TIN vs. singletons",
        ct[ct.comparison == "TIN median vs. singleton %"].iloc[0],
        "figure08c_TIN_vs_singletons.png")
