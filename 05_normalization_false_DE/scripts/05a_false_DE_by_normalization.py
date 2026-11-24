"""Python equivalent of 05a_false_DE_by_normalization.R (pydeseq2 instead of DESeq2).

Can normalization choice change false-DE counts? Run from the repository root.
Output: 05_normalization_false_DE/results/python/. pydeseq2 is not numerically identical to
DESeq2, so DE counts can differ slightly from the R results.
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
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")
WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
CPUS = int(os.environ.get("THREADS", "4"))
OUT = Path("05_normalization_false_DE/results/python")
OUT.mkdir(parents=True, exist_ok=True)


def read_count(run):
    d = pd.read_csv(WORK / "STAR_both/fC" / f"{run}.count", sep="\t", comment="#")
    return pd.Series(d.iloc[:, -1].values, index=d["Geneid"])


def median_of_ratios(m):
    """DESeq2 size factors (estimateSizeFactorsForMatrix); m: genes x samples."""
    with np.errstate(divide="ignore"):
        lg = np.log(m.values.astype(float))
    logmeans = lg.mean(axis=1)
    sf = np.array([np.exp(np.median((lg[:, j] - logmeans)[np.isfinite(logmeans) & (m.values[:, j] > 0)]))
                   for j in range(m.shape[1])])
    return sf / np.exp(np.mean(np.log(sf)))


def deseq_with_sf(counts, groups, sf):
    """DESeq2 pipeline (~group) with given size factors; returns the fitted dataset."""
    meta = pd.DataFrame({"group": groups}, index=counts.columns)
    dds = DeseqDataSet(counts=counts.T.astype(int), metadata=meta, design="~group", quiet=True, n_cpus=CPUS)
    dds.obs["size_factors"] = sf
    dds.layers["normed_counts"] = dds.X / sf[:, None]
    dds.var["_normed_means"] = dds.layers["normed_counts"].mean(axis=0)
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()
    if dds.refit_cooks:
        dds.refit()
    return dds


def results(dds, **kw):
    st = DeseqStats(dds, contrast=["group", "g2", "g1"], quiet=True, n_cpus=CPUS, **kw)
    st.summary()
    return st.results_df


# coding-gene counts (featureCounts -t gene has no rRNA rows) for the 14 ExpA/ExpB samples
smp = pd.read_csv("data/replicate_map.tsv", sep="\t")
smp = smp[smp.paper_experiment != "ExpC"]
pct = pd.read_csv("data/paper_TableS2C_rRNA.csv").set_index("replicate")["pct_rRNA"]
smp["p"] = smp.paper_replicate.map(pct).values / 100
coding = pd.DataFrame({f"r{rep}": read_count(r) for r, rep in zip(smp.run, smp.paper_replicate)})
coding = coding[coding.sum(axis=1) > 0]

# SYNTHETIC rRNA row (see topic 04), used only as an alternative source of size factors; DE is
# always tested on the same coding-only matrix.
p = pd.Series(smp.p.values, index=coding.columns)
with_rrna = pd.concat([coding, pd.DataFrame([np.round(coding.sum() * p / (1 - p))], index=["rRNA"])])

# The 8 splits drawn by the R script with set.seed(1): replicate 6 plus two others in group 2,
# three of the remaining four in group 1 (8 of the 60 possible comparisons).
SPLITS = [((4, 1, 5), (6, 2, 7)), ((2, 3, 4), (6, 1, 5)), ((3, 7, 1), (6, 2, 4)), ((3, 5, 7), (6, 1, 2)),
          ((2, 7, 5), (6, 1, 3)), ((2, 1, 3), (6, 4, 5)), ((2, 7, 4), (6, 3, 5)), ((2, 5, 7), (6, 1, 4))]


def n_de(norm_m, ids, total_norm):
    if total_norm:
        sf = norm_m[ids].sum().values.astype(float)
        sf = sf / sf.mean()
    else:
        sf = median_of_ratios(norm_m[ids])
    dds = deseq_with_sf(coding[ids], ["g1"] * 3 + ["g2"] * 3, sf)
    return int((results(dds).padj < 0.05).sum())


rows = []
for k, (g1, g2) in enumerate(SPLITS, start=1):
    ids = [f"r{i}" for i in g1 + g2]
    for norm, total in (("DESeq2 median-of-ratios", False), ("total-count scaling", True)):
        for rrna, m in (("included for size factors", with_rrna), ("ignored for size factors", coding)):
            rows.append(dict(split=k, normalisation=norm, rRNA=rrna, n=n_de(m, ids, total)))
            print(rows[-1], flush=True)
dB = pd.DataFrame(rows)
dB.to_csv(OUT / "table05a_false_DE_counts.csv", index=False)

# figure: one line per split, sqrt y axis, one panel per normalisation
fig, axes = plt.subplots(1, 2, figsize=(9, 6.5), sharey=True)
for ax, norm in zip(axes, ["DESeq2 median-of-ratios", "total-count scaling"]):
    d = dB[dB.normalisation == norm]
    for _, s in d.groupby("split"):
        y = [s.n[s.rRNA == "included for size factors"].iloc[0], s.n[s.rRNA == "ignored for size factors"].iloc[0]]
        ax.plot([0, 1], np.sqrt(y), color="grey", lw=1)
        ax.scatter([0, 1], np.sqrt(y), color="#2C7FB8", s=30, zorder=3)
    ax.set_xticks([0, 1], ["included", "ignored"])
    ax.set_xlim(-0.3, 1.3)
    ax.set_title(norm, fontsize=11)
    ax.set_xlabel("synthetic rRNA feature: source of size factors")
ticks = [0, 10, 100, 250, 500, 750]
axes[0].set_yticks(np.sqrt(ticks), [str(t) for t in ticks])
axes[0].set_ylabel("false-relative-to-label DE genes\n(padj < 0.05, sqrt scale)")
fig.suptitle("Same batch and genotype, no designed group difference:\n8 sampled 3-vs-3 splits from the 7 ExpA "
             "replicates (one line per split)\nReplicate 6 (23.7% rRNA, the highest-rRNA replicate retained in ExpA) "
             "is in every split", fontsize=10, x=0.02, ha="left")
fig.tight_layout()
fig.savefig(OUT / "figure05a_false_DE_by_normalization.png", dpi=300, facecolor="white")
print(dB.pivot_table(index="split", columns=["normalisation", "rRNA"], values="n"))
