"""Python equivalent of 06a_split_structure.R (pydeseq2 instead of DESeq2).

Why does one same-batch 3-vs-3 split give hundreds of DE genes? Run from the repository root.
Output: 06_split_structure/results/python/. pydeseq2 is not numerically identical to DESeq2, so
DE counts can differ slightly from the R results.
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
from scipy.stats import gaussian_kde, spearmanr

warnings.filterwarnings("ignore")
WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
CPUS = int(os.environ.get("THREADS", "4"))
OUT = Path("06_split_structure/results/python")
OUT.mkdir(parents=True, exist_ok=True)


def read_fc(run):
    return pd.read_csv(WORK / "STAR_both/fC" / f"{run}.count", sep="\t", comment="#")


def de(g1, g2):
    """DESeq2 defaults (median-of-ratios, ~group); returns the results table."""
    ids = list(g1) + list(g2)
    meta = pd.DataFrame({"group": ["g1"] * len(g1) + ["g2"] * len(g2)}, index=ids)
    dds = DeseqDataSet(counts=coding[ids].T, metadata=meta, design="~group", quiet=True, n_cpus=CPUS)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "g2", "g1"], quiet=True, n_cpus=CPUS)
    st.summary()
    return st.results_df


def n_de(res):
    return int((res.padj < 0.05).sum())


def vst(m):
    """Blind variance-stabilizing transformation; m: genes x samples."""
    m = m[m.sum(axis=1) > 0]
    meta = pd.DataFrame({"x": ["a"] * m.shape[1]}, index=m.columns)
    dds = DeseqDataSet(counts=m.T, metadata=meta, design="~1", quiet=True, n_cpus=CPUS)
    dds.vst(use_design=False)
    return pd.DataFrame(dds.layers["vst_counts"].T, index=m.index, columns=m.columns)


def median_of_ratios(m):
    with np.errstate(divide="ignore"):
        lg = np.log(m.values.astype(float))
    logmeans = lg.mean(axis=1)
    sf = np.array([np.exp(np.median((lg[:, j] - logmeans)[np.isfinite(logmeans) & (m.values[:, j] > 0)]))
                   for j in range(m.shape[1])])
    return sf / np.exp(np.mean(np.log(sf)))


# coding-gene counts (featureCounts -t gene) for the 14 ExpA/ExpB samples, filtered as in topic 05
smp = pd.read_csv("data/replicate_map.tsv", sep="\t")
smp = smp[smp.paper_experiment != "ExpC"]
pct = pd.read_csv("data/paper_TableS2C_rRNA.csv").set_index("replicate")["pct_rRNA"]
first = read_fc(smp.run.iloc[0])
chrom = pd.Series(first["Chr"].astype(str).values, index=first["Geneid"])
coding = pd.DataFrame({f"r{rep}": read_fc(r).set_index("Geneid").iloc[:, -1]
                       for r, rep in zip(smp.run, smp.paper_replicate)})
coding = coding[coding.sum(axis=1) > 0]
expA = [f"r{rep}" for rep in smp.paper_replicate[smp.paper_experiment == "ExpA"]]
pr = pd.Series(smp.paper_replicate.map(pct).values, index=coding.columns)

# split 6 of topic 05: replicates 1, 2, 3 (group 1) vs. 4, 5, 6 (group 2); replicate 7 unused
g1, g2 = ["r1", "r2", "r3"], ["r4", "r5", "r6"]
res = de(g1, g2)
sig = res[res.padj < 0.05]
c = chrom[sig.index]
print(f"split 6: {len(sig)} DE genes (padj < 0.05), {(sig.log2FoldChange > 0).sum()} higher in group 2; "
      f"{(~c.isin(['Pt', 'Mt'])).sum()} nuclear, {(c == 'Pt').sum()} Pt, {(c == 'Mt').sum()} Mt")
print(f"median |log2FC| of DE genes {sig.log2FoldChange.abs().median():.2f} (max {sig.log2FoldChange.abs().max():.2f}); "
      f"median baseMean DE {sig.baseMean.median():.0f} vs. all tested {res.baseMean[res.padj.notna()].median():.0f}")
sig_df = sig[["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]].copy()
sig_df.insert(0, "chr", c.values)
sig_df.insert(0, "gene", sig.index)
sig_df.sort_values("padj").to_csv(OUT / "table06a_split6_DE_genes.csv", index=False)

# 1. PCA of the 7 ExpA replicates (vst, 500 most variable genes)
v = vst(coding[expA])
top = v.loc[v.var(axis=1).sort_values(ascending=False).index[:500]].T
top = top - top.mean()
u, s, _ = np.linalg.svd(top.values, full_matrices=False)
pc1 = pd.Series(u[:, 0] * s[0], index=top.index)
ve = 100 * s ** 2 / (s ** 2).sum()
print(f"\nPCA, 500 most variable genes: PC1 {ve[0]:.0f}%, PC2 {ve[1]:.0f}% of variance")
print("PC1 per replicate (sign arbitrary):", ", ".join(f"{k} {x:.1f}" for k, x in pc1.sort_values().items()))

# 2. all 70 distinct 3-vs-3 splits of the 7 ExpA replicates (one left out)
splits = {}
for u_ in expA:
    rest = [x for x in expA if x != u_]
    for a in itertools.combinations(rest, 3):
        b = [x for x in rest if x not in a]
        key = " vs ".join(sorted(["+".join(a), "+".join(b)]))
        splits[key] = (list(a), b)
rows = []
for i, (k, (a, b)) in enumerate(splits.items(), start=1):
    rows.append(dict(split=k, n_de=n_de(de(a, b)),
                     log_rRNA_separation=abs(np.log(pr[a]).mean() - np.log(pr[b]).mean())))
    print(f"  {i}/70 {k}: {rows[-1]['n_de']}", flush=True)
all70 = pd.DataFrame(rows).sort_values("n_de", ascending=False, kind="stable")
all70.to_csv(OUT / "table06b_all_70_splits.csv", index=False)
rank6 = list(all70.split).index("r1+r2+r3 vs r4+r5+r6") + 1
print(f"\nall {len(all70)} splits: median {all70.n_de.median():g} DE genes, 90th percentile "
      f"{all70.n_de.quantile(0.9):g}, max {all70.n_de.iloc[0]} ({all70.split.iloc[0]}); split 6 ranks {rank6}")
print(f"Spearman rho, DE count vs. rRNA separation of the two groups: "
      f"{spearmanr(all70.n_de, all70.log_rRNA_separation)[0]:.2f}")

# 3. leave one out, and swap replicate 6 or 3 for replicate 7
print("\nvariants of split 6:")
for a, b in [(g1, ["r4", "r5"]), (g1, ["r4", "r6"]), (g1, ["r5", "r6"]), (["r2", "r3"], g2), (["r1", "r3"], g2),
             (["r1", "r2"], g2), (g1, ["r4", "r5", "r7"]), (["r1", "r2", "r7"], g2)]:
    print(f"  {','.join(a):<10} vs {','.join(b):<10} {n_de(de(a, b)):4d} DE genes")

# 4. per sample: mean signed log2 deviation from the group-1 mean on the DE genes
nc = coding[expA] / median_of_ratios(coding[expA])
l = np.log2(nc.loc[sig.index] + 1)
l = l.sub(l[g1].mean(axis=1), axis=0).mul(np.sign(sig.log2FoldChange), axis=0)
print("\nmean signed log2 deviation from the group-1 mean on the DE genes (positive = group-2 direction):")
print(l.mean().round(2).to_string())

# figure: A per-gene medians per group (densities), B rRNA per replicate
blue, orange = "#2C7FB8", "#D95F0E"
tested = res.index[res.padj.notna()]
nc6 = coding.loc[tested, g1 + g2] / median_of_ratios(coding[g1 + g2])
med1, med2 = np.log2(nc6[g1].median(axis=1) + 1), np.log2(nc6[g2].median(axis=1) + 1)
up, dn = sig.index[sig.log2FoldChange > 0], sig.index[sig.log2FoldChange < 0]
sets = [(f"all {len(tested):,} tested genes", tested), (f"{len(up)} DE genes higher in replicates 4, 5, 6", up),
        (f"{len(dn)} DE genes lower in replicates 4, 5, 6", dn)]
fig = plt.figure(figsize=(11, 7.5))
grid = np.linspace(2, 18, 400)
for i, (title, genes) in enumerate(sets):
    ax = fig.add_axes([0.07, 0.69 - i * 0.29, 0.52, 0.22])
    ax.plot(grid, gaussian_kde(med1[genes])(grid), color=blue, lw=1.5, label="replicates 1, 2, 3")
    ax.plot(grid, gaussian_kde(med2[genes])(grid), color=orange, lw=1.5, label="replicates 4, 5, 6")
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("density")
    if i == 0:
        ax.text(0, 1.18, "A. Gene medians per group", transform=ax.transAxes, fontweight="bold")
    if i == 2:
        ax.set_xlabel("per-gene median of normalized counts in the group (log2 + 1)")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, frameon=False)
ax = fig.add_axes([0.7, 0.12, 0.27, 0.74])
for x, grp, col in ((0, g1, blue), (1, g2, orange)):
    ax.bar(x, pr[grp].median(), width=0.6, color=col, alpha=0.25)
    ax.scatter([x] * 3, pr[grp], color=col, s=40, zorder=3)
    for r in grp:
        ax.annotate(r.replace("r", "rep "), (x, pr[r]), xytext=(12, 0), textcoords="offset points", va="center", fontsize=9)
ax.set_xticks([0, 1], ["replicates 1, 2, 3", "replicates 4, 5, 6"])
ax.set_ylabel("% rRNA (paper's Table S2C)")
ax.set_title("B. rRNA per replicate\nbar = group median", loc="left", fontweight="bold", fontsize=10)
fig.suptitle(f"Same genotype, same batch (ExpA), no designed difference: replicates 1, 2, 3 vs. 4, 5, 6 -> "
             f"{len(sig)} DE genes (padj < 0.05)", fontsize=11, fontweight="bold", x=0.02, ha="left")
fig.savefig(OUT / "figure06a_split_structure.png", dpi=300, facecolor="white")
