"""Python equivalent of 10a_mislabelling.R (pydeseq2 instead of DESeq2).

Mislabelling: the 3 ExpA and 3 ExpB replicates with the highest TIN (ExpA1-3, ExpB1-3 by TIN).
mislabelled: ExpA1, ExpA2, ExpB3 vs. ExpB1, ExpB2, ExpA3; correct: ExpA1-3 vs. ExpB1-3; batch ignored (~ group).
Sample QC on raw log2 counts of genes with > 2 reads in every sample: R-R plots of all 15 pairs, Pearson
correlation matrix, dendrogram of Euclidean distances (average linkage).
Run from the repository root. Output: 10_mislabelling/results/python/.
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
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

warnings.filterwarnings("ignore")
WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
CPUS = int(os.environ.get("THREADS", "4"))
OUT = Path("10_mislabelling/results/python")
OUT.mkdir(parents=True, exist_ok=True)
COLS = {"within ExpA": "#2C7FB8", "within ExpB": "#D95F0E", "ExpA vs. ExpB": "#636363"}

rmap = pd.read_csv("data/replicate_map.tsv", sep="\t")
rmap = rmap[rmap.paper_experiment.isin(["ExpA", "ExpB"])].copy()
rmap["TIN"] = [pd.read_csv(WORK / "TIN" / f"{r}_Aligned.sortedByCoord.out.summary.txt", sep="\t").iloc[0, 2]
               for r in rmap.run]
sel = rmap.sort_values(["paper_experiment", "TIN"], ascending=[True, False]).groupby("paper_experiment").head(3).copy()
sel["id"] = sel.paper_experiment + (sel.groupby("paper_experiment").cumcount() + 1).astype(str)
sel["label"] = sel.id + " (r" + sel.paper_replicate.astype(str) + ")"
sel["mislabelled_group"] = np.where(sel.id.isin(["ExpA1", "ExpA2", "ExpB3"]), "group 1", "group 2")
sel["correct_group"] = np.where(sel.paper_experiment == "ExpA", "group 1", "group 2")
sel = sel.set_index("id", drop=False)
sel[["id", "label", "run", "paper_replicate", "paper_experiment", "TIN", "mislabelled_group", "correct_group"]].to_csv(
    OUT / "table10a_samples.csv", index=False)
print(sel[["run", "paper_replicate", "TIN", "mislabelled_group", "correct_group"]])


def read_fc(run):
    d = pd.read_csv(WORK / "STAR_both/fC" / f"{run}.count", sep="\t", comment="#")
    return d.set_index("Geneid").iloc[:, -1]


counts = pd.DataFrame({i: read_fc(r) for i, r in zip(sel.id, sel.run)})


def de(groups):
    m = counts[counts.sum(axis=1) > 0]
    meta = pd.DataFrame({"group": np.where(groups == "group 1", "g1", "g2")}, index=m.columns)
    dds = DeseqDataSet(counts=m.T, metadata=meta, design="~group", quiet=True, n_cpus=CPUS)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "g2", "g1"], quiet=True, n_cpus=CPUS)
    st.summary()
    r = st.results_df
    return pd.DataFrame({"gene": r.index, "baseMean": r.baseMean, "log2FC": r.log2FoldChange, "padj": r.padj})


d_mis, d_cor = de(sel.mislabelled_group.values), de(sel.correct_group.values)


def sig(d, lfc=0):
    return set(d.gene[(d.padj < 0.05) & (d.log2FC.abs() > lfc)])


cnt = pd.DataFrame({
    "grouping": ["mislabelled (one sample swapped per group)", "correct (ExpA vs. ExpB)"],
    "group1": [", ".join(sel.id[sel.mislabelled_group == "group 1"]), ", ".join(sel.id[sel.correct_group == "group 1"])],
    "group2": [", ".join(sel.id[sel.mislabelled_group == "group 2"]), ", ".join(sel.id[sel.correct_group == "group 2"])],
    "DE_padj05": [len(sig(d_mis)), len(sig(d_cor))],
    "DE_padj05_lfc1": [len(sig(d_mis, 1)), len(sig(d_cor, 1))],
    "tested": [d_mis.padj.notna().sum(), d_cor.padj.notna().sum()],
    "mislabelled_DE_also_in_correct": [len(sig(d_mis) & sig(d_cor)), np.nan]})
cnt.to_csv(OUT / "table10b_DE_counts.csv", index=False)
print(cnt.to_string(index=False))
d_mis[d_mis.padj < 0.05].to_csv(OUT / "table10c_DE_genes_mislabelled.csv", index=False)
d_cor[d_cor.padj < 0.05].to_csv(OUT / "table10d_DE_genes_correct.csv", index=False)

# sample QC on raw counts
keep = (counts > 2).all(axis=1)
lg = np.log2(counts[keep])
print(f"\n{keep.sum()} genes with > 2 reads in all 6 samples")
pr = pd.DataFrame(list(itertools.combinations(sel.id, 2)), columns=["x", "y"])
pr["type"] = [f"within {a[:4]}" if a[:4] == b[:4] else "ExpA vs. ExpB" for a, b in zip(pr.x, pr.y)]
pr["r"] = [np.corrcoef(lg[a], lg[b])[0, 1] for a, b in zip(pr.x, pr.y)]
pr["type"] = pd.Categorical(pr.type, list(COLS))
pr = pr.sort_values("type", kind="stable").reset_index(drop=True)
pr.to_csv(OUT / "table10e_pairwise_correlation.csv", index=False)

lo, hi = lg.values.min(), lg.values.max()
fig, axes = plt.subplots(3, 5, figsize=(15, 9.5))
for ax, (_, p) in zip(axes.flat, pr.iterrows()):
    ax.scatter(lg[p.x], lg[p.y], s=0.3, alpha=0.15, color=COLS[p.type], rasterized=True)
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.5)
    ax.set(xlim=(lo, hi), ylim=(lo, hi), xlabel=sel.label[p.x], ylabel=sel.label[p.y]); ax.set_aspect("equal")
    ax.set_title(p.type, color=COLS[p.type], fontsize=8, fontweight="bold", loc="left")
    ax.text(lo + 0.5, hi - 0.5, f"r = {p.r:.4f}", va="top", fontsize=8)
    ax.tick_params(labelsize=6)
fig.suptitle(f"Replicate-replicate plots, raw counts (log2), genes with > 2 reads in all 6 samples ({keep.sum()} genes)",
             fontweight="bold", x=0.01, ha="left")
fig.tight_layout(); fig.savefig(OUT / "figure10a_RR_plots.png", dpi=300); plt.close(fig)

Z = linkage(pdist(lg.T.values, metric="euclidean"), method="average")
order = dendrogram(Z, no_plot=True, labels=list(lg.columns))["ivl"]
cm = lg.corr().loc[order, order]
cm.round(5).to_csv(OUT / "table10f_correlation_matrix.csv")
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(cm.values, cmap="OrRd")
for i in range(len(order)):
    for j in range(len(order)):
        ax.text(j, i, f"{cm.values[i, j]:.4f}", ha="center", va="center", fontsize=8)
ax.set_xticks(range(len(order)), order); ax.set_yticks(range(len(order)), order)
fig.colorbar(im, label="Pearson r")
ax.set_title("Correlation of raw log2 counts (order of the dendrogram)", loc="left", fontweight="bold", fontsize=10)
fig.tight_layout(); fig.savefig(OUT / "figure10b_correlation_matrix.png", dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 5.5))
dendrogram(Z, labels=list(lg.columns), ax=ax, color_threshold=0, above_threshold_color="black")
for t in ax.get_xticklabels():
    t.set_color("#2C7FB8" if t.get_text().startswith("ExpA") else "#D95F0E")
ax.set_ylabel("Euclidean distance (raw log2 counts)")
ax.set_title("Sample dendrogram (average linkage)\nshortest join = most similar samples", fontweight="bold")
fig.tight_layout(); fig.savefig(OUT / "figure10c_dendrogram.png", dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 5.5))
ax.bar(range(2), cnt.DE_padj05, color=["#D95F0E", "#2C7FB8"], width=0.6)
for i, rr in cnt.iterrows():
    ax.text(i, rr.DE_padj05, f"{rr.DE_padj05} DE genes\n({rr.DE_padj05_lfc1} with |log2FC| > 1)", ha="center", va="bottom")
ax.set_xticks(range(2), cnt.grouping, fontsize=9); ax.set_ylabel("DE genes (padj < 0.05)")
ax.set_ylim(0, cnt.DE_padj05.max() * 1.2 + 1)
ax.set_title(f"One swapped sample per group\nmislabelled: {cnt.group1[0]} vs. {cnt.group2[0]}\n"
             f"correct: {cnt.group1[1]} vs. {cnt.group2[1]}", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig(OUT / "figure10d_DE_counts.png", dpi=300); plt.close(fig)
