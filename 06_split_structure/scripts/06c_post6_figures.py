"""Python equivalent of 06c_post6_figures.R (pydeseq2 instead of DESeq2).

One single-panel figure per short post 6A-6D (volcano, rRNA swap, PCA, GO enrichment). Run from
the repository root after 06a_split_structure.py and 06b_volcano_GOEA.py (6D reads table06c).
Output: 06_split_structure/results/python/.
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
OUT = Path("06_split_structure/results/python")
blue, orange = "#2C7FB8", "#D95F0E"
plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})

smp = pd.read_csv("data/replicate_map.tsv", sep="\t")
smp = smp[smp.paper_experiment == "ExpA"]
pct = pd.read_csv("data/paper_TableS2C_rRNA.csv").set_index("replicate")["pct_rRNA"]
coding = pd.DataFrame({f"r{rep}": pd.read_csv(WORK / "STAR_both/fC" / f"{r}.count", sep="\t", comment="#")
                       .set_index("Geneid").iloc[:, -1] for r, rep in zip(smp.run, smp.paper_replicate)})
coding = coding[coding.sum(axis=1) > 0]
pr = pd.Series(smp.paper_replicate.map(pct).values, index=coding.columns)


def fit(g1, g2):
    ids = g1 + g2
    m = coding[ids]
    m = m[m.sum(axis=1) > 0]
    meta = pd.DataFrame({"group": ["g1"] * 3 + ["g2"] * 3}, index=ids)
    dds = DeseqDataSet(counts=m.T, metadata=meta, design="~group", quiet=True, n_cpus=CPUS)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "g2", "g1"], quiet=True, n_cpus=CPUS)
    st.summary()
    return st.results_df


def n_de(res):
    return int((res.padj < 0.05).sum())


def finish(fig, ax, title, subtitle, caption, name):
    fig.text(0.03, 0.97, title, fontsize=17, fontweight="bold", va="top")
    fig.text(0.03, 0.91, subtitle, fontsize=12, color="0.3", va="top")
    fig.text(0.03, 0.02, caption, fontsize=9, color="0.4")
    fig.savefig(OUT / name, dpi=200, facecolor="white")
    plt.close(fig)


def new_fig():
    fig = plt.figure(figsize=(7, 7))
    return fig, fig.add_axes([0.14, 0.12, 0.8, 0.68])


# 6A: volcano of split 6, DE = padj < 0.05 only (no fold-change cutoff)
res6 = fit(["r1", "r2", "r3"], ["r4", "r5", "r6"])
saved = pd.read_csv(OUT / "table06a_split6_DE_genes.csv").gene
assert set(res6.index[res6.padj < 0.05]) == set(saved), "DE genes differ from table06a (python)"
d = res6[res6.padj.notna()]
de = d[d.padj < 0.05]
small = int((de.log2FoldChange.abs() < 1).sum())
fig, ax = new_fig()
ax.scatter(d.log2FoldChange[d.padj >= 0.05], -np.log10(d.padj[d.padj >= 0.05]), s=3, color="0.8")
ax.scatter(de.log2FoldChange, -np.log10(de.padj), s=6, alpha=0.8,
           color=np.where(de.log2FoldChange > 0, orange, blue))
ax.axhline(-np.log10(0.05), ls="--", color="0.35")
ax.axvline(-1, ls=":", color="0.55")
ax.axvline(1, ls=":", color="0.55")
ax.text(-3.9, -np.log10(0.05), "padj = 0.05", va="bottom", fontsize=10, color="0.35")
ax.text(0, -np.log10(d.padj).max() * 0.97, f"{small} of {len(de)} DE genes\nchange less than 2-fold",
        ha="center", va="top", fontsize=11, color="0.2")
ax.set(xlim=(-4, 4), xlabel="log2 fold change (replicates 4, 5, 6 vs. 1, 2, 3)", ylabel="-log10 padj")
finish(fig, ax, f"Identical plants, {len(de)} DE genes",
       "Same genotype, same batch, no designed difference.\nDE = padj < 0.05 only, no fold-change cutoff.",
       "Froussios et al. 2019, ExpA; pydeseq2 defaults. Dotted lines: 2-fold.", "figure06c_post6A_volcano.png")

# 6B: swap the high-rRNA replicate 6 for low-rRNA replicate 7
res7 = fit(["r1", "r2", "r3"], ["r4", "r5", "r7"])
n = [n_de(res6), n_de(res7)]
fig, ax = new_fig()
ax.bar([0, 1], n, width=0.55, color=[orange, "0.55"])
for x, v in enumerate(n):
    ax.text(x, v, f"{v:,}", ha="center", va="bottom", fontsize=20, fontweight="bold")
ax.set_xticks([0, 1], [f"1, 2, 3 vs. 4, 5, 6\nreplicate 6: {pr['r6']:.1f}% rRNA",
                       f"1, 2, 3 vs. 4, 5, 7\nreplicate 7: {pr['r7']:.2f}% rRNA"])
ax.set_ylim(0, max(n) * 1.12)
ax.set_ylabel("DE genes (padj < 0.05)")
finish(fig, ax, "Is it the rRNA?", "Swap the highest-rRNA replicate for a low-rRNA one:\nDE genes go up, not down.",
       "Froussios et al. 2019, ExpA; rRNA % from the paper's Table S2C; pydeseq2 defaults.", "figure06d_post6B_rRNA_swap.png")

# 6C: PCA of the 7 replicates (vst, 500 most variable genes)
meta = pd.DataFrame({"x": ["a"] * coding.shape[1]}, index=coding.columns)
dds = DeseqDataSet(counts=coding.T, metadata=meta, design="~1", quiet=True, n_cpus=CPUS)
dds.vst(use_design=False)
v = pd.DataFrame(dds.layers["vst_counts"].T, index=coding.index, columns=coding.columns)
top = v.loc[v.var(axis=1).sort_values(ascending=False).index[:500]].T
top = top - top.mean()
u, s, _ = np.linalg.svd(top.values, full_matrices=False)
pc = pd.DataFrame(u[:, :2] * s[:2], index=top.index, columns=["PC1", "PC2"])
ve = np.round(100 * s ** 2 / (s ** 2).sum()).astype(int)
if pc.loc["r1", "PC1"] > 0:
    pc["PC1"] = -pc["PC1"]                                       # orient 1 -> 7 left to right
res_mix = fit(["r1", "r2", "r7"], ["r4", "r5", "r6"])
fig, ax = new_fig()
for r, row in pc.iterrows():
    col = blue if r in ("r1", "r2", "r3") else orange if r in ("r4", "r5", "r6") else "0.55"
    ax.scatter(row.PC1, row.PC2, s=900, color=col, edgecolors="white", lw=1)
    ax.text(row.PC1, row.PC2, r[1:], color="white", fontweight="bold", ha="center", va="center", fontsize=14)
y0 = pc.PC2.min() - 4
ax.annotate("", (pc.PC1.max(), y0), (pc.PC1.min(), y0), arrowprops=dict(arrowstyle="->", color="0.45"))
ax.text(pc.PC1.mean(), y0 - 1.5, f"1, 2, 3 vs. 4, 5, 7 (across the order): {n_de(res7):,} DE genes\n"
                                 f"1, 2, 7 vs. 4, 5, 6 (mixed): {n_de(res_mix)} DE genes",
        ha="center", va="top", fontsize=10, color="0.3")
ax.set_ylim(y0 - 8, pc.PC2.max() + 4)
ax.set(xlabel=f"PC1 ({ve[0]}% of variance)", ylabel=f"PC2 ({ve[1]}%)")
finish(fig, ax, "A hidden order",
       "The replicates line up almost by their number.\nBlue / orange: the two groups of split 6.",
       "Froussios et al. 2019, ExpA, 7 replicates; vst, 500 most variable genes.", "figure06e_post6C_PCA.png")

# 6D: GO enrichment, top terms with distinct gene sets; transcription factor term highlighted
goea = pd.read_csv(OUT / "table06c_GO_enrichment_split6.csv")
goea = goea[(goea.FDR < 0.05) & ~goea.DE_genes.duplicated()].head(8)
fig, ax = new_fig()
ax.set_position([0.5, 0.14, 0.45, 0.66])
ax.barh(range(len(goea))[::-1], -np.log10(goea.FDR), height=0.7,
        color=[orange if g == "GO:0003700" else "0.65" for g in goea.go])
ax.set_yticks(range(len(goea))[::-1], [f"{a}  ({b})" for a, b in zip(goea.name, goea.n_de)], fontsize=10)
ax.set_xlabel("-log10 FDR")
finish(fig, ax, "What are these genes?",
       "Transcription factors are the top GO term,\nin a comparison with no designed difference.",
       f"GO enrichment of the {len(saved)} DE genes (hypergeometric, BH). In brackets: DE genes in the term.",
       "figure06f_post6D_GO.png")

print(f"6A: {len(de)} DE ({small} < 2-fold); 6B: {n[0]} vs {n[1]}; 6C: PC1 {ve[0]}%, mix {n_de(res_mix)}; "
      "written figure06c-f")
