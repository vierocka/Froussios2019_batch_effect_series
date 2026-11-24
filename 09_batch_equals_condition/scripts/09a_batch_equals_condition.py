"""Python equivalent of 09a_batch_equals_condition.R (pydeseq2 + inmoose ComBat-seq instead of DESeq2 + sva).

When batch = condition: the 4 ExpA replicates with the highest TIN ("condition A") vs. the 4 ExpB replicates
with the highest TIN ("condition B"). Model 1 ignores the batch (~ condition); model 2 removes it with
ComBat-seq first (~ batch + condition is not estimable). Run from the repository root (needs inmoose:
pip install inmoose). Output: 09_batch_equals_condition/results/python/ (the Shiny apps read the R tables).
"""
import os
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from inmoose.pycombat import pycombat_seq
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import pearsonr

sys.path.insert(0, "lib")
from go_enrichment import gene_info, go_enrich  # noqa: E402

warnings.filterwarnings("ignore")
WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
CPUS = int(os.environ.get("THREADS", "4"))
OUT = Path("09_batch_equals_condition/results/python")
OUT.mkdir(parents=True, exist_ok=True)
PAL = {"ExpA": "#2C7FB8", "ExpB": "#D95F0E"}

# samples: 4 best TIN per experiment
rmap = pd.read_csv("data/replicate_map.tsv", sep="\t")
rmap = rmap[rmap.paper_experiment.isin(["ExpA", "ExpB"])].copy()
rmap["TIN"] = [pd.read_csv(WORK / "TIN" / f"{r}_Aligned.sortedByCoord.out.summary.txt", sep="\t").iloc[0, 2]
               for r in rmap.run]
sel = (rmap.sort_values(["paper_experiment", "TIN"], ascending=[True, False]).groupby("paper_experiment").head(4)
       .assign(sample=lambda d: d.paper_experiment + "_r" + d.paper_replicate.astype(str),
               batch=lambda d: d.paper_experiment,
               condition=lambda d: np.where(d.paper_experiment == "ExpA", "A", "B")).set_index("sample", drop=False))
sel[["sample", "run", "paper_replicate", "batch", "condition", "TIN"]].to_csv(OUT / "table09a_samples.csv", index=False)
print(sel[["run", "paper_replicate", "batch", "condition", "TIN"]])


def read_fc(run):
    d = pd.read_csv(WORK / "STAR_both/fC" / f"{run}.count", sep="\t", comment="#")
    return d.set_index("Geneid").iloc[:, -1]


counts = pd.DataFrame({s: read_fc(r) for s, r in zip(sel["sample"], sel.run)})
counts = counts[(counts >= 10).sum(axis=1) >= 4]
print(f"{len(counts)} genes with >= 10 reads in at least 4 samples")

# model 2: ~ batch + condition is not estimable (identical columns)
X = pd.get_dummies(sel[["batch", "condition"]], drop_first=True).astype(float)
X.insert(0, "intercept", 1.0)
rank = np.linalg.matrix_rank(X.values)
msg = f"design matrix has {X.shape[1]} columns but rank {rank}: not estimable"
print("~ batch + condition:", msg)
try:
    pycombat_seq(counts, sel.batch.tolist(), covar_mod=pd.DataFrame({"c": sel.condition.values}))
    msg_cs = "ran"
except Exception as e:      # noqa: BLE001
    msg_cs = str(e).strip().splitlines()[0]
print("ComBat-seq with condition as covariate:", msg_cs)
pd.DataFrame({"design": ["DESeq2 ~ batch + condition", "ComBat-seq, batch with condition as covariate"],
              "result": [msg, msg_cs]}).to_csv(OUT / "table09b_confounded_designs.csv", index=False)
combat = pycombat_seq(counts, sel.batch.tolist())
combat = pd.DataFrame(np.asarray(combat), index=counts.index, columns=counts.columns).round().astype(int)

info = gene_info().set_index("gene")
universe = list(counts.index)


def analyse(m, label, tag):
    meta = sel[["condition"]]
    dds = DeseqDataSet(counts=m.T, metadata=meta, design="~condition", quiet=True, n_cpus=CPUS)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "B", "A"], quiet=True, n_cpus=CPUS)
    st.summary()
    de = st.results_df.rename(columns={"log2FoldChange": "log2FC"})[["baseMean", "log2FC", "lfcSE", "pvalue", "padj"]]
    de.insert(0, "gene", de.index)
    de["symbol"] = de.gene.map(info.symbol)
    sig = de[de.padj < 0.05]
    strong = sig[sig.log2FC.abs() > 1]
    sig.sort_values("padj").to_csv(OUT / f"table09c_DE_genes_{tag}.csv", index=False)
    go_de = go_enrich(strong.gene, universe) if len(strong) >= 5 else go_enrich([], universe).iloc[0:0]
    go_de.to_csv(OUT / f"table09d_GO_DE_genes_{tag}.csv", index=False)

    dds.vst(use_design=False)
    vs = pd.DataFrame(dds.layers["vst_counts"], index=m.columns, columns=m.index)
    xc = vs - vs.mean(axis=0)
    u, s, vt = np.linalg.svd(xc.values, full_matrices=False)
    scores = u * s
    ve = 100 * s ** 2 / np.sum(s ** 2)
    sc = pd.DataFrame(scores[:, :3], index=m.columns, columns=["PC1", "PC2", "PC3"])
    sc["batch"] = sel.batch
    bnum = (sel.batch == "ExpB").astype(float).values
    pcb = []
    for k in range(3):
        r, p = pearsonr(scores[:, k], bnum)
        pcb.append(dict(model=label, PC=f"PC{k + 1}", pct_variance=ve[k], r_with_batch=r, R2=r ** 2, p=p))
    pcb = pd.DataFrame(pcb)
    contrib = pd.DataFrame(100 * vt[:3].T ** 2, index=m.index, columns=["contrib_PC1", "contrib_PC2", "contrib_PC3"])
    contrib["loading_PC1"] = vt[0]
    top20 = contrib.sort_values("contrib_PC1", ascending=False).head(20).copy()
    top20.insert(0, "gene", top20.index)
    top20.insert(1, "symbol", top20.gene.map(info.symbol))
    top20[["gene", "symbol", "contrib_PC1", "loading_PC1"]].to_csv(OUT / f"table09e_top20_PC1_genes_{tag}.csv", index=False)
    go_pc = go_enrich(top20.gene, universe)
    go_pc.to_csv(OUT / f"table09f_GO_top20_PC1_{tag}.csv", index=False)
    n_tested = de.padj.notna().sum()
    print(f"\n== {label}\nDE genes padj < 0.05: {len(sig)} ({100 * len(sig) / n_tested:.1f}% of {n_tested} tested); "
          f"with |log2FC| > 1: {len(strong)}")
    print(f"GO terms FDR < 0.05 for DE genes: {(go_de.FDR < 0.05).sum()}; for top 20 PC1 genes: {(go_pc.FDR < 0.05).sum()}")
    print(pcb.round(4).to_string(index=False))
    return dict(label=label, n_de=len(sig), n_strong=len(strong), sc=sc, ve=ve, pcb=pcb, top20=top20)


r1 = analyse(counts, "batch ignored (~ condition)", "batch_ignored")
r2 = analyse(combat, "batch removed (ComBat-seq, then ~ condition)", "batch_removed")
pd.concat([r1["pcb"], r2["pcb"]]).to_csv(OUT / "table09g_PC_batch_correlation.csv", index=False)
pd.DataFrame({"model": [r1["label"], r2["label"]], "DE_padj05": [r1["n_de"], r2["n_de"]],
              "DE_padj05_lfc1": [r1["n_strong"], r2["n_strong"]]}).to_csv(OUT / "table09h_DE_counts.csv", index=False)

# figures
fig, axes = plt.subplots(2, 2, figsize=(10, 9))
for row, r in zip(axes, (r1, r2)):
    for ax, (x, y) in zip(row, (("PC1", "PC2"), ("PC2", "PC3"))):
        for b, c in PAL.items():
            s = r["sc"][r["sc"].batch == b]
            ax.scatter(s[x], s[y], s=70, color=c, label=f"{b} = condition {'A' if b == 'ExpA' else 'B'}")
        for smp, p in r["sc"].iterrows():
            ax.annotate(smp.split("_")[1], (p[x], p[y]), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8)
        ax.set_xlabel(f"{x} ({r['ve'][int(x[2]) - 1]:.0f}%)"); ax.set_ylabel(f"{y} ({r['ve'][int(y[2]) - 1]:.0f}%)")
    row[0].set_title(f"{r['label']}\n{r['n_de']} DE genes (padj < 0.05)", loc="left", fontsize=10, fontweight="bold")
axes[1, 1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "figure09a_PCA.png", dpi=300); plt.close(fig)

pcb = pd.concat([r1["pcb"], r2["pcb"]])
fig, ax = plt.subplots(figsize=(7.5, 5))
for i, (lab, col) in enumerate(zip((r1["label"], r2["label"]), ("#636363", "#74C476"))):
    s = pcb[pcb.model == lab]
    xs = np.arange(3) + (i - 0.5) * 0.38
    ax.bar(xs, s.r_with_batch.abs(), width=0.36, color=col, label=lab)
    for xx, (_, rr) in zip(xs, s.iterrows()):
        ax.text(xx, abs(rr.r_with_batch) + 0.02, f"{rr.pct_variance:.0f}% var.\np = {rr.p:.2g}", ha="center", fontsize=7)
ax.set_xticks(range(3), ["PC1", "PC2", "PC3"]); ax.set_ylim(0, 1.2)
ax.set_ylabel("|Pearson r| between PC score and batch"); ax.legend(fontsize=8)
ax.set_title("How much of each PC is batch?", loc="left", fontweight="bold")
fig.tight_layout(); fig.savefig(OUT / "figure09b_PC_batch_correlation.png", dpi=300); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(13, 6.5))
for ax, r in zip(axes, (r1, r2)):
    t = r["top20"].iloc[::-1]
    lab = [g if (pd.isna(s) or s == g) else f"{g} ({s})" for g, s in zip(t.gene, t.symbol)]
    ax.barh(lab, t.contrib_PC1, color=["#D95F0E" if v > 0 else "#2C7FB8" for v in t.loading_PC1])
    pos = ", ".join(i.split("_")[1] for i in r["sc"].index[r["sc"].PC1 > 0])
    ax.set_title(f"{r['label']}\ntop 20 genes together: {t.contrib_PC1.sum():.1f}% of PC1; positive PC1 side: {pos}\n"
                 "orange: positive loading, blue: negative", loc="left", fontsize=9)
    ax.set_xlabel("% of PC1 variance explained by the gene"); ax.tick_params(axis="y", labelsize=7)
fig.tight_layout(); fig.savefig(OUT / "figure09c_top20_PC1_contributions.png", dpi=300); plt.close(fig)
