"""Python equivalent of 06b_volcano_GOEA.R (pydeseq2 instead of DESeq2).

Volcano plot and GO enrichment of split 6's DE genes, combined with panels A and B of
06a_split_structure into one 4-panel figure. Run from the repository root after
06a_split_structure.py. Output: 06_split_structure/results/python/.
"""
import gzip
import os
import urllib.request
import warnings
from functools import lru_cache
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import gaussian_kde, hypergeom
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
WORK = Path(os.environ.get("PIPELINE_DIR", "pipeline_data"))
CPUS = int(os.environ.get("THREADS", "4"))
OUT = Path("06_split_structure/results/python")
GO_DIR = WORK / "GO_annotation"
GO_DIR.mkdir(parents=True, exist_ok=True)

# GO annotation: TAIR GAF and go-basic.obo from the GO Consortium (downloaded once if missing)
gaf_file, obo_file = GO_DIR / "tair.gaf.gz", GO_DIR / "go-basic.obo"
for f, url in ((gaf_file, "http://current.geneontology.org/annotations/tair.gaf.gz"),
               (obo_file, "http://current.geneontology.org/ontology/go-basic.obo")):
    if not f.exists():
        urllib.request.urlretrieve(url, f)

# coding-gene counts and split 6, as in 06a_split_structure.py
smp = pd.read_csv("data/replicate_map.tsv", sep="\t")
smp = smp[smp.paper_experiment != "ExpC"]
pct = pd.read_csv("data/paper_TableS2C_rRNA.csv").set_index("replicate")["pct_rRNA"]
coding = pd.DataFrame({f"r{rep}": pd.read_csv(WORK / "STAR_both/fC" / f"{r}.count", sep="\t", comment="#")
                       .set_index("Geneid").iloc[:, -1] for r, rep in zip(smp.run, smp.paper_replicate)})
coding = coding[coding.sum(axis=1) > 0]
pr = pd.Series(smp.paper_replicate.map(pct).values, index=coding.columns)
g1, g2 = ["r1", "r2", "r3"], ["r4", "r5", "r6"]
meta = pd.DataFrame({"group": ["g1"] * 3 + ["g2"] * 3}, index=g1 + g2)
dds = DeseqDataSet(counts=coding[g1 + g2].T, metadata=meta, design="~group", quiet=True, n_cpus=CPUS)
dds.deseq2()
st = DeseqStats(dds, contrast=["group", "g2", "g1"], quiet=True, n_cpus=CPUS)
st.summary()
res = st.results_df
tested = res.index[res.padj.notna()]
sig = list(res.index[res.padj < 0.05])
saved = pd.read_csv(OUT / "table06a_split6_DE_genes.csv").gene
assert set(sig) == set(saved), "DE genes differ from table06a (python)"
print(f"split 6: {len(sig)} DE genes (padj < 0.05) out of {len(tested)} tested; matches table06a")

# GO: direct annotations (NOT qualifiers removed), propagated to all ancestors (is_a, part_of)
gaf = pd.read_csv(gaf_file, sep="\t", comment="!", header=None, dtype=str, quoting=3, compression="gzip")
gaf = gaf[~gaf[3].fillna("").str.contains("NOT")][[1, 4]]
gaf.columns = ["gene", "go"]
terms, cur = {}, None
with open(obo_file) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("["):
            cur = {"parents": []} if line == "[Term]" else None
            continue
        if cur is None or ": " not in line:
            continue
        k, val = line.split(": ", 1)
        if k == "id":
            cur["id"] = val
            terms[val] = cur
        elif k in ("name", "namespace"):
            cur[k] = val
        elif k == "is_a":
            cur["parents"].append(val.split()[0])
        elif k == "relationship" and val.startswith("part_of "):
            cur["parents"].append(val.split()[1])
        elif k == "is_obsolete" and val == "true":
            terms.pop(cur["id"], None)
            cur = None


@lru_cache(maxsize=None)
def ancestors(t):
    a = {t}
    for p in terms[t]["parents"]:
        if p in terms:
            a |= ancestors(p)
    return frozenset(a)


gaf = gaf[gaf.go.isin(terms.keys()) & gaf.gene.isin(set(tested))].drop_duplicates()
g2go = pd.DataFrame([(g, a) for g, t in zip(gaf.gene, gaf.go) for a in ancestors(t)],
                    columns=["gene", "go"]).drop_duplicates()
universe = set(g2go.gene)
sig_go = [g for g in sig if g in universe]
print(f"GO-annotated: {len(universe)} of {len(tested)} tested genes, {len(sig_go)} of {len(sig)} DE genes")

# hypergeometric test per term (>= 5 annotated tested genes, no upper limit), BH within each ontology
size = g2go.go.value_counts()
hits = g2go[g2go.gene.isin(set(sig_go))].go.value_counts()
keep = size.index[size >= 5]
N, n = len(universe), len(sig_go)
goea = pd.DataFrame({"go": keep, "n_term": size[keep].values, "n_de": hits.reindex(keep, fill_value=0).values})
goea["expected"] = goea.n_term * n / N
goea["fold"] = goea.n_de / goea.expected
goea["p"] = hypergeom.sf(goea.n_de - 1, N, goea.n_term, n)
goea["name"] = [terms[t]["name"] for t in goea.go]
goea["ontology"] = [terms[t]["namespace"] for t in goea.go]
goea["FDR"] = goea.groupby("ontology").p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
de_genes = g2go[g2go.gene.isin(set(sig_go))].groupby("go").gene.apply(lambda g: ";".join(sorted(g)))
goea["DE_genes"] = goea.go.map(de_genes).fillna("")
goea = goea.sort_values(["FDR", "p"], kind="stable")[
    ["ontology", "go", "name", "n_term", "n_de", "expected", "fold", "p", "FDR", "DE_genes"]]
goea.to_csv(OUT / "table06c_GO_enrichment_split6.csv", index=False)
print(f"GO terms tested: {len(goea)}; FDR < 0.05: {(goea.FDR < 0.05).sum()}")
print(goea[goea.FDR < 0.05][["ontology", "go", "name", "n_de", "n_term", "fold", "FDR"]].head(20).to_string(index=False))
tf_genes = set(g2go.gene[g2go.go == "GO:0003700"])     # DNA-binding transcription factor activity

# figure
blue, orange = "#2C7FB8", "#D95F0E"
sf = dds.obs["size_factors"].values
nc6 = coding.loc[tested, g1 + g2] / sf
med1, med2 = np.log2(nc6[g1].median(axis=1) + 1), np.log2(nc6[g2].median(axis=1) + 1)
up = [g for g in sig if res.log2FoldChange[g] > 0]
dn = [g for g in sig if res.log2FoldChange[g] < 0]
fig = plt.figure(figsize=(16, 12))
grid = np.linspace(2, 18, 400)
# A: per-gene medians per group
for i, (title, genes) in enumerate([(f"all {len(tested):,} tested genes", tested),
                                    (f"{len(up)} DE genes higher in replicates 4, 5, 6", up),
                                    (f"{len(dn)} DE genes lower in replicates 4, 5, 6", dn)]):
    ax = fig.add_axes([0.05, 0.80 - i * 0.135, 0.36, 0.1])
    ax.plot(grid, gaussian_kde(med1[genes])(grid), color=blue, label="replicates 1, 2, 3")
    ax.plot(grid, gaussian_kde(med2[genes])(grid), color=orange, label="replicates 4, 5, 6")
    ax.set_title(title, fontsize=9)
    ax.set_ylabel("density", fontsize=8)
    if i == 0:
        ax.text(0, 1.25, "A. Gene medians per group", transform=ax.transAxes, fontweight="bold", fontsize=12)
    if i == 2:
        ax.set_xlabel("per-gene median of normalized counts in the group (log2 + 1)")
        ax.legend(ncol=2, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.35))
# B: volcano (unshrunken log2FC, as used for the DE call)
ax = fig.add_axes([0.5, 0.53, 0.47, 0.37])
d = res.loc[tested]
nde = d.padj >= 0.05
ax.scatter(d.log2FoldChange[nde], -np.log10(d.padj[nde]), s=3, color="grey", alpha=0.4, label="not DE")
for mask, col, lab in ((~nde & (d.log2FoldChange > 0), orange, "DE, higher in replicates 4, 5, 6"),
                       (~nde & (d.log2FoldChange < 0), blue, "DE, lower in replicates 4, 5, 6")):
    ax.scatter(d.log2FoldChange[mask], -np.log10(d.padj[mask]), s=4, color=col, alpha=0.6, label=lab)
tf = d.index.isin(list(tf_genes)) & ~nde
ax.scatter(d.log2FoldChange[tf], -np.log10(d.padj[tf]), s=18, facecolors="none", edgecolors="black", lw=0.5)
ax.axhline(-np.log10(0.05), ls="--", color="grey")
for g, name in (("AT4G16780", "ATHB-2/HAT4"), ("AT4G17460", "HAT1"), ("AT2G44910", "ATHB4")):
    if g in d.index:
        ax.annotate(name, (d.log2FoldChange[g], -np.log10(d.padj[g])), xytext=(15, 0), textcoords="offset points",
                    fontsize=9, arrowprops=dict(arrowstyle="-", lw=0.5))
ax.set(xlabel="log2 fold change (replicates 4, 5, 6 vs. 1, 2, 3)", ylabel="-log10 padj")
ax.set_title(f"B. Volcano plot of split 6\ndashed line: padj = 0.05; black ring: DE gene annotated as a transcription "
             f"factor (GO:0003700, {tf.sum()} genes)", loc="left", fontsize=10)
ax.legend(loc="upper left", fontsize=8, markerscale=3)
# C: rRNA per replicate
ax = fig.add_axes([0.07, 0.07, 0.3, 0.33])
for x, grp, col in ((0, g1, blue), (1, g2, orange)):
    ax.bar(x, pr[grp].median(), width=0.6, color=col, alpha=0.25)
    ax.scatter([x] * 3, pr[grp], color=col, s=40, zorder=3)
    for r in grp:
        ax.annotate(r.replace("r", "rep "), (x, pr[r]), xytext=(12, 0), textcoords="offset points", va="center", fontsize=9)
ax.set_xticks([0, 1], ["replicates 1, 2, 3", "replicates 4, 5, 6"])
ax.set_ylabel("% rRNA (paper's Table S2C)")
ax.set_title("C. rRNA per replicate (bar = group median)", loc="left", fontweight="bold", fontsize=11)
# D: GO enrichment, 12 terms with the smallest FDR; GO:0003700 in orange
dD = goea[goea.FDR < 0.05].head(12)
if dD.empty:
    dD = goea.head(12)
abbr = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}
labels = [f"{r.name} ({abbr[r.ontology]}, {r.n_de}/{r.n_term})" for r in dD.itertuples()]
ax = fig.add_axes([0.72, 0.07, 0.25, 0.33])
ax.barh(range(len(dD))[::-1], -np.log10(dD.FDR), color=[orange if g == "GO:0003700" else "grey" for g in dD.go])
ax.set_yticks(range(len(dD))[::-1], labels, fontsize=8)
ax.axvline(-np.log10(0.05), ls="--", color="grey")
ax.set_xlabel("-log10 FDR (hypergeometric, BH within ontology)")
ax.set_title(f"D. GO enrichment of the {len(sig)} DE genes ({len(sig_go)} GO-annotated)\nbackground: {len(universe):,} "
             "GO-annotated tested genes; brackets: DE genes / genes in term", loc="right", fontsize=10)
fig.suptitle(f"Same genotype, same batch (ExpA), no designed difference: replicates 1, 2, 3 vs. 4, 5, 6 -> "
             f"{len(sig)} DE genes (padj < 0.05)", fontsize=12, fontweight="bold", x=0.02, ha="left")
fig.savefig(OUT / "figure06b_split_structure_4panels.png", dpi=300, facecolor="white")
print("written: figure06b, table06c")
