"""Python equivalent of 07a_effect_size.R (pydeseq2 instead of DESeq2 + apeglm).

Does normalization only move genes across padj < 0.05, or also shift effect sizes? Same 8 splits
and 4 size-factor treatments as topic 05; only the rule that calls a gene "DE" changes. Run from
the repository root. Output: 07_effect_size/results/python/. Shrinkage uses pydeseq2's
lfc_shrink (apeGLM-type Cauchy prior); counts can differ slightly from DESeq2/apeglm.
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
OUT = Path("07_effect_size/results/python")
OUT.mkdir(parents=True, exist_ok=True)


def read_count(run):
    d = pd.read_csv(WORK / "STAR_both/fC" / f"{run}.count", sep="\t", comment="#")
    return pd.Series(d.iloc[:, -1].values, index=d["Geneid"])


def median_of_ratios(m):
    with np.errstate(divide="ignore"):
        lg = np.log(m.values.astype(float))
    logmeans = lg.mean(axis=1)
    sf = np.array([np.exp(np.median((lg[:, j] - logmeans)[np.isfinite(logmeans) & (m.values[:, j] > 0)]))
                   for j in range(m.shape[1])])
    return sf / np.exp(np.mean(np.log(sf)))


def fit(norm_m, ids, total_norm):
    """DESeq2 pipeline on the coding-only matrix with size factors from norm_m."""
    if total_norm:
        sf = norm_m[ids].sum().values.astype(float)
        sf = sf / sf.mean()
    else:
        sf = median_of_ratios(norm_m[ids])
    meta = pd.DataFrame({"group": ["g1"] * 3 + ["g2"] * 3}, index=ids)
    dds = DeseqDataSet(counts=coding[ids].T.astype(int), metadata=meta, design="~group", quiet=True, n_cpus=CPUS)
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


def stats(dds, **kw):
    st = DeseqStats(dds, contrast=["group", "g2", "g1"], quiet=True, n_cpus=CPUS, **kw)
    st.summary()
    return st


def calls(dds):
    """Call rules: significance, raw effect size, shrunken effect size, and a threshold test."""
    r = stats(dds).results_df                                   # MLE log2FC, H0: log2FC = 0
    s = stats(dds)
    s.lfc_shrink(coeff="group[T.g2]")
    sh = s.results_df.log2FoldChange.reindex(r.index).abs()
    t1 = stats(dds, lfc_null=1, alt_hypothesis="greaterAbs").results_df   # H0: |log2FC| <= 1
    lfc, sig = r.log2FoldChange.abs(), r.padj < 0.05
    return {"padj < 0.05": sig.sum(),
            "|log2FC| > 1 (raw)": (lfc > 1).sum(), "|log2FC| > 2 (raw)": (lfc > 2).sum(),
            "|log2FC| > 1 (shrunken)": (sh > 1).sum(), "|log2FC| > 2 (shrunken)": (sh > 2).sum(),
            "padj < 0.05 & |log2FC| > 1 (raw)": (sig & (lfc > 1)).sum(),
            "padj < 0.05 & |log2FC| > 2 (raw)": (sig & (lfc > 2)).sum(),
            "padj < 0.05 & |log2FC| > 1 (shrunken)": (sig & (sh > 1)).sum(),
            "padj < 0.05 & |log2FC| > 2 (shrunken)": (sig & (sh > 2)).sum(),
            "test |log2FC| > 1, padj < 0.05": (t1.padj < 0.05).sum()}, r


# identical to topic 05: coding-gene counts, synthetic rRNA row, the 8 splits drawn by R (set.seed(1))
smp = pd.read_csv("data/replicate_map.tsv", sep="\t")
smp = smp[smp.paper_experiment != "ExpC"]
pct = pd.read_csv("data/paper_TableS2C_rRNA.csv").set_index("replicate")["pct_rRNA"]
coding = pd.DataFrame({f"r{rep}": read_count(r) for r, rep in zip(smp.run, smp.paper_replicate)})
coding = coding[coding.sum(axis=1) > 0]
p = pd.Series(smp.paper_replicate.map(pct).values / 100, index=coding.columns)
with_rrna = pd.concat([coding, pd.DataFrame([np.round(coding.sum() * p / (1 - p))], index=["rRNA"])])
SPLITS = [((4, 1, 5), (6, 2, 7)), ((2, 3, 4), (6, 1, 5)), ((3, 7, 1), (6, 2, 4)), ((3, 5, 7), (6, 1, 2)),
          ((2, 7, 5), (6, 1, 3)), ((2, 1, 3), (6, 4, 5)), ((2, 7, 4), (6, 3, 5)), ((2, 5, 7), (6, 1, 4))]
TREAT = [("DESeq2 median-of-ratios", "included for size factors", False, True),
         ("DESeq2 median-of-ratios", "ignored for size factors", False, False),
         ("total-count scaling", "included for size factors", True, True),
         ("total-count scaling", "ignored for size factors", True, False)]


def gsf(dds):
    l = np.log2(dds.obs["size_factors"].values)
    return l[3:].mean() - l[:3].mean()


res, shift, lowc = [], [], []
for k, (g1, g2) in enumerate(SPLITS, start=1):
    ids = [f"r{i}" for i in g1 + g2]
    fits = [fit(with_rrna if w else coding, ids, tot) for _, _, tot, w in TREAT]
    out = [calls(f) for f in fits]
    # who passes a raw log2FC cutoff in a split with no designed difference? (median-of-ratios, rRNA ignored)
    r0 = out[1][1]
    ok = r0.log2FoldChange.notna()
    for t in (1, 2):
        x = ok & (r0.log2FoldChange.abs() > t)
        lowc.append(dict(split=k, cutoff=t, n=int(x.sum()), median_baseMean=r0.baseMean[x].median(),
                         pct_baseMean_below_10=100 * (r0.baseMean[x] < 10).mean(),
                         median_baseMean_all=r0.baseMean[ok].median()))
    for (norm, rrna, _, _), (n, _) in zip(TREAT, out):
        res += [dict(split=k, normalisation=norm, rRNA=rrna, rule=rule, n=int(v)) for rule, v in n.items()]
    # how far does ignoring the rRNA feature move every gene's log2FC? (per normalisation)
    for m in (0, 2):
        d = out[m + 1][1].log2FoldChange - out[m][1].log2FoldChange
        shift.append(dict(split=k, normalisation=TREAT[m][0], median_shift=d.median(),
                          iqr_shift=d.quantile(0.75) - d.quantile(0.25), sf_shift=gsf(fits[m]) - gsf(fits[m + 1])))
    print(f"split {k} done", flush=True)
res, shift, lowc = pd.DataFrame(res), pd.DataFrame(shift), pd.DataFrame(lowc)
res.to_csv(OUT / "table07a_effect_size_counts.csv", index=False)
shift.to_csv(OUT / "table07b_log2FC_shift.csv", index=False)
lowc.to_csv(OUT / "table07c_raw_log2FC_expression.csv", index=False)

w = res.pivot_table(index=["split", "normalisation", "rule"], columns="rRNA", values="n").reset_index()
w["changed"] = w["included for size factors"] != w["ignored for size factors"]
print("splits (of 8) whose count changes when the rRNA feature is ignored:")
print(w.groupby(["rule", "normalisation"]).changed.sum().unstack().to_string())
print("\nrange of counts per rule and normalisation (min-max over splits and both rRNA treatments):")
print(res.groupby(["rule", "normalisation"]).n.agg(lambda x: f"{x.min()}-{x.max()}").unstack().to_string())
print("\nlog2FC shift when rRNA is ignored:")
print(shift.round(4).to_string(index=False))
print("\ngenes passing a raw |log2FC| cutoff (median-of-ratios, rRNA ignored): how well expressed are they?")
print(lowc.round(1).to_string(index=False))

# figure: median-of-ratios only (including or ignoring the rRNA feature gives the same size factors)
mor = res[(res.normalisation == "DESeq2 median-of-ratios") & (res.rRNA == "ignored for size factors")]
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), sharey=True)
panels = [(["padj < 0.05", "|log2FC| > 1 (raw)", "|log2FC| > 2 (raw)", "|log2FC| > 1 (shrunken)",
            "|log2FC| > 2 (shrunken)"], "A. Effect-size cutoffs alone",
           "8 same-batch 3-vs-3 splits, no designed difference; DESeq2 median-of-ratios"),
          (["padj < 0.05", "padj < 0.05 & |log2FC| > 1 (raw)", "padj < 0.05 & |log2FC| > 2 (raw)",
            "padj < 0.05 & |log2FC| > 1 (shrunken)", "padj < 0.05 & |log2FC| > 2 (shrunken)"],
           "B. Significance combined with effect size", "same splits; raw and shrunken log2FC")]
for ax, (rules, title, sub) in zip(axes, panels):
    for k, d in mor.groupby("split"):
        y = np.sqrt(d.set_index("rule").n.reindex(rules).values)
        ax.plot(range(len(rules)), y, color="0.75", lw=1)
        ax.scatter(range(len(rules)), y, s=30, zorder=3, color="#D95F0E" if k == 6 else "#2C7FB8",
                   label=("split 6 (topic 06)" if k == 6 else "other splits") if k in (1, 6) else None)
    ax.set_xticks(range(len(rules)), [r.replace(" (", "\n(").replace(" & ", "\n& ") for r in rules], fontsize=9)
    ax.set_title(f"{title}\n{sub}", loc="left", fontsize=11)
ticks = [0, 10, 100, 500, 1000, 2000]
axes[0].set_yticks(np.sqrt(ticks), [str(t) for t in ticks])
axes[0].set_ylabel("genes called (sqrt scale)")
axes[1].legend(loc="upper right", frameon=False)
fig.tight_layout()
fig.savefig(OUT / "figure07a_effect_size.png", dpi=300, facecolor="white")
