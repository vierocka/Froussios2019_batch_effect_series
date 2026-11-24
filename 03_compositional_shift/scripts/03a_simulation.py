"""Python equivalent of 03a_simulation.R.

Teaching example: 10 extreme genes (like rRNA) barely move the correlation of two log2-scale
expression vectors. Run from the repository root. Output: 03_compositional_shift/results/python/.
Random numbers differ from R (different generator), so exact values differ slightly.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

OUT = Path("03_compositional_shift/results/python")
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(1)
n, n_changed, n_loops, rho = 10000, 10, 1000, 0.9
lo, hi = np.log2(0.01), 12                       # log2 scale, range of typical rlog / log2 values

# two vectors with correlation ~0.9 that both span [lo, hi] (Gaussian copula)
z = rng.multivariate_normal([0, 0], [[1, rho], [rho, 1]], n)
x = lo + (hi - lo) * norm.cdf(z[:, 0])
y = lo + (hi - lo) * norm.cdf(z[:, 1])
r0 = np.corrcoef(x, y)[0, 1]


def share(v, idx):
    """Linear-scale share of the total taken by a set of genes."""
    return (2.0 ** v[idx]).sum() / (2.0 ** v).sum()


# scenario 1: 10 random genes jump to the top of the range
r_new, share_new = np.empty(n_loops), np.empty(n_loops)
for i in range(n_loops):
    idx = rng.choice(n, n_changed, replace=False)
    y2 = y.copy()
    y2[idx] = rng.uniform(hi - 1, hi, n_changed)
    r_new[i] = np.corrcoef(x, y2)[0, 1]
    share_new[i] = share(y2, idx)

# scenario 2: the 10 genes are set to hold 30% of the linear-scale total (like the 31% rRNA library)
r_30 = np.empty(n_loops)
for i in range(n_loops):
    idx = rng.choice(n, n_changed, replace=False)
    y2 = y.copy()
    rest = np.delete(y, idx)
    y2[idx] = np.log2(0.3 / 0.7 * (2.0 ** rest).sum() / n_changed) + rng.normal(0, 0.1, n_changed)
    r_30[i] = np.corrcoef(x, y2)[0, 1]

print(f"baseline r = {r0:.4f}")
print(f"scenario 1 (10 genes at 11-12, {100 * share_new.mean():.1f}% of linear total): "
      f"r {r_new.min():.4f} - {r_new.max():.4f} over {n_loops} loops")
print(f"scenario 2 (10 genes hold 30% of linear total):      r {r_30.min():.4f} - {r_30.max():.4f} over {n_loops} loops")

# figure
idx = rng.choice(n, n_changed, replace=False)
y2 = y.copy()
y2[idx] = rng.uniform(hi - 1, hi, n_changed)
changed = np.zeros(n, bool)
changed[idx] = True
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.5))
a.scatter(x[~changed], y2[~changed], s=1, alpha=0.15, color="grey")
a.scatter(x[changed], y2[changed], s=25, color="#2C7FB8")
a.set(xlabel="sample 1 (log2 scale)", ylabel="sample 2 (log2 scale)")
a.set_title(f"A. r = {np.corrcoef(x, y2)[0, 1]:.3f}; 10 genes (blue) pushed to the top of the range",
            fontsize=10, fontweight="bold", loc="left")
bins = np.arange(0.86, 0.91 + 0.0005, 0.0005)
b.hist(r_new, bins=bins, alpha=0.7, color="#2C7FB8", label="~1% of total signal")
b.hist(r_30, bins=bins, alpha=0.7, color="#D95F0E", label="30% of total signal")
b.axvline(r0, linestyle="--", color="black")
b.set(xlim=(0.86, 0.91), xlabel="correlation after the change", ylabel="loops (of 1000)")
b.set_title(f"B. {n_loops} random draws of 10 genes; dashed = original r", fontsize=10, fontweight="bold", loc="left")
b.legend(title="10 changed genes hold", loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False)
fig.tight_layout()
fig.savefig(OUT / "figure03a_compositional_shift.png", dpi=300, facecolor="white")
