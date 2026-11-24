"""Summary slide for post 5.5 (series summary of topics 02-05).

Run from the repository root. Output: 05_normalization_false_DE/results/figure05b_post05.5_summary_slide.png.
Numbers come from the topic READMEs and data/paper_TableS2C_rRNA.csv.
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch

OUT = Path("05_normalization_false_DE/results/figure05b_post05.5_summary_slide.png")

for f in ("NotoSans-Regular.ttf", "NotoSans-Bold.ttf"):
    p = Path("/usr/share/fonts/truetype/noto") / f
    if p.exists():
        font_manager.fontManager.addfont(str(p))
plt.rcParams["font.family"] = ["Noto Sans", "DejaVu Sans"]  # DejaVu supplies → and ≈

SURFACE = "#fcfcfb"
TILE = "#f1f0ec"
INK = "#0b0b0b"
INK2 = "#52514e"
ACCENT = "#2a78d6"

TILES = [
    ("Post 2 · detection",
     "<1% → >30%",
     "rRNA content among 17 replicates\ngrown under one condition",
     "Mean genome-wide correlation:\n≈ 0.997 for all 14 samples tested"),
    ("Post 3 · why",
     "few vs. thousands",
     "A handful of dominant features\nbarely moves a correlation",
     "computed over thousands of\ngenes that behave alike"),
    ("Post 4 · consequence",
     "+45%",
     "share of every coding gene when\n31% rRNA is left out of the total",
     "Raw counts changed:\nnot a single read"),
    ("Post 5 · normalization",
     "4 of 8",
     "splits changed DE calls under\ntotal-count scaling (23 → 8, 503 → 615)",
     "DESeq2 median-of-ratios:\nidentical counts in all 8"),
]

fig = plt.figure(figsize=(8, 10), dpi=150)
fig.patch.set_facecolor(SURFACE)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 125)
ax.axis("off")

# header
ax.text(6, 118, "Same genotype. Same condition.", fontsize=21, weight="bold", color=INK, va="top")
ax.text(6, 112.2, "What changed what we saw?", fontsize=21, weight="bold", color=ACCENT, va="top")
ax.text(6, 105.2, "Froussios et al. (2019, Bioinformatics) · 17 Arabidopsis Col-0 replicates, 3 experiments",
        fontsize=9.5, color=INK2, va="top")

# 2 x 2 tiles
W, H, GAP = 42, 30, 4
X0, Y_TOP = 6, 100
for i, (tag, hero, line1, line2) in enumerate(TILES):
    col, row = i % 2, i // 2
    x = X0 + col * (W + GAP)
    y = Y_TOP - (row + 1) * H - row * GAP
    ax.add_patch(FancyBboxPatch((x, y), W, H, boxstyle="round,pad=0,rounding_size=1.6",
                                fc=TILE, ec="none"))
    ax.text(x + 3, y + H - 3, tag.upper(), fontsize=8.5, weight="bold", color=INK2, va="top")
    ax.text(x + 3, y + H - 8, hero, fontsize=22 if len(hero) < 12 else 17,
            weight="bold", color=ACCENT, va="top")
    ax.text(x + 3, y + H - 15.5, line1, fontsize=9.5, color=INK, va="top", linespacing=1.35)
    ax.text(x + 3, y + 8, line2, fontsize=9.5, color=INK2, va="top", linespacing=1.35)

# split 6 strip: the open question
y = 20
ax.add_patch(FancyBboxPatch((6, y), 88, 11.5, boxstyle="round,pad=0,rounding_size=1.6",
                            fc="none", ec=ACCENT, lw=1.6))
ax.text(9, y + 8.6, "Open question  ·  Split 6", fontsize=9.5, weight="bold", color=ACCENT, va="top")
ax.text(9, y + 4.6, "> 500 DE genes under every normalization. The labels meant nothing biologically;\n"
        "the structure was already in the gene counts. Where does the contrast come from?",
        fontsize=9.5, color=INK, va="top", linespacing=1.35)

# take-home
ax.text(50, 13, "QC is not a report before the “real” analysis.", fontsize=13.5,
        weight="bold", color=INK, ha="center", va="top")
ax.text(50, 8.6, "It is part of scientific decision-making.", fontsize=13.5,
        weight="bold", color=INK, ha="center", va="top")
ax.text(94, 2.6, "#TransparentScience", fontsize=8, color=INK2, ha="right", va="bottom")

fig.savefig(OUT, facecolor=SURFACE)
print(OUT)
