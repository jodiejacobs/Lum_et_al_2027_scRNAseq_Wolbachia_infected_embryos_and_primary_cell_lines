#!/usr/bin/env python3
"""Figure 11: species timing of the EGFR loss (Pvr/Egfr ratio, Egfr, Star, rho by stage).
Run from the project root after pseudotime_rtk_switch."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = pd.read_csv("results/pseudotime/rtk_switch/switch_by_sample.csv")
d["lr"] = np.log2((d.cpm_Pvr + 1) / (d.cpm_Egfr + 1))
inf = d.condition == "Dsim6B-wMel"
base = d[~inf]                                      # cured Dsim6B represents the Dsim_w line stage
g = base.groupby(["lineage", "stage"])[["lr", "cpm_Egfr", "cpm_S", "cpm_rho"]].mean()
gi = d[inf][["lr", "cpm_Egfr", "cpm_S", "cpm_rho"]].mean()
ST = ["embryo", "primary", "line"]; X = np.arange(3)
BLUE, ORANGE = "#2a78d6", "#eb6834"
STY = {"Dmel_JupiterGFP": (BLUE, "-", "o"), "Dmel_Ubkhc": (BLUE, "--", "s"),
       "Dsim_Merrill23": (ORANGE, "-", "o"), "Dsim_w": (ORANGE, "--", "s")}
LAB = {"Dmel_JupiterGFP": "D. mel JupiterGFP", "Dmel_Ubkhc": "D. mel ubkhc",
       "Dsim_Merrill23": "D. sim Merrill23", "Dsim_w": "D. sim w−"}
panels = [("lr", "log2 (Pvr+1)/(Egfr+1)", "A  Pvr-to-Egfr ratio", False),
          ("cpm_Egfr", "CPM + 1 (log scale)", "B  Egfr", True),
          ("cpm_S", "CPM + 1 (log scale)", "C  Star (S)", True),
          ("cpm_rho", "CPM + 1 (log scale)", "D  rhomboid (rho)", True)]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#888",
                     "axes.labelcolor": "#333", "xtick.color": "#555", "ytick.color": "#555"})
fig, axes = plt.subplots(1, 4, figsize=(11, 3.1))
for ax, (col, yl, title, logy) in zip(axes, panels):
    for lin, (c, ls, m) in STY.items():
        y = g.loc[lin][col].reindex(ST).values
        ax.plot(X, y + (1 if logy else 0), ls=ls, color=c, lw=2, marker=m, ms=6, label=LAB[lin])
    yi = gi[col] + (1 if logy else 0)
    ax.plot(2.12, yi, marker="s", ms=7, mfc="white", mec=ORANGE, mew=1.8, ls="none",
            label="D. sim 6B-wMel (infected)")
    if logy:
        ax.set_yscale("log")
        from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
        ax.yaxis.set_major_locator(FixedLocator([1, 3, 10, 30, 100, 300]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks(X); ax.set_xticklabels(["embryo", "primary\ncell line", "cell line"])
    ax.set_xlim(-0.3, 2.4); ax.set_ylabel(yl); ax.set_title(title, loc="left", fontsize=10, color="#222")
    ax.grid(axis="y", color="#e6e6e6", lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0].axhline(0, color="#999", lw=0.8)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=5, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.1, 1, 1))
out = "results/pseudotime/rtk_switch/rtk_species_timing"
fig.savefig(out + ".png", dpi=250); fig.savefig(out + ".pdf")
print(g.round(1)); print(gi.round(1))
