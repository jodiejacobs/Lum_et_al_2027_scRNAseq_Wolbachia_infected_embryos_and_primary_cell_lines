#!/usr/bin/env python3
"""Overview of the computational pipeline for Lum et al. 2027 (SVG; PDF/PNG via cairosvg)."""
import math
import sys

W, H = 1600, 1100
FONT = "Arial, 'Liberation Sans', Helvetica, sans-serif"
EMB, PRI, LIN = "#D9822B", "#2A9D8F", "#3B6FB6"
WOL = "#C2185B"
INK, MUTED, LINE = "#222222", "#5F5F5F", "#9A9A9A"
STEP = "#F4F4F4"
ARM_A, ARM_B = "#EAF1FA", "#EAF6F3"
out = []
add = out.append


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, size=15, weight="normal", fill=INK, anchor="middle", style="normal"):
    parts = esc(s).split("*")
    body = "".join(f'<tspan font-style="italic">{p}</tspan>' if i % 2 else p for i, p in enumerate(parts))
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
        f'font-style="{style}" fill="{fill}" text-anchor="{anchor}">{body}</text>')


def lines(x, y, rows, size=13.5, gap=18, **kw):
    for i, r in enumerate(rows):
        text(x, y + i * gap, r, size=size, **kw)


def box(x, y, w, h, fill=STEP, stroke=LINE, r=10, sw=1.3):
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


def arrow(x1, y1, x2, y2, color=MUTED, sw=2.2, head=10):
    ang = math.atan2(y2 - y1, x2 - x1)
    bx, by = x2 - head * math.cos(ang), y2 - head * math.sin(ang)
    add(f'<line x1="{x1}" y1="{y1}" x2="{bx:.1f}" y2="{by:.1f}" stroke="{color}" stroke-width="{sw}"/>')
    p2 = (bx + head * .55 * math.sin(ang), by - head * .55 * math.cos(ang))
    p3 = (bx - head * .55 * math.sin(ang), by + head * .55 * math.cos(ang))
    add(f'<polygon points="{x2},{y2} {p2[0]:.1f},{p2[1]:.1f} {p3[0]:.1f},{p3[1]:.1f}" fill="{color}"/>')


def badge(x, y, n, color=INK):
    add(f'<circle cx="{x}" cy="{y}" r="15" fill="{color}"/>')
    text(x, y + 5.5, str(n), size=15, weight="bold", fill="white")


def step(n, x, y, w, h, title, rows, tool=None):
    box(x, y, w, h)
    badge(x + 26, y + 28, n)
    text(x + 50, y + 34, title, size=18, weight="bold", anchor="start")
    lines(x + 50, y + 60, rows, anchor="start")
    if tool:
        text(x + w - 16, y + 34, tool, size=13, fill=MUTED, anchor="end", style="italic")


add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
add(f'<rect width="{W}" height="{H}" fill="white"/>')
text(W / 2, 46, "Computational processing of *Wolbachia*-infected *Drosophila* scRNA-seq libraries", size=26, weight="bold")

# ---- inputs: four stocks x three stages
y0 = 82
text(W / 2, y0 + 8, "Input: 18 PIPseq libraries from four stocks (two *D. melanogaster* with *w*Mel, two *D. simulans* with *w*Ri)",
     size=15, fill=MUTED)
cols = [("Embryos", "4 libraries", EMB), ("Primary cell lines", "4 libraries, a few weeks in culture", PRI),
        ("Cell lines", "10 libraries (5 lines x 2 replicates)", LIN)]
bw, gap = 400, 40
x0 = (W - (3 * bw + 2 * gap)) / 2
for i, (t, s, c) in enumerate(cols):
    x = x0 + i * (bw + gap)
    box(x, y0 + 24, bw, 62, fill=c, stroke=c)
    text(x + bw / 2, y0 + 51, t, size=19, weight="bold", fill="white")
    text(x + bw / 2, y0 + 73, s, size=13.5, fill="white")
arrow(W / 2, y0 + 92, W / 2, y0 + 120)

# ---- shared per-sample steps
sx, sw_ = 60, 1480
step(1, sx, 212, sw_, 92, "Joint host + *Wolbachia* quantification",
     ["Pseudoalign each library to a combined reference: host transcriptome (*D. melanogaster* or *D. simulans*) + *Wolbachia* genome (*w*Mel or *w*Ri) + a 16S panel",
      "Per cell: host gene counts, *Wolbachia* gene counts, and 16S counts (to flag other bacteria)"],
     tool="kallisto | bustools (kb count)")
arrow(W / 2, 304, W / 2, 330)
step(2, sx, 332, sw_, 92, "Quality control and *Wolbachia* load per cell",
     ["Remove low-quality cells and doublets (Scrublet, score > 0.4); record UMIs, genes detected, % mitochondrial",
      "*Wolbachia* load = *Wolbachia* gene UMIs / all UMIs; rRNA-based titer kept for the atlas analyses"],
     tool="scanpy, scrublet")
arrow(W / 2, 424, W / 2, 450)
step(3, sx, 452, sw_, 74, "Common gene space across species",
     ["Map *D. simulans* genes to 1:1 *D. melanogaster* orthologs (reciprocal best hits), so every sample uses FlyBase *D. melanogaster* annotation"],
     tool="ortholog map")

# ---- two arms
ya = 560
arrow(W / 2, 526, W * 0.27, ya - 4)
arrow(W / 2, 526, W * 0.73, ya - 4)
aw, ah = 720, 470
ax, bx = 40, W - 40 - aw
box(ax, ya, aw, ah, fill=ARM_A, stroke="#B9CBE3", r=14)
box(bx, ya, aw, ah, fill=ARM_B, stroke="#B5D9D1", r=14)
text(ax + aw / 2, ya + 34, "A. Atlas projection: embryos vs cell lines", size=20, weight="bold")
text(ax + aw / 2, ya + 56, "13 libraries (embryos + cell lines)", size=13.5, fill=MUTED)
text(bx + aw / 2, ya + 34, "B. Stage comparison within each stock", size=20, weight="bold")
text(bx + aw / 2, ya + 56, "4 stocks x 3 stages (embryo, primary cell line, cell line), all 18 libraries", size=13.5, fill=MUTED)


def sub(x, y, w, h, n, title, rows, color):
    box(x, y, w, h, fill="white", stroke=color, r=8)
    badge(x + 22, y + 24, n, color)
    text(x + 44, y + 30, title, size=16, weight="bold", anchor="start")
    lines(x + 44, y + 52, rows, size=13, gap=17, anchor="start")


ca, cb = "#3B6FB6", "#2A9D8F"
sub(ax + 20, ya + 74, aw - 40, 92, "4a", "Project onto the Flysta3D-v2 embryo atlas",
    ["Fixed atlas PCA/UMAP; ingest maps each cell to its atlas neighbours",
     "(no joint re-clustering, so samples cannot shape each other's embedding)",
     "Transfer cell type, tissue, and germ layer with a confidence score"], ca)
arrow(ax + aw / 2, ya + 166, ax + aw / 2, ya + 186, color=ca)
sub(ax + 20, ya + 188, aw - 40, 110, "5a", "Compare composition and similarity",
    ["Cell-type composition and Shannon diversity per condition",
     "Pseudobulk Spearman correlations: conditions and atlas tissues",
     "*Wolbachia* titer by transferred cell type and condition",
     "Tests: exact Mann-Whitney; permutation over all embryo/line relabelings"], ca)
arrow(ax + aw / 2, ya + 298, ax + aw / 2, ya + 318, color=ca)
sub(ax + 20, ya + 320, aw - 40, 130, "6a", "Readouts",
    ["Cell lines converge on hemolymph-, gonad-, and CNS-like identities",
     "Cell lines resemble each other more than their own stock's embryos",
     "Embryos, not cell lines, match a single embryonic tissue",
     "Figures 1-4; Tables S9-S11 (results/atlas_stats/)"], ca)

colw = (aw - 60) / 2
lx, rx = bx + 20, bx + 40 + colw
sub(lx, ya + 74, aw - 40, 74, "4b", "Per-stock embedding and pseudobulk",
    ["Within-sample HVGs (no ribosomal, mito, heat-shock, stress genes); kNN stage mixing",
     "Pseudobulk counts per library; stocks are the biological replicates"], cb)
arrow(bx + aw / 2, ya + 148, bx + aw / 2, ya + 166, color=cb)
blocks = [
    ("5b", ["Differential", "expression"], ["PyDESeq2,", "~ stock + stage", "Depth-matched rerun", "GSEA (GO, KEGG)"]),
    ("5c", ["Marker-based", "cell states"], ["Module scores", "per cell", "Paired logit t-test", "Proliferating flag"]),
    ("5d", ["Polycomb", "marking"], ["Embryo H3K27me3", "S2, BG3 Pc states", "Fisher, logistic OR", "GO of off genes"]),
    ("5e", ["Primary-culture", "cells"], ["Rank test vs", "cell-line change", "*Wolbachia* load", "Species candidates"]),
]
bw2, bh2 = (aw - 40 - 30) / 4, 150
for i, (n, t, rows) in enumerate(blocks):
    x = bx + 20 + i * (bw2 + 10)
    box(x, ya + 168, bw2, bh2, fill="white", stroke=cb, r=8)
    badge(x + 20, ya + 190, n, cb)
    lines(x + 42, ya + 188, t, size=14, gap=17, anchor="start", weight="bold")
    lines(x + 12, ya + 240, rows, size=12.5, gap=17, anchor="start")
arrow(bx + aw / 2, ya + 318, bx + aw / 2, ya + 338, color=cb)
sub(bx + 20, ya + 340, aw - 40, 110, "6b", "Readouts",
    ["Shared stage differences across four stocks and two species",
     "Silenced developmental regulators are Polycomb targets",
     "Proliferating primary cells partly resemble cell lines",
     "Figures 5-9; Tables 1-2, S1-S8, S12-S13 (results/pseudotime/)"], cb)

# ---- footer
fy = ya + ah + 40
add(f'<line x1="60" y1="{fy - 18}" x2="{W - 60}" y2="{fy - 18}" stroke="#DDDDDD" stroke-width="1"/>')
text(W / 2, fy + 4, "Snakemake 9 workflow on a SLURM cluster; conda environments for scanpy, kallisto | bustools, and R. "
     "Statistics treat stocks or conditions, not cells, as replicates, except within-culture correlations; BH correction throughout.",
     size=13.5, fill=MUTED)
add('</svg>')

svg = "\n".join(out)
base = sys.argv[1] if len(sys.argv) > 1 else "pipeline_overview"
open(base + ".svg", "w").write(svg)
# PDF/PNG: rsvg-convert (librsvg) renders the italic tspans correctly; cairosvg misplaces them
import shutil
import subprocess
if shutil.which("rsvg-convert"):
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", base + ".pdf", base + ".svg"], check=True)
    subprocess.run(["rsvg-convert", "-z", "2", "-o", base + ".png", base + ".svg"], check=True)
else:
    print("rsvg-convert not found (brew install librsvg / mamba install librsvg); wrote SVG only")
print("wrote", base)
