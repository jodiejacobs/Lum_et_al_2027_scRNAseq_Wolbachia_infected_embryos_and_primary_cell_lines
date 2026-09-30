#!/usr/bin/env python3
"""Summary cartoon of the Lum et al. 2027 cell-line establishment model (SVG)."""
import math
import random

W, H = 1600, 1000
FONT = "Arial, 'Liberation Sans', Helvetica, sans-serif"

# palette
EMB, PRI, LIN = "#D9822B", "#2A9D8F", "#3B6FB6"
WOL = "#C2185B"
INK, MUTED, GRID = "#222222", "#6B6B6B", "#D0D0D0"
PC = "#8E3B46"         # Polycomb / H3K27me3
MAPK = "#E0A100"
CELLTYPES = {"neural": "#7F9CCB", "muscle": "#9CCB86", "epidermis": "#E8C07D",
             "gut": "#C9A0DC", "fatbody": "#F2D16B", "crystal": "#8FA3AD",
             "plasmatocyte": "#8CB4E0", "lamellocyte": "#B7D7C9"}

out = []
add = out.append
rnd = random.Random(7)


def text(x, y, s, size=16, weight="normal", fill=INK, anchor="middle", style="normal"):
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" font-style="{style}" fill="{fill}" text-anchor="{anchor}">{s}</text>')


def lines(x, y, rows, size=14, gap=18, **kw):
    for i, r in enumerate(rows):
        text(x, y + i * gap, r, size=size, **kw)


def circle(x, y, r, fill, stroke=INK, sw=1.2, dash=None, opacity=1):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="{fill}" stroke="{stroke}" '
        f'stroke-width="{sw}"{d} opacity="{opacity}"/>')


def wolb(x, y, ang, n=1, spread=6):
    for _ in range(n):
        dx, dy = rnd.uniform(-spread, spread), rnd.uniform(-spread, spread)
        a = ang + rnd.uniform(-60, 60)
        add(f'<rect x="{x + dx - 4:.1f}" y="{y + dy - 1.6:.1f}" width="8" height="3.2" rx="1.6" '
            f'fill="{WOL}" transform="rotate({a:.0f} {x + dx:.1f} {y + dy:.1f})"/>')


def arrow(x1, y1, x2, y2, color=INK, sw=3, dash=None, head=12):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    ang = math.atan2(y2 - y1, x2 - x1)
    bx, by = x2 - head * math.cos(ang), y2 - head * math.sin(ang)
    add(f'<line x1="{x1}" y1="{y1}" x2="{bx:.1f}" y2="{by:.1f}" stroke="{color}" '
        f'stroke-width="{sw}"{d}/>')
    p1 = (x2, y2)
    p2 = (bx + head * 0.55 * math.sin(ang), by - head * 0.55 * math.cos(ang))
    p3 = (bx - head * 0.55 * math.sin(ang), by + head * 0.55 * math.cos(ang))
    add(f'<polygon points="{p1[0]:.1f},{p1[1]:.1f} {p2[0]:.1f},{p2[1]:.1f} '
        f'{p3[0]:.1f},{p3[1]:.1f}" fill="{color}"/>')


def dividing(x, y, r, fill, ring=None):
    """Two daughter nuclei in a waisted cell (mitotic figure)."""
    if ring:
        circle(x, y, r + 7, "none", stroke=ring, sw=2.2)
    add(f'<ellipse cx="{x - r * 0.45:.1f}" cy="{y:.1f}" rx="{r * 0.75:.1f}" ry="{r * 0.85:.1f}" '
        f'fill="{fill}" stroke="{INK}" stroke-width="1.2"/>')
    add(f'<ellipse cx="{x + r * 0.45:.1f}" cy="{y:.1f}" rx="{r * 0.75:.1f}" ry="{r * 0.85:.1f}" '
        f'fill="{fill}" stroke="{INK}" stroke-width="1.2"/>')
    circle(x - r * 0.45, y, r * 0.3, "#3A3A3A", stroke="none")
    circle(x + r * 0.45, y, r * 0.3, "#3A3A3A", stroke="none")


def cross(x, y, s=9, color="#B03A2E", sw=3):
    add(f'<line x1="{x - s}" y1="{y - s}" x2="{x + s}" y2="{y + s}" stroke="{color}" stroke-width="{sw}"/>')
    add(f'<line x1="{x - s}" y1="{y + s}" x2="{x + s}" y2="{y - s}" stroke="{color}" stroke-width="{sw}"/>')


add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
add(f'<rect width="{W}" height="{H}" fill="white"/>')

# ───────────────────────── top band: three stages ─────────────────────────
text(40, 40, "A", size=26, weight="bold", anchor="start")
text(800, 42, "Cell lines arise by outgrowth of a proliferative subset, not from one parental tissue",
     size=19, weight="bold")
cx = [255, 800, 1345]
for x, lab, col in zip(cx, ["Embryo", "Primary cell line", "Cell line"], [EMB, PRI, LIN]):
    add(f'<rect x="{x - 150}" y="62" width="300" height="30" rx="15" fill="{col}"/>')
    text(x, 83, lab, size=17, weight="bold", fill="white")

# Embryo: ellipse with diverse tissues, graded infection
ex, ey = cx[0], 225
add(f'<ellipse cx="{ex}" cy="{ey}" rx="190" ry="95" fill="#FFF6EC" stroke="{EMB}" stroke-width="3"/>')
types = ["neural", "muscle", "epidermis", "gut", "neural", "muscle", "epidermis"]
for i in range(58):
    while True:
        px, py = rnd.uniform(-175, 175), rnd.uniform(-82, 82)
        if (px / 175) ** 2 + (py / 82) ** 2 < 0.92:
            break
    t = types[i % len(types)]
    circle(ex + px, ey + py, 10, CELLTYPES[t], sw=0.8)
    if rnd.random() < 0.9:
        wolb(ex + px, ey + py, rnd.uniform(0, 180), n=rnd.choice([1, 1, 2, 2, 3]), spread=4)
lines(ex, 345, ["Many tissues (neural, muscle, epidermis, gut)",
                "Graded Wolbachia titer",
                "Developmental TFs and Hox genes on;",
                "their loci already Polycomb targets"], size=14)

# Primary cell line: dish with mixed states, AMP bursts, proliferative subset
px0, py0 = cx[1], 225
add(f'<ellipse cx="{px0}" cy="{py0}" rx="195" ry="100" fill="#EEF8F6" stroke="{PRI}" stroke-width="3"/>')
spots = []
for i in range(34):
    while True:
        qx, qy = rnd.uniform(-175, 175), rnd.uniform(-82, 82)
        if (qx / 175) ** 2 + (qy / 82) ** 2 < 0.9 and all(
                (qx - a) ** 2 + (qy - b) ** 2 > 26 ** 2 for a, b in spots):
            break
    spots.append((qx, qy))
kinds = ["fatbody", "crystal", "plasmatocyte", "lamellocyte", "fatbody", "neural", "plasmatocyte"]
for i, (qx, qy) in enumerate(spots):
    x, y = px0 + qx, py0 + qy
    k = kinds[i % len(kinds)]
    if i in (3, 11, 19, 27):
        dividing(x, y, 11, CELLTYPES["plasmatocyte"], ring=LIN)
    elif k == "fatbody":
        circle(x, y, 13, CELLTYPES["fatbody"], sw=0.9)
        for j in range(3):
            circle(x + rnd.uniform(-6, 6), y + rnd.uniform(-6, 6), 2.6, "white", stroke="#C9A640", sw=0.6)
    elif k == "crystal":
        circle(x, y, 11, CELLTYPES["crystal"], sw=0.9)
        add(f'<rect x="{x - 4}" y="{y - 2}" width="8" height="4" fill="#5A6B73" '
            f'transform="rotate(35 {x} {y})"/>')
    elif k == "lamellocyte":
        add(f'<ellipse cx="{x}" cy="{y}" rx="16" ry="6" fill="{CELLTYPES["lamellocyte"]}" '
            f'stroke="{INK}" stroke-width="0.9"/>')
    else:
        circle(x, y, 10, CELLTYPES[k], sw=0.9)
    if rnd.random() < 0.75:
        wolb(x, y, rnd.uniform(0, 180), n=rnd.choice([1, 1, 2]), spread=4)
# AMP / injury bursts
for (bx, by) in [(-130, -55), (60, 62), (140, -20), (-40, -70)]:
    x, y = px0 + bx, py0 + by
    for a in range(0, 360, 45):
        r1, r2 = 5, 11
        add(f'<line x1="{x + r1 * math.cos(math.radians(a)):.1f}" y1="{y + r1 * math.sin(math.radians(a)):.1f}" '
            f'x2="{x + r2 * math.cos(math.radians(a)):.1f}" y2="{y + r2 * math.sin(math.radians(a)):.1f}" '
            f'stroke="#E4572E" stroke-width="2"/>')
lines(px0, 345, ["Injury/antimicrobial-peptide, fat-body and",
                 "crystal-cell programs switch on",
                 "Rare proliferating cells (ringed) already",
                 "carry part of the cell-line expression shift"], size=14)

# Cell line: uniform proliferative hemocyte-like cells; infection retained or lost
lx, ly = cx[2], 225
add(f'<ellipse cx="{lx}" cy="{ly}" rx="195" ry="100" fill="#EEF3FB" stroke="{LIN}" stroke-width="3"/>')
add(f'<line x1="{lx - 190}" y1="{ly}" x2="{lx + 190}" y2="{ly}" stroke="{GRID}" stroke-width="1.5" '
    f'stroke-dasharray="5,4"/>')
for half, infected in [(-1, True), (1, False)]:
    placed = []
    for i in range(16):
        for _ in range(200):
            qx = rnd.uniform(-165, 165)
            qy = half * rnd.uniform(30, 80)
            if (qx / 175) ** 2 + (qy / 88) ** 2 < 0.88 and all(
                    (qx - a) ** 2 + (qy - b) ** 2 > 25 ** 2 for a, b in placed):
                break
        placed.append((qx, qy))
        x, y = lx + qx, ly + qy
        if i % 3 == 0:
            dividing(x, y, 10, CELLTYPES["plasmatocyte"])
        else:
            circle(x, y, 10, CELLTYPES["plasmatocyte"], sw=0.9)
        if infected:
            wolb(x, y, rnd.uniform(0, 180), n=4, spread=5)
add(f'<rect x="{lx - 150}" y="{ly - 18}" width="300" height="36" rx="6" fill="white" '
    f'stroke="{GRID}" stroke-width="1"/>')
text(lx, ly - 4, "▲ infected: JW18wMel, Dsim6B-wMel", size=12, fill=WOL)
text(lx, ly + 12, "▼ infection lost: ubkhc, Dsim-Merrill23", size=12, fill=MUTED)
lines(lx, 345, ["One convergent state across 2 species and",
                "4 lineages: proliferative (27–66% of cells),",
                "hemocyte-like, injury program off",
                "Programs lost, none gained"], size=14)

# stage arrows
arrow(460, 225, 590, 225, color=PRI)
lines(525, 190, ["dissociation", "and culture"], size=13, fill=MUTED)
arrow(1005, 225, 1135, 225, color=LIN)
lines(1070, 180, ["outgrowth of", "proliferative", "cells"], size=13, fill=MUTED)
text(1070, 262, "(inferred)", size=12, fill=MUTED, style="italic")

# ───────────────────────── bottom band: shared cell-line program ─────────────────────────
add(f'<line x1="40" y1="430" x2="{W - 40}" y2="430" stroke="{GRID}" stroke-width="1.5"/>')
text(40, 470, "B", size=26, weight="bold", anchor="start")
text(80, 470, "Shared program of established cell lines", size=19, weight="bold", anchor="start")

# big cell
CX, CY, R = 640, 720, 200
circle(CX, CY, R, "#F4F8FD", stroke=LIN, sw=4)
# nucleus
NX, NY, NR = CX + 50, CY + 30, 100
circle(NX, NY, NR, "#FFFFFF", stroke=INK, sw=2)
text(NX + 25, NY + NR - 22, "nucleus", size=12, fill=MUTED, style="italic")
# Polycomb-silenced loci: chromatin fibre with nucleosomes + me3 flags
def locus(x0, y0, label, n=5):
    add(f'<path d="M{x0} {y0} q12 -10 24 0 t24 0 t24 0 t24 0" fill="none" stroke="#555" stroke-width="2"/>')
    for i in range(n):
        xx = x0 + 6 + i * 21
        circle(xx, y0, 6, "#C9C9C9", stroke="#555", sw=1)
        add(f'<line x1="{xx}" y1="{y0 - 6}" x2="{xx}" y2="{y0 - 15}" stroke="{PC}" stroke-width="1.8"/>')
        circle(xx, y0 - 17, 3.2, PC, stroke="none")
    text(x0 + 52, y0 + 22, label, size=12, style="italic")
locus(NX - 90, NY - 30, "Hox, svp, lz, fkh, dac")
locus(NX - 70, NY + 22, "grim, skl (H99)", n=4)
text(NX - 40, NY + 64, "● H3K27me3", size=12, fill=PC, weight="bold")
# MAPK targets in nucleus
add(f'<rect x="{NX + 5}" y="{NY - 88}" width="80" height="38" rx="6" fill="#FFF3CC" stroke="{MAPK}" stroke-width="1.5"/>')
lines(NX + 45, NY - 73, ["sty, pnt,", "kek1 ↑"], size=11, gap=13)

# PVR receptors on membrane (top-left arc) with Pvf2 autocrine loop
def receptor(angle_deg, color, crossed=False):
    a = math.radians(angle_deg)
    x, y = CX + R * math.cos(a), CY + R * math.sin(a)
    ox, oy = math.cos(a), math.sin(a)
    # stem crossing membrane
    add(f'<line x1="{x - 16 * ox:.1f}" y1="{y - 16 * oy:.1f}" x2="{x + 14 * ox:.1f}" y2="{y + 14 * oy:.1f}" '
        f'stroke="{color}" stroke-width="6" stroke-linecap="round"/>')
    # Y arms
    for s in (-1, 1):
        px, py = -oy * s, ox * s
        add(f'<line x1="{x + 14 * ox:.1f}" y1="{y + 14 * oy:.1f}" x2="{x + 26 * ox + 9 * px:.1f}" '
            f'y2="{y + 26 * oy + 9 * py:.1f}" stroke="{color}" stroke-width="5" stroke-linecap="round"/>')
    if crossed:
        cross(x + 20 * ox, y + 20 * oy, s=13)
    return x, y, ox, oy

pvr = [receptor(a, LIN) for a in (-140, -120, -100)]
egfr = receptor(-45, "#9A9A9A", crossed=True)
text(egfr[0] + 40, egfr[1] - 20, "EGFR", size=14, weight="bold", fill="#9A9A9A", anchor="start")
text(egfr[0] + 40, egfr[1] - 3, "and other RTKs lost", size=12, fill=MUTED, anchor="start")
text(pvr[2][0] + 20, pvr[2][1] - 8, "PVR ↑", size=15, weight="bold", fill=LIN, anchor="start")
# Pvf2 ligands secreted and binding (autocrine loop, dashed = weakly supported)
for (x, y, ox, oy) in pvr:
    circle(x + 36 * ox, y + 36 * oy, 5.5, MAPK, stroke=INK, sw=0.8)
sx, sy = CX + R * math.cos(math.radians(-162)), CY + R * math.sin(math.radians(-162))
add(f'<path d="M{CX - 110} {CY - 40} Q {sx - 10:.0f} {sy + 10:.0f}, {sx - 45:.0f} {sy - 30:.0f} '
    f'T {pvr[0][0] + 36 * pvr[0][2] - 8:.0f} {pvr[0][1] + 36 * pvr[0][3] + 6:.0f}" '
    f'fill="none" stroke="{MAPK}" stroke-width="2.5" stroke-dasharray="7,5"/>')
text(sx - 55, sy + 5, "Pvf2", size=14, weight="bold", fill="#9A6B00", anchor="end")
text(sx - 55, sy + 22, "(autocrine)", size=12, fill=MUTED, anchor="end")
# signalling to nucleus
arrow(pvr[1][0] + 15, pvr[1][1] + 30, NX + 5, NY - 70, color=MAPK, sw=3, dash="7,5")
text(pvr[1][0] + 8, pvr[1][1] + 72, "Ras/MAPK", size=13, weight="bold", fill="#9A6B00", anchor="end")

# cell cycle badge
add(f'<rect x="{CX - 180}" y="{CY + 20}" width="112" height="56" rx="8" fill="#EAF1FB" stroke="{LIN}" stroke-width="1.5"/>')
lines(CX - 124, CY + 40, ["stg, CycB ↑", "fzr ↓: mitosis,", "not endocycles"], size=11, gap=13)

# lost external signals (left)
text(60, 560, "Tissue signals lost", size=15, weight="bold", anchor="start")
lost = ["Dpp", "Activin (Actβ)", "Insulin-like (Ilp6)", "JAK/STAT ligand (upd1)", "Ecdysone response (Hr3, E93)"]
for i, s in enumerate(lost):
    y = 592 + i * 30
    cross(72, y - 5, s=7, sw=2.5)
    text(90, y, s, size=13, anchor="start")
lines(60, 760, ["Also lost: injury/AMP, fat-body and", "crystal-cell programs"], size=12,
      anchor="start", fill=MUTED)

# Wolbachia box (right)
BX, BY = 1000, 510
add(f'<rect x="{BX}" y="{BY}" width="270" height="200" rx="10" fill="#FCF0F5" stroke="{WOL}" stroke-width="1.5"/>')
text(BX + 135, BY + 26, "Wolbachia", size=15, weight="bold", fill=WOL, style="italic")
lines(BX + 15, BY + 52, ["2 of 3 untreated lineages lost", "infection after the primary stage",
                         "Hypothesis: host division outpaces",
                         "bacterial replication (dilution)",
                         "Infected Dsim6B-wMel re-expresses",
                         "rpr, hid vs. cured Dsim6B"],
      size=12, gap=19, anchor="start")
add(f'<line x1="{BX + 15}" y1="{BY + 81}" x2="{BX + 255}" y2="{BY + 81}" stroke="{GRID}"/>')
add(f'<line x1="{BX + 15}" y1="{BY + 138}" x2="{BX + 255}" y2="{BY + 138}" stroke="{GRID}"/>')

# evidence key (right)
KX, KY = 1300, 510
add(f'<rect x="{KX}" y="{KY}" width="260" height="200" rx="10" fill="#FAFAFA" stroke="{GRID}" stroke-width="1.5"/>')
text(KX + 130, KY + 26, "Strength of evidence", size=15, weight="bold")
rows = [("Polycomb silencing of", "developmental loci", "strong", "#1B7F3B"),
        ("Proliferative cells shift", "toward cell lines", "moderate", "#B7791F"),
        ("PVR/Pvf2/MAPK loop", "(dashed)", "weak", "#B03A2E"),
        ("Wolbachia dilution", "(hypothesis)", "untested", MUTED)]
for i, (a, b, c, col) in enumerate(rows):
    y = KY + 58 + i * 38
    text(KX + 15, y, a, size=12, anchor="start")
    text(KX + 15, y + 14, b, size=11, anchor="start", fill=MUTED)
    add(f'<rect x="{KX + 180}" y="{y - 12}" width="68" height="20" rx="10" fill="{col}"/>')
    text(KX + 214, y + 2, c, size=11, fill="white", weight="bold")

# footnote
text(40, 975, "Pseudobulk DE with 4 lineages as replicates; Polycomb maps: embryo H3K27me3, S2 and BG3 "
     "Pc-state (modENCODE). Dashed lines mark weakly supported or untested links.",
     size=12, fill=MUTED, anchor="start")
# Wolbachia legend
wolb(1470, 972, 20, n=1, spread=0)
text(1482, 976, "Wolbachia", size=12, fill=MUTED, anchor="start", style="italic")

add("</svg>")
open("model_cartoon.svg", "w").write("\n".join(out))
print("wrote model_cartoon.svg")
