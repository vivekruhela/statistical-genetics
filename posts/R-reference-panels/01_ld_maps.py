#!/usr/bin/env python
"""Post R, figure 1: the same stretch of DNA has a different LD map in different ancestries.

Left: r^2 between every pair of common variants in a 150 kb region, European people below the
diagonal and Yoruba people above it. Right: how fast LD fades with distance, averaged over all
ten teaching loci. Both groups have the same number of people (every k-th European is used),
because r^2 from a small sample is biased upwards by about 1/n.
  OUT/ld_decay.tsv, OUT/ld_map.tsv, OUT/fig1-ld-maps.html
Usage: python posts/R-reference-panels/01_ld_maps.py data results/post-R [L02]
R twin: 01_ld_maps.R (same numbers).
"""
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402
from figures import COLORS, SEQ, SURFACE, base_layout, save_figure  # noqa: E402

D, OUT = sys.argv[1:3]
SHOW = sys.argv[3] if len(sys.argv) > 3 else "L02"
BINS = np.array([0, 2, 5, 10, 20, 50, 100, 200]) * 1000   # distance bins for the decay curve
MAP_BP, MAP_MAX = 150000, 150                             # heatmap region and most variants shown

samples = pd.read_csv(f"{D}/cohort/samples.tsv", sep="\t")
loci = pd.read_csv(f"{D}/cohort/loci.tsv", sep="\t")
yri = samples.loc[samples.superpopulation == "AFR", "sample"].tolist()
eur_all = samples.loc[samples.superpopulation == "EUR", "sample"].tolist()
eur = [eur_all[i * len(eur_all) // len(yri)] for i in range(len(yri))]
groups = {"European": eur, "Yoruba": yri}


def maf(g):
    af = g.mean(1) / 2
    return np.minimum(af, 1 - af)


decay, ld_map = [], None
for lid in loci.locus:
    info, g_all, ids = sg.read_traw(f"{D}/cohort/genotypes/loci/{lid}.traw")
    col = {s: k for k, s in enumerate(ids)}
    for grp, people in groups.items():
        g = g_all[:, [col[s] for s in people]]
        keep = maf(g) >= 0.05
        g, pos = g[keep], info.pos.to_numpy()[keep]
        r2 = sg.ld_matrix(g) ** 2
        dist = np.abs(np.subtract.outer(pos, pos))
        up = np.triu(np.ones_like(r2, dtype=bool), 1) & (dist < BINS.max())
        b = np.searchsorted(BINS, dist[up], side="right")
        d = pd.DataFrame({"bin": b, "r2": r2[up]}).groupby("bin").r2.agg(["sum", "size"]).reset_index()
        d.insert(0, "group", grp)
        decay.append(d)
    if lid == SHOW:
        centre = info.pos[info.id == loci.anchor_lead_id[loci.locus == lid].iloc[0]].iloc[0]
        ge, gy = g_all[:, [col[s] for s in eur]], g_all[:, [col[s] for s in yri]]
        idx = np.flatnonzero((maf(ge) >= 0.05) & (maf(gy) >= 0.05) & (np.abs(info.pos - centre) <= MAP_BP / 2))
        idx = idx[np.unique(np.floor(np.linspace(0, len(idx) - 1, min(MAP_MAX, len(idx)))).astype(int))]
        m = sg.ld_matrix(ge[idx]) ** 2
        r_yri = sg.ld_matrix(gy[idx]) ** 2
        upper = np.triu_indices(len(idx), 1)
        m[upper] = r_yri[upper]
        k = len(idx)
        vid = info.id.to_numpy()[idx]
        # column-major, as R writes it: i (row) varies fastest
        ld_map = pd.DataFrame({"locus": lid, "i": np.tile(np.arange(1, k + 1), k), "j": np.repeat(np.arange(1, k + 1), k),
                               "id_i": np.tile(vid, k), "id_j": np.repeat(vid, k), "r2": m.T.reshape(-1)})

decay = pd.concat(decay).groupby(["group", "bin"]).agg(s=("sum", "sum"), n=("size", "sum")).reset_index()
decay = pd.DataFrame({"group": decay.group, "from_kb": BINS[decay.bin - 1] / 1000, "to_kb": BINS[decay.bin] / 1000,
                      "mean_r2": decay.s / decay.n, "n_pairs": decay.n}).sort_values(["group", "from_kb"])
sg.write_tsv(decay, f"{OUT}/ld_decay.tsv")
sg.write_tsv(ld_map, f"{OUT}/ld_map.tsv")
print(decay.to_string(index=False))

# ---- figure ----
k = ld_map.i.max()
z = ld_map.r2.to_numpy().reshape(k, k).T
vid = ld_map.id_i.to_numpy()[:k]
hover = [[f"{vid[i]}<br>{vid[j]}<br>{'Yoruba' if i < j else 'European'}: r² = {z[i, j]:.2f}" for j in range(k)]
         for i in range(k)]
fig = make_subplots(rows=1, cols=2, column_widths=[0.5, 0.5], horizontal_spacing=0.17)
fig.add_trace(go.Heatmap(x=np.arange(1, k + 1), y=np.arange(1, k + 1), z=z, zmin=0, zmax=1, text=hover,
                         hoverinfo="text", colorscale=SEQ, meta="seq",
                         colorbar=dict(title=dict(text="r²", side="top"), len=0.6, thickness=10, x=0.425)), 1, 1)
for grp, color in [("European", COLORS["aqua"]), ("Yoruba", COLORS["orange"])]:
    d = decay[decay.group == grp]
    fig.add_trace(go.Scatter(
        x=(d.from_kb + d.to_kb) / 2, y=d.mean_r2, mode="lines+markers", name=f"{grp} ({len(groups[grp])} people)",
        line=dict(width=2, color=color), marker=dict(size=8, line=dict(width=2, color=SURFACE)), hoverinfo="text",
        text=[f"{grp}, {a:g} to {b:g} kb apart<br>mean r² = {r:.3f} ({n:,} pairs)"
              for a, b, r, n in zip(d.from_kb, d.to_kb, d.mean_r2, d.n_pairs)]), 1, 2)
fig.update_layout(**base_layout(
    height=440, legend=dict(x=0.98, xanchor="right", y=0.98),
    xaxis=dict(title=f"Below the diagonal: European · above: Yoruba<br>({k} common variants across {MAP_BP // 1000} kb at {SHOW})", showticklabels=False,
               showgrid=False, fixedrange=True, scaleanchor="y", constrain="domain"),
    yaxis=dict(showticklabels=False, showgrid=False, fixedrange=True, autorange="reversed", constrain="domain"),
    xaxis2=dict(title="Distance between two variants (kb, log scale)", type="log", fixedrange=True),
    yaxis2=dict(title="Average r²", rangemode="tozero", fixedrange=True),
    annotations=[
        dict(x=0.02, y=1.06, xref="paper", yref="paper", showarrow=False, xanchor="left", font=dict(size=13),
             text="<b>One region, two LD maps</b>"),
        dict(x=0.59, y=1.06, xref="paper", yref="paper", showarrow=False, xanchor="left", font=dict(size=13),
             text="<b>LD fades faster in Yoruba genomes</b>")]))
save_figure(fig, f"{OUT}/fig1-ld-maps.html", min_width=760)
