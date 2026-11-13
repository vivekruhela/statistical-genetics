#!/usr/bin/env python
"""Post R, figure 2: how well each reference panel filled in the variants we hid.

For every hidden variant (not on the pretend chip from 02_imputation.sh) and every panel, the
squared correlation between the imputed dosage and the real genotype, separately in the
European and the Yoruba people of the cohort. That r^2 is the usual measure of imputation
accuracy: 1 is perfect, 0 means the panel knew nothing about this variant in these people.
  OUT/imputation_r2.tsv.gz       one row per variant, group and panel
  OUT/imputation_summary.tsv     average r^2 by group, panel and allele frequency
  OUT/fig2-imputation.html
Usage: python posts/R-reference-panels/03_imputation_accuracy.py data results/post-R/imputation results/post-R
R twin: 03_imputation_accuracy.R (same numbers).
"""
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402
from figures import COLORS, SURFACE, base_layout, save_figure  # noqa: E402

D, IMP, OUT = sys.argv[1:4]
PANELS = {"EUR": "European panel", "AFR": "African panel", "BOTH": "Both panels together"}
MAF_BINS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
BIN_LABELS = [f"{100 * a:g}-{100 * b:g}%" for a, b in zip(MAF_BINS[:-1], MAF_BINS[1:])]
GROUPS = {"EUR": "European", "AFR": "Yoruba"}

samples = pd.read_csv(f"{D}/cohort/samples.tsv", sep="\t")
loci = pd.read_csv(f"{D}/cohort/loci.tsv", sep="\t")


def row_r2(a, b):
    """Squared correlation between matching rows of a and b (0 where either row is constant)."""
    a = a - a.mean(1, keepdims=True)
    b = b - b.mean(1, keepdims=True)
    num, den = (a * b).sum(1), np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (num / den) ** 2, 0.0)


res = []
for lid in loci.locus:
    info, g, gids = sg.read_traw(f"{D}/cohort/genotypes/loci/{lid}.traw")
    gcol = {s: k for k, s in enumerate(gids)}
    grow = {v: k for k, v in enumerate(info.id)}
    chip = set(pd.read_csv(f"{IMP}/chip/{lid}.sites.txt", sep="\t", header=None)[2])
    for panel, plabel in PANELS.items():
        ds = pd.read_csv(f"{IMP}/imputed/{lid}.{panel}.ds.tsv.gz", sep="\t", index_col=0)
        ids = [v for v in ds.index if v not in chip and v in grow]
        for grp, glabel in GROUPS.items():
            people = samples.loc[samples.superpopulation == grp, "sample"].tolist()
            truth = g[np.ix_([grow[v] for v in ids], [gcol[s] for s in people])]
            af = truth.mean(1) / 2
            maf = np.minimum(af, 1 - af)
            keep = maf >= MAF_BINS[0]
            dmat = ds.loc[ids, people].to_numpy(float)
            res.append(pd.DataFrame({"locus": lid, "variant_id": np.array(ids)[keep], "group": glabel, "panel": plabel,
                                     "maf": maf[keep], "r2": row_r2(dmat[keep], truth[keep])}))
res = pd.concat(res, ignore_index=True)
sg.write_tsv(res, f"{OUT}/imputation_r2.tsv.gz")
res["bin"] = np.clip(np.searchsorted(MAF_BINS, res.maf, side="right") - 1, 0, len(BIN_LABELS) - 1)
summ = (res.groupby(["group", "panel", "bin"])
        .agg(n_variants=("r2", "size"), mean_r2=("r2", "mean"), share_r2_above_0_8=("r2", lambda r: (r >= 0.8).mean()))
        .reset_index().sort_values(["group", "panel", "bin"]))
summ.insert(2, "maf_bin", [BIN_LABELS[b] for b in summ.bin])
summ = summ.drop(columns="bin").rename(columns={"share_r2_above_0_8": "share_r2_above_0.8"})
sg.write_tsv(summ, f"{OUT}/imputation_summary.tsv")
print(summ.to_string(index=False))

# ---- figure: one panel per group of people, one line per reference panel ----
cols = [COLORS["aqua"], COLORS["orange"], COLORS["blue"]]
fig = make_subplots(rows=1, cols=2, shared_yaxes=True, horizontal_spacing=0.08)
for k, glabel in enumerate(GROUPS.values()):
    for j, (panel, plabel) in enumerate(PANELS.items()):
        d = summ[(summ.group == glabel) & (summ.panel == plabel)]
        fig.add_trace(go.Scatter(
            x=d.maf_bin, y=d.mean_r2, mode="lines+markers", name=plabel, legendgroup=plabel, showlegend=k == 0,
            line=dict(width=2, color=cols[j], dash="dot" if panel == "BOTH" else "solid"), marker=dict(size=8, line=dict(width=2, color=SURFACE)), hoverinfo="text",
            text=[f"{glabel} people, {plabel}<br>variants with MAF {b}: {n:,}<br>average r² = {r:.2f}"
                  f"<br>imputed well (r² ≥ 0.8): {100 * s:.0f}%"
                  for b, n, r, s in zip(d.maf_bin, d.n_variants, d.mean_r2, d["share_r2_above_0.8"])]), 1, k + 1)
xa = dict(title="Minor allele frequency of the hidden variant", categoryorder="array", categoryarray=BIN_LABELS,
          fixedrange=True)
fig.update_layout(**base_layout(
    height=430, legend=dict(orientation="h", x=0, y=-0.24), margin=dict(l=60, r=20, t=50, b=110),
    xaxis=xa, xaxis2=xa, yaxis=dict(title="Imputation accuracy (average r²)", range=[0, 1.02], fixedrange=True),
    yaxis2=dict(fixedrange=True),
    annotations=[dict(x=0, y=1.08, xref="paper", yref="paper", xanchor="left", showarrow=False,
                      text="<b>Filling in European genomes</b>", font=dict(size=13)),
                 dict(x=0.54, y=1.08, xref="paper", yref="paper", xanchor="left", showarrow=False,
                      text="<b>Filling in Yoruba genomes</b>", font=dict(size=13))]))
save_figure(fig, f"{OUT}/fig2-imputation.html", min_width=700)
