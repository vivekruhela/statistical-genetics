#!/usr/bin/env python
"""Post 0 concept figure: one causal variant lights up a whole tower of hits.

A toy locus, simulated from scratch so anyone can rerun it: 120 variants in strong LD, 1,500
people, and one causal variant (number 60) with a modest effect on a trait. The top panel is
what a GWAS sees (-log10 p for every variant, coloured by LD with the causal variant). The bottom
panel is what fine-mapping says (posterior inclusion probability, single causal variant, Wakefield
approximate Bayes factors) with the 95% credible set marked.

Usage: python posts/00-series-map/ld_tower.py ld-tower.html [toy_locus.tsv]
"""
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import norm

out = sys.argv[1]
SEED, N, M, RHO, CAUSAL, BETA = 178, 1500, 120, 0.985, 60, 0.22

# genotypes: two haplotypes per person from a correlated latent field, so neighbours are in LD
rng = np.random.default_rng(SEED)
L = np.linalg.cholesky(RHO ** np.abs(np.subtract.outer(np.arange(M), np.arange(M))))
thresh = rng.uniform(-0.6, 0.9, M)
hap = lambda: ((L @ rng.standard_normal((M, N))).T > thresh).astype(int)
G = (hap() + hap()).astype(float)
g = G[:, CAUSAL]
y = BETA * (g - g.mean()) / g.std() + rng.standard_normal(N)

# what a GWAS does: one regression per variant
Gc, yc = G - G.mean(0), y - y.mean()
b = Gc.T @ yc / (Gc ** 2).sum(0)
se = np.sqrt(((yc[:, None] - Gc * b) ** 2).sum(0) / (N - 2) / (Gc ** 2).sum(0))
z = b / se
mlogp = -np.log10(2 * norm.sf(np.abs(z)))
r2 = np.corrcoef(G.T)[CAUSAL] ** 2

# what fine-mapping does: Wakefield ABF, one causal variant, prior effect sd 0.15
V, W = se ** 2, 0.15 ** 2
lbf = 0.5 * np.log(V / (V + W)) + 0.5 * z ** 2 * W / (V + W)
pip = np.exp(lbf - lbf.max())
pip /= pip.sum()
order = np.argsort(-pip)
cs = set(order[: np.searchsorted(np.cumsum(pip[order]), 0.95) + 1])

d = pd.DataFrame({"variant": np.arange(1, M + 1), "kb": np.arange(M) * 10, "r2": r2, "mlog10p": mlogp,
                  "pip": pip, "in_cs": [i in cs for i in range(M)], "causal": np.arange(M) == CAUSAL})
if len(sys.argv) > 2:
    d.to_csv(sys.argv[2], sep="\t", index=False, float_format="%.5g")
top = int(np.argmax(mlogp))
print(f"{(mlogp > -np.log10(5e-8)).sum()} variants pass 5e-8; top hit is variant {top + 1}, "
      f"causal is {CAUSAL + 1}; 95% credible set has {len(cs)} variants: {sorted(int(i) + 1 for i in cs)}")

BLUE, ORANGE, INK = "#2a78d6", "#eb6834", ("#0b0b0b", "#52514e", "#c9c8c2", "#fcfcfb")
SEQ = [[0, "#dbe7f6"], [0.5, "#7fa9e0"], [1, "#1c4f8f"]]   # one hue, light to dark: r2 0 to 1
hover = [f"Variant {v} at {k} kb<br>LD with the causal variant: r² = {a:.2f}<br>-log10 p = {p:.1f}<br>PIP = {q:.2f}"
         for v, k, a, p, q in zip(d.variant, d.kb, d.r2, d.mlog10p, d.pip)]

fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.62, 0.38], vertical_spacing=0.09,
                    subplot_titles=("What a GWAS sees: every variant near the causal one looks guilty",
                                    "What fine-mapping says: probability that each variant is the causal one"))
fig.add_trace(go.Scatter(x=d.kb, y=d.mlog10p, mode="markers", hovertext=hover, hoverinfo="text", showlegend=False,
                         marker=dict(size=9, color=d.r2, colorscale=SEQ, cmin=0, cmax=1,
                                     line=dict(width=1, color=INK[3]),
                                     colorbar=dict(title=dict(text="LD (r²) with<br>causal variant", side="top"),
                                                   len=0.55, y=0.72, thickness=12, tickvals=[0, 0.5, 1]))),
              row=1, col=1)
c = d.iloc[CAUSAL]
fig.add_trace(go.Scatter(x=[c.kb], y=[c.mlog10p], mode="markers", hoverinfo="skip", showlegend=False,
                         marker=dict(size=17, symbol="circle-open", line=dict(width=3, color=ORANGE))), row=1, col=1)
fig.add_hline(y=-np.log10(5e-8), line=dict(width=1, dash="dot", color=INK[2]), row=1, col=1)
fig.add_annotation(x=c.kb, y=c.mlog10p, text="the real cause<br>(we planted it)", ax=-90, ay=-25, arrowcolor=ORANGE,
                   font=dict(color=ORANGE, size=12), row=1, col=1)
fig.add_annotation(x=d.kb[top], y=d.mlog10p[top], text="top hit", ax=55, ay=-30, row=1, col=1, font=dict(size=12))
fig.add_annotation(x=d.kb.max(), y=-np.log10(5e-8), text="genome-wide significance", showarrow=False, xanchor="right",
                   yanchor="bottom", font=dict(size=11), row=1, col=1)

fig.add_trace(go.Bar(x=d.kb, y=d.pip, hovertext=hover, hoverinfo="text", showlegend=False, width=7,
                     marker=dict(color=[BLUE if a else INK[2] for a in d.in_cs], line=dict(width=0))), row=2, col=1)
fig.add_annotation(x=d.kb[order[0]], y=d.pip.max(), text=f"95% credible set: {len(cs)} variants, one of them the real cause",
                   ax=70, ay=-10, xanchor="left", font=dict(size=12, color=BLUE), arrowcolor=BLUE, row=2, col=1)

fig.update_xaxes(title_text="Position in the region (kb)", row=2, col=1)
fig.update_xaxes(showgrid=False, zeroline=False, fixedrange=True, range=[-20, M * 10 + 10])
fig.update_yaxes(title_text="-log10 p", row=1, col=1, gridcolor="rgba(128,128,128,0.15)", zeroline=False, fixedrange=True)
fig.update_yaxes(title_text="PIP", range=[0, 1], row=2, col=1, gridcolor="rgba(128,128,128,0.15)", zeroline=False,
                 fixedrange=True)
fig.update_annotations(selector=dict(text=fig.layout.annotations[0].text), font=dict(size=13))
fig.update_layout(height=560, margin=dict(l=60, r=20, t=40, b=50), paper_bgcolor="rgba(0,0,0,0)",
                  plot_bgcolor="rgba(0,0,0,0)", bargap=0, dragmode=False, hoverlabel=dict(align="left"),
                  font=dict(family="system-ui, -apple-system, Segoe UI, sans-serif", size=12, color=INK[0]))

# dark mode follows the site's theme toggle when framed on the same site, else the OS setting;
# on a narrow screen the figure keeps its size and scrolls sideways
post_script = """
var gd = document.getElementById('{plot_id}');
function isDark() {
  try { var t = window.parent.document.documentElement.dataset.theme; if (t) return t === 'dark'; } catch (err) {}
  return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
}
function paint() {
  var d = isDark(), ink = d ? %s : %s;
  document.body.style.background = 'transparent';
  document.documentElement.style.colorScheme = d ? 'dark' : 'light';
  var ann = (gd.layout.annotations || []).map(function (a) {
    if (!a.font || !a.font.color || a.font.color === %s || a.font.color === %s) { a.font = Object.assign({}, a.font, {color: ink[0]}); }
    return a; });
  Plotly.relayout(gd, {'font.color': ink[0], annotations: ann, 'shapes[0].line.color': ink[2]});
  Plotly.restyle(gd, {'marker.line.color': ink[3]}, [0]);
  Plotly.restyle(gd, {'marker.color': [gd.data[2].marker.color.map(function (c) { return c === %s ? c : ink[2]; })]}, [2]);
}
if (window.innerWidth < 640) { document.body.style.overflowX = 'auto'; document.body.style.margin = '0'; Plotly.relayout(gd, {width: 640}); }
paint();
if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', paint);
try { new MutationObserver(paint).observe(window.parent.document.documentElement, {attributes: true, attributeFilter: ['data-theme']}); } catch (err) {}
""" % ('["#f0efec", "#c3c2b7", "#4a4a47", "#1a1a19"]', list(INK), '"#0b0b0b"', '"#f0efec"', f'"{BLUE}"')
fig.write_html(out, include_plotlyjs="cdn", full_html=True, post_script=post_script,
               config={"displayModeBar": False, "responsive": True})
print(f"wrote {out}")
