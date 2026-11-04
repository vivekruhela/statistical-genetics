#!/usr/bin/env python
"""Post 0 figure: the clickable series map (Python version; series_map.R is the R twin).

Reads posts/series.tsv and writes a standalone Plotly HTML page. Each post is a node; the
GWAS lane runs above the shared trunk and the expression lane below, and both meet in
fine-mapping. Published posts are filled and open when clicked; upcoming posts are rings.

Usage: python posts/00-series-map/series_map.py posts/series.tsv series-map.html [YYYY-MM-DD]
The optional date decides which posts count as published (default: today).
"""
import datetime as dt
import sys

import pandas as pd
import plotly.graph_objects as go

tsv, out = sys.argv[1:3]
today = sys.argv[3] if len(sys.argv) > 3 else dt.date.today().isoformat()
SITE = "https://vivekruhela.github.io/blog/"

# colours: the validated three-slot categorical set, light and dark steps
LANES = {  # lane: (label, light, dark)
    "shared": ("Shared trunk", "#2a78d6", "#3987e5"),
    "gwas": ("Genotype lane: GWAS and ROH", "#eb6834", "#d95926"),
    "expression": ("Expression lane", "#1baf7a", "#199e70"),
}
INK = {"light": ("#0b0b0b", "#52514e", "#c9c8c2", "#fcfcfb"), "dark": ("#f0efec", "#c3c2b7", "#4a4a47", "#1a1a19")}

p = pd.read_csv(tsv, sep="\t", dtype=str).fillna("")
p["x"], p["y"] = p.x.astype(float), p.y.astype(float)
p["published"] = p.publish <= today
p["url"] = [SITE + f"{d}-{s}/" if pub else "" for d, s, pub in zip(p.publish, p.slug, p.published)]
pos = p.set_index("id")[["x", "y"]]
# a label above a node would sit on an edge that leaves straight up; put those below
up = {a for r in p.itertuples() for a in r.after.split(";") if a and pos.x[a] == r.x and r.y > pos.y[a]}
p["textpos"] = ["bottom center" if i in up else "top center" for i in p.id]

fig = go.Figure()
# dependency edges, drawn first so the nodes sit on top
ex, ey = [], []
for r in p.itertuples():
    for a in filter(None, r.after.split(";")):
        ex += [pos.x[a], r.x, None]
        ey += [pos.y[a], r.y, None]
fig.add_trace(go.Scatter(x=ex, y=ey, mode="lines", line=dict(width=2, color=INK["light"][2]),
                         hoverinfo="skip", showlegend=False, name="edges"))

for lane, (label, col, _) in LANES.items():
    d = p[p.lane == lane]
    when = ["Published " if pub else "Coming " for pub in d.published]
    hover = [f"<b>{i}  {t}</b><br>{q}<br>{w}{dt.date.fromisoformat(day):%d %b %Y}" + ("<br><i>Click to open</i>" if pub else "")
             for i, t, q, w, day, pub in zip(d.id, d.title, d.question, when, d.publish, d.published)]
    fig.add_trace(go.Scatter(
        x=d.x, y=d.y, mode="markers+text", name=label, text=d.id, textposition=list(d.textpos),
        textfont=dict(size=13, color=INK["light"][0]),
        marker=dict(size=26, color=[col if pub else INK["light"][3] for pub in d.published],
                    line=dict(width=3, color=col)),
        customdata=d.url, hovertext=hover, hoverinfo="text"))

# part headers: above the trunk and genotype lane, below the expression lane
heads = p.groupby("part", sort=False).agg(x0=("x", "min"), x1=("x", "max"), lane=("lane", "first"))
for part, r in heads.iterrows():
    fig.add_annotation(x=(r.x0 + r.x1) / 2, y=-1.55 if r.lane == "expression" else 2.75,
                       text=f"<b>{part}</b>", showarrow=False, font=dict(size=12, color=INK["light"][1]))

fig.add_annotation(x=0, y=-0.11, xref="paper", yref="paper", xanchor="left", yanchor="top", showarrow=False,
                   text="Filled: published (click to open) · Ring: coming soon",
                   font=dict(size=12, color=INK["light"][1]))

fig.update_layout(
    template="plotly_white", height=470, margin=dict(l=10, r=10, t=10, b=75),
    paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
    font=dict(family="system-ui, -apple-system, Segoe UI, sans-serif", color=INK["light"][0]),
    xaxis=dict(visible=False, range=[-0.6, 13.6], fixedrange=True),
    yaxis=dict(visible=False, range=[-1.8, 3.0], fixedrange=True),
    legend=dict(orientation="h", x=0, y=-0.02, yanchor="top", font=dict(size=12), itemclick=False, itemdoubleclick=False),
    hoverlabel=dict(align="left"), dragmode=False,
)

# Click opens a published post in the page around the frame; dark mode follows the site's
# theme toggle when the frame is on the same site, else the OS setting.
post_script = """
var gd = document.getElementById('{plot_id}');
gd.on('plotly_click', function (e) {
  var u = e.points[0].customdata;
  if (u) { (window.top || window).location.href = u; }
});
function isDark() {
  try { var t = window.parent.document.documentElement.dataset.theme; if (t) return t === 'dark'; } catch (err) {}
  return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
}
function paint() {
  var d = isDark(), ink = d ? %s : %s, cols = d ? %s : %s;
  document.body.style.background = 'transparent';
  document.documentElement.style.colorScheme = d ? 'dark' : 'light';
  var n = gd.data.length;
  for (var k = 1; k < n; k++) {
    var pub = gd.data[k].customdata.map(function (u) { return !!u; });
    var c = cols[k - 1];
    Plotly.restyle(gd, {'marker.color': [pub.map(function (p) { return p ? c : ink[3]; })],
                        'marker.line.color': c, 'textfont.color': ink[0]}, [k]);
  }
  Plotly.restyle(gd, {'line.color': ink[2]}, [0]);
  var ann = (gd.layout.annotations || []).map(function (a) { a.font.color = ink[1]; return a; });
  Plotly.relayout(gd, {'font.color': ink[0], 'legend.font.color': ink[0], annotations: ann});
}
paint();
if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', paint);
try { new MutationObserver(paint).observe(window.parent.document.documentElement, {attributes: true, attributeFilter: ['data-theme']}); } catch (err) {}
""" % (list(INK["dark"]), list(INK["light"]), [c[2] for c in LANES.values()], [c[1] for c in LANES.values()])

fig.write_html(out, include_plotlyjs="cdn", full_html=True, post_script=post_script,
               config={"displayModeBar": False, "responsive": True})
print(f"wrote {out}: {len(p)} posts, {int(p.published.sum())} published as of {today}")
