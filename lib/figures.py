"""Figure helpers for the posts (Python track); lib/figures.R is the R twin.

Figures are drawn in the light theme; lib/figure_theme.js switches them to dark where the
reader's site or system is dark. Pages load plotly.js from a CDN, so each file stays small.
"""
import os

COLORS = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
          "magenta": "#e87ba4", "green": "#008300"}
INK, INK2, LINE, SURFACE = "#0b0b0b", "#52514e", "#c9c8c2", "#fcfcfb"
SEQ = [[0, SURFACE], [0.15, "#cde2fb"], [0.5, "#6da7ec"], [1, "#104281"]]
FONT = dict(family="system-ui, -apple-system, Segoe UI, sans-serif", size=12, color=INK)
_JS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figure_theme.js")


def base_layout(**kw):
    """Layout shared by every figure: transparent background, quiet grid, no zoom."""
    out = dict(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", font=FONT, dragmode=False,
               hoverlabel=dict(align="left"), margin=dict(l=60, r=20, t=40, b=50))
    out.update(kw)
    return out


def save_figure(fig, path, min_width=640):
    """Write a standalone HTML page that follows the site's light or dark theme."""
    with open(_JS) as fh:
        js = fh.read()
    post = js + "\nstatgenTheme(document.getElementById('{plot_id}'), %d);\n" % min_width
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.write_html(path, include_plotlyjs="cdn", full_html=True, post_script=post,
                   config={"displayModeBar": False, "responsive": True})
    print(f"wrote {path}")
