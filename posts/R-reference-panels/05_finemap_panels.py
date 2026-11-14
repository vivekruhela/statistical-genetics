#!/usr/bin/env python
"""Post R, analysis 3b: fine-mapping the same GWAS with LD from different reference panels.

Python twin of 05_finemap_panels.R, with lib/susie_rss.py in place of susieR. For each window
from 04_gwas_ld, SuSiE-RSS fine-maps the European GWAS z-scores with LD from each of these, and
the answer key in data/truth scores the result:
  insample  the GWAS people themselves (the ideal, and rarely available)
  eur       the European people of our cohort (a matched panel, like 1000 Genomes EUR)
  eur50     50 of them (a matched but small panel)
  afr       the African-ancestry template people (a mismatched panel)
  mixed     European cohort people plus the African templates (a panel that mixes ancestries)
  flip      in-sample LD, but one variant's alleles swapped in the summary statistics
A variant that does not vary in a panel is dropped for that panel, as it would be in practice.
Each fit is also checked with ld_check() (lib/statgen.py): how well the panel's LD predicts each
z-score from its neighbours. The mismatch score s is near 0 when panel and GWAS agree.
  OUT/finemap_variants.tsv.gz, OUT/finemap_sets.tsv, OUT/finemap_summary.tsv
  OUT/fig3-finemap-panels.html, OUT/fig4-z-check.html, OUT/fig5-panel-scorecard.html
Usage: python posts/R-reference-panels/05_finemap_panels.py data results/post-R [L05_trait_q]
"""
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import norm

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402
from figures import COLORS, INK2, LINE, base_layout, save_figure  # noqa: E402
from susie_rss import SusieError, susie_rss  # noqa: E402

D, OUT = sys.argv[1:3]
SHOW = sys.argv[3] if len(sys.argv) > 3 else None
MIN_AF, FLAG_SD = 0.01, 4      # a z-score more than 4 SD from its prediction is flagged
PANELS = {"insample": "GWAS people (in-sample)", "eur": "European panel", "eur50": "Small European panel",
          "afr": "African panel", "mixed": "Mixed panel", "flip": "In-sample, one allele flipped"}

windows = pd.read_csv(f"{OUT}/windows.tsv", sep="\t")
samples = pd.read_csv(f"{D}/cohort/samples.tsv", sep="\t")
loci = pd.read_csv(f"{D}/cohort/loci.tsv", sep="\t")
truth = pd.read_csv(f"{D}/truth/truth.tsv", sep="\t")
eur = samples.loc[samples.superpopulation == "EUR", "sample"].tolist()
eur50 = [eur[i * len(eur) // 50] for i in range(50)]


def panel_genotypes(lid, ids):
    info, g, gids = sg.read_traw(f"{D}/cohort/genotypes/loci/{lid}.traw")
    g = g[[{v: k for k, v in enumerate(info.id)}[v] for v in ids]]
    col = {s: k for k, s in enumerate(gids)}
    h = pd.read_csv(f"{D}/work/gwas_templates/AFR/{lid}.hap.gz", sep="\t", index_col=0).reindex(ids)
    H = h.to_numpy(float)
    afr = np.nan_to_num(H[:, 0::2] + H[:, 1::2])   # a variant missing from the panel: treat it as not varying
    ge = g[:, [col[s] for s in eur]]
    return {"eur": ge, "eur50": g[:, [col[s] for s in eur50]], "afr": afr, "mixed": np.hstack([ge, afr])}


vars_, sets, summ = [], [], []
for w in windows.itertuples():
    ss = pd.read_csv(f"{D}/gwas/EUR/{w.trait}.sumstats.tsv.gz", sep="\t")
    ss = ss[(ss.locus == w.locus) & (ss.af >= MIN_AF) & (ss.af <= 1 - MIN_AF) & (ss.pos >= w.start)
            & (ss.pos <= w.end)].reset_index(drop=True)
    R_in, ld_ids = sg.read_ld(f"{OUT}/ld/insample/{w.locus}_{w.trait}")
    assert ld_ids == ss.variant_id.tolist()
    ep = "disease_liability" if w.trait == "disease" else "trait_q"
    causal = set(truth.variant_id[(truth.locus == w.locus) & (truth.endpoint == ep)]) & set(ss.variant_id)
    pg = panel_genotypes(w.locus, ss.variant_id.tolist())
    flip_at = np.argsort(-ss.z.abs().to_numpy(), kind="stable")[9]   # the variant with the 10th strongest signal
    for pn, plabel in PANELS.items():
        R = R_in if pn in ("insample", "flip") else sg.ld_matrix(pg[pn])
        z = ss.z.to_numpy().copy()
        if pn == "flip":
            z[flip_at] = -z[flip_at]
        n_panel = w.n_gwas if pn in ("insample", "flip") else pg[pn].shape[1]
        ok = ~np.isnan(np.diag(R))
        R, z, ids = R[np.ix_(ok, ok)], z[ok], ss.variant_id.to_numpy()[ok]
        try:
            fit = susie_rss(z, R, w.n_gwas)
            status = "converged" if fit["converged"] else "did not converge"
            pip, cs = fit["pip"], fit["sets"]
        except SusieError as e:
            status = "stopped: z and LD disagree" if "prior variance" in str(e) else "stopped: error"
            pip, cs = np.full(len(z), np.nan), []
        chk = sg.ld_check(z, R, w.n_gwas)
        # number the credible sets from left to right, by the position of their top variant
        cs = sorted(cs, key=lambda c: c["variants"][int(np.argmax(pip[c["variants"]]))])
        cs_of = np.full(len(z), np.nan)
        is_causal = np.isin(ids, list(causal))
        for j, c in enumerate(cs, start=1):
            idx = np.array(c["variants"])
            cs_of[idx] = j
            sets.append(dict(locus=w.locus, trait=w.trait, panel=plabel, cs=j, size=len(idx), purity=c["purity"],
                             coverage=c["coverage"], holds_causal=bool(is_causal[idx].any()),
                             lead_variant=ids[idx][np.argmax(pip[idx])]))
        vars_.append(pd.DataFrame(dict(
            locus=w.locus, trait=w.trait, panel=plabel, variant_id=ids, pos=ss.pos.to_numpy()[ok], z=z, pip=pip,
            cs=pd.array(cs_of, dtype="Int64"), causal=is_causal,
            flipped=(pn == "flip") & (ids == ss.variant_id[flip_at]),
            expected_z=chk.expected.to_numpy(), z_std_diff=chk.z_std_diff.to_numpy())))
        true_cs = [bool(is_causal[c["variants"]].any()) for c in cs]
        in_cs = set(ids[[v for c in cs for v in c["variants"]]])
        summ.append(dict(
            locus=w.locus, trait=w.trait, panel=plabel, panel_people=n_panel, n_variants=int(ok.sum()),
            n_dropped=int((~ok).sum()), status=status, n_cs=len(cs), n_cs_true=sum(true_cs),
            n_cs_false=len(cs) - sum(true_cs), n_causal=len(causal), n_causal_found=len(causal & in_cs),
            causal_pip_sum=np.nan if np.isnan(pip).all() else pip[is_causal].sum(), mismatch_s=chk.s.iloc[0],
            n_flagged=int((chk.z_std_diff.abs() > FLAG_SD).sum()), max_abs_z_std_diff=chk.z_std_diff.abs().max()))
    print(f"{w.locus} {w.trait} done")
vars_ = pd.concat(vars_, ignore_index=True)
sets = pd.DataFrame(sets)
summ = pd.DataFrame(summ)
sg.write_tsv(vars_, f"{OUT}/finemap_variants.tsv.gz")
sg.write_tsv(sets, f"{OUT}/finemap_sets.tsv")
sg.write_tsv(summ, f"{OUT}/finemap_summary.tsv")
print(summ[["locus", "panel", "n_dropped", "status", "n_cs", "n_cs_true", "n_cs_false", "n_causal", "n_causal_found",
            "mismatch_s", "n_flagged"]].round(3).to_string(index=False))

# ---- figure 3: one locus, one row per panel ----------------------------------------
SHOW = SHOW or f"{windows.locus[0]}_{windows.trait[0]}"
show = windows[windows.locus + "_" + windows.trait == SHOW].iloc[0]
fig_panels = list(PANELS)[:5]
cols = [COLORS["blue"], COLORS["aqua"], COLORS["yellow"], COLORS["orange"], COLORS["magenta"]]
v = vars_[(vars_.locus == show.locus) & (vars_.trait == show.trait)]
ring = dict(size=13, color="rgba(0,0,0,0)", line=dict(width=2.5, color=COLORS["orange"]))
fig3 = make_subplots(rows=6, cols=1, shared_xaxes=True, row_heights=[0.2] + [0.16] * 5, vertical_spacing=0.022)
gw = v[v.panel == PANELS["insample"]]
mlogp = -(norm.logsf(gw.z.abs()) + np.log(2)) / np.log(10)
fig3.add_trace(go.Scatter(x=gw.pos / 1e6, y=mlogp, mode="markers", marker=dict(size=6, color=LINE), showlegend=False,
                          hoverinfo="text", text=[f"{a}<br>z = {b:.1f}" for a, b in zip(gw.variant_id, gw.z)]), 1, 1)
fig3.add_trace(go.Scatter(x=gw.pos[gw.causal] / 1e6, y=mlogp[gw.causal], mode="markers", marker=ring,
                          hoverinfo="skip", showlegend=False), 1, 1)
fig3.update_yaxes(title_text="GWAS<br>-log10 p", row=1, col=1)
annotations = []
for j, pn in enumerate(fig_panels):
    d = v[v.panel == PANELS[pn]]
    sm = summ[(summ.locus == show.locus) & (summ.trait == show.trait) & (summ.panel == PANELS[pn])].iloc[0]
    fig3.add_trace(go.Scatter(
        x=d.pos / 1e6, y=d.pip, mode="markers", showlegend=False, hoverinfo="text",
        marker=dict(size=7, color=[LINE if pd.isna(c) else cols[j] for c in d.cs]),
        text=[f"{a}<br>PIP = {p:.2f}" + ("" if pd.isna(c) else f"<br>credible set {c}") + ("<br><b>a real cause</b>" if k else "")
              for a, p, c, k in zip(d.variant_id, d.pip, d.cs, d.causal)]), j + 2, 1)
    fig3.add_trace(go.Scatter(x=d.pos[d.causal] / 1e6, y=d.pip[d.causal], mode="markers", marker=ring,
                              hoverinfo="skip", showlegend=False), j + 2, 1)
    fig3.update_yaxes(title_text="PIP", range=[-0.05, 1.1], row=j + 2, col=1)
    yref = "y domain" if j + 2 == 1 else f"y{j + 2} domain"
    annotations.append(dict(x=0, y=1, xref="x domain", yref=yref, xanchor="left", yanchor="bottom", showarrow=False,
                            font=dict(size=12),
                            text=f"<b>{PANELS[pn]}</b> ({sm.panel_people:,} people): {sm.n_cs} credible "
                                 f"set{'' if sm.n_cs == 1 else 's'}, {sm.n_cs_false} false"))
chrom = loci.chr[loci.locus == show.locus].iloc[0]
fig3.update_xaxes(fixedrange=True)
fig3.update_yaxes(fixedrange=True)
fig3.update_xaxes(title_text=f"Position on {chrom} (Mb)", row=6, col=1)
fig3.update_layout(**base_layout(height=860, margin=dict(l=70, r=20, t=30, b=50), annotations=annotations))
save_figure(fig3, f"{OUT}/fig3-finemap-panels.html", min_width=640)

# ---- figure 4: does the panel's LD fit the z-scores? -------------------------------------
chk_panels = ["insample", "afr", "flip"]
fig4 = make_subplots(rows=1, cols=3, horizontal_spacing=0.07)
annotations = []
for j, pn in enumerate(chk_panels):
    d = v[v.panel == PANELS[pn]]
    flag = d.z_std_diff.abs() > FLAG_SD
    lim = max(d.z.abs().max(), d.expected_z.abs().max())
    fig4.add_trace(go.Scatter(
        x=d.expected_z, y=d.z, mode="markers", showlegend=False, hoverinfo="text",
        marker=dict(size=6, color=[COLORS["orange"] if f else LINE for f in flag]),
        text=[f"{a}<br>observed z = {b:.1f}<br>predicted from neighbours = {c:.1f}<br>{abs(e):.1f} SD apart"
              + ("<br><b>the flipped variant</b>" if fl else "")
              for a, b, c, e, fl in zip(d.variant_id, d.z, d.expected_z, d.z_std_diff, d.flipped)]), 1, j + 1)
    fig4.add_trace(go.Scatter(x=[-lim, lim], y=[-lim, lim], mode="lines", hoverinfo="skip", showlegend=False,
                              line=dict(width=1, color=INK2, dash="dot")), 1, j + 1)
    fig4.update_xaxes(title_text="z predicted from neighbours", row=1, col=j + 1)
    s = summ[(summ.locus == show.locus) & (summ.trait == show.trait) & (summ.panel == PANELS[pn])].mismatch_s.iloc[0]
    annotations.append(dict(x=0.5, y=1.02, xref=f"x{'' if j == 0 else j + 1} domain", yref="paper", yanchor="bottom",
                            showarrow=False, font=dict(size=12),
                            text=f"<b>{PANELS[pn]}</b><br>s = {s:.3f}, {int(flag.sum())} flagged"))
fig4.update_yaxes(title_text="Observed z", row=1, col=1)
fig4.update_xaxes(fixedrange=True)
fig4.update_yaxes(fixedrange=True)
fig4.update_layout(**base_layout(height=400, margin=dict(l=60, r=20, t=70, b=50), annotations=annotations))
save_figure(fig4, f"{OUT}/fig4-z-check.html", min_width=760)

# ---- figure 5: the scorecard over all windows --------------------------------------
sc = summ.groupby("panel", sort=False).agg(found=("n_causal_found", "sum"), causal=("n_causal", "sum"),
                                           false_cs=("n_cs_false", "sum"),
                                           stopped=("status", lambda s: s.str.startswith("stopped").sum()))
sc = sc.loc[list(PANELS.values())[::-1]].reset_index()
lab = [f"{p}<br>real causes caught: {f} of {c}<br>false credible sets: {x}"
       + (f"<br>SuSiE refused {s} window{'' if s == 1 else 's'}" if s else "")
       for p, f, c, x, s in zip(sc.panel, sc.found, sc.causal, sc.false_cs, sc.stopped)]
fig5 = make_subplots(rows=1, cols=2, shared_yaxes=True, horizontal_spacing=0.06)
fig5.add_trace(go.Bar(y=sc.panel, x=sc.found, orientation="h", marker=dict(color=COLORS["blue"]), hovertext=lab,
                      hoverinfo="text", text=sc.found, textposition="outside", cliponaxis=False, showlegend=False), 1, 1)
fig5.add_trace(go.Bar(y=sc.panel, x=sc.false_cs, orientation="h", marker=dict(color=COLORS["orange"]), hovertext=lab,
                      hoverinfo="text", text=sc.false_cs, textposition="outside", cliponaxis=False, showlegend=False), 1, 2)
fig5.update_xaxes(title_text=f"Real causal variants caught (of {sc.causal.iloc[0]})", range=[0, sc.causal.iloc[0]],
                  fixedrange=True, row=1, col=1)
fig5.update_xaxes(title_text="False credible sets", rangemode="tozero", fixedrange=True, row=1, col=2)
fig5.update_yaxes(fixedrange=True)
fig5.update_layout(**base_layout(height=330, margin=dict(l=190, r=20, t=20, b=50), bargap=0.35))
save_figure(fig5, f"{OUT}/fig5-panel-scorecard.html", min_width=640)
