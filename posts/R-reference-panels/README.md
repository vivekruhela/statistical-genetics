# Post R: reference panels

Code for the post *Reference panels: which one, when, how*. Every number and figure in the post
comes from these scripts, run on the series data (`simulate/run_all.sh`). Each analysis has an R
version (the main track) and a Python twin that gives the same numbers.

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `01_ld_maps.{R,py}` | LD between common variants, European vs Yoruba people of the cohort (same number of each) | `ld_decay.tsv`, `ld_map.tsv`, `fig1-ld-maps.html` |
| 2 | `02_imputation.sh` | hides 9 in 10 common variants, then imputes them with Beagle 5 from a European, an African and a combined panel (people not in the cohort) | `imputation/` |
| 3 | `03_imputation_accuracy.{R,py}` | imputed vs real genotypes: r² by ancestry, panel and allele frequency | `imputation_summary.tsv`, `imputation_r2.tsv.gz`, `fig2-imputation.html` |
| 4 | `04_gwas_ld.{R,py}` | fine-mapping windows (±250 kb around each GWAS hit with p < 5e-8, capstone locus left out), and the in-sample LD of the 10,000 European GWAS people, rebuilt exactly as step 9 made them | `windows.tsv`, `ld/insample/` |
| 5 | `05_finemap_panels.{R,py}` | SuSiE-RSS on the same GWAS with LD from six sources, scored against the answer key, plus the z-score check of each panel | `finemap_*.tsv`, `fig3-…`, `fig4-…`, `fig5-…` |
| | `decision-chart.html` | the interactive "which panel?" chart; plain HTML, no code to run | |

The R track uses susieR; the Python track uses `lib/susie_rss.py`, a port of `susieR::susie_rss`
that gives the same PIPs (checked in `tests/test_lib.py`). `lib/statgen.{R,py}` has the shared LD
helpers (`ld_matrix`, `read_ld`/`write_ld`, and `ld_check`, the z-score check), and
`lib/figures.{R,py}` with `lib/figure_theme.js` gives every figure the site's light and dark
themes.

## Run it

From the repository root, in the `statgen-r` environment (which has R, PLINK 2, BCFtools and
Beagle):

```bash
bash posts/R-reference-panels/run_all.sh data results/post-R               # R track
TRACK=python bash posts/R-reference-panels/run_all.sh data results/post-R  # Python track
bash posts/R-reference-panels/check.sh data results/post-R-check           # both, and compare
```

On a 4-core machine the whole post takes about 15 to 30 minutes: a few minutes of Beagle, a few
for the in-sample LD, and most of the rest for fine-mapping 6 panels × the windows. Set
`THREADS` for Beagle and `SHOW` (default `L05_trait_q`) for the locus drawn in figures 3 and 4.

Step 4 has to use the language that built the GWAS cohort (the released data are the R build),
because R and Python draw different random genomes; `run_all.sh` picks the right one by itself.
For a quick trial on the chromosome 22 test build, add `STATGEN_CONFIG=tests/config_test.json`.

## On the site

The post is `2026-11-21-finemap-r-reference-panels.md` (goes in `src/content/blog/`). The figures
go in `public/finemap/R/`, using the Python pages, which load plotly.js from a CDN and stay small:

| Page | From |
|---|---|
| `ld-maps.html` | `fig1-ld-maps.html` |
| `imputation.html` | `fig2-imputation.html` |
| `finemap-panels.html` | `fig3-finemap-panels.html` |
| `z-check.html` | `fig4-z-check.html` |
| `scorecard.html` | `fig5-panel-scorecard.html` |
| `decision-chart.html` | `decision-chart.html` |
