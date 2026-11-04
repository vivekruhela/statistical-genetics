# Post 0: series map and setup

- `2026-11-04-finemap-00-series-map-and-setup.md`: the post, in the site's format (`draft: true`).
  It goes in the site's `src/content/blog/`.
- `series_map.R` and `series_map.py`: the clickable series map, from `posts/series.tsv`.
  The page goes in the site's `public/finemap/00/series-map.html`, which the post embeds.

```bash
Rscript posts/00-series-map/series_map.R posts/series.tsv series-map.html 2026-11-04
python  posts/00-series-map/series_map.py posts/series.tsv series-map.html 2026-11-04
```

The date decides which posts are drawn as published. Re-run both on each publish day so the
newly published node becomes clickable. Without pandoc, the R page keeps its JavaScript in a
`series-map_files/` folder next to it; the Python page loads Plotly from a CDN.

`ld_tower.py` makes the concept figure (a toy locus: one causal variant, a tower of GWAS hits,
and the credible set fine-mapping gives). It simulates its own data, so it needs no download:

```bash
python posts/00-series-map/ld_tower.py ld-tower.html
```

The page goes in the site's `public/finemap/00/ld-tower.html`.
