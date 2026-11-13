#!/usr/bin/env Rscript
# Post R, figure 1: the same stretch of DNA has a different LD map in different ancestries.
#
# Left: r^2 between every pair of common variants in a 150 kb region, European people below the
# diagonal and Yoruba people above it. Right: how fast LD fades with distance, averaged over all
# ten teaching loci. Both groups have the same number of people (every k-th European is used),
# because r^2 from a small sample is biased upwards by about 1/n.
#   OUT/ld_decay.tsv, OUT/ld_map.tsv, OUT/fig1-ld-maps.html
# Usage: Rscript posts/R-reference-panels/01_ld_maps.R data results/post-R [L02]
# Python twin: 01_ld_maps.py (same numbers).
source("lib/statgen.R"); source("lib/figures.R")
args <- commandArgs(trailingOnly = TRUE)
D <- args[1]; OUT <- args[2]; SHOW <- if (length(args) > 2) args[3] else "L02"
BINS <- c(0, 2, 5, 10, 20, 50, 100, 200) * 1000   # distance bins for the decay curve
MAP_BP <- 150000; MAP_MAX <- 150                 # heatmap region and most variants shown

samples <- fread(file.path(D, "cohort/samples.tsv"))
loci <- fread(file.path(D, "cohort/loci.tsv"))
yri <- samples[superpopulation == "AFR", sample]
eur_all <- samples[superpopulation == "EUR", sample]
eur <- eur_all[floor((seq_along(yri) - 1) * length(eur_all) / length(yri)) + 1]
groups <- list(European = eur, Yoruba = yri)

maf <- function(g) { af <- rowMeans(g) / 2; pmin(af, 1 - af) }
decay <- list(); map <- NULL
for (lid in loci$locus) {
  geno <- read_traw(file.path(D, "cohort/genotypes/loci", paste0(lid, ".traw")))
  for (grp in names(groups)) {
    g <- geno$g[, groups[[grp]]]
    keep <- maf(g) >= 0.05
    g <- g[keep, ]; pos <- geno$info$pos[keep]
    r2 <- ld_matrix(g)^2
    dist <- abs(outer(pos, pos, "-"))
    up <- upper.tri(r2) & dist < max(BINS)
    bin <- findInterval(dist[up], BINS)
    decay[[length(decay) + 1]] <- data.table(locus = lid, group = grp, bin = bin, r2 = r2[up])[
      , .(sum_r2 = sum(r2), n_pairs = .N), by = .(locus, group, bin)]
  }
  if (lid == SHOW) {
    centre <- geno$info$pos[geno$info$id == loci[locus == lid, anchor_lead_id]]
    ok <- maf(geno$g[, eur]) >= 0.05 & maf(geno$g[, yri]) >= 0.05 &
      abs(geno$info$pos - centre) <= MAP_BP / 2
    idx <- which(ok)
    idx <- idx[unique(floor(seq(0, length(idx) - 1, length.out = min(MAP_MAX, length(idx)))) + 1)]
    r_eur <- ld_matrix(geno$g[idx, eur])^2
    r_yri <- ld_matrix(geno$g[idx, yri])^2
    m <- r_eur; m[upper.tri(m)] <- r_yri[upper.tri(r_yri)]
    map <- data.table(locus = lid, i = rep(seq_along(idx), times = length(idx)) , j = rep(seq_along(idx), each = length(idx)),
                      id_i = rep(geno$info$id[idx], times = length(idx)), id_j = rep(geno$info$id[idx], each = length(idx)),
                      r2 = as.vector(m))
  }
}
decay <- rbindlist(decay)[, .(mean_r2 = sum(sum_r2) / sum(n_pairs), n_pairs = sum(n_pairs)), by = .(group, bin)][order(group, bin)]
decay[, `:=`(from_kb = BINS[bin] / 1000, to_kb = BINS[bin + 1] / 1000)]
decay <- decay[, .(group, from_kb, to_kb, mean_r2, n_pairs)]
write_tsv(decay, file.path(OUT, "ld_decay.tsv"))
write_tsv(map, file.path(OUT, "ld_map.tsv"))
print(decay)

# ---- figure ----
k <- max(map$i)
z <- matrix(map$r2, k, k)
ids <- map$id_i[seq_len(k)]
kb <- (as.integer(sub("^[^:]*:([0-9]+):.*", "\\1", ids)) - min(as.integer(sub("^[^:]*:([0-9]+):.*", "\\1", ids)))) / 1000
hover <- matrix(sprintf("%s<br>%s<br>%s: r² = %.2f", rep(ids, times = k), rep(ids, each = k),
                        ifelse(outer(seq_len(k), seq_len(k), "<"), "Yoruba", "European"), as.vector(z)), k, k)
hm <- plot_ly(x = seq_len(k), y = seq_len(k), z = z, type = "heatmap", zmin = 0, zmax = 1, text = hover,
              hoverinfo = "text", colorscale = SEQ, meta = "seq",
              colorbar = list(title = list(text = "r²", side = "top"), len = 0.6, thickness = 10, x = 0.425))
mids <- (BINS[-1] + BINS[-length(BINS)]) / 2000
dc <- plot_ly()
for (grp in names(groups)) {
  d <- decay[group == grp]
  dc <- add_trace(dc, x = (d$from_kb + d$to_kb) / 2, y = d$mean_r2, type = "scatter", mode = "lines+markers",
                  name = sprintf("%s (%d people)", grp, length(groups[[grp]])),
                  line = list(width = 2, color = if (grp == "European") COLORS[["aqua"]] else COLORS[["orange"]]),
                  marker = list(size = 8, line = list(width = 2, color = SURFACE)),
                  text = sprintf("%s, %g to %g kb apart<br>mean r² = %.3f (%s pairs)", grp, d$from_kb, d$to_kb,
                                 d$mean_r2, format(d$n_pairs, big.mark = ",")), hoverinfo = "text")
}
fig <- subplot(hm, dc, widths = c(0.5, 0.5), margin = 0.085, titleX = TRUE, titleY = TRUE) |>
  base_layout(height = 440, legend = list(x = 0.98, xanchor = "right", y = 0.98),
              xaxis = list(title = sprintf("Below the diagonal: European \u00b7 above: Yoruba<br>(%d common variants across %d kb at %s)", k, MAP_BP / 1000, SHOW),
                           showticklabels = FALSE, showgrid = FALSE, fixedrange = TRUE, scaleanchor = "y", constrain = "domain"),
              yaxis = list(showticklabels = FALSE, showgrid = FALSE, fixedrange = TRUE, autorange = "reversed", constrain = "domain"),
              xaxis2 = list(title = "Distance between two variants (kb, log scale)", type = "log", fixedrange = TRUE),
              yaxis2 = list(title = "Average r²", rangemode = "tozero", fixedrange = TRUE),
              annotations = list(
                list(x = 0.02, y = 1.06, xref = "paper", yref = "paper", showarrow = FALSE, xanchor = "left",
                     text = "<b>One region, two LD maps</b>", font = list(size = 13)),
                list(x = 0.59, y = 1.06, xref = "paper", yref = "paper", showarrow = FALSE, xanchor = "left",
                     text = "<b>LD fades faster in Yoruba genomes</b>", font = list(size = 13))))
save_figure(fig, file.path(OUT, "fig1-ld-maps.html"), "LD maps by ancestry", min_width = 760)
