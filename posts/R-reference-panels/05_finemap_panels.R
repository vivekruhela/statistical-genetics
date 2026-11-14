#!/usr/bin/env Rscript
# Post R, analysis 3b: fine-mapping the same GWAS with LD from different reference panels.
#
# For each window from 04_gwas_ld.R, SuSiE-RSS (susieR) fine-maps the European GWAS z-scores
# with LD from each of these, and the answer key in data/truth scores the result:
#   insample  the GWAS people themselves (the ideal, and rarely available)
#   eur       the European people of our cohort (a matched panel, like 1000 Genomes EUR)
#   eur50     50 of them (a matched but small panel)
#   afr       the African-ancestry template people (a mismatched panel)
#   mixed     European cohort people plus the African templates (a panel that mixes ancestries)
#   flip      in-sample LD, but one variant's alleles swapped in the summary statistics
# A variant that does not vary in a panel is dropped for that panel, as it would be in practice.
# Each fit is also checked with ld_check() (lib/statgen.R): how well the panel's LD predicts each
# z-score from its neighbours. The mismatch score s is near 0 when panel and GWAS agree.
#   OUT/finemap_variants.tsv.gz   one row per variant and panel: z, PIP, credible set, z check
#   OUT/finemap_sets.tsv          one row per credible set: size, purity, does it hold a causal variant
#   OUT/finemap_summary.tsv       one row per window and panel
#   OUT/fig3-finemap-panels.html, OUT/fig4-z-check.html, OUT/fig5-panel-scorecard.html
# Usage: Rscript posts/R-reference-panels/05_finemap_panels.R data results/post-R [L05_trait_q]
# Python twin: 05_finemap_panels.py (same numbers, with lib/susie_rss.py in place of susieR).
suppressPackageStartupMessages(library(susieR))
source("lib/statgen.R"); source("lib/figures.R")
args <- commandArgs(trailingOnly = TRUE)
D <- args[1]; OUT <- args[2]; SHOW <- if (length(args) > 2) args[3] else NA
MIN_AF <- 0.01; FLAG_SD <- 4      # a z-score more than 4 SD from its prediction is flagged
PANELS <- c(insample = "GWAS people (in-sample)", eur = "European panel", eur50 = "Small European panel",
            afr = "African panel", mixed = "Mixed panel", flip = "In-sample, one allele flipped")

windows <- fread(file.path(OUT, "windows.tsv"))
samples <- fread(file.path(D, "cohort/samples.tsv"))
loci <- fread(file.path(D, "cohort/loci.tsv"))
truth <- fread(file.path(D, "truth/truth.tsv"))
eur <- samples[superpopulation == "EUR", sample]
eur50 <- eur[floor((0:49) * length(eur) / 50) + 1]

panel_genotypes <- function(lid, ids) {
  geno <- read_traw(file.path(D, "cohort/genotypes/loci", paste0(lid, ".traw")))
  g <- geno$g[match(ids, geno$info$id), , drop = FALSE]
  H <- read_haplotypes(file.path(D, "work/gwas_templates/AFR", paste0(lid, ".hap.gz")), ids)
  afr <- H[, seq(1, ncol(H), 2), drop = FALSE] + H[, seq(2, ncol(H), 2), drop = FALSE]
  afr[is.na(afr)] <- 0L   # a variant missing from the panel tells us nothing: treat it as not varying
  rownames(g) <- rownames(afr) <- ids
  list(eur = g[, eur], eur50 = g[, eur50], afr = afr, mixed = cbind(g[, eur], afr))
}
fit_status <- function(fit) {
  if (inherits(fit, "error")) {
    if (grepl("prior variance is unreasonably large", conditionMessage(fit))) "stopped: z and LD disagree"
    else "stopped: error"
  } else if (isTRUE(fit$converged)) "converged" else "did not converge"
}

vars <- list(); sets <- list(); summ <- list()
for (k in seq_len(nrow(windows))) {
  w <- windows[k]
  ss <- fread(file.path(D, "gwas/EUR", paste0(w$trait, ".sumstats.tsv.gz")))
  ss <- ss[locus == w$locus & af >= MIN_AF & af <= 1 - MIN_AF & pos >= w$start & pos <= w$end]
  R_in <- read_ld(file.path(OUT, "ld/insample", paste0(w$locus, "_", w$trait)))
  stopifnot(identical(rownames(R_in), ss$variant_id))
  ep <- if (w$trait == "disease") "disease_liability" else "trait_q"
  causal <- unique(truth[locus == w$locus & endpoint == ep, variant_id])
  causal <- intersect(causal, ss$variant_id)
  pg <- panel_genotypes(w$locus, ss$variant_id)
  flip_at <- order(-abs(ss$z))[10]   # swap the alleles of the variant with the 10th strongest signal
  for (pn in names(PANELS)) {
    R <- switch(pn, insample = , flip = R_in, ld_matrix(pg[[pn]]))
    z <- ss$z
    if (pn == "flip") z[flip_at] <- -z[flip_at]
    n_panel <- switch(pn, insample = , flip = w$n_gwas, ncol(pg[[pn]]))
    ok <- !is.na(diag(R))
    R <- R[ok, ok]; z <- z[ok]; ids <- ss$variant_id[ok]
    fit <- tryCatch(suppressWarnings(susie_rss(z, R, n = w$n_gwas)), error = function(e) e)
    status <- fit_status(fit)
    chk <- ld_check(z, R, w$n_gwas)
    pip <- if (inherits(fit, "error")) rep(NA_real_, length(z)) else fit$pip
    cs_of <- rep(NA_integer_, length(z))
    cs <- if (inherits(fit, "error") || is.null(fit$sets$cs)) list() else fit$sets$cs
    purity <- if (length(cs)) fit$sets$purity$min.abs.corr else numeric(0)
    coverage <- if (length(cs)) fit$sets$coverage else numeric(0)
    # number the credible sets from left to right, by the position of their top variant
    lead <- vapply(cs, function(x) x[which.max(pip[x])], integer(1))
    o <- order(lead); cs <- cs[o]; purity <- purity[o]; coverage <- coverage[o]
    for (j in seq_along(cs)) {
      cs_of[cs[[j]]] <- j
      sets[[length(sets) + 1]] <- data.table(locus = w$locus, trait = w$trait, panel = PANELS[[pn]], cs = j,
        size = length(cs[[j]]), purity = purity[j], coverage = coverage[j],
        holds_causal = any(ids[cs[[j]]] %in% causal), lead_variant = ids[cs[[j]]][which.max(pip[cs[[j]]])])
    }
    vars[[length(vars) + 1]] <- data.table(locus = w$locus, trait = w$trait, panel = PANELS[[pn]], variant_id = ids,
      pos = ss$pos[ok], z = z, pip = pip, cs = cs_of, causal = ids %in% causal, flipped = pn == "flip" & ids == ss$variant_id[flip_at],
      expected_z = chk$expected, z_std_diff = chk$z_std_diff)
    true_cs <- vapply(cs, function(x) any(ids[x] %in% causal), logical(1))
    found <- causal %in% ids[unlist(cs)]
    summ[[length(summ) + 1]] <- data.table(locus = w$locus, trait = w$trait, panel = PANELS[[pn]], panel_people = n_panel,
      n_variants = sum(ok), n_dropped = sum(!ok), status = status, n_cs = length(cs), n_cs_true = sum(true_cs),
      n_cs_false = sum(!true_cs), n_causal = length(causal), n_causal_found = sum(found),
      causal_pip_sum = if (inherits(fit, "error")) NA_real_ else sum(pip[ids %in% causal]),
      mismatch_s = chk$s[1], n_flagged = sum(abs(chk$z_std_diff) > FLAG_SD), max_abs_z_std_diff = max(abs(chk$z_std_diff)))
  }
  cat(sprintf("%s %s done\n", w$locus, w$trait))
}
vars <- rbindlist(vars); sets <- rbindlist(sets); summ <- rbindlist(summ)
write_tsv(vars, file.path(OUT, "finemap_variants.tsv.gz"))
write_tsv(sets, file.path(OUT, "finemap_sets.tsv"))
write_tsv(summ, file.path(OUT, "finemap_summary.tsv"))
print(summ[, .(locus, panel, n_dropped, status, n_cs, n_cs_true, n_cs_false, n_causal, n_causal_found,
               s = round(mismatch_s, 3), n_flagged)], nrows = 200)

# ---- figure 3: one locus, one row per panel ----------------------------------------
if (is.na(SHOW)) SHOW <- paste0(windows$locus[1], "_", windows$trait[1])
show <- windows[paste0(locus, "_", trait) == SHOW]
fig_panels <- names(PANELS)[1:5]
cols <- c(COLORS[["blue"]], COLORS[["aqua"]], COLORS[["yellow"]], COLORS[["orange"]], COLORS[["magenta"]])
v <- vars[locus == show$locus & trait == show$trait]
gw <- v[panel == PANELS[["insample"]]]
rows <- list(plot_ly(gw, x = ~pos / 1e6, y = ~-(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10), type = "scatter",
                     mode = "markers", marker = list(size = 6, color = LINE), hoverinfo = "text",
                     text = ~sprintf("%s<br>z = %.1f", variant_id, z), showlegend = FALSE) |>
               add_trace(data = gw[causal == TRUE], x = ~pos / 1e6, y = ~-(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10),
                         marker = list(size = 13, color = "rgba(0,0,0,0)", line = list(width = 2.5, color = COLORS[["orange"]])),
                         hoverinfo = "skip", showlegend = FALSE) |>
               layout(yaxis = list(title = "GWAS<br>-log10 p", fixedrange = TRUE)))
for (j in seq_along(fig_panels)) {
  d <- v[panel == PANELS[[fig_panels[j]]]]
  sm <- summ[locus == show$locus & trait == show$trait & panel == PANELS[[fig_panels[j]]]]
  d[, col := ifelse(is.na(cs), LINE, cols[j])]
  rows[[j + 1]] <- plot_ly(d, x = ~pos / 1e6, y = ~pip, type = "scatter", mode = "markers",
                           marker = list(size = 7, color = ~I(col)), hoverinfo = "text", showlegend = FALSE,
                           text = ~sprintf("%s<br>PIP = %.2f%s%s", variant_id, pip,
                                           ifelse(is.na(cs), "", sprintf("<br>credible set %d", cs)),
                                           ifelse(causal, "<br><b>a real cause</b>", ""))) |>
    add_trace(data = d[causal == TRUE], x = ~pos / 1e6, y = ~pip,
              marker = list(size = 13, color = "rgba(0,0,0,0)", line = list(width = 2.5, color = COLORS[["orange"]])),
              hoverinfo = "skip") |>
    layout(yaxis = list(title = "PIP", range = c(-0.05, 1.1), fixedrange = TRUE),
           annotations = list(list(x = 0, y = 1, xref = "paper", yref = "paper", xanchor = "left", yanchor = "bottom",
             showarrow = FALSE, font = list(size = 12),
             text = sprintf("<b>%s</b> (%s people): %d credible set%s, %d false", PANELS[[fig_panels[j]]],
                            format(sm$panel_people, big.mark = ","), sm$n_cs, ifelse(sm$n_cs == 1, "", "s"), sm$n_cs_false))))
}
fig3 <- subplot(rows, nrows = length(rows), shareX = TRUE, titleY = TRUE, heights = c(0.2, rep(0.16, 5)), margin = c(0, 0, 0.022, 0.022)) |>
  base_layout(height = 860, margin = list(l = 70, r = 20, t = 30, b = 50),
              xaxis = list(title = sprintf("Position on %s (Mb)", loci[locus == show$locus, chr]), fixedrange = TRUE))
save_figure(fig3, file.path(OUT, "fig3-finemap-panels.html"), "Fine-mapping with different LD panels", min_width = 640)

# ---- figure 4: does the panel's LD fit the z-scores? -------------------------------------
chk_panels <- c("insample", "afr", "flip")
figs <- lapply(chk_panels, function(pn) {
  d <- v[panel == PANELS[[pn]]]
  d[, flag := abs(z_std_diff) > FLAG_SD]
  lim <- max(abs(c(d$z, d$expected_z)))
  plot_ly(d, x = ~expected_z, y = ~z, type = "scatter", mode = "markers", showlegend = FALSE,
          marker = list(size = 6, color = ~I(ifelse(flag, COLORS[["orange"]], LINE))), hoverinfo = "text",
          text = ~sprintf("%s<br>observed z = %.1f<br>predicted from neighbours = %.1f<br>%.1f SD apart%s",
                          variant_id, z, expected_z, abs(z_std_diff), ifelse(flipped, "<br><b>the flipped variant</b>", ""))) |>
    add_trace(x = c(-lim, lim), y = c(-lim, lim), mode = "lines", line = list(width = 1, color = INK2, dash = "dot"),
              marker = NULL, hoverinfo = "skip", inherit = FALSE, type = "scatter", showlegend = FALSE) |>
    layout(xaxis = list(title = "z predicted from neighbours", fixedrange = TRUE),
           yaxis = list(title = "Observed z", fixedrange = TRUE),
           annotations = list(list(x = 0.5, y = 1.02, xref = "paper", yref = "paper", yanchor = "bottom", showarrow = FALSE,
             font = list(size = 12), text = sprintf("<b>%s</b><br>s = %.3f, %d flagged", PANELS[[pn]],
             summ[locus == show$locus & trait == show$trait & panel == PANELS[[pn]], mismatch_s],
             sum(d$flag)))))
})
fig4 <- subplot(figs, nrows = 1, titleX = TRUE, titleY = TRUE, margin = 0.035) |>
  base_layout(height = 400, margin = list(l = 60, r = 20, t = 70, b = 50))
save_figure(fig4, file.path(OUT, "fig4-z-check.html"), "Does the LD fit the z-scores?", min_width = 760)

# ---- figure 5: the scorecard over all windows --------------------------------------
sc <- summ[, .(found = sum(n_causal_found), causal = sum(n_causal), false_cs = sum(n_cs_false),
               stopped = sum(startsWith(status, "stopped"))), by = panel]
sc[, panel := factor(panel, levels = rev(PANELS))]
lab <- sprintf("%s<br>real causes caught: %d of %d<br>false credible sets: %d%s", sc$panel, sc$found, sc$causal,
               sc$false_cs, ifelse(sc$stopped > 0, sprintf("<br>SuSiE refused %d window%s", sc$stopped, ifelse(sc$stopped == 1, "", "s")), ""))
f1 <- plot_ly(sc, y = ~panel, x = ~found, type = "bar", orientation = "h", marker = list(color = COLORS[["blue"]]),
              hovertext = lab, hoverinfo = "text", text = ~found, textposition = "outside", cliponaxis = FALSE,
              showlegend = FALSE) |>
  layout(xaxis = list(title = sprintf("Real causal variants caught (of %d)", sc$causal[1]), range = c(0, sc$causal[1]), fixedrange = TRUE))
f2 <- plot_ly(sc, y = ~panel, x = ~false_cs, type = "bar", orientation = "h", marker = list(color = COLORS[["orange"]]),
              hovertext = lab, hoverinfo = "text", text = ~false_cs, textposition = "outside", cliponaxis = FALSE,
              showlegend = FALSE) |>
  layout(xaxis = list(title = "False credible sets", rangemode = "tozero", fixedrange = TRUE))
fig5 <- subplot(f1, f2, nrows = 1, shareY = TRUE, titleX = TRUE, margin = 0.03) |>
  base_layout(height = 330, margin = list(l = 190, r = 20, t = 20, b = 50), bargap = 0.35,
              yaxis = list(fixedrange = TRUE))
save_figure(fig5, file.path(OUT, "fig5-panel-scorecard.html"), "Reference panel scorecard", min_width = 640)
