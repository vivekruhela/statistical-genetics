#!/usr/bin/env Rscript
# Post R, figure 2: how well each reference panel filled in the variants we hid.
#
# For every hidden variant (not on the pretend chip from 02_imputation.sh) and every panel, the
# squared correlation between the imputed dosage and the real genotype, separately in the
# European and the Yoruba people of the cohort. That r^2 is the usual measure of imputation
# accuracy: 1 is perfect, 0 means the panel knew nothing about this variant in these people.
#   OUT/imputation_r2.tsv.gz       one row per variant, group and panel
#   OUT/imputation_summary.tsv     average r^2 by group, panel and allele frequency
#   OUT/fig2-imputation.html
# Usage: Rscript posts/R-reference-panels/03_imputation_accuracy.R data results/post-R/imputation results/post-R
# Python twin: 03_imputation_accuracy.py (same numbers).
source("lib/statgen.R"); source("lib/figures.R")
args <- commandArgs(trailingOnly = TRUE)
D <- args[1]; IMP <- args[2]; OUT <- args[3]
PANELS <- c(EUR = "European panel", AFR = "African panel", BOTH = "Both panels together")
MAF_BINS <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
GROUPS <- c(EUR = "European", AFR = "Yoruba")

samples <- fread(file.path(D, "cohort/samples.tsv"))
loci <- fread(file.path(D, "cohort/loci.tsv"))
row_r2 <- function(a, b) {   # squared correlation between matching rows of a and b
  a <- a - rowMeans(a); b <- b - rowMeans(b)
  num <- rowSums(a * b); den <- sqrt(rowSums(a^2) * rowSums(b^2))
  ifelse(den > 0, (num / den)^2, 0)
}
res <- list()
for (lid in loci$locus) {
  geno <- read_traw(file.path(D, "cohort/genotypes/loci", paste0(lid, ".traw")))
  chip <- fread(file.path(IMP, "chip", paste0(lid, ".sites.txt")), header = FALSE)$V3
  for (panel in names(PANELS)) {
    ds <- fread(file.path(IMP, "imputed", paste0(lid, ".", panel, ".ds.tsv.gz")))
    ids <- intersect(setdiff(ds$id, chip), geno$info$id)
    dmat <- as.matrix(ds[match(ids, ds$id), -1])
    for (grp in names(GROUPS)) {
      people <- samples[superpopulation == grp, sample]
      truth <- geno$g[ids, people, drop = FALSE]
      af <- rowMeans(truth) / 2
      maf <- pmin(af, 1 - af)
      keep <- maf >= min(MAF_BINS)
      res[[length(res) + 1]] <- data.table(locus = lid, variant_id = ids[keep], group = GROUPS[[grp]],
                                           panel = PANELS[[panel]], maf = maf[keep],
                                           r2 = row_r2(dmat[keep, people, drop = FALSE], truth[keep, , drop = FALSE]))
    }
  }
}
res <- rbindlist(res)
res[, maf_bin := cut(maf, MAF_BINS, right = FALSE, include.lowest = TRUE,
                     labels = sprintf("%g-%g%%", 100 * head(MAF_BINS, -1), 100 * MAF_BINS[-1]))]
write_tsv(res[, -"maf_bin"], file.path(OUT, "imputation_r2.tsv.gz"))
summ <- res[, .(n_variants = .N, mean_r2 = mean(r2), share_r2_above_0.8 = mean(r2 >= 0.8)),
            by = .(group, panel, maf_bin)][order(group, panel, maf_bin)]
summ[, maf_bin := as.character(maf_bin)]
write_tsv(summ, file.path(OUT, "imputation_summary.tsv"))
print(summ)

# ---- figure: one panel per group of people, one line per reference panel ----
cols <- c(COLORS[["aqua"]], COLORS[["orange"]], COLORS[["blue"]])
bins <- unique(summ$maf_bin)
figs <- lapply(seq_along(GROUPS), function(k) {
  f <- plot_ly()
  for (j in seq_along(PANELS)) {
    d <- summ[group == GROUPS[[k]] & panel == PANELS[[j]]]
    f <- add_trace(f, x = d$maf_bin, y = d$mean_r2, type = "scatter", mode = "lines+markers", name = PANELS[[j]],
                   legendgroup = PANELS[[j]], showlegend = k == 1,
                   line = list(width = 2, color = cols[j], dash = if (j == 3) "dot" else "solid"), marker = list(size = 8, line = list(width = 2, color = SURFACE)),
                   text = sprintf("%s people, %s<br>variants with MAF %s: %s<br>average r² = %.2f<br>imputed well (r² ≥ 0.8): %.0f%%",
                                  GROUPS[[k]], PANELS[[j]], d$maf_bin, format(d$n_variants, big.mark = ","),
                                  d$mean_r2, 100 * d$share_r2_above_0.8), hoverinfo = "text")
  }
  f
})
fig <- subplot(figs, nrows = 1, shareY = TRUE, margin = 0.04, titleX = TRUE) |>
  base_layout(height = 430, legend = list(orientation = "h", x = 0, y = -0.24),
              margin = list(l = 60, r = 20, t = 50, b = 110),
              xaxis = list(title = "Minor allele frequency of the hidden variant", categoryorder = "array",
                           categoryarray = bins, fixedrange = TRUE),
              xaxis2 = list(title = "Minor allele frequency of the hidden variant", categoryorder = "array",
                            categoryarray = bins, fixedrange = TRUE),
              yaxis = list(title = "Imputation accuracy (average r²)", range = c(0, 1.02), fixedrange = TRUE),
              annotations = list(
                list(x = 0, y = 1.08, xref = "paper", yref = "paper", xanchor = "left", showarrow = FALSE,
                     text = "<b>Filling in European genomes</b>", font = list(size = 13)),
                list(x = 0.52, y = 1.08, xref = "paper", yref = "paper", xanchor = "left", showarrow = FALSE,
                     text = "<b>Filling in Yoruba genomes</b>", font = list(size = 13))))
save_figure(fig, file.path(OUT, "fig2-imputation.html"), "Imputation by reference panel", min_width = 700)
