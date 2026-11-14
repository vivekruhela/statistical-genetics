#!/usr/bin/env Rscript
# Post R, analysis 3a: the fine-mapping windows, and the in-sample LD of the GWAS people.
#
# A window is 250 kb either side of the strongest GWAS variant at a teaching locus, for each
# trait (trait_q, disease) whose strongest variant passes p < 5e-8 in the European GWAS. The
# capstone locus (L10) is left out on purpose. For each window this rebuilds the 10,000 European
# GWAS genomes exactly as step 9 made them (same seed, same code), and writes their LD matrix:
# the "in-sample" LD that a summary-statistics method would ideally use, and almost never has.
#   OUT/windows.tsv
#   OUT/ld/insample/<locus>_<trait>.ld.{bin,ids}
# The R build must be read with this script and a Python build with 04_gwas_ld.py, because the
# two languages draw different random genomes; the script checks gwas/EUR/summary.tsv.
# Usage: Rscript posts/R-reference-panels/04_gwas_ld.R data results/post-R
source("lib/statgen.R")
args <- commandArgs(trailingOnly = TRUE)
D <- args[1]; OUT <- args[2]
HALF_WINDOW <- 250000; MIN_MLOG10P <- -log10(5e-8); MIN_AF <- 0.01

cfg <- read_config(); gc_cfg <- cfg$gwas_cohort
loci <- fread(file.path(D, "cohort/loci.tsv"))
built <- fread(file.path(D, "gwas/EUR/summary.tsv"))
if (built$built_with != "R") stop("these GWAS genomes were built with ", built$built_with, "; run 04_gwas_ld.py instead")
if (built$n != gc_cfg$n[["EUR"]]) stop("the GWAS has ", built$n, " people but the config says ", gc_cfg$n[["EUR"]],
                                       ": set STATGEN_CONFIG to the config the data were built with")

windows <- list()
for (trait in c("trait_q", "disease")) {
  ss <- fread(file.path(D, "gwas/EUR", paste0(trait, ".sumstats.tsv.gz")))
  for (lid in setdiff(loci$locus, loci[scenario == "capstone", locus])) {
    x <- ss[locus == lid & af >= MIN_AF & af <= 1 - MIN_AF]
    lead <- x[which.max(abs(z))]
    if (lead$mlog10p < MIN_MLOG10P) next
    w <- x[abs(pos - lead$pos) <= HALF_WINDOW]
    windows[[length(windows) + 1]] <- data.table(locus = lid, trait = trait, lead_id = lead$variant_id,
      lead_pos = lead$pos, lead_mlog10p = lead$mlog10p, start = min(w$pos), end = max(w$pos),
      n_variants = nrow(w), n_gwas = lead$n)
  }
}
windows <- rbindlist(windows)[order(locus, trait)]
write_tsv(windows, file.path(OUT, "windows.tsv"))
print(windows)

# the GWAS genomes at one locus, exactly as in simulate/09_gwas_cohort.R
locus_variants <- function(lid) {
  pv <- fread(cmd = paste("grep -v '^##'", shQuote(file.path(D, "cohort/genotypes/loci", paste0(lid, ".pvar")))))
  setnames(pv, 1:5, c("chr", "pos", "id", "ref", "alt"))
  pv[, 1:5]
}
build_locus <- function(pop_i, i, n) {
  set.seed(cfg$seed + 1000L * pop_i + i)
  pv <- locus_variants(loci$locus[i])
  H <- read_haplotypes(file.path(D, "work/gwas_templates", names(gc_cfg$n)[pop_i], paste0(loci$locus[i], ".hap.gz")), pv$id)
  list(pv = pv, G = mosaic_dosage(H, pv$pos, n, gc_cfg$segment_kb * 1000, gc_cfg$miscopy_rate))
}

pop_i <- match("EUR", names(gc_cfg$n))
for (k in seq_len(nrow(windows))) {
  w <- windows[k]
  prefix <- file.path(OUT, "ld/insample", paste0(w$locus, "_", w$trait))
  if (file.exists(paste0(prefix, ".ld.bin"))) { cat(w$locus, w$trait, "already done\n"); next }
  b <- build_locus(pop_i, match(w$locus, loci$locus), gc_cfg$n[["EUR"]])
  ss <- fread(file.path(D, "gwas/EUR", paste0(w$trait, ".sumstats.tsv.gz")))
  ids <- ss[locus == w$locus & af >= MIN_AF & af <= 1 - MIN_AF & pos >= w$start & pos <= w$end, variant_id]
  g <- b$G[match(ids, b$pv$id), , drop = FALSE]
  rownames(g) <- ids
  # check that these are the genomes behind the summary statistics
  if (max(abs(rowMeans(g) / 2 - ss[match(ids, variant_id), af])) > 1e-9)
    stop("rebuilt genomes do not match the summary statistics: was the data built with another config or seed?")
  write_ld(ld_matrix(g), ids, prefix)
  cat(sprintf("%s %s: in-sample LD for %d variants in %d people\n", w$locus, w$trait, length(ids), ncol(g)))
}
