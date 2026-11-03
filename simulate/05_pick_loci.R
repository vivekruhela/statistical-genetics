#!/usr/bin/env Rscript
# Step 5 (R): genome-wide cis-eQTL scan on the real expression, then pick the teaching loci.
# Each locus is a window around a gene with a strong real eQTL, so post 3.3 has a real showpiece,
# and each window later gets one simulated scenario (config "scenarios", in order).
# Usage: Rscript simulate/05_pick_loci.R data/work/expression_real data/work/plink data/work/loci
# Python twin: simulate/05_pick_loci.py (same outputs, same numbers).
source("lib/statgen.R")
args <- commandArgs(trailingOnly = TRUE)
expr_dir <- args[1]; plink_dir <- args[2]; out <- args[3]
cfg <- read_config()

samples <- fread(file.path(expr_dir, "samples.tsv"))
genes <- fread(file.path(expr_dir, "genes.tsv"))
expr <- fread(file.path(expr_dir, "expression_int.tsv.gz"))
emat <- as.matrix(expr[, -1]); rownames(emat) <- expr$gene_id
emat <- emat[, samples$sample]
epcs <- fread(file.path(expr_dir, "expression_pcs.tsv"))
gpcs <- fread(file.path(plink_dir, "cohort_pca.eigenvec"))
setnames(gpcs, 1, "sample")
no_geno <- setdiff(samples$sample, gpcs$sample)
if (length(no_geno)) stop("no genotypes (so no genotype PCs) for: ", paste(no_geno, collapse = ", "),
                          "\nRerun from step 0 so step 3 leaves them out.")
covar <- cbind(
  sex = as.numeric(samples$sex == "male"),
  as.matrix(gpcs[match(samples$sample, sample), -1])[, seq_len(cfg$genotype$n_pcs), drop = FALSE],
  as.matrix(epcs[match(samples$sample, sample), -1])[, seq_len(cfg$expression$n_expression_pcs), drop = FALSE]
)
eres <- residualise(emat, covar)
rownames(eres) <- rownames(emat)

hw <- cfg$loci$half_window_bp
scan <- list()
for (f in sort(list.files(file.path(plink_dir, "scan"), pattern = "\\.traw$", full.names = TRUE))) {
  tr <- read_traw(f)
  stopifnot(all(samples$sample %in% colnames(tr$g)))
  tr$g <- tr$g[, samples$sample, drop = FALSE]
  chr_name <- tr$info$chr[1]
  gi <- genes[chr == chr_name]
  if (!nrow(gi)) next
  gi <- gi[order(gene_id)]
  gres <- residualise(tr$g, covar)
  rownames(gres) <- rownames(tr$g)
  res <- cis_scan(eres[gi$gene_id, , drop = FALSE], gi$tss, gres, tr$info$pos, hw, ncol(covar))
  cat(sprintf("%s: %d genes scanned\n", chr_name, nrow(res)))
  scan[[chr_name]] <- res
}
scan <- rbindlist(scan)
scan <- merge(genes[, .(gene_id, gene_name, chr, tss)], scan, by = "gene_id")
setorder(scan, -lead_mlog10p, gene_id)
write_tsv(scan, file.path(out, "cis_scan_lead.tsv"))

# ---- choose the loci: strongest eQTL genes, spread out, away from awkward regions ----
lc <- cfg$loci
cand <- scan[lead_mlog10p >= lc$min_lead_mlog10p & n_cis >= lc$min_variants & n_cis <= lc$max_variants]
cand[, `:=`(start = pmax(1, tss - hw), end = tss + hw)]
excl <- rbindlist(lapply(lc$exclude_regions, as.data.table))
bad <- vapply(seq_len(nrow(cand)), function(i)
  any(excl$chr == cand$chr[i] & excl$start <= cand$end[i] & excl$end >= cand$start[i]), logical(1))
cand <- cand[!bad]
chosen <- cand[0]
for (i in seq_len(nrow(cand))) {
  if (nrow(chosen) == lc$n_loci) break
  x <- cand[i]
  same <- chosen[chr == x$chr]
  if (nrow(same) >= lc$max_loci_per_chr) next
  if (any(same$start - lc$min_gap_bp <= x$end & same$end + lc$min_gap_bp >= x$start)) next
  chosen <- rbind(chosen, x)
}
if (nrow(chosen) < lc$n_loci) stop(sprintf("only %d loci passed the filters", nrow(chosen)))
scen <- vapply(cfg$scenarios, `[[`, "", "name")
loci <- chosen[, .(locus = sprintf("L%02d", seq_len(.N)), scenario = scen[seq_len(.N)],
                   chr, start, end, anchor_gene_id = gene_id, anchor_gene = gene_name, anchor_tss = tss,
                   anchor_lead_id = lead_id, anchor_lead_mlog10p = lead_mlog10p,
                   anchor_second_mlog10p = second_mlog10p, n_scan_variants = n_cis)]
write_tsv(loci, file.path(out, "loci.tsv"))
print(loci[, .(locus, scenario, chr, start, anchor_gene, round(anchor_lead_mlog10p, 1))])
