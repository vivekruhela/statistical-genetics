#!/usr/bin/env Rscript
# Step 3 (R): match GEUVADIS libraries to genotyped people, filter and normalise expression.
# Usage: Rscript simulate/03_prepare_expression.R data/geuvadis data/meta data/work/expression_real
# Python twin: simulate/03_prepare_expression.py (same outputs, same numbers).
source("lib/statgen.R")
args <- commandArgs(trailingOnly = TRUE)
geuv_dir <- args[1]; meta_dir <- args[2]; out <- args[3]
cfg <- read_config()

counts <- read_geuvadis_counts(file.path(geuv_dir, "GD660.GeneQuantCount.txt.gz"))
cohort <- fread(file.path(meta_dir, "cohort_samples.tsv"))
keep <- sort(intersect(colnames(counts), cohort$sample))
cat(sprintf("GEUVADIS individuals: %d; with genotypes and metadata: %d\n", ncol(counts), length(keep)))
counts <- counts[, keep, drop = FALSE]
samples <- cohort[match(keep, sample)]
print(samples[, .N, by = .(population, sex)][order(population, sex)])

genes <- read_genes(file.path(meta_dir, "gencode.v39.annotation.gtf.gz"))
norm <- normalise_expression(counts, genes, cfg)
pcs <- sample_pcs(norm$int, cfg$expression$n_expression_pcs)
cat(sprintf("genes kept after filtering: %d of %d\n", nrow(norm$int), nrow(counts)))

dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv(samples, file.path(out, "samples.tsv"))
writeLines(samples$sample, file.path(out, "samples.txt"))
write_tsv(matrix_table(counts[rownames(counts) %in% genes$gene_id, ], "gene_id"),
          file.path(out, "counts_all.tsv.gz"))
write_tsv(genes[match(rownames(norm$int), gene_id)], file.path(out, "genes.tsv"))
write_tsv(data.table(sample = names(norm$size_factors), size_factor = norm$size_factors),
          file.path(out, "size_factors.tsv"))
write_tsv(matrix_table(norm$int, "gene_id"), file.path(out, "expression_int.tsv.gz"))
write_tsv(matrix_table(pcs, "sample"), file.path(out, "expression_pcs.tsv"))
