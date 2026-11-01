#!/usr/bin/env Rscript
# Step 9 (R): a large GWAS cohort at the teaching loci, and its summary statistics.
#
# 445 people cannot power a GWAS, and in real studies the GWAS and the eQTL data come from
# different people anyway (which is what coloc and SMR assume). So for each ancestry in
# config "gwas_cohort" we build new people as mosaics of the step 8 template haplotypes
# (real LD, none of the cohort's people), run them through the same causal graph as the
# small cohort (data/truth/edges.tsv and nodes.tsv), and test every locus variant for
# the quantitative trait and the disease, adjusting for age and sex.
#
# Two passes per locus keep memory low: pass 1 builds each locus and keeps only the causal
# variants; after the phenotypes are drawn, pass 2 rebuilds each locus from the same seed
# and computes the statistics.
#   OUT/POP/trait_q.sumstats.tsv.gz, OUT/POP/disease.sumstats.tsv.gz, OUT/POP/summary.tsv
# Usage: Rscript simulate/09_gwas_cohort.R data/cohort data/truth data/work/gwas_templates data/gwas
# Python twin: simulate/09_gwas_cohort.py (same model, NumPy random numbers).
source("lib/statgen.R")
args <- commandArgs(trailingOnly = TRUE)
cohort_dir <- args[1]; truth_dir <- args[2]; tmpl_dir <- args[3]; out <- args[4]
cfg <- read_config(); gc_cfg <- cfg$gwas_cohort
loci <- fread(file.path(cohort_dir, "loci.tsv"))
edges <- fread(file.path(truth_dir, "edges.tsv"))
nodes <- fread(file.path(truth_dir, "nodes.tsv"))
model <- fromJSON(file.path(truth_dir, "model.json"))
causal <- unique(sub("^var:", "", edges[startsWith(from, "var:")]$from))

locus_variants <- function(lid) {
  pv <- fread(cmd = paste("grep -v '^##'", shQuote(file.path(cohort_dir, "genotypes", "loci", paste0(lid, ".pvar")))))
  setnames(pv, 1:5, c("chr", "pos", "id", "ref", "alt"))
  pv[, 1:5]
}
build_locus <- function(pop_i, i, n) {
  set.seed(cfg$seed + 1000L * pop_i + i)
  pv <- locus_variants(loci$locus[i])
  H <- read_haplotypes(file.path(tmpl_dir, names(gc_cfg$n)[pop_i], paste0(loci$locus[i], ".hap.gz")), pv$id)
  if (anyNA(H)) stop(loci$locus[i], ": template haplotypes miss some cohort variants")
  list(pv = pv, G = mosaic_dosage(H, pv$pos, n, gc_cfg$segment_kb * 1000, gc_cfg$miscopy_rate))
}

for (pop_i in seq_along(gc_cfg$n)) {
  pop <- names(gc_cfg$n)[pop_i]; n <- gc_cfg$n[[pop]]
  cat(sprintf("%s: building %d people\n", pop, n))
  inputs <- list()
  for (i in seq_len(nrow(loci))) {
    b <- build_locus(pop_i, i, n)
    k <- which(b$pv$id %in% causal)
    for (j in k) inputs[[paste0("var:", b$pv$id[j])]] <- b$G[j, ]
  }
  set.seed(cfg$seed + 1000L * pop_i)
  covs <- simulate_covariates(n)
  inputs <- c(inputs, standard_covariates(covs$age, covs$male, covs$bmi))
  vals <- simulate_dag(nodes, edges, inputs, n)
  trait_q <- model$trait_q$offset + model$trait_q$scale * vals$trait_q
  disease <- as.integer(vals$disease_liability > model$disease_threshold)
  X <- cbind(1, covs$age, covs$male)

  res <- list()
  for (i in seq_len(nrow(loci))) {
    b <- build_locus(pop_i, i, n)
    af <- rowMeans(b$G) / 2
    for (s in split(seq_len(nrow(b$pv)), ceiling(seq_len(nrow(b$pv)) / 1000))) {
      G <- b$G[s, , drop = FALSE]
      for (tr in c("trait_q", "disease")) {
        st <- assoc_stats(G, if (tr == "trait_q") trait_q else disease, X, binary = tr == "disease")
        res[[length(res) + 1]] <- cbind(data.table(trait = tr, locus = loci$locus[i]),
                                        b$pv[s, .(variant_id = id, chr, pos, ref, effect_allele = alt)],
                                        data.table(af = af[s], n = n), st)
      }
    }
    cat(sprintf("  %s %s done\n", pop, loci$locus[i]))
  }
  res <- rbindlist(res)
  res <- res[af > 0 & af < 1]
  for (tr in c("trait_q", "disease")) {
    x <- res[trait == tr, -"trait"]
    if (tr == "disease") x[, `:=`(n_cases = sum(disease), n_controls = n - sum(disease))]
    write_tsv(x, file.path(out, pop, paste0(tr, ".sumstats.tsv.gz")))
  }
  summ <- data.table(population = pop, n = n, n_cases = sum(disease), trait_mean = mean(trait_q),
                     trait_sd = sd(trait_q), mean_age = mean(covs$age), male_fraction = mean(covs$male),
                     built_with = "R")
  write_tsv(summ, file.path(out, pop, "summary.tsv"))
  print(summ)
}
