#!/usr/bin/env Rscript
# Step 7 (R): plant the simulated layers on the real cohort and write the truth.
#
# For each teaching locus this picks causal variants from the real genotypes, adds a teaching
# gene (SIMGxx) next to the real anchor gene, and wires variants, CpGs, expression, a
# quantitative trait and a disease into one causal graph following the locus scenario. Then it:
#   - simulates age, BMI, the trait and the disease (sex is real)
#   - simulates CpG M-values (promoter clusters at each teaching gene plus null CpGs)
#   - adds teaching-gene counts to the real GEUVADIS counts
#   - spikes a disease effect into some real genes, for the differential expression post
#   - renormalises expression and writes everything under data/cohort and data/truth
#
# The released dataset is the output of this R script. simulate/07_simulate_layers.py builds
# the same model with Python's random numbers (same design, different draws).
#
# Usage: Rscript simulate/07_simulate_layers.R data/work data/cohort data/truth data/meta
source("lib/statgen.R")
args <- commandArgs(trailingOnly = TRUE)
work <- args[1]; cohort_dir <- args[2]; truth_dir <- args[3]; meta_dir <- args[4]
cfg <- read_config()
ef <- cfg$effects; cv <- cfg$causal_variants; me <- cfg$methylation
set.seed(cfg$seed)

samples <- fread(file.path(work, "expression_real", "samples.tsv"))
n <- nrow(samples)
is_eur <- samples$region == "EUR"
loci <- fread(file.path(work, "loci", "loci.tsv"))

# ---- covariates (sex is real) -------------------------------------------------------
covs <- simulate_covariates(n, male = as.integer(samples$sex == "male"))
inputs <- standard_covariates(covs$age, covs$male, covs$bmi)

# ---- causal variants -------------------------------------------------------------------
maf <- function(p) pmin(p, 1 - p)
pick1 <- function(idx, what, locus) {
  if (!length(idx)) stop(sprintf("%s: no variant fits the rule for %s", locus, what))
  idx[sample.int(length(idx), 1)]
}
edges <- list(); nodes <- list(); genes_sim <- list(); cpgs <- list(); variants <- list()
add_edge <- function(locus, scenario, from, to, r2, sign, center = 0, var_ref = 1, role = "") {
  edges[[length(edges) + 1]] <<- data.table(
    locus = locus, scenario = scenario, from = from, to = to, role = role, r2 = r2, sign = sign,
    center = center, var_ref = var_ref, weight = sign * sqrt(r2 / var_ref))
}
rsign <- function() sample(c(-1, 1), 1)

for (i in seq_len(nrow(loci))) {
  L <- loci[i]; lid <- L$locus; sc <- L$scenario
  tr <- read_traw(file.path(cohort_dir, "genotypes", "loci", paste0(lid, ".traw")))
  g <- tr$g[, samples$sample, drop = FALSE]
  info <- tr$info
  info[, `:=`(af = rowMeans(g) / 2, af_eur = rowMeans(g[, is_eur]) / 2, af_yri = rowMeans(g[, !is_eur]) / 2)]
  tss <- L$anchor_tss + ef$teaching_gene_offset_bp
  if (tss > L$end) tss <- L$anchor_tss - ef$teaching_gene_offset_bp
  snv <- nchar(info$ref) == 1 & nchar(info$alt) == 1
  ok <- snv & maf(info$af_eur) >= cv$min_maf_eur & maf(info$af_yri) >= cv$min_maf_yri
  near <- abs(info$pos - tss) <= ef$cis_window_bp
  r2_eur <- function(a) suppressWarnings(as.vector(cor(g[a, is_eur], t(g[, is_eur]))^2))

  A <- if (sc == "ancestry") {
    pick1(which(snv & near & maf(info$af_eur) >= cv$ancestry_maf_eur[1] &
                  maf(info$af_eur) <= cv$ancestry_maf_eur[2] &
                  maf(info$af_yri) >= cv$ancestry_maf_yri_min), "A (ancestry)", lid)
  } else if (sc == "gwas_only") {
    pick1(which(ok), "A", lid)
  } else {
    pick1(which(ok & near), "A", lid)
  }
  rA <- r2_eur(A); rA[is.na(rA)] <- 1
  B <- C <- NA_integer_
  if (sc %in% c("two_signals", "capstone")) B <- pick1(which(ok & near & rA < cv$low_ld_r2_max), "B", lid)
  if (sc == "gwas_two") B <- pick1(which(ok & rA < cv$low_ld_r2_max), "B", lid)
  if (sc == "linkage") B <- pick1(which(ok & abs(info$pos - info$pos[A]) > 1000 &
                                          rA >= cv$linkage_r2[1] & rA <= cv$linkage_r2[2]), "B (linkage)", lid)
  if (sc == "capstone") {
    rB <- r2_eur(B); rB[is.na(rB)] <- 1
    C <- pick1(which(ok & rA < cv$low_ld_r2_max & rB < cv$low_ld_r2_max), "C", lid)
  }
  roles <- c(A = A, B = B, C = C); roles <- roles[!is.na(roles)]
  vnode <- function(k) paste0("var:", info$id[k])
  vcenter <- function(k) 2 * info$af[k]
  vvar <- function(k) 2 * info$af[k] * (1 - info$af[k])
  for (r in names(roles)) {
    k <- roles[[r]]
    inputs[[vnode(k)]] <- g[k, ]
    variants[[length(variants) + 1]] <- cbind(data.table(locus = lid, scenario = sc, role = r),
                                              info[k, .(id, chr, pos, ref, alt, af, af_eur, af_yri)])
  }

  # teaching gene, its promoter CpG cluster and the cluster's shared factor
  gene <- sprintf("SIMG%02d", i); expr <- paste0("expr:", gene); fac <- paste0("factor:", lid)
  genes_sim[[i]] <- data.table(gene_id = gene, gene_name = gene, gene_type = "teaching",
                               chr = L$chr, start = tss, end = tss + 20000, strand = "+", tss = tss)
  cpg_ids <- sprintf("simcg_%s_%02d", lid, seq_len(me$cpgs_per_locus))
  cpgs[[length(cpgs) + 1]] <- data.table(cpg_id = cpg_ids, chr = L$chr,
                                         pos = sort(tss + sample(-2000:2000, me$cpgs_per_locus)),
                                         locus = lid, cluster = gene)
  nodes[[length(nodes) + 1]] <- data.table(node = fac, type = "factor", order = 1L)
  nodes[[length(nodes) + 1]] <- data.table(node = expr, type = "expression", order = 3L)
  for (cg in cpg_ids) add_edge(lid, sc, fac, paste0("cpg:", cg), me$regional_factor_r2, 1, role = "regional")
  cpg1 <- paste0("cpg:", cpg_ids[1])

  # scenario wiring; r2 targets are per causal variant (see config "_r2_note")
  eq <- function(k, r2, role) add_edge(lid, sc, vnode(k), expr, r2, rsign(), vcenter(k), vvar(k), role)
  direct <- function(k, r2, to, role) add_edge(lid, sc, vnode(k), to, r2, rsign(), vcenter(k), vvar(k), role)
  via_expr <- function(target_r2, chain_r2, to) add_edge(lid, sc, expr, to, target_r2 / chain_r2, rsign(), role = "expr_to_trait")
  mediated_chain <- function(k) {
    add_edge(lid, sc, vnode(k), cpg1, ef$mqtl_r2, rsign(), vcenter(k), vvar(k), "mqtl")
    add_edge(lid, sc, cpg1, expr, ef$cpg_to_expr_r2, -1, role = "cpg_to_expr")
    ef$mqtl_r2 * ef$cpg_to_expr_r2
  }
  switch(sc,
    mediated    = { ch <- mediated_chain(A); via_expr(ef$trait_r2, ch, "trait_q") },
    two_signals = { eq(A, ef$eqtl_r2, "eqtl"); eq(B, ef$eqtl_second_r2, "eqtl"); via_expr(ef$trait_r2, ef$eqtl_r2, "trait_q") },
    linkage     = { eq(A, ef$eqtl_r2, "eqtl"); direct(B, ef$trait_r2, "trait_q", "direct") },
    pleiotropy  = { eq(A, ef$eqtl_r2, "eqtl"); direct(A, ef$trait_r2, "trait_q", "direct") },
    gwas_only   = { direct(A, ef$trait_strong_r2, "trait_q", "direct") },
    eqtl_only   = { eq(A, ef$eqtl_r2, "eqtl") },
    gwas_two    = { eq(A, ef$eqtl_r2, "eqtl"); via_expr(ef$trait_r2, ef$eqtl_r2, "trait_q"); direct(B, ef$trait_r2, "trait_q", "direct") },
    binary      = { eq(A, ef$eqtl_r2, "eqtl"); via_expr(ef$disease_r2, ef$eqtl_r2, "disease_liability") },
    ancestry    = { eq(A, ef$eqtl_r2, "eqtl"); via_expr(ef$trait_r2, ef$eqtl_r2, "trait_q") },
    capstone    = { ch <- mediated_chain(A); eq(B, ef$eqtl_second_r2, "eqtl"); via_expr(ef$trait_r2, ch, "trait_q"); direct(C, ef$trait_r2, "trait_q", "direct") },
    stop("unknown scenario ", sc)
  )
}

# null CpGs, spread over the locus windows, and age effects on a share of all CpGs
nnull <- me$null_cpgs
li <- sample(seq_len(nrow(loci)), nnull, replace = TRUE)
cpgs[[length(cpgs) + 1]] <- data.table(
  cpg_id = sprintf("simcg_null_%04d", seq_len(nnull)), chr = loci$chr[li],
  pos = round(loci$start[li] + runif(nnull) * (loci$end[li] - loci$start[li])),
  locus = loci$locus[li], cluster = "none")
cpgs <- rbindlist(cpgs)
cpgs[, `:=`(m_mean = ifelse(cluster == "none", runif(.N, -3, 3), runif(.N, -4, -1)),
            m_scale = runif(.N, 0.3, 0.9))]
setorder(cpgs, chr, pos)
age_cpg <- cpgs$cpg_id[runif(nrow(cpgs)) < me$age_cpg_fraction]
for (cg in age_cpg) add_edge(cpgs[cpg_id == cg]$locus, "methylation_age", "cov:age", paste0("cpg:", cg), me$age_r2, rsign(), role = "age")
nodes[[length(nodes) + 1]] <- data.table(node = paste0("cpg:", cpgs$cpg_id), type = "cpg", order = 2L)

# covariates on the trait and the disease
for (to in c("trait_q", "disease_liability")) {
  for (cn in names(ef$covariate_r2)) add_edge("global", "covariates", paste0("cov:", cn), to, ef$covariate_r2[[cn]], 1, role = "covariate")
}
nodes[[length(nodes) + 1]] <- data.table(node = c("trait_q", "disease_liability"), type = c("trait", "liability"), order = 4L)

edges <- rbindlist(edges); nodes <- rbindlist(nodes)
inflow <- edges[, .(r2_in = sum(r2)), by = to]
nodes <- merge(nodes, inflow, by.x = "node", by.y = "to", all.x = TRUE)
nodes[is.na(r2_in), r2_in := 0]
if (any(nodes$r2_in >= 1)) stop("variance explained >= 1 for: ", paste(nodes[r2_in >= 1]$node, collapse = ", "))
nodes[, resid_sd := sqrt(1 - r2_in)]
setorder(nodes, order, node)

# ---- simulate the graph for the cohort ------------------------------------------------
vals <- simulate_dag(nodes, edges, inputs, n)
threshold <- qnorm(1 - ef$disease_prevalence)
clinical <- data.table(sample = samples$sample, population = samples$population,
                       superpopulation = samples$region, sex = samples$sex,
                       age = covs$age, bmi = covs$bmi,
                       trait_q = round(50 + 10 * vals$trait_q, 2),
                       disease = as.integer(vals$disease_liability > threshold))
cat(sprintf("disease cases: %d of %d\n", sum(clinical$disease), n))

meth <- t(vapply(seq_len(nrow(cpgs)), function(k) cpgs$m_mean[k] + cpgs$m_scale[k] * vals[[paste0("cpg:", cpgs$cpg_id[k])]], numeric(n)))
dimnames(meth) <- list(cpgs$cpg_id, samples$sample)

# ---- expression: teaching genes and the disease spike -----------------------------------
counts <- fread(file.path(work, "expression_real", "counts_all.tsv.gz"))
cmat <- as.matrix(counts[, -1]); rownames(cmat) <- counts$gene_id
cmat <- cmat[, samples$sample]
sf <- fread(file.path(work, "expression_real", "size_factors.tsv"))
sf <- sf$size_factor[match(samples$sample, sf$sample)]
genes_sim <- rbindlist(genes_sim)
sim_counts <- t(vapply(genes_sim$gene_id, function(gn) {
  mu <- 2^(ef$teaching_gene_log2_mean + ef$teaching_gene_log2_sd * vals[[paste0("expr:", gn)]]) * sf
  rnbinom(n, mu = mu, size = 1 / ef$teaching_gene_dispersion)
}, numeric(n)))
colnames(sim_counts) <- samples$sample

genes <- read_genes(file.path(meta_dir, "gencode.v39.annotation.gtf.gz"))
genes_real <- fread(file.path(work, "expression_real", "genes.tsv"))
in_loci <- vapply(seq_len(nrow(genes_real)), function(k)
  any(loci$chr == genes_real$chr[k] & loci$start - 1e6 <= genes_real$tss[k] & loci$end + 1e6 >= genes_real$tss[k]), logical(1))
well_expressed <- genes_real$gene_id %in% rownames(cmat)[apply(cmat, 1, median) >= 30]
pool <- sort(genes_real$gene_id[!in_loci & well_expressed])
de <- cfg$de
de_ids <- sample(pool, de$n_up + de$n_down)
de_tab <- data.table(gene_id = de_ids, direction = rep(c("up", "down"), c(de$n_up, de$n_down)))
de_tab[, log2fc := ifelse(direction == "up", 1, -1) * pmax(0.2, rnorm(.N, de$log2fc_mean, de$log2fc_sd))]
cases <- which(clinical$disease == 1)
for (k in seq_len(nrow(de_tab))) {
  x <- cmat[de_tab$gene_id[k], cases]
  fc <- 2^de_tab$log2fc[k]
  cmat[de_tab$gene_id[k], cases] <- if (fc > 1) x + rpois(length(x), x * (fc - 1)) else rbinom(length(x), x, fc)
}
de_tab <- merge(de_tab, genes_real[, .(gene_id, gene_name, chr, tss)], by = "gene_id")
setorder(de_tab, direction, gene_id)

all_counts <- rbind(cmat, sim_counts)
all_genes <- rbind(genes, genes_sim)
cfg$expression$gene_types <- c(cfg$expression$gene_types, "teaching")
norm <- normalise_expression(all_counts, all_genes, cfg)
pcs <- sample_pcs(norm$int, cfg$expression$n_expression_pcs)

# ---- write -------------------------------------------------------------------------------
cat("writing outputs\n")
write_tsv(clinical, file.path(cohort_dir, "clinical", "clinical.tsv"))
write_tsv(samples[, .(sample, population, superpopulation = region, sex)], file.path(cohort_dir, "samples.tsv"))
write_tsv(matrix_table(round(meth, 4), "cpg_id"), file.path(cohort_dir, "methylation", "methylation_m.tsv.gz"))
write_tsv(cpgs[, .(cpg_id, chr, pos, promoter_of = ifelse(cluster == "none", "", cluster))], file.path(cohort_dir, "methylation", "cpgs.tsv"))
ex <- file.path(cohort_dir, "expression")
write_tsv(matrix_table(all_counts[rownames(all_counts) %in% all_genes$gene_id, ], "gene_id"), file.path(ex, "counts.tsv.gz"))
write_tsv(all_genes[match(rownames(norm$int), gene_id)], file.path(ex, "genes.tsv"))
write_tsv(matrix_table(norm$int, "gene_id"), file.path(ex, "expression_int.tsv.gz"))
write_tsv(matrix_table(round(norm$log2norm, 6), "gene_id"), file.path(ex, "expression_log2norm.tsv.gz"))
write_tsv(data.table(sample = names(norm$size_factors), size_factor = norm$size_factors), file.path(ex, "size_factors.tsv"))
write_tsv(matrix_table(pcs, "sample"), file.path(ex, "expression_pcs.tsv"))
write_tsv(loci, file.path(cohort_dir, "loci.tsv"))

variants <- rbindlist(variants)
write_tsv(edges, file.path(truth_dir, "edges.tsv"))
write_tsv(nodes, file.path(truth_dir, "nodes.tsv"))
write_tsv(variants, file.path(truth_dir, "causal_variants.tsv"))
write_tsv(de_tab, file.path(truth_dir, "de_genes.tsv"))
write_tsv(cpgs, file.path(truth_dir, "cpg_parameters.tsv"))
write_tsv(genes_sim, file.path(truth_dir, "teaching_genes.tsv"))
write_tsv(data.table(sample = samples$sample, as.data.table(vals[grep("^(expr|trait|disease)", names(vals))])),
          file.path(truth_dir, "latent_small_cohort.tsv.gz"))
writeLines(toJSON(list(seed = cfg$seed, disease_threshold = threshold,
                       disease_prevalence = ef$disease_prevalence,
                       trait_q = list(offset = 50, scale = 10),
                       teaching_gene = list(log2_mean = ef$teaching_gene_log2_mean, log2_sd = ef$teaching_gene_log2_sd,
                                            dispersion = ef$teaching_gene_dispersion),
                       built_with = "R"), auto_unbox = TRUE, pretty = TRUE, digits = NA),
           file.path(truth_dir, "model.json"))

# truth.tsv: one readable row per causal variant effect, the answer key for the posts
v_edges <- edges[startsWith(from, "var:")]
v_edges[, id := sub("^var:", "", from)]
truth <- merge(v_edges, unique(variants[, .(id, role, chr, pos, ref, alt, af, af_eur, af_yri)]), by = "id")
downstream <- edges[role %in% c("cpg_to_expr", "expr_to_trait")]
# follow a node down the cpg -> expression -> trait chain; returns the path and the product
# of the edge weights, so a variant's effect on the endpoint is per_allele_effect * gain
follow <- function(to) {
  out <- to; cur <- to; gain <- 1
  repeat {
    nx <- downstream[from == cur]
    if (!nrow(nx)) break
    gain <- gain * nx$weight[1]; cur <- nx$to[1]; out <- paste(out, "->", cur)
  }
  list(path = out, endpoint = cur, gain = gain)
}
truth[, path := vapply(to, function(x) follow(x)$path, "")]
truth[, endpoint := vapply(to, function(x) follow(x)$endpoint, "")]
truth[, endpoint_effect := weight * vapply(to, function(x) follow(x)$gain, 0)]
truth <- truth[, .(locus, scenario, variant_id = id, variant_role = role.y, chr, pos, ref, effect_allele = alt,
                   af, af_eur, af_yri, acts_on = to, effect_type = role.x, path, r2_target = r2,
                   per_allele_effect = weight, endpoint, endpoint_effect)]
setorder(truth, locus, variant_role, acts_on)
write_tsv(truth, file.path(truth_dir, "truth.tsv"))
print(truth[, .(locus, scenario, variant_role, variant_id, effect_type, path)])
