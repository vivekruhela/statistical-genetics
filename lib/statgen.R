# Shared helpers for the data build (R track).
# lib/statgen.py is the Python twin: the same functions, the same numbers.
# Every function here is deterministic; all randomness lives in the simulate/ scripts.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# A config may name a "_base" config and override only some of its values.
read_config <- function(path = Sys.getenv("STATGEN_CONFIG", "config/cohort.json")) {
  cfg <- fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  if (!is.null(cfg[["_base"]])) {
    base <- read_config(cfg[["_base"]])
    cfg[["_base"]] <- NULL
    cfg <- modifyList(base, cfg)
  }
  cfg
}

# ---- genes -------------------------------------------------------------------

# One row per gene from a GENCODE GTF: id without version, name, type, position, TSS.
read_genes <- function(gtf) {
  g <- fread(cmd = paste("zcat", shQuote(gtf), "| grep -v '^#' | awk -F'\\t' '$3==\"gene\"'"),
             sep = "\t", header = FALSE, quote = "",
             col.names = c("chr", "src", "feature", "start", "end", "score", "strand", "frame", "attr"))
  attr_value <- function(key) sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", g$attr)
  out <- data.table(
    gene_id   = sub("\\..*", "", attr_value("gene_id")),
    gene_name = attr_value("gene_name"),
    gene_type = attr_value("gene_type"),
    chr = g$chr, start = g$start, end = g$end, strand = g$strand
  )
  out[, tss := ifelse(strand == "+", start, end)]
  out <- out[!grepl("_PAR_Y", attr_value("gene_id"))]
  unique(out, by = "gene_id")
}

# ---- expression --------------------------------------------------------------

# GEUVADIS count table -> integer matrix, genes x individuals.
# Library columns are named like HG00096.1.M_111124_6; the individual is the part before the
# first dot. Where one person has several libraries we keep the deepest one.
read_geuvadis_counts <- function(path) {
  d <- fread(path)
  id_col <- names(d)[1]
  info_cols <- intersect(names(d), c("TargetID", "Gene_Symbol", "Chr", "Coord"))
  libs <- setdiff(names(d), info_cols)
  m <- as.matrix(d[, ..libs])
  rownames(m) <- sub("\\..*", "", d[[id_col]])
  m <- m[!duplicated(rownames(m)), , drop = FALSE]
  ind <- sub("\\..*", "", libs)
  depth <- colSums(m)
  ord <- order(ind, -depth)
  keep <- ord[!duplicated(ind[ord])]
  m <- m[, keep, drop = FALSE]
  colnames(m) <- ind[keep]
  m <- m[, order(colnames(m)), drop = FALSE]
  storage.mode(m) <- "double"
  round(m)
}

# Rank-based inverse normal transform (Blom offsets), ties get average ranks.
inverse_normal <- function(x) {
  r <- rank(x, ties.method = "average")
  qnorm((r - 3 / 8) / (length(x) + 1 / 4))
}

# DESeq2-style median-of-ratios size factors, using genes with no zero count.
size_factors <- function(counts) {
  pos <- rowSums(counts == 0) == 0
  logc <- log(counts[pos, , drop = FALSE])
  log_geo <- rowMeans(logc)
  exp(apply(logc - log_geo, 2, median))
}

# Filter genes, normalise and transform. Returns a list of matrices (genes x samples).
normalise_expression <- function(counts, genes, cfg) {
  ecfg <- cfg$expression
  autosomes <- paste0("chr", 1:22)
  ann <- genes[gene_type %in% ecfg$gene_types & chr %in% autosomes]
  counts <- counts[rownames(counts) %in% ann$gene_id, , drop = FALSE]
  cpm <- t(t(counts) / colSums(counts)) * 1e6
  ok <- (counts >= ecfg$min_count) & (cpm >= ecfg$min_cpm)
  counts <- counts[rowMeans(ok) >= ecfg$min_sample_fraction, , drop = FALSE]
  counts <- counts[order(rownames(counts)), , drop = FALSE]
  sf <- size_factors(counts)
  log2norm <- log2(t(t(counts) / sf) + 1)
  int <- t(apply(log2norm, 1, inverse_normal))
  dimnames(int) <- dimnames(log2norm)
  list(counts = counts, size_factors = sf, log2norm = log2norm, int = int)
}

# ---- linear algebra -----------------------------------------------------------

# Top k principal component scores of the samples (columns of mat).
# Sign rule shared with Python: the score with the largest magnitude is positive.
sample_pcs <- function(mat, k) {
  x <- scale(t(mat), center = TRUE, scale = FALSE)
  s <- svd(x, nu = k, nv = 0)
  scores <- s$u[, seq_len(k), drop = FALSE] %*% diag(s$d[seq_len(k)], k, k)
  for (j in seq_len(k)) {
    if (scores[which.max(abs(scores[, j])), j] < 0) scores[, j] <- -scores[, j]
  }
  dimnames(scores) <- list(colnames(mat), paste0("PC", seq_len(k)))
  scores
}

# Residuals of each row of y (features x samples) after regressing on covariates
# (samples x q, an intercept is added).
residualise <- function(y, covar) {
  q <- qr(cbind(1, as.matrix(covar)))
  t(qr.resid(q, t(y)))
}

# Centre each row and scale it to unit length; constant rows become all zero.
unit_rows <- function(x) {
  x <- x - rowMeans(x)
  nrm <- sqrt(rowSums(x^2))
  nrm[nrm == 0] <- Inf
  x / nrm
}

# -log10 two-sided p-value of a correlation r with df residual degrees of freedom.
cor_mlog10p <- function(r, df) {
  r <- pmin(abs(r), 1 - 1e-15)
  tstat <- r * sqrt(df / (1 - r^2))
  -(pt(tstat, df, lower.tail = FALSE, log.p = TRUE) + log(2)) / log(10)
}

# ---- genotypes ----------------------------------------------------------------

# Read a PLINK 2 .traw (from --export A-transpose) into a dosage matrix of ALT allele counts,
# variants x samples, plus a variant table. Variant IDs are chr:pos:ref:alt.
read_traw <- function(path) {
  d <- fread(path)
  setnames(d, 1:6, c("chr", "id", "cm", "pos", "counted", "other"))
  sample_cols <- names(d)[-(1:6)]
  g <- as.matrix(d[, ..sample_cols])
  alt <- sub(".*:", "", d$id)
  flip <- d$counted != alt
  g[flip, ] <- 2 - g[flip, ]
  rownames(g) <- d$id
  colnames(g) <- sub("^[^_]*_", "", sample_cols)
  parts <- tstrsplit(d$id, ":", fixed = TRUE)
  info <- data.table(id = d$id, chr = parts[[1]], pos = as.integer(parts[[2]]),
                     ref = parts[[3]], alt = parts[[4]])
  list(info = info, g = g)
}

# ---- cis-eQTL scan ---------------------------------------------------------------

# For each gene: the lead cis variant, and the best variant after conditioning on the lead.
# expr_res and geno_res must already be residualised on the same covariates; n_covar is the
# number of covariate columns used (for degrees of freedom).
cis_scan <- function(expr_res, gene_tss, geno_res, snp_pos, half_window, n_covar) {
  n <- ncol(expr_res)
  df <- n - 2 - n_covar
  e <- unit_rows(expr_res)
  g <- unit_rows(geno_res)
  ord <- order(snp_pos)
  stopifnot(identical(ord, seq_along(snp_pos)))
  lo <- findInterval(gene_tss - half_window - 1, snp_pos) + 1
  hi <- findInterval(gene_tss + half_window, snp_pos)
  res <- vector("list", nrow(e))
  for (i in seq_len(nrow(e))) {
    if (hi[i] < lo[i]) next
    idx <- lo[i]:hi[i]
    r <- as.vector(g[idx, , drop = FALSE] %*% e[i, ])
    best <- which.max(abs(r))
    lead <- idx[best]
    # conditional: remove the lead variant from expression and from the other variants
    gl <- g[lead, ]
    e2 <- e[i, ] - sum(e[i, ] * gl) * gl
    g2 <- g[idx, , drop = FALSE] - (g[idx, , drop = FALSE] %*% gl) %*% t(gl)
    r2 <- as.vector(unit_rows(g2) %*% (e2 / sqrt(sum(e2^2))))
    r2[best] <- 0
    best2 <- which.max(abs(r2))
    res[[i]] <- data.table(
      gene_id = rownames(e)[i], n_cis = length(idx),
      lead_id = rownames(g)[lead], lead_r = r[best],
      lead_mlog10p = cor_mlog10p(r[best], df),
      second_id = rownames(g)[idx[best2]],
      second_mlog10p = cor_mlog10p(r2[best2], df - 1)
    )
  }
  rbindlist(res)
}

# ---- simulation helpers ------------------------------------------------------------

# Squared correlation between rows of a dosage matrix.
ld_r2 <- function(g) {
  u <- unit_rows(g)
  (u %*% t(u))^2
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, path, sep = "\t", quote = FALSE)
}

# Matrix (rows = features) to a table with the row names as the first column.
matrix_table <- function(m, id_name) {
  out <- data.table(rownames(m), m)
  setnames(out, 1, id_name)
  out
}

# ---- the causal model (shared by steps 7 and 9) ----------------------------------

# Simulate every latent node of the causal graph for n people.
#   nodes:  node, order, resid_sd (one row per latent node: factor, cpg, expr, trait, liability)
#   edges:  from, to, weight, center   value(to) += weight * (value(from) - center)
#   inputs: named list of observed vectors (variant dosages "var:<id>", covariates "cov:<name>")
# Nodes are drawn in `order`, one rnorm(n) each, so a run is reproducible from the seed.
simulate_dag <- function(nodes, edges, inputs, n) {
  vals <- inputs
  nodes <- nodes[order(order, node)]
  for (k in seq_len(nrow(nodes))) {
    nd <- nodes$node[k]
    ein <- edges[to == nd]
    v <- nodes$resid_sd[k] * rnorm(n)
    for (j in seq_len(nrow(ein))) {
      v <- v + ein$weight[j] * (vals[[ein$from[j]]] - ein$center[j])
    }
    vals[[nd]] <- v
  }
  vals[nodes$node]
}

# Covariates on the standard scale used by the causal model (fixed constants, not sample
# statistics, so the small cohort and the GWAS cohort share one model).
standard_covariates <- function(age, male, bmi) {
  list(`cov:age` = (age - 47.5) / (55 / sqrt(12)),
       `cov:sex` = (male - 0.5) / 0.5,
       `cov:bmi` = (log(bmi) - log(26)) / 0.15)
}

simulate_covariates <- function(n, male = NULL) {
  age <- round(runif(n, 20, 75))
  if (is.null(male)) male <- rbinom(n, 1, 0.5)
  bmi <- round(exp(rnorm(n, log(26) + 0.004 * (age - 47.5), 0.14)), 1)
  data.table(age = age, male = male, bmi = bmi)
}

# ---- the GWAS cohort (step 9) -------------------------------------------------------

# Read a template haplotype table from step 8, rows aligned to the variant IDs `ids`.
read_haplotypes <- function(path, ids) {
  h <- fread(path)
  m <- as.matrix(h[, -1]); storage.mode(m) <- "integer"
  m[match(ids, h$id), , drop = FALSE]
}

# New people as mosaics of template haplotypes (a simple copying model, as in HAPGEN):
# each new haplotype copies one template, switches to a random template at breakpoints
# that fall on average every segment_bp, and miscopies an allele with probability miscopy.
# Returns ALT dosages, variants x n_people (integer).
mosaic_dosage <- function(H, pos, n_people, segment_bp, miscopy, batch = 500) {
  m <- nrow(H); nt <- ncol(H)
  span <- max(pos) - min(pos) + 1
  out <- matrix(0L, m, n_people)
  for (b0 in seq(1, n_people, by = batch)) {
    people <- b0:min(n_people, b0 + batch - 1)
    nh <- 2 * length(people)
    idx <- matrix(0L, m, nh)
    for (h in seq_len(nh)) {
      k <- rpois(1, span / segment_bp)
      brk <- sort(runif(k, min(pos), max(pos)))
      tmpl <- sample.int(nt, k + 1, replace = TRUE)
      idx[, h] <- tmpl[findInterval(pos, brk) + 1]
    }
    hap <- matrix(H[cbind(rep(seq_len(m), nh), as.vector(idx))], m, nh)
    nflip <- rbinom(1, length(hap), miscopy)
    if (nflip > 0) {
      f <- sample.int(length(hap), nflip)
      hap[f] <- 1L - hap[f]
    }
    out[, people] <- hap[, seq(1, nh, 2), drop = FALSE] + hap[, seq(2, nh, 2), drop = FALSE]
  }
  out
}

# Null logistic regression by IRLS; returns fitted probabilities.
logistic_null <- function(y, X, iter = 50) {
  beta <- rep(0, ncol(X))
  for (it in seq_len(iter)) {
    mu <- 1 / (1 + exp(-as.vector(X %*% beta)))
    w <- mu * (1 - mu)
    step <- solve(crossprod(X, X * w), crossprod(X, y - mu))
    beta <- beta + as.vector(step)
    if (max(abs(step)) < 1e-10) break
  }
  1 / (1 + exp(-as.vector(X %*% beta)))
}

# Per-variant association statistics, adjusting for covariates X (n x q, with intercept).
# Quantitative y: least squares (Frisch-Waugh). Binary y: score test from the null logistic
# fit, with beta = U / I, the usual one-step approximation to the log odds ratio.
# G is variants x people. Returns a table of beta, se, z, mlog10p.
assoc_stats <- function(G, y, X, binary = FALSE) {
  G <- G * 1.0
  if (!binary) {
    q <- qr(X)
    yr <- qr.resid(q, y)
    gy <- as.vector(G %*% yr)
    CG <- t(qr.Q(q)) %*% t(G)
    gg <- rowSums(G^2) - colSums(CG^2)
    gg[gg <= 1e-8] <- NA  # monomorphic in this cohort
    beta <- gy / gg
    df <- length(y) - ncol(X) - 1
    s2 <- (sum(yr^2) - beta * gy) / df
    se <- sqrt(s2 / gg)
  } else {
    mu <- logistic_null(y, X)
    w <- mu * (1 - mu)
    U <- as.vector(G %*% (y - mu))
    XtWX_inv <- solve(crossprod(X, X * w))
    CWG <- crossprod(X * w, t(G))
    I <- as.vector((G^2) %*% w) - colSums(CWG * (XtWX_inv %*% CWG))
    I[I <= 1e-8] <- NA
    beta <- U / I
    se <- 1 / sqrt(I)
  }
  z <- beta / se
  data.table(beta = beta, se = se, z = z,
             mlog10p = -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10))
}
