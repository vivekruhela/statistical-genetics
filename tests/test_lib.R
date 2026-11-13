# Unit tests for lib/statgen.R. Run: Rscript -e 'testthat::test_file("tests/test_lib.R")'
library(testthat)
source(if (file.exists("lib/statgen.R")) "lib/statgen.R" else "../lib/statgen.R")

test_that("inverse normal transform is rank based and symmetric", {
  x <- c(5, 1, 3, 3, 10)
  z <- inverse_normal(x)
  expect_equal(order(z), order(rank(x, ties.method = "first")))
  expect_equal(z[3], z[4])
  expect_equal(sum(inverse_normal(1:7)), 0, tolerance = 1e-12)
})

test_that("size factors match the median-of-ratios definition", {
  m <- matrix(c(10, 20, 30, 20, 40, 60, 5, 10, 15), 3)
  expect_equal(unname(size_factors(m)), c(1, 2, 0.5) / exp(mean(log(c(1, 2, 0.5)))))
})

test_that("cis scan agrees with lm", {
  set.seed(1)
  n <- 80
  g <- matrix(rbinom(5 * n, 2, 0.4), 5, n, dimnames = list(paste0("v", 1:5), NULL))
  cv <- matrix(rnorm(2 * n), n)
  e <- matrix(0.5 * g[3, ] + cv[, 1] + rnorm(n), 1, dimnames = list("gene", NULL))
  res <- cis_scan(residualise(e, cv), 300, residualise(g, cv), c(100, 200, 300, 400, 500), 1000, 2)
  fit <- summary(lm(e[1, ] ~ g[3, ] + cv))
  expect_equal(res$lead_id, "v3")
  expect_equal(res$lead_mlog10p, -log10(fit$coefficients[2, 4]), tolerance = 1e-8)
})

test_that("linear association statistics match lm", {
  set.seed(2)
  n <- 200
  G <- matrix(rbinom(3 * n, 2, 0.3), 3)
  X <- cbind(1, rnorm(n), rbinom(n, 1, 0.5))
  y <- 0.3 * G[2, ] + X[, 2] + rnorm(n)
  st <- assoc_stats(G, y, X)
  fit <- summary(lm(y ~ G[2, ] + X[, -1]))
  expect_equal(st$beta[2], fit$coefficients[2, 1], tolerance = 1e-10)
  expect_equal(st$se[2], fit$coefficients[2, 2], tolerance = 1e-10)
})

test_that("logistic score test matches glm's score test", {
  set.seed(3)
  n <- 500
  G <- matrix(rbinom(n, 2, 0.3), 1)
  X <- cbind(1, rnorm(n))
  y <- rbinom(n, 1, plogis(-0.5 + 0.4 * G[1, ] + 0.3 * X[, 2]))
  st <- assoc_stats(G, y, X, binary = TRUE)
  null <- glm(y ~ X[, 2], family = binomial)
  rao <- anova(null, glm(y ~ X[, 2] + G[1, ], family = binomial), test = "Rao")
  expect_equal(st$z^2, rao$Rao[2], tolerance = 1e-4)  # glm stops at its own convergence tolerance
})

test_that("mosaic people carry template alleles and the requested size", {
  set.seed(4)
  H <- matrix(rbinom(50 * 20, 1, 0.5), 50)
  d <- mosaic_dosage(H, seq(1000, 50000, by = 1000), 30, 10000, 0)
  expect_equal(dim(d), c(50L, 30L))
  expect_true(all(d %in% 0:2))
})

test_that("the causal graph gives unit-variance nodes", {
  set.seed(5)
  n <- 20000
  nodes <- data.table(node = c("a", "b"), order = 1:2, resid_sd = c(1, sqrt(1 - 0.3)))
  edges <- data.table(from = "a", to = "b", weight = sqrt(0.3), center = 0)
  v <- simulate_dag(nodes, edges, list(), n)
  expect_equal(var(v$b), 1, tolerance = 0.03)
  expect_equal(cor(v$a, v$b)^2, 0.3, tolerance = 0.03)
})

# ---- LD and reference panels (post R) ----
toy_ld <- function() {
  p <- 30; n <- 1000
  R <- 0.9^abs(outer(1:p, 1:p, "-"))
  b <- rep(0, p); b[8] <- 0.15; b[22] <- -0.12
  list(R = R, n = n, z = as.vector(sqrt(n) * R %*% b + 0.8 * sin(1:p)))
}

test_that("ld_matrix gives correlations and NA for a variant that does not vary", {
  g <- rbind(c(0, 1, 2, 1, 0), c(0, 1, 2, 2, 0), c(1, 1, 1, 1, 1))
  r <- ld_matrix(g)
  expect_equal(r[1, 2], cor(g[1, ], g[2, ]))
  expect_true(all(is.na(r[3, ])))
})

test_that("an LD matrix survives a round trip to disk", {
  t <- toy_ld(); f <- tempfile()
  write_ld(t$R, paste0("v", 1:30), f)
  r <- read_ld(f)
  expect_equal(unname(r), t$R, tolerance = 1e-7)
  expect_equal(rownames(r), paste0("v", 1:30))
})

test_that("ld_check agrees with susieR's estimate_s_rss and kriging_rss", {
  skip_if_not_installed("susieR")
  t <- toy_ld()
  R2 <- 0.7^abs(outer(1:30, 1:30, "-"))
  ch <- ld_check(t$z, R2, t$n)
  k <- suppressWarnings(susieR::kriging_rss(t$z, R2, n = t$n))$conditional_dist
  expect_equal(ch$s[1], susieR::estimate_s_rss(t$z, R2, t$n), tolerance = 1e-6)
  expect_equal(ch$expected, k$condmean, tolerance = 1e-6)
  expect_equal(ch$z_std_diff, k$z_std_diff, tolerance = 1e-6)
})
