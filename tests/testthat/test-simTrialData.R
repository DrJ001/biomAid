# tests/testthat/test-simTrialData.R
# Tests for simTrialData() — pure R, no ASReml needed.
# Note: n_fa must be < ngroup (nsite for MET-only, ntreat*nsite for multi-treatment).
# Small tests use n_fa=1 with nsite=3, or n_fa=2 with nsite>=4.

# ---------------------------------------------------------------------------
# 1. Basic output structure
# ---------------------------------------------------------------------------
test_that("simTrialData returns a named list with 'data' and 'params'", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      seed = 1L, verbose = FALSE)
  expect_named(out, c("data", "params"))
  expect_s3_class(out$data, "data.frame")
  expect_true(is.list(out$params))
})

# ---------------------------------------------------------------------------
# 2. Data dimensions — MET-only
# ---------------------------------------------------------------------------
test_that("MET-only: total plots = nvar * nsite * nrep", {
  nvar <- 4L; nsite <- 3L; nrep <- 2L
  out <- simTrialData(nvar = nvar, nsite = nsite, n_fa = 1L, nrep = nrep,
                      seed = 2L, verbose = FALSE)
  expect_equal(nrow(out$data), nvar * nsite * nrep)
})

test_that("MET-only: data has required columns (no Treatment / TSite)", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      seed = 3L, verbose = FALSE)
  expect_true(all(c("Site", "Variety", "Rep", "Row", "Column", "yield")
                  %in% names(out$data)))
  expect_false("Treatment" %in% names(out$data))
  expect_false("TSite"     %in% names(out$data))
})

# ---------------------------------------------------------------------------
# 3. Data dimensions — multi-treatment
# ---------------------------------------------------------------------------
test_that("multi-treatment: total plots = nvar * nsite * ntreat * nrep", {
  nvar <- 4L; nsite <- 3L; ntreat <- 2L; nrep <- 2L
  out <- simTrialData(nvar = nvar, nsite = nsite,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = nrep,
                      seed = 4L, verbose = FALSE)
  expect_equal(nrow(out$data), nvar * nsite * ntreat * nrep)
})

test_that("multi-treatment: data has Treatment and TSite columns", {
  out <- simTrialData(nvar = 4L, nsite = 3L,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = 2L,
                      seed = 5L, verbose = FALSE)
  expect_true(all(c("Site", "Treatment", "TSite", "Variety",
                    "Rep", "Row", "Column", "yield") %in% names(out$data)))
})

# ---------------------------------------------------------------------------
# 4. Factor columns
# ---------------------------------------------------------------------------
test_that("MET-only: Site, Variety, Rep are factors", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      seed = 6L, verbose = FALSE)
  d <- out$data
  expect_s3_class(d$Site,    "factor")
  expect_s3_class(d$Variety, "factor")
  expect_s3_class(d$Rep,     "factor")
})

test_that("multi-treatment: Treatment and TSite are factors", {
  out <- simTrialData(nvar = 4L, nsite = 3L,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = 2L,
                      seed = 7L, verbose = FALSE)
  expect_s3_class(out$data$Treatment, "factor")
  expect_s3_class(out$data$TSite,     "factor")
})

# ---------------------------------------------------------------------------
# 5. Factor level counts
# ---------------------------------------------------------------------------
test_that("number of Site levels equals nsite", {
  out <- simTrialData(nvar = 4L, nsite = 5L, n_fa = 2L, nrep = 2L,
                      seed = 8L, verbose = FALSE)
  expect_equal(nlevels(out$data$Site), 5L)
})

test_that("number of Variety levels equals nvar", {
  out <- simTrialData(nvar = 9L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      seed = 9L, verbose = FALSE)
  expect_equal(nlevels(out$data$Variety), 9L)
})

test_that("TSite has nsite * ntreat levels", {
  nsite <- 3L; ntreat <- 3L
  out <- simTrialData(nvar = 4L, nsite = nsite,
                      treatments = c("T0", "T1", "T2"), n_fa = 2L, nrep = 2L,
                      seed = 10L, verbose = FALSE)
  expect_equal(nlevels(out$data$TSite), nsite * ntreat)
})

test_that("Treatment factor levels match 'treatments' argument", {
  trt <- c("Dry", "Irr", "Wet")
  out <- simTrialData(nvar = 6L, nsite = 3L, treatments = trt,
                      n_fa = 2L, nrep = 2L, seed = 11L, verbose = FALSE)
  expect_equal(levels(out$data$Treatment), trt)
})

test_that("Rep factor has nrep levels", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 3L,
                      seed = 12L, verbose = FALSE)
  expect_equal(nlevels(out$data$Rep), 3L)
})

# ---------------------------------------------------------------------------
# 6. TSite separator
# ---------------------------------------------------------------------------
test_that("TSite labels contain sep character", {
  out <- simTrialData(nvar = 4L, nsite = 3L,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = 2L,
                      sim.options = list(sep = "_"),
                      seed = 13L, verbose = FALSE)
  expect_true(all(grepl("_", levels(out$data$TSite), fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 7. Reproducibility
# ---------------------------------------------------------------------------
test_that("same seed produces identical data", {
  out1 <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                       seed = 42L, verbose = FALSE)
  out2 <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                       seed = 42L, verbose = FALSE)
  expect_identical(out1$data$yield, out2$data$yield)
})

test_that("different seeds produce different data", {
  out1 <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                       seed = 1L, verbose = FALSE)
  out2 <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                       seed = 2L, verbose = FALSE)
  expect_false(identical(out1$data$yield, out2$data$yield))
})

# ---------------------------------------------------------------------------
# 8. params structure — MET-only
# ---------------------------------------------------------------------------
test_that("MET-only params has G, Lambda, Psi, site_means, n_fa", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      seed = 14L, verbose = FALSE)
  expect_true(all(c("G", "Lambda", "Psi", "site_means", "n_fa")
                  %in% names(out$params)))
  expect_false("treat_effects" %in% names(out$params))
  expect_false("g_arr"         %in% names(out$params))
})

test_that("MET-only: G is nsite x nsite", {
  nsite <- 4L
  out <- simTrialData(nvar = 4L, nsite = nsite, n_fa = 2L, nrep = 2L,
                      seed = 15L, verbose = FALSE)
  expect_equal(dim(out$params$G), c(nsite, nsite))
})

test_that("MET-only: Lambda is nsite x n_fa", {
  nsite <- 4L; n_fa <- 2L
  out <- simTrialData(nvar = 4L, nsite = nsite, n_fa = n_fa, nrep = 2L,
                      seed = 16L, verbose = FALSE)
  expect_equal(dim(out$params$Lambda), c(nsite, n_fa))
})

test_that("MET-only: Psi has length nsite", {
  nsite <- 4L
  out <- simTrialData(nvar = 4L, nsite = nsite, n_fa = 2L, nrep = 2L,
                      seed = 17L, verbose = FALSE)
  expect_length(out$params$Psi, nsite)
})

test_that("MET-only: n_fa stored correctly in params", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      seed = 18L, verbose = FALSE)
  expect_equal(out$params$n_fa, 1L)
})

# ---------------------------------------------------------------------------
# 9. params structure — multi-treatment
# ---------------------------------------------------------------------------
test_that("multi-treatment params has treat_effects and g_arr", {
  out <- simTrialData(nvar = 4L, nsite = 3L,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = 2L,
                      seed = 19L, verbose = FALSE)
  expect_true("treat_effects" %in% names(out$params))
  expect_true("g_arr"         %in% names(out$params))
})

test_that("multi-treatment: G is (ntreat*nsite) x (ntreat*nsite)", {
  nsite <- 3L; ntreat <- 2L
  out <- simTrialData(nvar = 4L, nsite = nsite,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = 2L,
                      seed = 20L, verbose = FALSE)
  ngroup <- nsite * ntreat
  expect_equal(dim(out$params$G), c(ngroup, ngroup))
})

test_that("multi-treatment: g_arr is nvar x ngroup matrix", {
  nvar <- 6L; nsite <- 3L; ntreat <- 2L
  out <- simTrialData(nvar = nvar, nsite = nsite,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = 2L,
                      seed = 21L, verbose = FALSE)
  expect_equal(dim(out$params$g_arr), c(nvar, nsite * ntreat))
})

test_that("treat_effects names match treatments argument", {
  trt <- c("N0", "N1", "N2")
  out <- simTrialData(nvar = 4L, nsite = 3L,
                      treatments = trt, n_fa = 2L, nrep = 2L,
                      seed = 22L, verbose = FALSE)
  expect_equal(names(out$params$treat_effects), trt)
})

# ---------------------------------------------------------------------------
# 10. G matrix properties
# ---------------------------------------------------------------------------
test_that("G is symmetric positive definite", {
  out <- simTrialData(nvar = 6L, nsite = 4L, n_fa = 2L, nrep = 2L,
                      seed = 23L, verbose = FALSE)
  G <- out$params$G
  expect_true(isSymmetric(G, tol = 1e-10))
  expect_true(all(eigen(G, only.values = TRUE)$values > 0))
})

test_that("G = Lambda Lambda' + diag(Psi)", {
  out <- simTrialData(nvar = 4L, nsite = 4L, n_fa = 2L, nrep = 2L,
                      seed = 24L, verbose = FALSE)
  G_reconstructed <- out$params$Lambda %*% t(out$params$Lambda) +
                     diag(out$params$Psi)
  expect_equal(out$params$G, G_reconstructed, tolerance = 1e-10)
})

test_that("genetic correlations are in (-1, 1)", {
  out <- simTrialData(nvar = 4L, nsite = 4L, n_fa = 2L, nrep = 2L,
                      seed = 25L, verbose = FALSE)
  R_off <- cov2cor(out$params$G)[upper.tri(out$params$G)]
  expect_true(all(R_off > -1 & R_off < 1))
})

# ---------------------------------------------------------------------------
# 11. Ordering and yield quality
# ---------------------------------------------------------------------------
test_that("data is sorted by Site then Row then Column", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      seed = 26L, verbose = FALSE)
  d   <- out$data
  expect_equal(seq_len(nrow(d)), order(d$Site, d$Row, d$Column))
})

test_that("yield is numeric, finite, no NAs", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      seed = 27L, verbose = FALSE)
  expect_true(is.numeric(out$data$yield))
  expect_false(anyNA(out$data$yield))
  expect_true(all(is.finite(out$data$yield)))
})

# ---------------------------------------------------------------------------
# 12. sim.options — labels
# ---------------------------------------------------------------------------
test_that("variety_prefix applied to Variety labels", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      sim.options = list(variety_prefix = "Genotype"),
                      seed = 28L, verbose = FALSE)
  expect_true(all(grepl("^Genotype", levels(out$data$Variety))))
})

test_that("site_prefix applied to Site labels", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      sim.options = list(site_prefix = "Trial"),
                      seed = 29L, verbose = FALSE)
  expect_true(all(grepl("^Trial", levels(out$data$Site))))
})

# ---------------------------------------------------------------------------
# 13. Layout — auto-computed per site based on observed variety count
# ---------------------------------------------------------------------------
test_that("balanced layout: total plots = nvar * nsite * nrep (multi-treatment)", {
  out <- simTrialData(nvar = 6L, nsite = 3L,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = 2L,
                      incidence = "balanced",
                      seed = 30L, verbose = FALSE)
  expect_equal(nrow(out$data), 6L * 3L * 2L * 2L)
})

test_that("unbalanced layout: total plots < nvar * nsite * nrep * ntreat", {
  out <- simTrialData(nvar = 10L, nsite = 4L,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = 2L,
                      incidence = "unbalanced",
                      seed = 31L, verbose = FALSE)
  expect_lt(nrow(out$data), 10L * 4L * 2L * 2L)
})

# ---------------------------------------------------------------------------
# 14. sim.options — treat_effects shift means
# ---------------------------------------------------------------------------
test_that("positive treat_effects increases mean yield for that treatment", {
  out <- simTrialData(nvar = 10L, nsite = 5L,
                      treatments = c("T0", "T1"), n_fa = 2L,
                      sim.options = list(treat_effects = c(0, 2000),
                                         sigma_genetic  = 10,
                                         error_sd       = 50,
                                         rep_sd         = 20,
                                         row_sd         = 10,
                                         col_sd         = 10),
                      seed = 31L, verbose = FALSE)
  means <- tapply(out$data$yield, out$data$Treatment, mean)
  expect_gt(means["T1"], means["T0"])
})

test_that("treat_effects wrong length errors", {
  expect_error(
    simTrialData(nsite = 3L, treatments = c("T0", "T1", "T2"),
                 n_fa = 2L,
                 sim.options = list(treat_effects = c(0, 100)),
                 verbose = FALSE),
    "'treat_effects' must have length"
  )
})

# ---------------------------------------------------------------------------
# 15. sim.options — outfile
# ---------------------------------------------------------------------------
test_that("outfile writes CSV", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
               sim.options = list(outfile = tmp),
               seed = 32L, verbose = FALSE)
  expect_true(file.exists(tmp))
  d <- read.csv(tmp, stringsAsFactors = FALSE)
  expect_true("yield" %in% names(d))
})

# ---------------------------------------------------------------------------
# 16. sim.options — unknown key warns
# ---------------------------------------------------------------------------
test_that("unknown sim.options key produces a warning", {
  expect_warning(
    simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                 sim.options = list(not_a_real_option = 99),
                 seed = 33L, verbose = FALSE),
    "Unknown sim.options"
  )
})

# ---------------------------------------------------------------------------
# 17. Error handling
# ---------------------------------------------------------------------------
test_that("single treatment errors", {
  expect_error(
    simTrialData(nsite = 3L, treatments = "T0", verbose = FALSE),
    "'treatments' must have at least 2"
  )
})

test_that("sep in treatment label errors", {
  expect_error(
    simTrialData(nsite = 3L, treatments = c("T-0", "T1"), verbose = FALSE),
    "'sep'.*must not appear"
  )
})

test_that("n_fa >= ngroup errors (MET-only: ngroup = nsite)", {
  expect_error(
    simTrialData(nvar = 4L, nsite = 2L, n_fa = 2L, nrep = 2L,
                 verbose = FALSE),
    "'n_fa'.*must be less than"
  )
})

test_that("n_fa < 1 errors", {
  expect_error(
    simTrialData(nvar = 4L, nsite = 4L, n_fa = 0L, nrep = 2L,
                 verbose = FALSE),
    "'n_fa' must be at least 1"
  )
})

# ---------------------------------------------------------------------------
# 18. Prime nvar produces 1 x nvar strip layout
# ---------------------------------------------------------------------------
test_that("prime nvar produces rows_per_block = 1 layout", {
  # nvar=7 is prime -> .best_dims returns rows=1, cols=7
  # max Row per site = nrep * 1 = nrep
  out <- simTrialData(nvar = 7L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      seed = 34L, verbose = FALSE)
  site_max_row <- tapply(out$data$Row, out$data$Site, max)
  expect_true(all(site_max_row == 2L))   # nrep=2, rows_per_block=1
})

# ---------------------------------------------------------------------------
# 19. Row values in valid range (MET-only, balanced — fixed layout per site)
# ---------------------------------------------------------------------------
test_that("Row values within expected range per site (balanced)", {
  nvar <- 4L; nrep <- 2L
  out <- simTrialData(nvar = nvar, nsite = 3L, n_fa = 1L, nrep = nrep,
                      incidence = "balanced",
                      seed = 35L, verbose = FALSE)
  d    <- out$data
  dims <- c(rows = as.integer(floor(sqrt(nvar))),
            cols = as.integer(nvar %/% floor(sqrt(nvar))))
  nrow_e <- nrep * dims["rows"]
  expect_true(all(d$Row >= 1L & d$Row <= nrow_e))
})

# ---------------------------------------------------------------------------
# 21. incidence = "balanced" (default) — all varieties in all sites
# ---------------------------------------------------------------------------
test_that("balanced incidence: all nvar varieties in every site", {
  nvar <- 6L; nsite <- 3L
  out <- simTrialData(nvar = nvar, nsite = nsite, n_fa = 1L, nrep = 2L,
                      incidence = "balanced", seed = 40L, verbose = FALSE)
  expect_equal(unname(colSums(out$params$incidence)), rep(nvar, nsite))
  expect_equal(nrow(out$data), nvar * nsite * 2L)
})

# ---------------------------------------------------------------------------
# 22. incidence = "unbalanced" — fewer plots than balanced
# ---------------------------------------------------------------------------
test_that("unbalanced: total plots < nvar * nsite * nrep", {
  out <- simTrialData(nvar = 20L, nsite = 6L, n_fa = 2L, nrep = 2L,
                      incidence = "unbalanced", seed = 41L, verbose = FALSE)
  expect_lt(nrow(out$data), 20L * 6L * 2L)
})

test_that("unbalanced: every variety in >= 2 sites", {
  out <- simTrialData(nvar = 20L, nsite = 6L, n_fa = 2L, nrep = 2L,
                      incidence = "unbalanced", seed = 42L, verbose = FALSE)
  sites_per_var <- rowSums(out$params$incidence)
  expect_true(all(sites_per_var >= 2L))
})

test_that("unbalanced: every site has >= 2 varieties", {
  out <- simTrialData(nvar = 20L, nsite = 6L, n_fa = 2L, nrep = 2L,
                      incidence = "unbalanced", seed = 43L, verbose = FALSE)
  vars_per_site <- colSums(out$params$incidence)
  expect_true(all(vars_per_site >= 2L))
})

test_that("unbalanced: incidence matrix stored in params", {
  out <- simTrialData(nvar = 10L, nsite = 4L, n_fa = 1L, nrep = 2L,
                      incidence = "unbalanced", seed = 44L, verbose = FALSE)
  inc <- out$params$incidence
  expect_true(is.matrix(inc))
  expect_equal(dim(inc), c(10L, 4L))
  expect_true(all(inc %in% c(0L, 1L)))
})

test_that("unbalanced: plot count matches incidence colSums * nrep", {
  out <- simTrialData(nvar = 15L, nsite = 5L, n_fa = 2L, nrep = 2L,
                      incidence = "unbalanced", seed = 45L, verbose = FALSE)
  expected <- sum(out$params$incidence) * 2L   # each obs appears nrep times
  expect_equal(nrow(out$data), expected)
})

# ---------------------------------------------------------------------------
# 23. User-supplied incidence matrix
# ---------------------------------------------------------------------------
test_that("user matrix: correct plot count", {
  inc <- matrix(1L, nrow = 8L, ncol = 4L)
  inc[1:2, 1:2] <- 0L   # first 2 varieties absent from first 2 sites
  out <- simTrialData(nvar = 8L, nsite = 4L, n_fa = 1L, nrep = 2L,
                      incidence = inc, seed = 50L, verbose = FALSE)
  expected <- sum(inc) * 2L
  expect_equal(nrow(out$data), expected)
})

test_that("user matrix: incidence stored correctly in params", {
  inc <- matrix(c(1,1,0, 1,1,1, 1,0,1, 1,1,1), nrow = 4L, ncol = 3L)
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      incidence = inc, seed = 51L, verbose = FALSE)
  expect_equal(out$params$incidence, matrix(as.integer(inc), nrow=4L, ncol=3L,
                                             dimnames=dimnames(out$params$incidence)))
})

test_that("user logical matrix accepted", {
  inc <- matrix(TRUE, nrow = 6L, ncol = 3L)
  inc[1, 1] <- FALSE
  expect_no_error(
    simTrialData(nvar = 6L, nsite = 3L, n_fa = 1L, nrep = 2L,
                 incidence = inc, seed = 52L, verbose = FALSE)
  )
})

test_that("user matrix wrong dimensions errors", {
  inc <- matrix(1L, nrow = 5L, ncol = 3L)   # wrong nvar
  expect_error(
    simTrialData(nvar = 6L, nsite = 3L, n_fa = 1L,
                 incidence = inc, verbose = FALSE),
    "must be 6 x 3"
  )
})

test_that("user matrix non-01 values errors", {
  inc <- matrix(c(0,1,2,1,1,0,1,1,1), nrow = 3L, ncol = 3L)
  expect_error(
    simTrialData(nvar = 3L, nsite = 3L, n_fa = 1L,
                 incidence = inc, verbose = FALSE),
    "0/1 or TRUE/FALSE"
  )
})

test_that("user matrix all-zero row errors", {
  inc <- matrix(1L, nrow = 4L, ncol = 3L)
  inc[2, ] <- 0L   # variety 2 in no sites
  expect_error(
    simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L,
                 incidence = inc, verbose = FALSE),
    "appear in no sites"
  )
})

test_that("user matrix site with <2 varieties errors", {
  inc <- matrix(1L, nrow = 4L, ncol = 3L)
  inc[2:4, 1] <- 0L   # only 1 variety in site 1
  expect_error(
    simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L,
                 incidence = inc, verbose = FALSE),
    "fewer than 2 varieties"
  )
})

# ---------------------------------------------------------------------------
# 24. Bad incidence string errors informatively
# ---------------------------------------------------------------------------
test_that("unknown incidence string errors", {
  expect_error(
    simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L,
                 incidence = "random", verbose = FALSE),
    "\"balanced\", \"unbalanced\""
  )
})

# ---------------------------------------------------------------------------
# 20. G rownames/colnames match group labels
# ---------------------------------------------------------------------------
test_that("MET-only: G rownames match site labels", {
  out <- simTrialData(nvar = 4L, nsite = 3L, n_fa = 1L, nrep = 2L,
                      sim.options = list(site_prefix = "S"),
                      seed = 36L, verbose = FALSE)
  expect_equal(rownames(out$params$G), levels(out$data$Site))
})

test_that("multi-treatment: G rownames match TSite levels", {
  out <- simTrialData(nvar = 4L, nsite = 3L,
                      treatments = c("T0", "T1"), n_fa = 1L, nrep = 2L,
                      seed = 37L, verbose = FALSE)
  expect_equal(rownames(out$params$G), levels(out$data$TSite))
})

# ---------------------------------------------------------------------------
# 25. sim.options — unbalanced incidence controls
# ---------------------------------------------------------------------------

test_that("core_pct=0.50 increases proportion of core varieties", {
  # With a larger core_pct, more varieties get the dense coverage pattern.
  # We can verify the structure computes without error and incidence is valid.
  out <- simTrialData(nvar = 20L, nsite = 6L, n_fa = 2L, nrep = 2L,
                      incidence   = "unbalanced",
                      sim.options = list(core_pct = 0.50),
                      seed = 70L, verbose = FALSE)
  inc <- out$params$incidence
  expect_equal(dim(inc), c(20L, 6L))
  expect_true(all(inc %in% c(0L, 1L)))
  expect_true(all(rowSums(inc) >= 1L))
  expect_true(all(colSums(inc) >= 2L))
})

test_that("reg_max_sites_pct=1.0 keeps all varieties fully connected", {
  # Setting both min and max to 1.0 forces every variety into every site
  out <- simTrialData(nvar = 10L, nsite = 5L, n_fa = 2L, nrep = 2L,
                      incidence   = "unbalanced",
                      sim.options = list(core_pct           = 0.20,
                                         reg_min_sites_pct  = 1.00,
                                         reg_max_sites_pct  = 1.00),
                      seed = 71L, verbose = FALSE)
  # All regular varieties must appear in all sites
  inc <- out$params$incidence
  n_core <- max(1L, ceiling(10L * 0.20))
  expect_true(all(rowSums(inc[(n_core + 1L):10L, ]) == 5L))
})

test_that("min_vars_pct=0.90 ensures sites have >= 90% of nvar varieties", {
  nvar <- 20L; nsite <- 6L
  out <- simTrialData(nvar = nvar, nsite = nsite, n_fa = 2L, nrep = 2L,
                      incidence   = "unbalanced",
                      sim.options = list(min_vars_pct = 0.90),
                      seed = 72L, verbose = FALSE)
  min_expected <- ceiling(nvar * 0.90)
  expect_true(all(colSums(out$params$incidence) >= min_expected))
})

test_that("unbalanced incidence options ignored when incidence = 'balanced'", {
  # Passing unbalanced options with balanced incidence should not warn or error
  expect_no_warning(
    simTrialData(nvar = 6L, nsite = 3L, n_fa = 1L, nrep = 2L,
                 incidence   = "balanced",
                 sim.options = list(core_pct = 0.50),
                 seed = 73L, verbose = FALSE)
  )
})

test_that("invalid core_pct errors informatively", {
  expect_error(
    simTrialData(nvar = 10L, nsite = 4L, n_fa = 1L,
                 incidence   = "unbalanced",
                 sim.options = list(core_pct = 1.5),
                 verbose = FALSE),
    "core_pct"
  )
})

test_that("reg_max_sites_pct < reg_min_sites_pct errors", {
  expect_error(
    simTrialData(nvar = 10L, nsite = 4L, n_fa = 1L,
                 incidence   = "unbalanced",
                 sim.options = list(reg_min_sites_pct = 0.70,
                                    reg_max_sites_pct = 0.50),
                 verbose = FALSE),
    "reg_max_sites_pct.*>=.*reg_min_sites_pct"
  )
})

test_that("invalid min_vars_pct errors informatively", {
  expect_error(
    simTrialData(nvar = 10L, nsite = 4L, n_fa = 1L,
                 incidence   = "unbalanced",
                 sim.options = list(min_vars_pct = 0),
                 verbose = FALSE),
    "min_vars_pct"
  )
})
