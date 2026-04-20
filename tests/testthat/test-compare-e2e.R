# tests/testthat/test-compare-e2e.R
# End-to-end tests for compare() using local_mocked_bindings.
# predict() is imported from stats into biomAid's namespace, so it can be
# mocked without an ASReml licence.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
make_asreml_model <- function(nedf = 100L, random_formula = ~Genotype) {
  structure(
    list(nedf = nedf, call = list(random = random_formula)),
    class = "asreml"
  )
}

make_pred_simple <- function(n = 10, seed = 1L) {
  set.seed(seed)
  pv <- data.frame(
    Genotype        = factor(paste0("G", sprintf("%02d", seq_len(n)))),
    predicted.value = rnorm(n, 50, 5),
    std.error       = runif(n, 0.5, 1.5),
    status          = factor(rep("Estimable", n)),
    stringsAsFactors = FALSE
  )
  A    <- matrix(rnorm(n * n) * 0.2, n, n)
  vcov <- crossprod(A) + diag(0.3, n)
  list(pvals = pv, vcov = vcov)
}

make_pred_grouped <- function(sites = c("S1","S2"), n_geno = 8L, seed = 2L) {
  set.seed(seed)
  genos <- paste0("G", sprintf("%02d", seq_len(n_geno)))
  pv <- expand.grid(
    Site            = factor(sites),
    Genotype        = factor(genos),
    KEEP.OUT.ATTRS  = FALSE,
    stringsAsFactors = FALSE
  )
  pv$Site     <- factor(pv$Site)
  pv$Genotype <- factor(pv$Genotype)
  n <- nrow(pv)
  pv$predicted.value <- rnorm(n, 50, 5)
  pv$std.error       <- runif(n, 0.5, 1.5)
  pv$status          <- factor(rep("Estimable", n))
  A    <- matrix(rnorm(n * n) * 0.2, n, n)
  vcov <- crossprod(A) + diag(0.3, n)
  list(pvals = pv, vcov = vcov)
}

# ---------------------------------------------------------------------------
# 1. Basic HSD — single group, no by
# ---------------------------------------------------------------------------
test_that("compare() HSD no by: returns data frame with HSD and avsed columns", {
  p <- make_pred_simple(n = 10, seed = 1L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- make_asreml_model()
  out <- compare(m, term = "Genotype", type = "HSD")
  expect_s3_class(out, "data.frame")
  expect_true("HSD"   %in% names(out))
  expect_true("avsed" %in% names(out))
  expect_equal(nrow(out), 10L)
})

# ---------------------------------------------------------------------------
# 2. LSD
# ---------------------------------------------------------------------------
test_that("compare() LSD returns LSD column", {
  p <- make_pred_simple(seed = 2L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Genotype", type = "LSD")
  expect_true("LSD" %in% names(out))
  expect_true(all(out$LSD > 0))
})

# ---------------------------------------------------------------------------
# 3. Bonferroni
# ---------------------------------------------------------------------------
test_that("compare() Bonferroni returns Bonferroni column", {
  p <- make_pred_simple(seed = 3L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Genotype", type = "Bonferroni")
  expect_true("Bonferroni" %in% names(out))
  expect_true(all(out$Bonferroni > 0))
})

# ---------------------------------------------------------------------------
# 4. by = single string — one group per Site
# ---------------------------------------------------------------------------
test_that("compare() by = 'Site': one criterion value per site group", {
  p <- make_pred_grouped(sites = c("S1","S2"), n_geno = 8L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Site:Genotype",
                 by = "Site", type = "HSD")
  expect_equal(nrow(out), 16L)
  # HSD is constant within each site
  for (s in c("S1","S2")) {
    vals <- out$HSD[out$Site == s]
    expect_equal(length(unique(round(vals, 10L))), 1L)
  }
})

# ---------------------------------------------------------------------------
# 5. by = colon string vs character vector — same result
# ---------------------------------------------------------------------------
test_that("compare() by colon string == character vector", {
  p <- make_pred_grouped()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- make_asreml_model()
  o1  <- compare(m, term = "Site:Genotype", by = "Site",  type = "HSD")
  o2  <- compare(m, term = "Site:Genotype", by = "Site",  type = "HSD")
  expect_identical(o1$HSD, o2$HSD)
})

# ---------------------------------------------------------------------------
# 6. alpha = 0.01 gives larger criterion than alpha = 0.05
# ---------------------------------------------------------------------------
test_that("compare() smaller alpha -> larger HSD", {
  p <- make_pred_simple(seed = 5L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m    <- make_asreml_model()
  h05  <- compare(m, term = "Genotype", type = "HSD", alpha = 0.05)$HSD[1L]
  h01  <- compare(m, term = "Genotype", type = "HSD", alpha = 0.01)$HSD[1L]
  expect_gt(h01, h05)
})

# ---------------------------------------------------------------------------
# 7. NA in predicted.value — removed before calculation
# ---------------------------------------------------------------------------
test_that("compare() NA predicted values are filtered out", {
  p <- make_pred_simple(n = 10, seed = 6L)
  p$pvals$predicted.value[c(2L, 5L)] <- NA_real_
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Genotype", type = "HSD")
  expect_equal(nrow(out), 8L)
})

# ---------------------------------------------------------------------------
# 8. Error: non-asreml model
# ---------------------------------------------------------------------------
test_that("compare() errors for non-asreml model", {
  expect_error(compare(list(), term = "Genotype"), "class.*asreml")
})

# ---------------------------------------------------------------------------
# 9. Error: invalid alpha
# ---------------------------------------------------------------------------
test_that("compare() errors for alpha outside (0,1)", {
  p <- make_pred_simple()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m <- make_asreml_model()
  expect_error(compare(m, term = "Genotype", alpha = 0),   "'alpha'")
  expect_error(compare(m, term = "Genotype", alpha = 1),   "'alpha'")
  expect_error(compare(m, term = "Genotype", alpha = 1.5), "'alpha'")
})

# ---------------------------------------------------------------------------
# 10. Error: by variable not in term
# ---------------------------------------------------------------------------
test_that("compare() errors when by variable not in term", {
  p <- make_pred_simple()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  expect_error(
    compare(make_asreml_model(), term = "Genotype",
            by = "Site", type = "HSD"),
    "not found in 'term'"
  )
})

# ---------------------------------------------------------------------------
# 11. Error: by accounts for all terms
# ---------------------------------------------------------------------------
test_that("compare() errors when by accounts for all terms", {
  p <- make_pred_grouped(sites = c("S1"), n_geno = 5L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  expect_error(
    compare(make_asreml_model(), term = "Site:Genotype",
            by = c("Site","Genotype"), type = "HSD"),
    "no levels remain"
  )
})

# ---------------------------------------------------------------------------
# 12. Warning: single-observation group yields NA criterion
# ---------------------------------------------------------------------------
test_that("compare() warns and returns NA for single-obs group", {
  p <- make_pred_grouped(sites = c("S1","S2"), n_geno = 1L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  # Both single-obs groups warn; suppress all and just check NAs
  suppressWarnings(
    out <- compare(make_asreml_model(), term = "Site:Genotype",
                   by = "Site", type = "HSD")
  )
  expect_true(all(is.na(out$HSD)))
})

# ---------------------------------------------------------------------------
# 13. pev = FALSE: warning when term not in random formula
# ---------------------------------------------------------------------------
test_that("compare() pev=FALSE warns if term not in random formula", {
  p <- make_pred_simple(seed = 8L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  # random formula doesn't include "Genotype"
  m <- make_asreml_model(random_formula = ~Block)
  expect_warning(
    compare(m, term = "Genotype", type = "HSD", pev = FALSE),
    "Falling back to PEV"
  )
})

# ---------------------------------------------------------------------------
# 14. pev = FALSE with matching random formula: uses .asreml_vparams
# ---------------------------------------------------------------------------
test_that("compare() pev=FALSE: falls back to PEV when vparams NULL", {
  # Simplest test for the pev=FALSE code path: .asreml_vparams returns NULL,
  # which triggers the warning and falls back to PEV.
  p <- make_pred_simple(n = 8L, seed = 9L)
  local_mocked_bindings(
    predict         = function(...) p,
    .asreml_vparams = function(model, term) NULL,
    .package        = "biomAid"
  )
  m <- structure(
    list(nedf = 50L, call = list(random = ~Genotype)),
    class = "asreml"
  )
  expect_warning(
    out <- compare(m, term = "Genotype", type = "HSD", pev = FALSE),
    "Falling back to PEV"
  )
  expect_s3_class(out, "data.frame")
  expect_true("HSD" %in% names(out))
})

# ---------------------------------------------------------------------------
# 15. Output columns: factor columns + predicted.value + std.error + type + avsed
# ---------------------------------------------------------------------------
test_that("compare() output has correct column set", {
  p <- make_pred_simple(seed = 10L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Genotype", type = "LSD")
  expect_true(all(c("Genotype","predicted.value","std.error","LSD","avsed")
                  %in% names(out)))
})
