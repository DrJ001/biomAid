# tests/testthat/test-waldTest-e2e.R
# End-to-end coverage for waldTest.asreml() and print.waldTest().
# predict() is imported from stats so can be mocked.

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
make_wt_pred <- function(treatments = c("N0","N1","N2"), seed = 1L) {
  set.seed(seed)
  n  <- length(treatments)
  pv <- data.frame(
    Treatment       = factor(treatments, levels = treatments),
    predicted.value = rnorm(n, 50, 5),
    std.error       = runif(n, 0.5, 1.5),
    status          = factor(rep("Estimable", n)),
    stringsAsFactors = FALSE
  )
  A    <- matrix(rnorm(n * n) * 0.3, n, n)
  vcov <- crossprod(A) + diag(0.3, n)
  list(pvals = pv, vcov = vcov)
}

# ---------------------------------------------------------------------------
# 1. waldTest.asreml: F-test uses model$nedf
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() F-test uses model$nedf", {
  p <- make_wt_pred()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- structure(list(nedf = 80L), class = "asreml")
  res <- waldTest.asreml(m,
                          classify = "Treatment",
                          cc   = list(list(coef  = c("N0","N1"),
                                           type  = "con",
                                           comp  = c(-1, 1))),
                          test = "F")
  expect_false(is.null(res$Contrasts))
  expect_true("F.Statistic" %in% names(res$Contrasts))
  expect_equal(res$test, "F")
})

# ---------------------------------------------------------------------------
# 2. waldTest.asreml: Wald test (no df_error needed)
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() Wald test works", {
  p <- make_wt_pred(seed = 2L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- structure(list(nedf = 50L), class = "asreml")
  res <- waldTest.asreml(m,
                          classify = "Treatment",
                          cc   = list(list(coef  = c("N0","N1","N2"),
                                           type  = "con",
                                           comp  = "pairwise")),
                          test = "Wald")
  expect_true("Wald.Statistic" %in% names(res$Contrasts))
  expect_equal(nrow(res$Contrasts), 3L)
})

# ---------------------------------------------------------------------------
# 3. waldTest.asreml: zero test
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() zero test works", {
  p <- make_wt_pred(seed = 3L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- structure(list(nedf = 60L), class = "asreml")
  res <- waldTest.asreml(m,
                          classify = "Treatment",
                          cc   = list(list(coef  = c("N1","N2"),
                                           type  = "zero",
                                           group = "joint")),
                          test = "Wald")
  expect_false(is.null(res$Zero))
  expect_equal(res$Zero$df, 2L)
})

# ---------------------------------------------------------------------------
# 4. waldTest.asreml: errors for non-asreml object
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() errors for non-asreml object", {
  expect_error(
    waldTest.asreml(list(), classify = "Treatment",
                    cc = list(list(coef = "N0", type = "con", comp = 1L))),
    "class.*asreml"
  )
})

# ---------------------------------------------------------------------------
# 5. waldTest.asreml: errors if nedf NULL and test = "F"
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() errors when nedf is NULL and test='F'", {
  p <- make_wt_pred(seed = 4L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m <- structure(list(nedf = NULL), class = "asreml")
  expect_error(
    waldTest.asreml(m, classify = "Treatment",
                    cc   = list(list(coef = c("N0","N1"), type = "con", comp = c(-1,1))),
                    test = "F"),
    "nedf.*NULL"
  )
})

# ---------------------------------------------------------------------------
# 6. print.waldTest: prints Contrasts section
# ---------------------------------------------------------------------------
test_that("print.waldTest() prints contrast header", {
  p <- make_wt_pred(seed = 5L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))))
  # print.waldTest uses cat() with "Contrast Tests"
  expect_output(print(res), "Contrast")
})

# ---------------------------------------------------------------------------
# 7. print.waldTest: prints Zero section
# ---------------------------------------------------------------------------
test_that("print.waldTest() prints zero test header", {
  p <- make_wt_pred(seed = 6L)
  res <- waldTest(p,
                  cc = list(list(coef  = c("N0","N1","N2"),
                                 type  = "zero",
                                 group = "all")))
  expect_output(print(res), "Zero")
})

# ---------------------------------------------------------------------------
# 8. print.waldTest: shows adjustment method when not "none"
# ---------------------------------------------------------------------------
test_that("print.waldTest() shows p-value adjustment method", {
  p <- make_wt_pred(seed = 7L)
  res <- waldTest(p,
                  cc     = list(list(coef = c("N0","N1","N2"),
                                     type = "con", comp = "pairwise")),
                  adjust = "bonferroni")
  expect_output(print(res), "bonferroni")
})
