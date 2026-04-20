# tests/testthat/test-fast-e2e.R
# End-to-end tests for fast() using local_mocked_bindings.
# .fa_asreml() is an internal biomAid wrapper, so it can be mocked directly.

# ---------------------------------------------------------------------------
# Fixture (same as test-fast.R)
# ---------------------------------------------------------------------------
.loads <- matrix(
  c( 2.0,  1.0,
     1.5, -0.8,
     1.8,  0.6,
     1.2, -1.2),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("E1","E2","E3","E4"), NULL)
)
.scores <- matrix(
  c( 1.0,  0.5,
     0.2, -0.3,
    -0.5,  0.8),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("G1","G2","G3"), NULL)
)
.spec   <- c(E1 = 0.10, E2 = 0.20, E3 = 0.15, E4 = 0.25)
.term   <- "fa(Site, 2):Genotype"
.sc_long <- data.frame(
  Site  = rep(paste0("Comp", 1:2), each = 3L),
  blupr = c(.scores[, 1L], .scores[, 2L]),
  stringsAsFactors = FALSE
)
.mock_sfa <- list(
  gammas = setNames(list(list("rotated loads" = .loads, "specific var" = .spec)), .term),
  blups  = setNames(list(list(scores = .sc_long)), .term)
)

make_fast_model <- function() {
  # fast() calls eval(model$call$data); store the data frame directly so
  # eval() can find it regardless of which environment fast() runs in.
  fast_mock_data <<- expand.grid(
    Site     = factor(c("E1","E2","E3","E4")),
    Genotype = factor(c("G1","G2","G3")),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  m <- list(call = list(data = quote(fast_mock_data)))
  class(m) <- "asreml"
  m
}

run_fast_e2e <- function(type = "all", ic.num = 2L) {
  local_mocked_bindings(.fa_asreml = function(...) .mock_sfa, .package = "biomAid")
  fast(make_fast_model(), term = .term, type = type, ic.num = ic.num)
}

# ---------------------------------------------------------------------------
# 1. type = "all" returns correct dimensions and columns
# ---------------------------------------------------------------------------
test_that("fast() type='all': correct rows and columns", {
  out <- run_fast_e2e("all")
  expect_equal(nrow(out), 4L * 3L)
  expected <- c("Site","Genotype","loads1","loads2","spec.var",
                "score1","score2","CVE","VE","fitted1","fitted2",
                "OP","dev","stab","iclass","iClassOP","iClassRMSD")
  expect_true(all(expected %in% names(out)))
})

# ---------------------------------------------------------------------------
# 2. CVE = fitted1 + fitted2
# ---------------------------------------------------------------------------
test_that("fast() CVE == fitted1 + fitted2", {
  out <- run_fast_e2e("all")
  expect_equal(out$CVE, out$fitted1 + out$fitted2, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3. VE = CVE + spec.var
# ---------------------------------------------------------------------------
test_that("fast() VE == CVE + spec.var", {
  out <- run_fast_e2e("all")
  expect_equal(out$VE, out$CVE + out$spec.var, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 4. OP = mean(loads1) * score1 — constant within genotype
# ---------------------------------------------------------------------------
test_that("fast() OP is constant within genotype", {
  out <- run_fast_e2e("FAST")
  for (g in c("G1","G2","G3")) {
    ops <- out$OP[out$Genotype == g]
    expect_equal(length(unique(round(ops, 12L))), 1L)
  }
})

# ---------------------------------------------------------------------------
# 5. stab = RMSD of dev across all environments
# ---------------------------------------------------------------------------
test_that("fast() stab is RMSD of dev", {
  out <- run_fast_e2e("FAST")
  for (g in c("G1","G2","G3")) {
    devs <- out$dev[out$Genotype == g]
    expect_equal(unique(round(out$stab[out$Genotype == g], 10L)),
                 round(sqrt(mean(devs^2)), 10L))
  }
})

# ---------------------------------------------------------------------------
# 6. type = "FAST" has no iClass columns
# ---------------------------------------------------------------------------
test_that("fast() type='FAST': no iClass columns", {
  out <- run_fast_e2e("FAST")
  expect_false("iclass"    %in% names(out))
  expect_false("iClassOP"  %in% names(out))
  expect_false("iClassRMSD" %in% names(out))
})

# ---------------------------------------------------------------------------
# 7. type = "iClass" has no FAST-only columns
# ---------------------------------------------------------------------------
test_that("fast() type='iClass': no FAST-only columns", {
  out <- run_fast_e2e("iClass")
  expect_false("OP"   %in% names(out))
  expect_false("stab" %in% names(out))
  expect_false("dev"  %in% names(out))
  expect_true("iclass" %in% names(out))
})

# ---------------------------------------------------------------------------
# 8. iClassRMSD = 0 when ic.num = k
# ---------------------------------------------------------------------------
test_that("fast() iClassRMSD = 0 when ic.num == k", {
  out <- run_fast_e2e("iClass", ic.num = 2L)
  expect_equal(max(out$iClassRMSD), 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 9. iClass sign patterns: E1 and E3 -> "pp", E2 and E4 -> "pn"
# ---------------------------------------------------------------------------
test_that("fast() iclass labels match loading sign patterns", {
  out <- run_fast_e2e("iClass")
  expect_equal(as.character(out$iclass[out$Site == "E1"][1L]), "pp")
  expect_equal(as.character(out$iclass[out$Site == "E2"][1L]), "pn")
})

# ---------------------------------------------------------------------------
# 10. Error: term without fa() component
# ---------------------------------------------------------------------------
test_that("fast() errors if term has no fa() component", {
  local_mocked_bindings(.fa_asreml = function(...) .mock_sfa, .package = "biomAid")
  expect_error(fast(make_fast_model(), term = "Site:Genotype"),
               "must contain a 'fa\\(")
})

# ---------------------------------------------------------------------------
# 11. Error: ic.num out of range
# ---------------------------------------------------------------------------
test_that("fast() errors if ic.num > k", {
  local_mocked_bindings(.fa_asreml = function(...) .mock_sfa, .package = "biomAid")
  expect_error(fast(make_fast_model(), term = .term, type = "iClass", ic.num = 5L),
               "ic.num.*must be between")
})
