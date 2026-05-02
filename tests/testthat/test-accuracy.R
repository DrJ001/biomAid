# tests/testthat/test-accuracy.R
# Comprehensive tests for accuracy() and its private helpers.
#
# Sections:
#   A. Private helper unit tests  — .split_top, .top_colon, .inside,
#                                   .parse_met_formula, .extract_G_diag
#   B. Error conditions
#   C. End-to-end — FA structure (group-level)
#   D. End-to-end — diag structure (group-level)
#   E. End-to-end — single-environment id() structure
#   F. by_variety = TRUE path
#   G. present = TRUE vs FALSE
#   H. metric argument — column presence / absence
#   I. Edge cases
#
# Mocking strategy:
#   summary()  -> .package = "base"    (not imported explicitly)
#   predict()  -> .package = "biomAid" (importFrom(stats, predict) in NAMESPACE)

# ===========================================================================
# Shared fixtures
# ===========================================================================

make_acc_model <- function(formula_str, data = NULL) {
  frm <- stats::as.formula(formula_str)
  # Store the actual data frame (not a quoted symbol) so .get_model_data()
  # can retrieve it via eval(model$call$data) without needing the symbol in scope.
  m   <- list(
    call     = list(random = frm, data = data),
    formulae = list(random = frm)
  )
  class(m) <- "asreml"
  m
}

make_varcomp_fa <- function(grp_names, nfa = 2L, seed = 1L) {
  set.seed(seed)
  rows <- list()
  for (g in grp_names)
    rows[[paste0("!", g, "!var")]] <- abs(rnorm(1, 0.5, 0.1))
  for (f in seq_len(nfa))
    for (g in grp_names)
      rows[[paste0("!", g, "!fa", f)]] <- rnorm(1, 1, 0.3)
  rows[["units!R"]] <- 1.0
  data.frame(component = unlist(rows), row.names = names(rows),
             stringsAsFactors = FALSE)
}

make_varcomp_diag <- function(grp_names, seed = 2L) {
  set.seed(seed)
  rows <- list()
  for (g in grp_names)
    rows[[paste0("Env:Variety!Env_", g)]] <- abs(rnorm(1, 0.8, 0.2))
  rows[["units!R"]] <- 1.0
  data.frame(component = unlist(rows), row.names = names(rows),
             stringsAsFactors = FALSE)
}

make_pvals_acc <- function(grp_names, var_names, include_sed = FALSE,
                            seed = 3L) {
  set.seed(seed)
  pv <- expand.grid(Env = factor(grp_names), Variety = factor(var_names),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  pv$predicted.value <- rnorm(nrow(pv))
  pv$std.error       <- runif(nrow(pv), 0.15, 0.50)
  pv$status          <- "Estimable"
  sed <- NULL
  if (include_sed) {
    n   <- nrow(pv)
    raw <- matrix(runif(n * n, 0.1, 0.8), n, n)
    sed <- (raw + t(raw)) / 2
    diag(sed) <- 0
  }
  list(pvals = pv, sed = sed)
}

make_acc_dat <- function(grp_names, var_names) {
  expand.grid(Env = factor(grp_names), Variety = factor(var_names),
              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

# Shared group / variety names used across sections C-I
fa_grps <- c("E1", "E2", "E3", "E4")
fa_vars <- paste0("Var", sprintf("%02d", 1:10))

make_fa_e2e <- function(grps = fa_grps, vars = fa_vars,
                         nfa = 2L, include_sed = TRUE, seed = 10L) {
  dat <- make_acc_dat(grps, vars)
  m   <- make_acc_model("~ fa(Env, 2):id(Variety)", data = dat)
  vc  <- make_varcomp_fa(grps, nfa = nfa, seed = seed)
  pv  <- make_pvals_acc(grps, vars, include_sed = include_sed, seed = seed)
  list(dat = dat, model = m, vc = vc, pv = pv)
}

diag_grps <- c("S1", "S2", "S3")
diag_vars <- paste0("G", sprintf("%02d", 1:8))

make_diag_e2e <- function(grps = diag_grps, vars = diag_vars,
                            include_sed = TRUE, seed = 20L) {
  dat <- make_acc_dat(grps, vars)
  m   <- make_acc_model("~ diag(Env):id(Variety)", data = dat)
  vc  <- make_varcomp_diag(grps, seed = seed)
  pv  <- make_pvals_acc(grps, vars, include_sed = include_sed, seed = seed)
  list(dat = dat, model = m, vc = vc, pv = pv)
}

# ===========================================================================
# SECTION A: Private helper unit tests
# ===========================================================================

test_that(".split_top splits on + at depth 0 only", {
  s   <- "fa(Env,2):id(Var)+at(Site):Rep"
  res <- biomAid:::.split_top(s, "+")
  expect_length(res, 2L)
  expect_equal(res[1L], "fa(Env,2):id(Var)")
  expect_equal(res[2L], "at(Site):Rep")
})

test_that(".split_top does not split inside parentheses", {
  s <- "fa(A+B,2):id(C)"
  expect_length(biomAid:::.split_top(s, "+"), 1L)
})

test_that(".split_top returns whole string when sep absent", {
  s <- "diag(Env):id(Var)"
  expect_equal(biomAid:::.split_top(s, "+"), s)
})

test_that(".top_colon finds colon at depth 0", {
  s   <- "fa(Env,2):id(Var)"
  pos <- biomAid:::.top_colon(s)
  expect_equal(substr(s, pos, pos), ":")
  expect_equal(substr(s, 1L, pos - 1L), "fa(Env,2)")
})

test_that(".top_colon returns NA when no top-level colon", {
  expect_true(is.na(biomAid:::.top_colon("diag(Env)")))
})

test_that(".top_colon ignores colons inside parens", {
  s   <- "fa(E:nv,2):id(Var)"
  pos <- biomAid:::.top_colon(s)
  expect_equal(substr(s, pos, pos), ":")
  expect_gt(pos, regexpr(")", s, fixed = TRUE)[1L])
})

test_that(".inside extracts content inside outer parens", {
  expect_equal(biomAid:::.inside("id(Variety)"), "Variety")
  expect_equal(biomAid:::.inside("fa(Env, 2)"),  "Env, 2")
  expect_equal(biomAid:::.inside("diag(Site)"),  "Site")
})

test_that(".inside returns string unchanged when no parens", {
  expect_equal(biomAid:::.inside("Variety"), "Variety")
})

test_that(".parse_met_formula: FA formula sets correct fields", {
  m <- make_acc_model("~ fa(Env, 2):id(Variety) + at(Site):Rep")
  p <- biomAid:::.parse_met_formula(m)
  expect_equal(p$group_var,  "Env")
  expect_equal(p$by_var,     "Variety")
  expect_equal(p$n_fa,       2L)
  expect_equal(p$only_term,  "fa(Env, 2):Variety")
  expect_equal(p$classify,   "Env:Variety")
})

test_that(".parse_met_formula: diag formula has NULL n_fa", {
  m <- make_acc_model("~ diag(Env):id(Variety)")
  p <- biomAid:::.parse_met_formula(m)
  expect_equal(p$group_var, "Env")
  expect_equal(p$by_var,    "Variety")
  expect_null(p$n_fa)
  expect_equal(p$only_term, "Env:Variety")
})

test_that(".parse_met_formula: single-env id() has NULL group_var", {
  m <- make_acc_model("~ id(Variety)")
  p <- biomAid:::.parse_met_formula(m)
  expect_null(p$group_var)
  expect_equal(p$by_var,    "Variety")
  expect_equal(p$only_term, "Variety")
  expect_equal(p$classify,  "Variety")
})

test_that(".parse_met_formula: corgh formula detected correctly", {
  m <- make_acc_model("~ corgh(Env):id(Variety)")
  p <- biomAid:::.parse_met_formula(m)
  expect_equal(p$group_var, "Env")
  expect_equal(p$only_term, "Env:Variety")
})

test_that(".extract_G_diag: FA structure returns named positive vector", {
  grps <- c("E1", "E2", "E3")
  vc   <- make_varcomp_fa(grps, nfa = 2L)
  m    <- make_acc_model("~ fa(Env, 2):id(Variety)")
  local_mocked_bindings(
    summary = function(model, ...) list(varcomp = vc),
    .package = "base"
  )
  G <- biomAid:::.extract_G_diag(m, grps)
  expect_named(G, grps)
  expect_true(all(G > 0))
})

test_that(".extract_G_diag: diag structure returns named positive vector", {
  grps <- c("E1", "E2", "E3")
  vc   <- make_varcomp_diag(grps)
  m    <- make_acc_model("~ diag(Env):id(Variety)")
  local_mocked_bindings(
    summary = function(model, ...) list(varcomp = vc),
    .package = "base"
  )
  G <- biomAid:::.extract_G_diag(m, grps)
  expect_named(G, grps)
  expect_true(all(G > 0))
})

test_that(".extract_G_diag: no matching rows returns NULL", {
  vc <- data.frame(component = 1.0, row.names = "units!R",
                   stringsAsFactors = FALSE)
  m  <- make_acc_model("~ diag(Env):id(Variety)")
  local_mocked_bindings(
    summary = function(model, ...) list(varcomp = vc),
    .package = "base"
  )
  expect_null(biomAid:::.extract_G_diag(m, c("E1", "E2")))
})

test_that(".extract_G_diag: FA G_jj = sum(lambda^2) + psi", {
  grps <- c("E1", "E2")
  psi  <- c(E1 = 0.3,  E2 = 0.4)
  lam  <- c(E1 = 1.2,  E2 = 0.8)
  rows <- list()
  for (g in grps) rows[[paste0("!", g, "!var")]] <- psi[g]
  for (g in grps) rows[[paste0("!", g, "!fa1")]] <- lam[g]
  vc <- data.frame(component = unlist(rows), row.names = names(rows),
                   stringsAsFactors = FALSE)
  m  <- make_acc_model("~ fa(Env, 1):id(Variety)")
  local_mocked_bindings(
    summary = function(model, ...) list(varcomp = vc),
    .package = "base"
  )
  G <- biomAid:::.extract_G_diag(m, grps)
  expect_equal(unname(G["E1"]), unname(lam["E1"]^2 + psi["E1"]), tolerance = 1e-10)
  expect_equal(unname(G["E2"]), unname(lam["E2"]^2 + psi["E2"]), tolerance = 1e-10)
})

# ===========================================================================
# SECTION B: Error conditions
# ===========================================================================

test_that("accuracy() errors when model has no data", {
  m <- make_acc_model("~ diag(Env):id(Variety)")
  # dat doesn't exist in this scope so .get_model_data fails
  expect_error(accuracy(m), "Cannot retrieve data")
})

test_that("accuracy() errors when G extraction returns NULL", {
  grps <- c("E1", "E2")
  vars <- paste0("V", 1:5)
  dat  <- make_acc_dat(grps, vars)
  m    <- make_acc_model("~ diag(Env):id(Variety)", data = dat)
  vc_empty <- data.frame(component = 1.0, row.names = "units!R",
                         stringsAsFactors = FALSE)
  pv <- make_pvals_acc(grps, vars)
  local_mocked_bindings(
    summary = function(model, ...) list(varcomp = vc_empty),
    .package = "base"
  )
  local_mocked_bindings(
    predict = function(...) pv,
    .package = "biomAid"
  )
  expect_error(accuracy(m), "Cannot extract G_jj")
})

test_that("metric match.arg rejects unknown string", {
  m <- make_acc_model("~ fa(Env, 2):id(Variety)")
  expect_error(accuracy(m, metric = "cullis"), "should be one of")
})

# ===========================================================================
# SECTION C: End-to-end — FA structure, group-level
# ===========================================================================

test_that("FA group-level: returns data.frame with correct columns (both metrics)", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = c("accuracy", "gen.H2"))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("group","n_vars","G_jj","mean_acc","sd_acc","gen.H2")
                  %in% names(res)))
})

test_that("FA group-level: one row per group", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  expect_equal(nrow(accuracy(e$model)), length(fa_grps))
})

test_that("FA group-level: group column matches group names", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  expect_equal(sort(accuracy(e$model)$group), sort(fa_grps))
})

test_that("FA group-level: mean_acc in [0, 1]", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "accuracy")
  expect_true(all(res$mean_acc >= 0 & res$mean_acc <= 1, na.rm = TRUE))
})

test_that("FA group-level: gen.H2 in [0, 1]", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "gen.H2")
  expect_true(all(res$gen.H2 >= 0 & res$gen.H2 <= 1, na.rm = TRUE))
})

test_that("FA group-level: G_jj is positive", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  expect_true(all(accuracy(e$model)$G_jj > 0, na.rm = TRUE))
})

test_that("FA group-level: sd_acc is non-negative", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "accuracy")
  expect_true(all(res$sd_acc >= 0, na.rm = TRUE))
})

test_that("FA group-level: n_vars equals number of varieties per group", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "accuracy")
  expect_true(all(res$n_vars == length(fa_vars)))
})

test_that("FA group-level: accuracy only — no gen.H2 column", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "accuracy")
  expect_false("gen.H2"   %in% names(res))
  expect_true("mean_acc"  %in% names(res))
})

test_that("FA group-level: gen.H2 only — no mean_acc or sd_acc columns", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "gen.H2")
  expect_false("mean_acc" %in% names(res))
  expect_false("sd_acc"   %in% names(res))
  expect_true("gen.H2"    %in% names(res))
})

test_that("FA group-level: classify override accepted without error", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  expect_no_error(accuracy(e$model, classify = "Env:Variety",
                            metric = "accuracy"))
})

# ===========================================================================
# SECTION D: End-to-end — diag structure
# ===========================================================================

test_that("diag group-level: returns data.frame with correct structure", {
  e <- make_diag_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), length(diag_grps))
  expect_true(all(c("group","n_vars","G_jj","mean_acc","sd_acc","gen.H2")
                  %in% names(res)))
})

test_that("diag group-level: G_jj matches varcomp component", {
  e <- make_diag_e2e(seed = 21L)
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "accuracy")
  for (g in diag_grps) {
    vc_row <- e$vc[endsWith(rownames(e$vc), paste0("_", g)), , drop = FALSE]
    if (nrow(vc_row) == 1L)
      expect_equal(res$G_jj[res$group == g], vc_row$component,
                   tolerance = 1e-10)
  }
})

test_that("diag group-level: mean_acc in [0, 1]", {
  e <- make_diag_e2e(seed = 22L)
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "accuracy")
  expect_true(all(res$mean_acc >= 0 & res$mean_acc <= 1, na.rm = TRUE))
})

# ===========================================================================
# SECTION E: End-to-end — single-environment id() structure
# ===========================================================================

make_single_env_e2e <- function(vars = paste0("V", 1:12), seed = 30L) {
  dat <- data.frame(Variety = factor(vars), stringsAsFactors = FALSE)
  m   <- make_acc_model("~ id(Variety)", data = dat)
  set.seed(seed)
  vc  <- data.frame(component = c(abs(rnorm(1, 1.0, 0.2)), 1.0),
                    row.names  = c("Variety", "units!R"),
                    stringsAsFactors = FALSE)
  set.seed(seed)
  pv_df <- data.frame(Variety         = factor(vars),
                      predicted.value = rnorm(length(vars)),
                      std.error       = runif(length(vars), 0.15, 0.40),
                      status          = "Estimable",
                      stringsAsFactors = FALSE)
  n   <- nrow(pv_df)
  raw <- matrix(runif(n * n, 0.1, 0.6), n, n)
  sed <- (raw + t(raw)) / 2; diag(sed) <- 0
  list(dat = dat, model = m, vc = vc, pv = list(pvals = pv_df, sed = sed))
}

test_that("single-env: returns one row labelled 'all'", {
  e <- make_single_env_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model)
  expect_equal(nrow(res), 1L)
  expect_equal(res$group, "all")
})

test_that("single-env: mean_acc in [0, 1]", {
  e <- make_single_env_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "accuracy")
  expect_true(res$mean_acc >= 0 && res$mean_acc <= 1)
})

test_that("single-env: gen.H2 in [0, 1]", {
  e <- make_single_env_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "gen.H2")
  expect_true(res$gen.H2 >= 0 && res$gen.H2 <= 1)
})

# ===========================================================================
# SECTION F: by_variety = TRUE
# ===========================================================================

test_that("by_variety: returns one row per variety x group", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, by_variety = TRUE, metric = "accuracy")
  expect_equal(nrow(res), length(fa_grps) * length(fa_vars))
})

test_that("by_variety: has group, variety, G_jj, pev, accuracy columns", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, by_variety = TRUE, metric = "accuracy")
  expect_true(all(c("group","variety","G_jj","pev","accuracy") %in% names(res)))
})

test_that("by_variety: no n_vars or sd_acc columns", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, by_variety = TRUE, metric = "accuracy")
  expect_false("n_vars" %in% names(res))
  expect_false("sd_acc" %in% names(res))
})

test_that("by_variety: accuracy values in [0, 1]", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, by_variety = TRUE, metric = "accuracy")
  expect_true(all(res$accuracy >= 0 & res$accuracy <= 1, na.rm = TRUE))
})

test_that("by_variety: pev values are positive", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, by_variety = TRUE, metric = "accuracy")
  expect_true(all(res$pev > 0, na.rm = TRUE))
})

test_that("by_variety + gen.H2: gen.H2 is constant within each group", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, by_variety = TRUE,
                  metric = c("accuracy", "gen.H2"))
  for (g in unique(res$group))
    expect_equal(length(unique(res$gen.H2[res$group == g])), 1L)
})

test_that("by_variety: variety column is character", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, by_variety = TRUE, metric = "accuracy")
  expect_type(res$variety, "character")
})

# ===========================================================================
# SECTION G: present = TRUE vs FALSE
# ===========================================================================

test_that("present=FALSE returns >= rows as present=TRUE (by_variety)", {
  grps <- c("E1", "E2", "E3")
  vars <- paste0("Var", 1:8)

  # Build an unbalanced observed data frame (only a subset of combos)
  obs_dat <- make_acc_dat(grps[1:2], vars[1:5])   # fewer combos than full grid
  m_obs   <- make_acc_model("~ fa(Env, 2):id(Variety)", data = obs_dat)
  vc      <- make_varcomp_fa(grps, nfa = 2L, seed = 40L)

  # pvals covers full grid; some std.error = 0 (mark as unobserved)
  pv <- make_pvals_acc(grps, vars, include_sed = FALSE, seed = 40L)
  pv$pvals$std.error[pv$pvals$Env == "E3"] <- 0   # E3 all invalid

  local_mocked_bindings(summary = function(model, ...) list(varcomp = vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) pv, .package = "biomAid")
  res_T <- accuracy(m_obs, by_variety = TRUE, metric = "accuracy",
                    present = TRUE)
  res_F <- accuracy(m_obs, by_variety = TRUE, metric = "accuracy",
                    present = FALSE)
  expect_lte(nrow(res_T), nrow(res_F))
})

# ===========================================================================
# SECTION H: metric argument — column presence / absence
# ===========================================================================

test_that("default metric returns both mean_acc and gen.H2 columns", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model)
  expect_true("mean_acc" %in% names(res))
  expect_true("gen.H2"   %in% names(res))
})

test_that("metric='accuracy' — mean_acc present, gen.H2 absent", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "accuracy")
  expect_true("mean_acc" %in% names(res))
  expect_false("gen.H2"  %in% names(res))
})

test_that("metric='gen.H2' — gen.H2 present, mean_acc and sd_acc absent", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  res <- accuracy(e$model, metric = "gen.H2")
  expect_true("gen.H2"    %in% names(res))
  expect_false("mean_acc" %in% names(res))
  expect_false("sd_acc"   %in% names(res))
})

test_that("group and G_jj columns always present for any metric", {
  e <- make_fa_e2e()
  for (m_arg in list("accuracy", "gen.H2", c("accuracy", "gen.H2"))) {
    local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                          .package = "base")
    local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
    res <- accuracy(e$model, metric = m_arg)
    expect_true("group" %in% names(res))
    expect_true("G_jj"  %in% names(res))
  }
})

# ===========================================================================
# SECTION I: Edge cases
# ===========================================================================

test_that("group with < 2 valid rows returns NA accuracy and 0L n_vars", {
  grps <- c("E1", "E2")
  vars <- paste0("V", 1:5)
  dat  <- make_acc_dat(grps, vars)
  m    <- make_acc_model("~ diag(Env):id(Variety)", data = dat)
  vc   <- make_varcomp_diag(grps, seed = 50L)
  set.seed(50L)
  pv <- expand.grid(Env = factor(grps), Variety = factor(vars),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  pv$predicted.value <- rnorm(nrow(pv))
  pv$std.error <- runif(nrow(pv), 0.1, 0.4)
  pv$status    <- "Estimable"
  # Only 1 valid (std.error > 0) row for E2; use present=FALSE so validity
  # comes from std.error filter rather than data lookup.
  pv$std.error[pv$Env == "E2"][-1L] <- 0

  local_mocked_bindings(summary = function(model, ...) list(varcomp = vc),
                        .package = "base")
  local_mocked_bindings(
    predict = function(...) list(pvals = pv, sed = NULL),
    .package = "biomAid"
  )
  res    <- accuracy(m, metric = "accuracy", present = FALSE)
  e2_row <- res[res$group == "E2", ]
  expect_equal(e2_row$n_vars, 0L)
  expect_true(is.na(e2_row$mean_acc))
})

test_that("all groups returned as data.frame even when one group is empty", {
  grps <- c("E1", "E2", "E3")
  vars <- paste0("V", 1:4)
  dat  <- make_acc_dat(grps, vars)
  m    <- make_acc_model("~ diag(Env):id(Variety)", data = dat)
  vc   <- make_varcomp_diag(grps, seed = 51L)
  pv   <- make_pvals_acc(grps, vars, include_sed = FALSE, seed = 51L)
  pv$pvals$std.error[pv$pvals$Env == "E3"] <- 0   # E3 all invalid

  local_mocked_bindings(summary = function(model, ...) list(varcomp = vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) pv, .package = "biomAid")
  res <- accuracy(m, metric = "accuracy")
  expect_equal(nrow(res), 3L)
  expect_equal(sort(res$group), sort(grps))
})

test_that("pworkspace argument accepted without error", {
  e <- make_fa_e2e()
  local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                        .package = "base")
  local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
  expect_no_error(
    accuracy(e$model, metric = "accuracy", pworkspace = "4gb")
  )
})

test_that("results are reproducible — same inputs give identical output", {
  e <- make_fa_e2e(seed = 60L)
  results <- lapply(1:2, function(i) {
    local_mocked_bindings(summary = function(model, ...) list(varcomp = e$vc),
                          .package = "base")
    local_mocked_bindings(predict = function(...) e$pv, .package = "biomAid")
    accuracy(e$model, metric = "accuracy")
  })
  expect_equal(results[[1L]], results[[2L]])
})
