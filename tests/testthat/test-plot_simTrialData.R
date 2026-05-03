# =============================================================================
# test-plot_simTrialData.R
# =============================================================================
#
# Tests for plot_simTrialData() and the simTrialData() change that now always
# returns params$g_arr.
#
# Shared simulation objects (created once per session for speed):
#   sim_bal   -- balanced MET-only, 10 varieties, 3 sites
#   sim_unbal -- unbalanced MET-only, 15 varieties, 5 sites
#   sim_multi -- multi-treatment, 10 varieties, 3 sites, 2 treatments
# =============================================================================

sim_bal   <- simTrialData(nvar = 10, nsite = 3, seed = 1,  verbose = FALSE)
sim_unbal <- simTrialData(nvar = 15, nsite = 5,
                          incidence = "unbalanced",
                          seed = 2, verbose = FALSE)
sim_multi <- simTrialData(nvar = 10, nsite = 3,
                          treatments = c("T0", "T1"),
                          seed = 3, verbose = FALSE)


# =============================================================================
# simTrialData() — g_arr always returned
# =============================================================================

test_that("simTrialData always returns g_arr for MET-only", {
  expect_true("g_arr" %in% names(sim_bal$params))
})

test_that("simTrialData g_arr has correct dimensions for MET-only", {
  g <- sim_bal$params$g_arr
  expect_equal(nrow(g), 10L)
  expect_equal(ncol(g), 3L)
})

test_that("simTrialData g_arr colnames match site names for MET-only", {
  g     <- sim_bal$params$g_arr
  sites <- levels(sim_bal$data$Site)
  expect_equal(colnames(g), sites)
})

test_that("simTrialData always returns g_arr for multi-treatment", {
  expect_true("g_arr" %in% names(sim_multi$params))
})

test_that("simTrialData g_arr dimensions correct for multi-treatment", {
  g <- sim_multi$params$g_arr
  expect_equal(nrow(g), 10L)
  expect_equal(ncol(g), 6L)   # 2 treatments x 3 sites
})

test_that("simTrialData g_arr colnames are TSite labels for multi-treatment", {
  g       <- sim_multi$params$g_arr
  tsites  <- levels(sim_multi$data$TSite)
  expect_equal(sort(colnames(g)), sort(tsites))
})

test_that("simTrialData g_arr rownames are variety labels", {
  g <- sim_bal$params$g_arr
  expect_equal(rownames(g), levels(sim_bal$data$Variety))
})

test_that("simTrialData g_arr is numeric", {
  expect_true(is.numeric(sim_bal$params$g_arr))
})


# =============================================================================
# plot_simTrialData() — input validation
# =============================================================================

test_that("error on non-list res", {
  expect_error(plot_simTrialData(42),
               "'res' must be the list returned by simTrialData()")
})

test_that("error on res missing data element", {
  bad <- list(params = sim_bal$params)
  expect_error(plot_simTrialData(bad),
               "'res' must be the list returned by simTrialData()")
})

test_that("error on res missing params element", {
  bad <- list(data = sim_bal$data)
  expect_error(plot_simTrialData(bad),
               "'res' must be the list returned by simTrialData()")
})

test_that("error when params missing g_arr", {
  bad <- sim_bal
  bad$params$g_arr <- NULL
  expect_error(plot_simTrialData(bad), "missing")
})

test_that("error when params missing G", {
  bad <- sim_bal
  bad$params$G <- NULL
  expect_error(plot_simTrialData(bad), "missing")
})

test_that("error when params missing incidence", {
  bad <- sim_bal
  bad$params$incidence <- NULL
  expect_error(plot_simTrialData(bad), "missing")
})

test_that("error on non-logical sort", {
  expect_error(plot_simTrialData(sim_bal, type = "incidence", sort = "yes"),
               "'sort' must be a single logical")
})

test_that("error on non-logical return_data", {
  expect_error(plot_simTrialData(sim_bal, return_data = "yes"),
               "'return_data' must be a single logical")
})

test_that("invalid type is caught by match.arg", {
  expect_error(plot_simTrialData(sim_bal, type = "banana"))
})

test_that("fill ignored with warning for non-trial types", {
  expect_warning(plot_simTrialData(sim_bal, type = "incidence", fill = "yield"),
                 "fill.*ignored")
})

test_that("label ignored with warning for non-trial types", {
  expect_warning(plot_simTrialData(sim_bal, type = "blup", label = "Variety"),
                 "label.*ignored")
})

test_that("ncol ignored with warning for non-trial types", {
  expect_warning(plot_simTrialData(sim_bal, type = "correlation", ncol = 2L),
                 "ncol.*ignored")
})

test_that("sites ignored with warning for incidence type", {
  expect_warning(
    plot_simTrialData(sim_bal, type = "incidence",
                     sites = levels(sim_bal$data$Site)[1L]),
    "sites.*ignored"
  )
})


# =============================================================================
# plot_simTrialData() — all four types return ggplot
# =============================================================================

test_that("type = 'trial' returns a ggplot", {
  p <- plot_simTrialData(sim_bal)
  expect_s3_class(p, "ggplot")
})

test_that("type = 'incidence' returns a ggplot", {
  p <- plot_simTrialData(sim_bal, type = "incidence")
  expect_s3_class(p, "ggplot")
})

test_that("type = 'correlation' returns a ggplot", {
  p <- plot_simTrialData(sim_bal, type = "correlation")
  expect_s3_class(p, "ggplot")
})

test_that("type = 'blup' returns a ggplot", {
  p <- plot_simTrialData(sim_bal, type = "blup")
  expect_s3_class(p, "ggplot")
})

test_that("all four types render without error (ggplot_build)", {
  for (tp in c("trial", "incidence", "correlation", "blup")) {
    expect_no_error(ggplot2::ggplot_build(
      plot_simTrialData(sim_bal, type = tp)
    ))
  }
})


# =============================================================================
# return_data = TRUE
# =============================================================================

test_that("return_data=TRUE returns a list with $plot and $data", {
  out <- plot_simTrialData(sim_bal, return_data = TRUE)
  expect_type(out, "list")
  expect_named(out, c("plot", "data"))
  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$data, "data.frame")
})

test_that("return_data=TRUE $data has rows for trial type", {
  out <- plot_simTrialData(sim_bal, return_data = TRUE)
  expect_gt(nrow(out$data), 0L)
})

test_that("return_data=TRUE for incidence has Variety, Site, Present columns", {
  out <- plot_simTrialData(sim_bal, type = "incidence", return_data = TRUE)
  expect_true(all(c("Variety", "Site", "Present") %in% names(out$data)))
})

test_that("return_data=TRUE for correlation has Group1, Group2, Correlation cols", {
  out <- plot_simTrialData(sim_bal, type = "correlation", return_data = TRUE)
  expect_true(all(c("Group1", "Group2", "Correlation") %in% names(out$data)))
})

test_that("return_data=TRUE for blup has Variety, Group, BLUP columns", {
  out <- plot_simTrialData(sim_bal, type = "blup", return_data = TRUE)
  expect_true(all(c("Variety", "Group", "BLUP") %in% names(out$data)))
})


# =============================================================================
# "trial" type — fill argument
# =============================================================================

test_that("trial: fill=NULL defaults to Rep colouring (no error)", {
  expect_no_error(plot_simTrialData(sim_bal, type = "trial", fill = NULL))
})

test_that("trial: fill='yield' produces continuous fill scale", {
  out <- plot_simTrialData(sim_bal, type = "trial",
                           fill = "yield", return_data = TRUE)
  expect_true(is.numeric(out$data$yield))
})

test_that("trial: fill='Variety' produces discrete fill", {
  p <- plot_simTrialData(sim_bal, type = "trial", fill = "Variety")
  expect_s3_class(p, "ggplot")
})

test_that("trial: fill='Rep' works explicitly", {
  p <- plot_simTrialData(sim_bal, type = "trial", fill = "Rep")
  expect_s3_class(p, "ggplot")
})

test_that("trial: fill='Treatment' works for multi-treatment data", {
  p <- plot_simTrialData(sim_multi, type = "trial", fill = "Treatment")
  expect_s3_class(p, "ggplot")
})

test_that("trial: error on unknown fill column", {
  expect_error(plot_simTrialData(sim_bal, fill = "NonExistent"),
               "fill column.*not found")
})


# =============================================================================
# "trial" type — label argument
# =============================================================================

test_that("trial: label='Variety' overlays text without error", {
  p <- plot_simTrialData(sim_bal, type = "trial", label = "Variety")
  expect_s3_class(p, "ggplot")
})

test_that("trial: label='Rep' works", {
  p <- plot_simTrialData(sim_bal, type = "trial", label = "Rep")
  expect_s3_class(p, "ggplot")
})

test_that("trial: label='yield' works (numeric column)", {
  p <- plot_simTrialData(sim_bal, type = "trial",
                         fill = "yield", label = "yield")
  expect_s3_class(p, "ggplot")
})

test_that("trial: error on unknown label column", {
  expect_error(plot_simTrialData(sim_bal, label = "NoSuchCol"),
               "label column.*not found")
})

test_that("trial: fill and label can both be specified", {
  p <- plot_simTrialData(sim_bal, fill = "yield", label = "Variety")
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# "trial" type — sites argument
# =============================================================================

test_that("trial: sites subsets to requested sites", {
  site1 <- levels(sim_bal$data$Site)[1L]
  out   <- plot_simTrialData(sim_bal, sites = site1, return_data = TRUE)
  expect_equal(nlevels(out$data$Site), 1L)
})

test_that("trial: sites produces warning for unknown sites", {
  expect_warning(
    plot_simTrialData(sim_bal, sites = c(levels(sim_bal$data$Site)[1L], "FakeEnv")),
    "not found.*ignored"
  )
})

test_that("trial: error when all sites are unknown", {
  expect_error(
    suppressWarnings(plot_simTrialData(sim_bal, sites = "ZZZ")),
    "No data remaining"
  )
})

test_that("trial: ncol argument is accepted without error", {
  p <- plot_simTrialData(sim_bal, ncol = 2L)
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# "incidence" type
# =============================================================================

test_that("incidence: data has nvar x nsite rows", {
  out <- plot_simTrialData(sim_bal, type = "incidence", return_data = TRUE)
  expected_rows <- nrow(sim_bal$params$incidence) * ncol(sim_bal$params$incidence)
  expect_equal(nrow(out$data), expected_rows)
})

test_that("incidence: Present column is 0 or 1", {
  out <- plot_simTrialData(sim_bal, type = "incidence", return_data = TRUE)
  expect_true(all(out$data$Present %in% c(0L, 1L)))
})

test_that("incidence: sort=TRUE reorders variety levels by site count", {
  out      <- plot_simTrialData(sim_unbal, type = "incidence",
                                sort = TRUE, return_data = TRUE)
  inc      <- sim_unbal$params$incidence
  var_levs <- levels(out$data$Variety)         # bottom-to-top
  spc      <- rowSums(inc)[var_levs]
  expect_true(all(diff(spc) >= 0))             # increasing bottom-to-top
})

test_that("incidence: sort=FALSE does not reorder by count", {
  out_sorted   <- plot_simTrialData(sim_unbal, type = "incidence",
                                    sort = TRUE,  return_data = TRUE)
  out_unsorted <- plot_simTrialData(sim_unbal, type = "incidence",
                                    sort = FALSE, return_data = TRUE)
  expect_false(identical(levels(out_sorted$data$Variety),
                         levels(out_unsorted$data$Variety)))
})

test_that("incidence: works for unbalanced design", {
  p <- plot_simTrialData(sim_unbal, type = "incidence")
  expect_s3_class(p, "ggplot")
})

test_that("incidence: balanced design has all-1 Present", {
  out <- plot_simTrialData(sim_bal, type = "incidence", return_data = TRUE)
  expect_true(all(out$data$Present == 1L))
})


# =============================================================================
# "correlation" type
# =============================================================================

test_that("correlation: data has ngroup^2 rows", {
  out   <- plot_simTrialData(sim_bal, type = "correlation", return_data = TRUE)
  nG    <- nrow(sim_bal$params$G)
  expect_equal(nrow(out$data), nG^2L)
})

test_that("correlation: Correlation values in [-1, 1]", {
  out <- plot_simTrialData(sim_bal, type = "correlation", return_data = TRUE)
  expect_true(all(out$data$Correlation >= -1 - 1e-9))
  expect_true(all(out$data$Correlation <=  1 + 1e-9))
})

test_that("correlation: diagonal values are 1", {
  out <- plot_simTrialData(sim_bal, type = "correlation", return_data = TRUE)
  diag_rows <- as.character(out$data$Group1) == as.character(out$data$Group2)
  expect_true(all(abs(out$data$Correlation[diag_rows] - 1) < 1e-9))
})

test_that("correlation: works for multi-treatment (larger G)", {
  p <- plot_simTrialData(sim_multi, type = "correlation")
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# "blup" type
# =============================================================================

test_that("blup: data has nvar x ngroup rows", {
  out <- plot_simTrialData(sim_bal, type = "blup", return_data = TRUE)
  nv  <- nrow(sim_bal$params$g_arr)
  ng  <- ncol(sim_bal$params$g_arr)
  expect_equal(nrow(out$data), nv * ng)
})

test_that("blup: sort=TRUE orders varieties by increasing mean BLUP (bottom to top)", {
  out      <- plot_simTrialData(sim_bal, type = "blup",
                                sort = TRUE, return_data = TRUE)
  g        <- sim_bal$params$g_arr
  var_levs <- levels(out$data$Variety)          # bottom-to-top
  means    <- rowMeans(g)[var_levs]
  expect_true(all(diff(means) >= -1e-9))        # non-decreasing
})

test_that("blup: sort=FALSE produces different ordering than sort=TRUE", {
  out_s <- plot_simTrialData(sim_bal, type = "blup",
                             sort = TRUE,  return_data = TRUE)
  out_u <- plot_simTrialData(sim_bal, type = "blup",
                             sort = FALSE, return_data = TRUE)
  # Levels will differ when there is any variance in mean BLUPs
  expect_false(identical(levels(out_s$data$Variety),
                         levels(out_u$data$Variety)))
})

test_that("blup: sites subsetting works for MET-only", {
  site1 <- colnames(sim_bal$params$g_arr)[1L]
  out   <- plot_simTrialData(sim_bal, type = "blup",
                             sites = site1, return_data = TRUE)
  expect_equal(nlevels(out$data$Group), 1L)
})

test_that("blup: sites substring matching works for multi-treatment", {
  site1 <- levels(sim_multi$data$Site)[1L]
  expect_message(
    out <- plot_simTrialData(sim_multi, type = "blup",
                             sites = site1, return_data = TRUE),
    "matched as substrings"
  )
  # Should have ntreat groups for that site
  expect_equal(nlevels(out$data$Group), 2L)
})

test_that("blup: error when no sites match", {
  expect_error(
    plot_simTrialData(sim_bal, type = "blup", sites = "ZZZnosuchsite"),
    "None of the 'sites' values matched"
  )
})

test_that("blup: works for multi-treatment (TSite groups)", {
  p <- plot_simTrialData(sim_multi, type = "blup")
  expect_s3_class(p, "ggplot")
})

test_that("blup: renders without error for unbalanced design", {
  p <- plot_simTrialData(sim_unbal, type = "blup")
  expect_s3_class(p, "ggplot")
})
