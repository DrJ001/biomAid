# tests/testthat/test-plot_randomRegress.R
# Tests for plot_randomRegress() — no ASReml needed.

library(ggplot2)

# ---- Shared mock result ------------------------------------------------
#
# Shaped like a real randomRegress() result: the first two columns of $blups
# are the section label and the unit factor, and $section / $unit name them.
# `legacy = TRUE` drops those two elements and uses the pre-rename column
# names, which is how the positional fallback is exercised.
make_mock_res <- function(ns = 3L, nvar = 8L,
                          levs = c("T0", "T1", "T2"),
                          ucol = "Variety", legacy = FALSE) {
  set.seed(42)
  sites <- paste0("E", seq_len(ns))
  scol  <- if (legacy) "Site" else "Section"

  blups <- do.call(rbind, lapply(sites, function(s) {
    u0 <- rnorm(nvar, 0, 200)
    u1 <- 0.9  * u0 + rnorm(nvar, 0, 80)
    u2 <- 1.2  * u0 + rnorm(nvar, 0, 120)
    d <- data.frame(
      section = s,
      unit    = paste0("V", seq_len(nvar)),
      T0 = u0, T1 = u1, T2 = u2,
      adj.T1 = u1 - 0.9  * u0,
      adj.T2 = u2 - 1.2  * u0,
      HSD.T1  = NA_real_,
      HSD.T2  = NA_real_,
      stringsAsFactors = FALSE
    )
    names(d)[1:2] <- c(scol, ucol)
    d
  }))

  beta <- list(
    T1 = matrix(runif(ns, 0.7, 1.1), ns, 1L,
                dimnames = list(sites, "T0")),
    T2 = matrix(runif(ns, 1.0, 1.4), ns, 1L,
                dimnames = list(sites, "T0"))
  )

  nts  <- length(levs) * ns
  nms  <- paste0(rep(levs, ns), "-", rep(sites, each = length(levs)))
  G    <- crossprod(matrix(rnorm(nts * nts), nts, nts))
  dimnames(G) <- list(nms, nms)

  out <- list(
    sep       = "-",
    blups     = blups,
    beta      = beta,
    Gmat      = G,
    TGmat     = G,
    sigmat    = matrix(runif(ns * 2L, 100, 400), ns, 2L,
                       dimnames = list(sites, c("T1", "T2"))),
    tmat      = diag(nts),
    cond_list = list(T0 = NULL, T1 = "T0", T2 = "T0"),
    type      = "baseline"
  )
  if (!legacy) { out$section <- scol; out$unit <- ucol }
  out
}

# ---- Input validation --------------------------------------------------

test_that("invalid type gives informative error", {
  expect_error(plot_randomRegress(make_mock_res(), type = "circles"),
               regexp = "arg")
})

test_that("missing res elements give informative error", {
  bad <- make_mock_res()
  bad$Gmat <- NULL
  expect_error(plot_randomRegress(bad), regexp = "Gmat")
})

test_that("non-theme theme argument gives informative error", {
  expect_error(
    plot_randomRegress(make_mock_res(), theme = "bw"),
    regexp = "theme"
  )
})

test_that("treatments with no match gives informative error", {
  expect_error(
    plot_randomRegress(make_mock_res(), type = "regress",
                       treatments = "N99"),
    regexp = "No conditioned treatments"
  )
})

# ---- Return types -------------------------------------------------------

test_that("all three types return a ggplot object", {
  res <- make_mock_res()
  for (tp in c("regress", "quadrant", "gmat"))
    expect_s3_class(plot_randomRegress(res, type = tp), "ggplot")
})

test_that("highlight = NULL suppresses highlighting without error", {
  res <- make_mock_res()
  for (tp in c("regress", "quadrant"))
    expect_s3_class(
      plot_randomRegress(res, type = tp, highlight = NULL), "ggplot"
    )
})

test_that("user-specified highlight returns ggplot", {
  res  <- make_mock_res()
  vars <- unique(res$blups$Variety)[1:3]
  expect_s3_class(
    plot_randomRegress(res, type = "quadrant", highlight = vars),
    "ggplot"
  )
})

test_that("non-existent highlight variety raises a warning", {
  res <- make_mock_res()
  expect_warning(
    plot_randomRegress(res, type = "quadrant",
                       highlight = c("V1", "DoesNotExist")),
    regexp = "not found"
  )
})

test_that("default highlights return data frame with Variety and group cols", {
  res   <- make_mock_res(ns = 4L, nvar = 12L)
  qdata <- plot_randomRegress(res, type = "quadrant", return_data = TRUE)
  hl    <- biomAid:::.rreg_default_highlights(qdata)
  expect_s3_class(hl, "data.frame")
  expect_true(all(c("Variety", "group") %in% names(hl)))
  expect_true(all(hl$group %in% c("tr", "bl")))
})

test_that("default highlights contain at most 6 varieties", {
  res   <- make_mock_res(ns = 4L, nvar = 12L)
  qdata <- plot_randomRegress(res, type = "quadrant", return_data = TRUE)
  hl    <- biomAid:::.rreg_default_highlights(qdata)
  expect_lte(nrow(hl), 6L)
  expect_gte(nrow(hl), 1L)
})

test_that("return_data = TRUE returns a data.frame for all types", {
  res <- make_mock_res()
  for (tp in c("regress", "quadrant", "gmat")) {
    df <- plot_randomRegress(res, type = tp, return_data = TRUE)
    expect_s3_class(df, "data.frame")
    expect_gt(nrow(df), 0L)
  }
})

# ---- Data frame column checks ------------------------------------------

test_that("regress return_data has required columns", {
  df <- plot_randomRegress(make_mock_res(), type = "regress",
                           return_data = TRUE)
  expect_true(all(c("Section", "Variety", "x", "y",
                    "pair_label", "beta") %in% names(df)))
})

test_that("quadrant return_data has required columns", {
  df <- plot_randomRegress(make_mock_res(), type = "quadrant",
                           return_data = TRUE)
  expect_true(all(c("Section", "Variety", "x", "y",
                    "pair_label") %in% names(df)))
})

test_that("gmat return_data has row_var, col_var, corr", {
  df <- plot_randomRegress(make_mock_res(), type = "gmat",
                           return_data = TRUE)
  expect_true(all(c("row_var", "col_var", "corr") %in% names(df)))
})

# ---- Correctness checks ------------------------------------------------

test_that("gmat correlations are within [-1, 1]", {
  df <- plot_randomRegress(make_mock_res(), type = "gmat",
                           return_data = TRUE)
  expect_true(all(df$corr >= -1 - 1e-9 & df$corr <= 1 + 1e-9))
})

test_that("gmat diagonal correlations equal 1", {
  df <- plot_randomRegress(make_mock_res(), type = "gmat",
                           return_data = TRUE)
  diag_rows <- as.character(df$row_var) == as.character(df$col_var)
  expect_equal(df$corr[diag_rows], rep(1, sum(diag_rows)), tolerance = 1e-9)
})

test_that("regress data has one row per variety x site x conditioned treatment", {
  res <- make_mock_res(ns = 3L, nvar = 8L)
  df  <- plot_randomRegress(res, type = "regress", return_data = TRUE)
  # 2 conditioned treatments x 3 sites x 8 varieties = 48
  expect_equal(nrow(df), 2L * 3L * 8L)
})

test_that("treatments filter reduces rows in regress output", {
  res  <- make_mock_res(ns = 3L, nvar = 8L)
  df_all <- plot_randomRegress(res, type = "regress",
                               return_data = TRUE)
  df_sub <- plot_randomRegress(res, type = "regress",
                               treatments  = "T1",
                               return_data = TRUE)
  expect_lt(nrow(df_sub), nrow(df_all))
  expect_true(all(df_sub$pair_label == "T1 | T0"))
})

test_that("regress return_data contains beta column with finite values", {
  df <- plot_randomRegress(make_mock_res(), type = "regress",
                           return_data = TRUE)
  expect_true("beta" %in% names(df))
  expect_true(all(is.finite(df$beta)))
})

test_that("custom theme is applied without error", {
  res <- make_mock_res()
  expect_s3_class(
    plot_randomRegress(res, type = "quadrant",
                       theme = ggplot2::theme_classic()),
    "ggplot"
  )
})

test_that("centre = TRUE returns ggplot for regress and quadrant", {
  res <- make_mock_res()
  for (tp in c("regress", "quadrant"))
    expect_s3_class(
      plot_randomRegress(res, type = tp, centre = TRUE), "ggplot"
    )
})

test_that("centre = TRUE shifts x values by adding section mean", {
  res    <- make_mock_res(ns = 3L, nvar = 8L)
  df_raw <- plot_randomRegress(res, type = "quadrant",
                               centre = FALSE, return_data = TRUE)
  df_cen <- plot_randomRegress(res, type = "quadrant",
                               centre = TRUE,  return_data = TRUE)
  # After adding the section mean, x should differ from raw by that mean
  sec_shifts <- ave(df_raw$x, df_raw$Section, FUN = mean)
  expect_equal(df_cen$x, df_raw$x + sec_shifts, tolerance = 1e-9)
})

test_that("invalid centre value gives informative error", {
  expect_error(
    plot_randomRegress(make_mock_res(), centre = "yes"),
    regexp = "logical"
  )
})

# ---- cond_x tests -------------------------------------------------------

# Mock result with partial conditioning (all 3 treatments conditioned)
make_mock_res_partial <- function(ns = 2L, nvar = 8L,
                                  levs = c("T0","T1","T2")) {
  set.seed(7L)
  sites <- paste0("E", seq_len(ns))
  blups <- do.call(rbind, lapply(sites, function(s) {
    u0 <- rnorm(nvar, 0, 200); u1 <- 0.8*u0 + rnorm(nvar,0,80)
    u2 <- 0.7*u0 + rnorm(nvar,0,100)
    data.frame(Section=s, Variety=paste0("V",seq_len(nvar)),
               T0=u0, T1=u1, T2=u2,
               adj.T0 = u0 - 0.5*u1 - 0.3*u2,
               adj.T1 = u1 - 0.4*u0 - 0.2*u2,
               adj.T2 = u2 - 0.3*u0 - 0.2*u1,
               HSD.T0=NA_real_, HSD.T1=NA_real_, HSD.T2=NA_real_,
               stringsAsFactors=FALSE)
  }))
  beta <- list(
    T0 = matrix(c(0.5,0.3, 0.5,0.3), ns, 2L,
                dimnames = list(sites, c("T1","T2"))),
    T1 = matrix(c(0.4,0.2, 0.4,0.2), ns, 2L,
                dimnames = list(sites, c("T0","T2"))),
    T2 = matrix(c(0.3,0.2, 0.3,0.2), ns, 2L,
                dimnames = list(sites, c("T0","T1")))
  )
  nts <- length(levs)*ns
  nms <- paste0(rep(levs,ns),"-",rep(sites,each=length(levs)))
  G   <- crossprod(matrix(rnorm(nts*nts),nts,nts)); dimnames(G)<-list(nms,nms)
  list(
    sep       = "-",
    blups     = blups,
    beta      = beta,
    Gmat      = G,
    TGmat     = G,
    sigmat    = matrix(runif(ns*3,10,50), ns, 3L,
                       dimnames=list(sites, levs)),
    tmat      = diag(nts),
    cond_list = list(T0=c("T1","T2"), T1=c("T0","T2"), T2=c("T0","T1")),
    type      = "partial",
    section   = "Section",
    unit      = "Variety"
  )
}

test_that("cond_x = 1L uses A_j[1L] for all panels (default behaviour)", {
  res <- make_mock_res()   # baseline: A_j = list(T1={T0}, T2={T0})
  df_default <- plot_randomRegress(res, type = "regress",
                                   cond_x = 1L, return_data = TRUE)
  # x should be T0 BLUPs for all rows
  t0_blups <- res$blups$T0
  expect_equal(df_default$x, rep(t0_blups, 2L), tolerance = 1e-9)
})

test_that("cond_x = 2L uses AVP: x is T2 partial residual (T2 - gamma*T1)", {
  res <- make_mock_res_partial()
  # T0 conditioned on {T1, T2}: cond_x=2 -> k=T2, A_rest={T1}
  # x = T2 - (G_ss[T1,T2]/G_ss[T1,T1]) * T1
  df_cx <- plot_randomRegress(res, type = "regress",
                               cond_x = 2L, return_data = TRUE)
  t0_rows_e1 <- df_cx[df_cx$pair_label == "T0 | T1, T2" &
                       df_cx$Section == "E1", ]

  # Extract G_ss for site E1: columns 1:3 of Gmat (T0-E1, T1-E1, T2-E1)
  nms    <- colnames(res$Gmat)
  e1_idx <- which(grepl("E1", nms) & (grepl("^T0-", nms) |
                                       grepl("^T1-", nms) |
                                       grepl("^T2-", nms)))
  G_ss   <- res$Gmat[e1_idx, e1_idx]
  rownames(G_ss) <- colnames(G_ss) <- c("T0","T1","T2")

  gamma_k   <- G_ss["T1","T2"] / G_ss["T1","T1"]
  blup_e1   <- res$blups[res$blups$Section == "E1", ]
  x_expected <- blup_e1$T2 - gamma_k * blup_e1$T1

  expect_equal(t0_rows_e1$x, x_expected, tolerance = 1e-9)
})

test_that("cond_x = 2L AVP: avp_active attribute is TRUE", {
  res <- make_mock_res_partial()
  df  <- plot_randomRegress(res, type = "regress",
                             cond_x = 2L, return_data = TRUE)
  expect_true(isTRUE(attr(df, "avp_active")))
})

test_that("cond_x = 1L (baseline): avp_active is FALSE", {
  res <- make_mock_res()   # baseline, |A_j| = 1
  df  <- plot_randomRegress(res, type = "regress",
                             cond_x = 1L, return_data = TRUE)
  expect_false(isTRUE(attr(df, "avp_active")))
})

test_that("cond_x vector: T0 panel uses AVP, T2 panel (A_j size=2, cx=1) also AVP", {
  res <- make_mock_res_partial()
  # All panels have |A_j| = 2, so AVP is always active regardless of cond_x
  df_cx <- plot_randomRegress(res, type = "regress",
                               cond_x = c(2L, 1L, 1L), return_data = TRUE)
  expect_true(isTRUE(attr(df_cx, "avp_active")))
  # T1 panel: k = A_T1[1] = T0, A_rest = {T2}
  # x = T0 - (G_ss[T2,T0]/G_ss[T2,T2]) * T2
  t1_rows_e1 <- df_cx[df_cx$pair_label == "T1 | T0, T2" &
                       df_cx$Section == "E1", ]
  nms    <- colnames(res$Gmat)
  e1_idx <- which(grepl("E1", nms) & (grepl("^T0-", nms) |
                                       grepl("^T1-", nms) |
                                       grepl("^T2-", nms)))
  G_ss   <- res$Gmat[e1_idx, e1_idx]
  rownames(G_ss) <- colnames(G_ss) <- c("T0","T1","T2")
  gamma_k    <- G_ss["T2","T0"] / G_ss["T2","T2"]
  blup_e1    <- res$blups[res$blups$Section == "E1", ]
  x_expected <- blup_e1$T0 - gamma_k * blup_e1$T2
  expect_equal(t1_rows_e1$x, x_expected, tolerance = 1e-9)
})

test_that("cond_x out-of-range index warns and falls back to 1", {
  res <- make_mock_res()   # baseline: all A_j size = 1
  expect_warning(
    plot_randomRegress(res, type = "regress", cond_x = 5L),
    regexp = "out of range"
  )
})

test_that("cond_x scalar recycled to all panels produces ggplot", {
  res <- make_mock_res_partial()
  expect_s3_class(
    plot_randomRegress(res, type = "regress", cond_x = 2L),
    "ggplot"
  )
})

test_that("cond_x invalid type gives informative error", {
  expect_error(
    plot_randomRegress(make_mock_res(), cond_x = "T0"),
    regexp = "positive integer"
  )
})

test_that("cond_x attr on data frame is consistent treatment label", {
  res <- make_mock_res()   # baseline: all A_j = {T0}, cond_lv always T0
  df  <- plot_randomRegress(res, type = "regress",
                             cond_x = 1L, return_data = TRUE)
  expect_equal(attr(df, "cond_lv"), "T0")
})

test_that("cond_x attr is NULL when x-axis treatment varies across panels", {
  res <- make_mock_res_partial()
  # cond_x = c(2,1,1): T0 uses A_T0[2]=T2; T1 uses A_T1[1]=T0; T2 uses A_T2[1]=T0
  # cond_lv_used = c("T2","T0","T0") -> not all same -> attr NULL
  df  <- plot_randomRegress(res, type = "regress",
                             cond_x = c(2L, 1L, 1L), return_data = TRUE)
  expect_null(attr(df, "cond_lv"))
})

test_that("cond_x = 2L scalar: all panels use A_j[2], attr is NULL (T2, T2, T1)", {
  res <- make_mock_res_partial()
  # T0: A_j[2]=T2; T1: A_j[2]=T2; T2: A_j[2]=T1 -> c("T2","T2","T1") -> not all same
  df  <- plot_randomRegress(res, type = "regress",
                             cond_x = 2L, return_data = TRUE)
  expect_null(attr(df, "cond_lv"))
})

# ---- Identifying column names ------------------------------------------
#
# The section and unit columns are named by the analysis, not by this
# function, so the plot frames must follow whatever randomRegress() recorded
# in $section / $unit.

test_that("return_data column names match those of res$blups", {
  res <- make_mock_res()
  for (tp in c("regress", "quadrant")) {
    df <- plot_randomRegress(res, type = tp, return_data = TRUE)
    expect_equal(names(df)[1:2], names(res$blups)[1:2])
    expect_equal(names(df)[1:2], c(res$section, res$unit))
  }
})

test_that("a unit factor called Genotype is carried through to the plot data", {
  res <- make_mock_res(ucol = "Genotype")
  df  <- plot_randomRegress(res, type = "regress", return_data = TRUE)
  expect_equal(names(df)[1:2], c("Section", "Genotype"))
  expect_false("Variety" %in% names(df))
  expect_setequal(unique(df$Genotype), unique(res$blups$Genotype))
  # the names are also what the builders tag on for the plot layers
  expect_equal(attr(df, "unit_col"),    "Genotype")
  expect_equal(attr(df, "section_col"), "Section")
})

test_that("Genotype-named units plot, highlight and label without error", {
  res <- make_mock_res(ucol = "Genotype")
  for (tp in c("regress", "quadrant")) {
    expect_s3_class(plot_randomRegress(res, type = tp), "ggplot")       # default hl
    expect_s3_class(
      plot_randomRegress(res, type = tp, highlight = c("V1", "V2")), "ggplot")
    expect_s3_class(
      plot_randomRegress(res, type = tp, highlight = NULL), "ggplot")
  }
  # The highlight frame is keyed by the unit column, so merge() must match
  hl <- biomAid:::.rreg_default_highlights(
    plot_randomRegress(res, type = "quadrant", return_data = TRUE))
  expect_named(hl, c("Genotype", "group"))
  expect_gt(nrow(hl), 0L)
})

test_that("the plots build (not merely construct) for a renamed unit factor", {
  # ggplot() is lazy: a bad aes or facet specification only fails on build,
  # so constructing the object is not evidence that it renders.
  res <- make_mock_res(ucol = "Genotype")
  for (tp in c("regress", "quadrant", "gmat"))
    expect_s3_class(
      ggplot2::ggplot_build(plot_randomRegress(res, type = tp)),
      "ggplot_built"
    )
})

test_that("legacy result without $section/$unit falls back to column position", {
  # Results produced before the rename carry Site / Variety and neither
  # element; the first two columns of $blups still identify the row.
  res <- make_mock_res(legacy = TRUE)
  expect_null(res$section)
  expect_null(res$unit)
  expect_equal(biomAid:::.rreg_cols(res),
               list(section = "Site", unit = "Variety"))

  df <- plot_randomRegress(res, type = "regress", return_data = TRUE)
  expect_equal(names(df)[1:2], c("Site", "Variety"))
  for (tp in c("regress", "quadrant"))
    expect_s3_class(
      ggplot2::ggplot_build(plot_randomRegress(res, type = tp)),
      "ggplot_built"
    )
})

test_that("stored names that are not columns of blups fall back to position", {
  res <- make_mock_res()
  res$section <- "NoSuchColumn"
  res$unit    <- NA_character_
  expect_equal(biomAid:::.rreg_cols(res),
               list(section = "Section", unit = "Variety"))
})

test_that(".rreg_cols() rejects a blups frame that cannot identify a row", {
  expect_error(
    biomAid:::.rreg_cols(list(blups = data.frame(x = 1))),
    regexp = "section column"
  )
})
