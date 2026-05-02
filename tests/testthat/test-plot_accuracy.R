# tests/testthat/test-plot_accuracy.R
# Comprehensive tests for plot_accuracy() — no ASReml needed.
#
# Sections:
#   A. Shared fixtures
#   B. Input validation errors
#   C. All 6 types return ggplot objects
#   D. return_data = TRUE — structure checks for all 6 types
#   E. Single-metric calls (accuracy only / gen.H2 only) — facet behaviour
#   F. Both-metric calls — facet presence
#   G. Violin: gen.H2 suppression
#   H. Comparison plots — data structure (dumbbell, scatter, diff)
#   I. Group ordering
#   J. Custom theme
#   K. Edge cases

library(ggplot2)

# ===========================================================================
# SECTION A: Shared fixtures
# ===========================================================================

# Group-level accuracy() output (as returned by accuracy(..., by_variety=FALSE))
make_mock_grp <- function(groups = paste0("E", 1:5), seed = 1L) {
  set.seed(seed)
  data.frame(
    group    = groups,
    n_vars   = rep(20L, length(groups)),
    G_jj     = abs(rnorm(length(groups), 0.8, 0.15)),
    mean_acc = runif(length(groups), 0.60, 0.92),
    sd_acc   = runif(length(groups), 0.03, 0.10),
    gen.H2   = runif(length(groups), 0.45, 0.82),
    stringsAsFactors = FALSE
  )
}

# By-variety accuracy() output (as returned by accuracy(..., by_variety=TRUE))
make_mock_bv <- function(groups = paste0("E", 1:4),
                          vars   = paste0("Var", sprintf("%02d", 1:15)),
                          seed   = 2L) {
  set.seed(seed)
  d <- expand.grid(group = groups, variety = vars,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  d$G_jj     <- rep(abs(rnorm(length(groups), 0.8, 0.15)), each = length(vars))
  d$pev      <- runif(nrow(d), 0.05, 0.50)
  d$accuracy <- sqrt(pmax(0, 1 - d$pev / d$G_jj))
  # gen.H2 is group-level constant — same value for all varieties in a group
  h2_vals    <- runif(length(groups), 0.45, 0.82)
  d$gen.H2   <- rep(h2_vals, each = length(vars))
  d
}

# Accuracy-only versions (no gen.H2 column)
make_mock_grp_acc_only <- function(groups = paste0("E", 1:4), seed = 3L) {
  d <- make_mock_grp(groups, seed)
  d[, setdiff(names(d), "gen.H2")]
}

make_mock_bv_acc_only <- function(groups = paste0("E", 1:4),
                                   vars   = paste0("V", 1:10), seed = 4L) {
  d <- make_mock_bv(groups, vars, seed)
  d[, setdiff(names(d), "gen.H2")]
}

# gen.H2-only versions (no mean_acc / sd_acc)
make_mock_grp_h2_only <- function(groups = paste0("E", 1:4), seed = 5L) {
  d <- make_mock_grp(groups, seed)
  d[, setdiff(names(d), c("mean_acc", "sd_acc"))]
}

# ===========================================================================
# SECTION B: Input validation errors
# ===========================================================================

test_that("unknown type errors", {
  expect_error(
    plot_accuracy(make_mock_grp(), type = "boxplot"),
    regexp = "should be one of|arg"
  )
})

test_that("comparison type without res2 errors", {
  for (tp in c("dumbbell", "scatter", "diff"))
    expect_error(
      plot_accuracy(make_mock_grp(), type = tp),
      "requires res2"
    )
})

test_that("violin without by_variety data errors", {
  expect_error(
    plot_accuracy(make_mock_grp(), type = "violin"),
    "requires by_variety"
  )
})

test_that("heatmap without by_variety data errors", {
  expect_error(
    plot_accuracy(make_mock_grp(), type = "heatmap"),
    "requires by_variety"
  )
})

test_that("scatter without by_variety data in res errors", {
  expect_error(
    plot_accuracy(make_mock_grp(), make_mock_grp(), type = "scatter"),
    "requires by_variety"
  )
})

test_that("scatter without by_variety data in res2 errors", {
  expect_error(
    plot_accuracy(make_mock_bv(), make_mock_grp(), type = "scatter"),
    "requires by_variety.*res2"
  )
})

test_that("metric='accuracy' when no accuracy column in res errors", {
  d <- make_mock_grp_h2_only()
  expect_error(
    plot_accuracy(d, type = "lollipop", metric = "accuracy"),
    "'accuracy' not in"
  )
})

test_that("metric='gen.H2' when no gen.H2 column in res errors", {
  d <- make_mock_grp_acc_only()
  expect_error(
    plot_accuracy(d, type = "lollipop", metric = "gen.H2"),
    "'gen.H2' not in"
  )
})

test_that("single-model type with res2 supplied warns (not errors)", {
  for (tp in c("lollipop", "violin", "heatmap"))
    expect_warning(
      plot_accuracy(
        if (tp == "lollipop") make_mock_grp() else make_mock_bv(),
        res2 = make_mock_grp(),
        type = tp, metric = "accuracy"
      ),
      "ignores res2"
    )
})

# ===========================================================================
# SECTION C: All 6 types return ggplot objects
# ===========================================================================

test_that("lollipop returns ggplot", {
  expect_s3_class(
    plot_accuracy(make_mock_grp(), type = "lollipop", metric = "accuracy"),
    "ggplot"
  )
})

test_that("violin returns ggplot", {
  expect_s3_class(
    plot_accuracy(make_mock_bv(), type = "violin", metric = "accuracy"),
    "ggplot"
  )
})

test_that("heatmap returns ggplot", {
  expect_s3_class(
    plot_accuracy(make_mock_bv(), type = "heatmap", metric = "accuracy"),
    "ggplot"
  )
})

test_that("dumbbell returns ggplot", {
  expect_s3_class(
    plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                  type = "dumbbell", metric = "accuracy"),
    "ggplot"
  )
})

test_that("scatter returns ggplot", {
  expect_s3_class(
    plot_accuracy(make_mock_bv(seed=1L), make_mock_bv(seed=2L),
                  type = "scatter", metric = "accuracy"),
    "ggplot"
  )
})

test_that("diff returns ggplot", {
  expect_s3_class(
    plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                  type = "diff", metric = "accuracy"),
    "ggplot"
  )
})

# ===========================================================================
# SECTION D: return_data = TRUE — structure checks for all 6 types
# ===========================================================================

test_that("lollipop return_data is list with $plot and $data", {
  rd <- plot_accuracy(make_mock_grp(), type = "lollipop",
                      metric = "accuracy", return_data = TRUE)
  expect_named(rd, c("plot", "data"))
  expect_s3_class(rd$plot, "ggplot")
  expect_s3_class(rd$data, "data.frame")
})

test_that("lollipop return_data$data has group, value, metric_label columns", {
  rd <- plot_accuracy(make_mock_grp(), type = "lollipop",
                      metric = "accuracy", return_data = TRUE)
  expect_true(all(c("group", "value", "metric_label") %in% names(rd$data)))
})

test_that("lollipop return_data: one row per group x metric", {
  n_grp <- 5L
  rd <- plot_accuracy(make_mock_grp(paste0("E", 1:n_grp)), type = "lollipop",
                      metric = "accuracy", return_data = TRUE)
  expect_equal(nrow(rd$data), n_grp)
})

test_that("violin return_data is list with $plot and $data", {
  rd <- plot_accuracy(make_mock_bv(), type = "violin",
                      metric = "accuracy", return_data = TRUE)
  expect_named(rd, c("plot", "data"))
  expect_s3_class(rd$data, "data.frame")
})

test_that("violin return_data$data has group, variety, value columns", {
  rd <- plot_accuracy(make_mock_bv(), type = "violin",
                      metric = "accuracy", return_data = TRUE)
  expect_true(all(c("group", "variety", "value") %in% names(rd$data)))
})

test_that("heatmap return_data$data has group, variety, value columns", {
  rd <- plot_accuracy(make_mock_bv(), type = "heatmap",
                      metric = "accuracy", return_data = TRUE)
  expect_true(all(c("group", "variety", "value") %in% names(rd$data)))
})

test_that("dumbbell return_data is list with $plot, $data (list with model1/model2/segments)", {
  rd <- plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                      type = "dumbbell", metric = "accuracy", return_data = TRUE)
  expect_named(rd, c("plot", "data"))
  expect_true(is.list(rd$data))
  expect_true(all(c("model1", "model2", "segments") %in% names(rd$data)))
})

test_that("dumbbell return_data$data$segments has gain_dir column", {
  rd <- plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                      type = "dumbbell", metric = "accuracy", return_data = TRUE)
  expect_true("gain_dir" %in% names(rd$data$segments))
  expect_true(all(rd$data$segments$gain_dir %in% c("Improved", "Declined")))
})

test_that("scatter return_data$data has group, variety, .value_x, .value_y", {
  rd <- plot_accuracy(make_mock_bv(seed=1L), make_mock_bv(seed=2L),
                      type = "scatter", metric = "accuracy", return_data = TRUE)
  expect_true(all(c("group", "variety", ".value_x", ".value_y")
                  %in% names(rd$data)))
})

test_that("diff return_data$data has group, diff_val, gain_dir columns", {
  rd <- plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                      type = "diff", metric = "accuracy", return_data = TRUE)
  expect_true(all(c("group", "diff_val", "gain_dir") %in% names(rd$data)))
  expect_true(all(rd$data$gain_dir %in% c("Improved", "Declined")))
})

test_that("diff return_data: diff_val = model2 - model1", {
  g1 <- make_mock_grp(seed = 10L)
  g2 <- make_mock_grp(seed = 11L)
  rd <- plot_accuracy(g1, g2, type = "diff", metric = "accuracy",
                      return_data = TRUE)
  merged <- merge(g1[, c("group","mean_acc")], g2[, c("group","mean_acc")],
                  by = "group", suffixes = c("_1","_2"))
  merged$expected_diff <- merged$mean_acc_2 - merged$mean_acc_1
  for (grp in rd$data$group) {
    row_d  <- rd$data[rd$data$group == grp, ]
    row_m  <- merged[merged$group == grp, ]
    expect_equal(row_d$diff_val, row_m$expected_diff, tolerance = 1e-10)
  }
})

# ===========================================================================
# SECTION E: Single-metric calls — facet behaviour
# ===========================================================================

test_that("lollipop single metric: no facet_wrap in plot", {
  p  <- plot_accuracy(make_mock_grp_acc_only(), type = "lollipop",
                      metric = "accuracy")
  # ggplot build should not throw; check no FacetWrap
  pb <- ggplot2::ggplot_build(p)
  expect_false(inherits(pb$layout$facet, "FacetWrap"))
})

test_that("lollipop both metrics: facet_wrap present in plot", {
  p  <- plot_accuracy(make_mock_grp(), type = "lollipop",
                      metric = c("accuracy", "gen.H2"))
  pb <- ggplot2::ggplot_build(p)
  expect_true(inherits(pb$layout$facet, "FacetWrap"))
})

test_that("dumbbell both metrics: facet_wrap present", {
  p  <- plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                      type = "dumbbell",
                      metric = c("accuracy", "gen.H2"))
  pb <- ggplot2::ggplot_build(p)
  expect_true(inherits(pb$layout$facet, "FacetWrap"))
})

test_that("diff both metrics: facet_wrap present", {
  p  <- plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                      type = "diff",
                      metric = c("accuracy", "gen.H2"))
  pb <- ggplot2::ggplot_build(p)
  expect_true(inherits(pb$layout$facet, "FacetWrap"))
})

test_that("heatmap both metrics: facet_wrap present", {
  p  <- plot_accuracy(make_mock_bv(), type = "heatmap",
                      metric = c("accuracy", "gen.H2"))
  pb <- ggplot2::ggplot_build(p)
  expect_true(inherits(pb$layout$facet, "FacetWrap"))
})

# ===========================================================================
# SECTION F: gen.H2 only — lollipop works cleanly
# ===========================================================================

test_that("lollipop gen.H2 only returns ggplot without error", {
  expect_s3_class(
    plot_accuracy(make_mock_grp_h2_only(), type = "lollipop",
                  metric = "gen.H2"),
    "ggplot"
  )
})

test_that("lollipop gen.H2 only return_data has metric_label == 'Gen. H\u00b2'", {
  rd <- plot_accuracy(make_mock_grp_h2_only(), type = "lollipop",
                      metric = "gen.H2", return_data = TRUE)
  expect_true(all(rd$data$metric_label == "Gen. H\u00b2"))
})

# ===========================================================================
# SECTION G: Violin — gen.H2 suppression
# ===========================================================================

test_that("violin both metrics: gen.H2 suppressed with message", {
  expect_message(
    plot_accuracy(make_mock_bv(), type = "violin",
                  metric = c("accuracy", "gen.H2")),
    "gen.H2.*suppressed|Gen\\..*suppressed|Cullis.*suppressed"
  )
})

test_that("violin both metrics after suppression: returns ggplot (accuracy only)", {
  p <- suppressMessages(
    plot_accuracy(make_mock_bv(), type = "violin",
                  metric = c("accuracy", "gen.H2"))
  )
  expect_s3_class(p, "ggplot")
})

test_that("violin gen.H2 only (constant within group): returns ggplot", {
  bv_h2 <- make_mock_bv()
  bv_h2 <- bv_h2[, setdiff(names(bv_h2), "accuracy")]
  expect_s3_class(
    plot_accuracy(bv_h2, type = "violin", metric = "gen.H2"),
    "ggplot"
  )
})

test_that("violin accuracy only: no suppression message", {
  expect_no_message(
    plot_accuracy(make_mock_bv_acc_only(), type = "violin",
                  metric = "accuracy")
  )
})

# ===========================================================================
# SECTION H: Comparison plots — data structure
# ===========================================================================

test_that("dumbbell: segments have one row per group x metric", {
  n_grp <- 5L; n_metric <- 1L
  rd <- plot_accuracy(make_mock_grp(paste0("E", 1:n_grp), seed=1L),
                      make_mock_grp(paste0("E", 1:n_grp), seed=2L),
                      type = "dumbbell", metric = "accuracy",
                      return_data = TRUE)
  expect_equal(nrow(rd$data$segments), n_grp * n_metric)
})

test_that("dumbbell both metrics: segments have one row per group x metric", {
  n_grp <- 4L; n_metric <- 2L
  rd <- plot_accuracy(make_mock_grp(paste0("E", 1:n_grp), seed=1L),
                      make_mock_grp(paste0("E", 1:n_grp), seed=2L),
                      type = "dumbbell", metric = c("accuracy", "gen.H2"),
                      return_data = TRUE)
  expect_equal(nrow(rd$data$segments), n_grp * n_metric)
})

test_that("scatter: row count = n_groups x n_varieties", {
  groups <- paste0("E", 1:3)
  vars   <- paste0("V", 1:10)
  rd <- plot_accuracy(make_mock_bv(groups, vars, seed=1L),
                      make_mock_bv(groups, vars, seed=2L),
                      type = "scatter", metric = "accuracy",
                      return_data = TRUE)
  expect_equal(nrow(rd$data), length(groups) * length(vars))
})

test_that("diff: diff_val is finite", {
  rd <- plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                      type = "diff", metric = "accuracy", return_data = TRUE)
  expect_true(all(is.finite(rd$data$diff_val)))
})

test_that("label1 and label2 appear in dumbbell data", {
  rd <- plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                      type        = "dumbbell",
                      metric      = "accuracy",
                      label1      = "Diag",
                      label2      = "FA(2)",
                      return_data = TRUE)
  labels <- unique(c(rd$data$model1$model_label, rd$data$model2$model_label))
  expect_true("Diag"  %in% labels)
  expect_true("FA(2)" %in% labels)
})

# ===========================================================================
# SECTION I: Group ordering
# ===========================================================================

test_that("lollipop: groups ordered by ascending mean accuracy", {
  grp <- make_mock_grp(seed = 20L)
  rd  <- plot_accuracy(grp, type = "lollipop", metric = "accuracy",
                       return_data = TRUE)
  # group factor levels should be in ascending mean accuracy order
  ord      <- order(grp$mean_acc)
  expected <- grp$group[ord]
  expect_equal(levels(rd$data$group), expected)
})

test_that("diff: groups ordered by ascending diff_val", {
  g1 <- make_mock_grp(seed = 21L)
  g2 <- make_mock_grp(seed = 22L)
  rd <- plot_accuracy(g1, g2, type = "diff", metric = "accuracy",
                      return_data = TRUE)
  expect_equal(levels(rd$data$group),
               as.character(rd$data$group[order(rd$data$diff_val)])
  )
})

# ===========================================================================
# SECTION J: Custom theme
# ===========================================================================

test_that("custom theme applied to lollipop without error", {
  expect_s3_class(
    plot_accuracy(make_mock_grp(), type = "lollipop", metric = "accuracy",
                  theme = ggplot2::theme_classic()),
    "ggplot"
  )
})

test_that("custom theme applied to violin without error", {
  expect_s3_class(
    plot_accuracy(make_mock_bv(), type = "violin", metric = "accuracy",
                  theme = ggplot2::theme_minimal()),
    "ggplot"
  )
})

test_that("custom theme applied to dumbbell without error", {
  expect_s3_class(
    plot_accuracy(make_mock_grp(seed=1L), make_mock_grp(seed=2L),
                  type  = "dumbbell", metric = "accuracy",
                  theme = ggplot2::theme_light()),
    "ggplot"
  )
})

# ===========================================================================
# SECTION K: Edge cases
# ===========================================================================

test_that("single group works for lollipop", {
  d <- make_mock_grp(groups = "E1")
  expect_s3_class(
    plot_accuracy(d, type = "lollipop", metric = "accuracy"),
    "ggplot"
  )
})

test_that("single group and single variety works for violin", {
  d <- make_mock_bv(groups = "E1", vars = paste0("V", 1:5))
  expect_s3_class(
    plot_accuracy(d, type = "violin", metric = "accuracy"),
    "ggplot"
  )
})

test_that("groups with no sd_acc (NA) still render lollipop", {
  d         <- make_mock_grp(seed = 30L)
  d$sd_acc  <- NA_real_
  expect_s3_class(
    plot_accuracy(d, type = "lollipop", metric = "accuracy"),
    "ggplot"
  )
})

test_that("return_data = FALSE returns ggplot (not a list)", {
  p <- plot_accuracy(make_mock_grp(), type = "lollipop",
                     metric = "accuracy", return_data = FALSE)
  expect_s3_class(p, "ggplot")
  expect_false(is.list(p) && "plot" %in% names(p))
})

test_that("all 6 types build without ggplot2 error", {
  for (tp in c("lollipop", "violin", "heatmap")) {
    d <- if (tp == "lollipop") make_mock_grp() else make_mock_bv()
    expect_no_error(ggplot2::ggplot_build(
      plot_accuracy(d, type = tp, metric = "accuracy")
    ))
  }
  for (tp in c("dumbbell", "scatter", "diff")) {
    d1 <- if (tp == "scatter") make_mock_bv(seed=1L) else make_mock_grp(seed=1L)
    d2 <- if (tp == "scatter") make_mock_bv(seed=2L) else make_mock_grp(seed=2L)
    expect_no_error(ggplot2::ggplot_build(
      plot_accuracy(d1, d2, type = tp, metric = "accuracy")
    ))
  }
})
