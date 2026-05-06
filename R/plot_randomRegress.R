# ============================================================
#  plot_randomRegress.R
#  ggplot2 visualisation for randomRegress() output.
# ============================================================

utils::globalVariables(c(
  "x", "y", "Variety", "group", "beta_label", "col_var", "row_var",
  "corr", "cov2cor", "ave", "aggregate", "quantile"
))

#' @importFrom stats aggregate ave cov2cor quantile
NULL

# ---- Data preparation helpers ------------------------------------------

# ---- AVP helpers ---------------------------------------------------------

#' Extract the T x T within-site G-matrix block from the full Gmat
#'
#' @param Gmat       Full Gmat with column names combining treatment and site
#' @param site       Site label (e.g. "Env1")
#' @param treatments Character vector of treatment labels needed
#' @param sep        Separator used to join treatment and site labels
#' @return T x T matrix with rownames/colnames = treatments, or NULL on failure
#' @noRd
.rreg_site_Gmat <- function(Gmat, site, treatments, sep) {

  tsnams <- colnames(Gmat)

  if (!any(grepl(sep, tsnams, fixed = TRUE))) {
    # No separator: tsnams ARE the treatment labels (single group, no site)
    ok <- treatments[treatments %in% tsnams]
    if (length(ok) == 0L) return(NULL)
    return(Gmat[ok, ok, drop = FALSE])
  }

  st    <- strsplit(tsnams, split = sep, fixed = TRUE)
  part1 <- vapply(st, `[`, character(1L), 1L)
  part2 <- vapply(st, `[`, character(1L), 2L)

  # Determine which part carries treatment vs site labels
  if (all(treatments %in% part1)) {
    tnam <- part1; snam <- part2
  } else {
    tnam <- part2; snam <- part1
  }

  site_cols <- which(snam == site & tnam %in% treatments)
  if (length(site_cols) == 0L) return(NULL)

  G_ss <- Gmat[site_cols, site_cols, drop = FALSE]
  dimnames(G_ss) <- list(tnam[site_cols], tnam[site_cols])
  G_ss
}

#' Compute added-variable-plot residuals for one conditioning panel
#'
#' Projects both the x-axis treatment (k) and the conditioned treatment (j)
#' onto the orthogonal complement of A_rest using the G-matrix, so that the
#' slope of the resulting scatter equals the stored partial regression
#' coefficient exactly (in population).
#'
#' @param blups_j    Numeric vector — raw BLUPs for conditioned treatment j
#' @param blups_k    Numeric vector — raw BLUPs for x-axis treatment k
#' @param blups_rest Data frame or matrix — raw BLUPs for A_rest treatments
#' @param G_ss       T x T within-site G-matrix (rownames = treatment labels)
#' @param j,k        Treatment label strings
#' @param A_rest     Character vector of remaining conditioning treatments
#' @return List with elements x and y (partial residuals)
#' @noRd
.rreg_avp_xy <- function(blups_j, blups_k, blups_rest, G_ss, j, k, A_rest) {

  if (length(A_rest) == 0L)
    return(list(x = blups_k, y = blups_j))

  G_rr <- G_ss[A_rest, A_rest, drop = FALSE]
  G_rk <- G_ss[A_rest, k,      drop = FALSE]
  G_rj <- G_ss[A_rest, j,      drop = FALSE]

  solve_safe <- function(A, b) {
    tryCatch(
      if (nrow(A) == 1L) drop(b) / A[1L, 1L] else drop(solve(A, b)),
      error = function(e) NULL
    )
  }

  gamma_k <- solve_safe(G_rr, G_rk)
  gamma_j <- solve_safe(G_rr, G_rj)

  # Fallback to raw BLUPs if G is singular or A_rest has NAs
  if (is.null(gamma_k) || is.null(gamma_j) || anyNA(blups_rest))
    return(list(x = blups_k, y = blups_j))

  u_rest <- as.matrix(blups_rest)
  list(
    x = blups_k - drop(u_rest %*% gamma_k),
    y = blups_j - drop(u_rest %*% gamma_j)
  )
}

# ---- Regression data builder ---------------------------------------------

#' @noRd
.rreg_regress_data <- function(res, treatments, centre, cond_x = 1L) {

  blups     <- res$blups
  cond_list <- res$cond_list
  Gmat      <- res$Gmat
  sep       <- if (!is.null(res$sep)) res$sep else "-"  # stored by randomRegress()
  conditioned <- names(Filter(Negate(is.null), cond_list))

  if (!is.null(treatments))
    conditioned <- intersect(conditioned, treatments)
  if (length(conditioned) == 0L)
    stop("No conditioned treatments to plot. Check 'treatments'.")

  n_panels <- length(conditioned)
  cx       <- rep_len(as.integer(cond_x), n_panels)

  # Validate each index against its panel's A_j length
  for (ci in seq_len(n_panels)) {
    lv_j <- conditioned[ci]
    aj_n <- length(cond_list[[lv_j]])
    if (cx[ci] < 1L || cx[ci] > aj_n) {
      warning("cond_x[", ci, "] = ", cx[ci], " is out of range for '", lv_j,
              "' (conditioning set has ", aj_n, " member(s)). Using 1.")
      cx[ci] <- 1L
    }
  }

  cond_lv_used <- character(n_panels)
  avp_active   <- FALSE

  rows <- lapply(seq_len(n_panels), function(ci) {

    lv_j    <- conditioned[ci]
    A_j     <- cond_list[[lv_j]]
    cond_lv <- A_j[cx[ci]]
    A_rest  <- A_j[-cx[ci]]          # all A_j members EXCEPT the selected one
    cond_lv_used[ci] <<- cond_lv

    beta_j    <- res$beta[[lv_j]]    # ns x |A_j|, rownames = sites
    site_beta <- setNames(beta_j[, cond_lv], rownames(beta_j))

    if (length(A_rest) == 0L) {
      # |A_j| = 1: no projection needed — use raw BLUPs
      df <- data.frame(
        Site       = blups$Site,
        Variety    = blups$Variety,
        x          = blups[[cond_lv]],
        y          = blups[[lv_j]],
        pair_label = paste0(lv_j, " | ", paste(A_j, collapse = ", ")),
        beta       = site_beta[as.character(blups$Site)],
        stringsAsFactors = FALSE
      )
    } else {
      # |A_j| >= 2: added variable plot — project per site
      avp_active <<- TRUE
      site_rows <- lapply(unique(blups$Site), function(s) {
        idx <- blups$Site == s
        bs  <- blups[idx, , drop = FALSE]

        G_ss <- .rreg_site_Gmat(Gmat, s,
                                unique(c(lv_j, A_j)), sep)
        avp  <- .rreg_avp_xy(
          blups_j    = bs[[lv_j]],
          blups_k    = bs[[cond_lv]],
          blups_rest = bs[, A_rest, drop = FALSE],
          G_ss       = G_ss,
          j          = lv_j,
          k          = cond_lv,
          A_rest     = A_rest
        )
        data.frame(
          Site       = s,
          Variety    = bs$Variety,
          x          = avp$x,
          y          = avp$y,
          pair_label = paste0(lv_j, " | ", paste(A_j, collapse = ", ")),
          beta       = site_beta[s],
          stringsAsFactors = FALSE
        )
      })
      df <- do.call(rbind, site_rows)
    }

    if (centre) {
      df$x <- df$x + ave(df$x, df$Site, FUN = mean)
      df$y <- df$y + ave(df$y, df$Site, FUN = mean)
    }
    df
  })

  out  <- do.call(rbind, rows)
  x_lv <- if (length(unique(cond_lv_used)) == 1L) cond_lv_used[1L] else NULL
  attr(out, "cond_lv")    <- x_lv
  attr(out, "avp_active") <- avp_active
  out
}

#' @noRd
.rreg_quadrant_data <- function(res, treatments, centre) {

  blups       <- res$blups
  cond_list   <- res$cond_list
  conditioned <- names(Filter(Negate(is.null), cond_list))

  if (!is.null(treatments))
    conditioned <- intersect(conditioned, treatments)
  if (length(conditioned) == 0L)
    stop("No conditioned treatments to plot. Check 'treatments'.")

  # Under "partial" conditioning every treatment has a non-NULL conditioning
  # set, so there is no unconditional (efficiency) treatment. Fall back to
  # the raw BLUP of the first treatment in levs as the x-axis reference.
  uncond <- names(Filter(is.null, cond_list))
  eff_lv <- if (length(uncond) > 0L) uncond[1L] else names(cond_list)[1L]

  rows <- lapply(conditioned, function(lv_j) {
    resp_col <- paste0("resp.", lv_j)
    if (!(resp_col %in% names(blups)))
      stop("Column '", resp_col, "' not found in res$blups.")
    cond_lv <- cond_list[[lv_j]][1L]
    df <- data.frame(
      Site       = blups$Site,
      Variety    = blups$Variety,
      x          = blups[[eff_lv]],
      y          = blups[[resp_col]],
      pair_label = paste0(lv_j, " | ", paste(cond_list[[lv_j]], collapse = ", ")),
      stringsAsFactors = FALSE
    )
    if (centre)
      df$x <- df$x + ave(df$x, df$Site, FUN = mean)
    df
  })
  out <- do.call(rbind, rows)
  attr(out, "eff_lv") <- eff_lv   # carry forward for axis labelling
  out
}

#' @noRd
.rreg_gmat_data <- function(res) {

  Gmat <- res$Gmat
  if (is.null(Gmat)) stop("res$Gmat is NULL.")

  cormat <- cov2cor(Gmat)
  nms    <- rownames(cormat)
  n      <- length(nms)

  rows <- vector("list", n * n)
  idx  <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      rows[[idx]] <- data.frame(row_var = nms[i], col_var = nms[j],
                                corr = cormat[i, j],
                                stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
  df         <- do.call(rbind, rows)
  df$row_var <- factor(df$row_var, levels = nms)
  df$col_var <- factor(df$col_var, levels = rev(nms))
  df
}

# ---- Default highlight selection ---------------------------------------

# Selects up to 3 varieties per quadrant to annotate, choosing the most
# extreme varieties (by distance from origin) that exceed the within-quadrant
# median distance.  The top-right quadrant is shown in orange and the
# bottom-left in blue.
# @param qdata  Data frame from .rreg_quadrant_data()
# @return Data frame with columns Variety and group ("tr" or "bl")
#' @noRd
.rreg_default_highlights <- function(qdata) {

  # Collapse to one representative point per variety
  var_means <- aggregate(cbind(x, y) ~ Variety, data = qdata, FUN = mean)

  .pick_quad <- function(df, sx, sy, grp) {
    in_q <- df[df$x * sx > 0 & df$y * sy > 0, , drop = FALSE]
    if (nrow(in_q) == 0L)
      return(data.frame(Variety = character(0L), group = character(0L),
                        stringsAsFactors = FALSE))

    in_q$dist <- sqrt(in_q$x^2 + in_q$y^2)

    # Keep varieties beyond the within-quadrant median distance
    cands <- in_q[in_q$dist > quantile(in_q$dist, 0.5), , drop = FALSE]

    if (nrow(cands) == 0L)
      return(data.frame(Variety = character(0L), group = character(0L),
                        stringsAsFactors = FALSE))

    # Return up to 3, ordered by decreasing distance (most extreme first)
    cands  <- cands[order(cands$dist, decreasing = TRUE), , drop = FALSE]
    chosen <- head(cands$Variety, 3L)

    data.frame(Variety = chosen, group = grp, stringsAsFactors = FALSE)
  }

  rbind(
    .pick_quad(var_means,  1,  1, "tr"),
    .pick_quad(var_means, -1, -1, "bl")
  )
}

# ---- Shared highlight layer builder ------------------------------------

# Returns a list of ggplot2 layers (geom_point + geom_text) for highlighted
# varieties.  Returns an empty list when hl is NULL.
#' @noRd
.rreg_highlight_layers <- function(df, hl, label_col = "Variety") {

  if (is.null(hl) || nrow(hl) == 0L) return(list())

  hl_data <- merge(df, hl, by = "Variety")

  # Colour palette: top-right warm, bottom-left cool, custom green
  col_map <- c(tr = "#E07B39", bl = "#4E79A7", custom = "#59A14F")

  list(
    ggplot2::geom_point(
      data        = hl_data[hl_data$group == "tr", ],
      ggplot2::aes(x = x, y = y),
      colour      = col_map["tr"],
      size        = 2.8,
      inherit.aes = FALSE
    ),
    ggplot2::geom_point(
      data        = hl_data[hl_data$group == "bl", ],
      ggplot2::aes(x = x, y = y),
      colour      = col_map["bl"],
      size        = 2.8,
      inherit.aes = FALSE
    ),
    ggplot2::geom_point(
      data        = hl_data[hl_data$group == "custom", ],
      ggplot2::aes(x = x, y = y),
      colour      = col_map["custom"],
      size        = 2.8,
      inherit.aes = FALSE
    ),
    ggplot2::geom_text(
      data         = hl_data,
      ggplot2::aes(x = x, y = y, label = Variety,
                   colour = group),
      size         = 2.5,
      vjust        = -0.7,
      fontface     = "bold",
      inherit.aes  = FALSE,
      show.legend  = FALSE
    ),
    ggplot2::scale_colour_manual(
      values = col_map,
      guide  = "none"
    )
  )
}

# ---- Plot builders -------------------------------------------------------

#' @noRd
.rreg_plot_regress <- function(df, hl, centre, theme, ...) {

  panel_df <- unique(df[, c("pair_label", "Site", "beta")])
  panel_df$beta_label <- paste0("\u03b2 = ", round(panel_df$beta, 2L))

  base_col <- if (is.null(hl)) "#4E79A7" else "grey50"

  # x-axis label and caption reflect whether AVP projection was applied.
  cond_lv    <- attr(df, "cond_lv")
  avp_active <- isTRUE(attr(df, "avp_active"))

  x_label <- if (!is.null(cond_lv)) {
    base <- if (avp_active) paste0(cond_lv, " BLUP (partial residual)")
            else paste0(cond_lv, " BLUP")
    if (centre) paste0(base, " + site mean") else base
  } else {
    if (centre) "Conditioning BLUP (+ site mean)" else "Conditioning treatment BLUP"
  }

  plot_caption <- if (avp_active)
    "Added variable plot \u2014 x and y are partial residuals; dotted line: partial regression coefficient"
  else
    "Dotted line: site-specific random regression"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
    ggplot2::geom_point(size = 1.6, alpha = 0.7, colour = base_col, ...) +
    .rreg_highlight_layers(df, hl) +
    ggplot2::geom_abline(
      data     = panel_df,
      ggplot2::aes(slope = beta, intercept = 0),
      linetype = "dotted", linewidth = 0.7, colour = "grey20"
    ) +
    ggplot2::geom_text(
      data        = panel_df,
      ggplot2::aes(x = -Inf, y = Inf, label = beta_label),
      hjust       = -0.15, vjust = 1.5,
      size        = 3.2, colour = "grey20", fontface = "italic",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(pair_label ~ Site, scales = "free") +
    ggplot2::labs(
      x       = x_label,
      y       = if (centre) "Conditioned BLUP (+ site mean)" else
                             "Conditioned treatment BLUP",
      caption = plot_caption
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.spacing    = ggplot2::unit(0.4, "lines")
    )
  p
}

#' @noRd
.rreg_plot_quadrant <- function(df, hl, centre, theme, ...) {

  base_col <- if (is.null(hl)) "#E15759" else "grey50"

  # eff_lv is attached by .rreg_quadrant_data() so axis label reflects
  # whether the x-axis is an unconditional efficiency BLUP or the first
  # treatment in levs (partial conditioning fallback).
  eff_lv <- attr(df, "eff_lv")
  x_label <- if (centre)
    paste0(eff_lv, " BLUP (+ site mean)")
  else
    paste0(eff_lv, " BLUP")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                        linewidth  = 0.7, colour = "grey30") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted",
                        linewidth  = 0.7, colour = "grey30") +
    ggplot2::geom_point(size = 1.6, alpha = 0.7, colour = base_col, ...) +
    .rreg_highlight_layers(df, hl) +
    ggplot2::facet_grid(pair_label ~ Site, scales = "free") +
    ggplot2::labs(
      x       = x_label,
      y       = "Responsiveness (conditional BLUP)",
      caption = "Dotted lines at zero divide each panel into four quadrants"
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.spacing    = ggplot2::unit(0.4, "lines")
    )
  p
}

#' @noRd
.rreg_plot_gmat <- function(df, theme, ...) {

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col_var, y = row_var,
                                         fill = corr)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(abs(corr) > 0.01, round(corr, 2L), "")),
      size = 2.5, colour = "grey20"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#4575B4", mid = "white", high = "#D73027",
      midpoint = 0, limits = c(-1, 1), name = "Correlation"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7)
    )
  p
}


# ---- Public function -----------------------------------------------------

#' Plot Output from randomRegress()
#'
#' @description
#' Produces one of three ggplot2 visualisations from a [randomRegress()] result
#' object.  All plot types return a `ggplot` object that can be further
#' customised with the standard `+` operator.
#'
#' @details
#' The three `type` options are:
#' \describe{
#'   \item{`"regress"`}{Grid of scatter plots faceted by BLUP pair (rows) and
#'     site (columns).  Each panel plots the raw conditioned-treatment BLUPs
#'     (y) against the conditioning-treatment BLUPs (x) for one site x one
#'     treatment pair.  A dotted random regression line with the site-specific
#'     beta slope passes through the origin.  The site-specific
#'     \eqn{\hat{\beta}} is annotated in the top-left corner of each panel.}
#'   \item{`"quadrant"`}{Grid of scatter plots faceted by BLUP pair (rows) and
#'     site (columns).  Each panel plots responsiveness BLUPs (y) against the
#'     conditioning treatment BLUPs (x = efficiency) for one site x one
#'     treatment pair.  Dotted zero reference lines on both axes divide each
#'     panel into four quadrants.}
#'   \item{`"gmat"`}{Heatmap of the G-matrix converted to a correlation
#'     matrix via [stats::cov2cor()].  Fill uses a diverging palette centred
#'     at zero.}
#' }
#'
 #' **Variety highlighting** (`type = "regress"` and `"quadrant"` only):
#' By default, up to six varieties are identified and annotated across all
#' panels — up to three from the top-right quadrant of the efficiency x
#' responsiveness space (above average on both axes, shown in orange) and up
#' to three from the bottom-left quadrant (below average on both axes, shown
#' in blue).  Within each quadrant only varieties whose distance from the
#' origin exceeds the within-quadrant median are considered, and the final
#' selection is the most extreme of those candidates ordered by decreasing
#' distance.  The same varieties are consistently annotated across all site
#' and treatment-pair panels.
#'
#' @param res         A list returned by [randomRegress()].
#' @param type        Character string selecting the plot type. One of
#'   `"regress"` (default), `"quadrant"`, or `"gmat"`.
#' @param treatments  Character vector restricting which conditioned treatments
#'   are included in the plot.  `NULL` (default) includes all conditioned
#'   treatments.  Ignored for `type = "gmat"`.
#' @param centre      Logical.  If `FALSE` (default), BLUPs are plotted on
#'   their natural scale (already centred near zero by the mixed model).
#'   If `TRUE`, the within-site mean of the unconditional treatment is added
#'   back to the x-axis values, placing BLUPs on an approximate absolute
#'   yield scale.  For true ASReml BLUPs the site mean is effectively zero
#'   so the change is minimal; this option is mainly useful when BLUPs have
#'   been computed from treatment means (e.g. in the demo).
#'   Ignored for `type = "gmat"`.
#' @param cond_x      Positive integer or integer vector (default `1L`).
#'   Selects which member of the conditioning set \eqn{A_j} is placed on
#'   the x-axis of the `"regress"` plot for each conditioned-treatment
#'   panel.  A scalar is recycled across all panels; a vector of the same
#'   length as the number of conditioned treatments plotted sets each panel
#'   independently.  For example, with three conditioned treatments under
#'   partial conditioning:
#'   \itemize{
#'     \item `cond_x = 1L` (default) — first member of \eqn{A_j} for every
#'       panel, e.g. `c(1, 1, 1)`.
#'     \item `cond_x = 2L` — second member for every panel.
#'     \item `cond_x = c(2, 1, 2)` — panel-specific selection.
#'   }
#'   An index that exceeds \eqn{|A_j|} for a given panel triggers a warning
#'   and falls back to 1.  Has no effect when \eqn{|A_j| = 1} (e.g.
#'   `type = "baseline"`).  Ignored for `type = "quadrant"` and `"gmat"`.
#' @param highlight   Controls variety annotation for `"regress"` and
#'   `"quadrant"` plots. One of:
#'   \describe{
#'     \item{`"default"`}{Automatically selects 6 varieties (3 top-right,
#'       3 bottom-left) using the polar-angle algorithm described in Details.}
#'     \item{character vector}{Highlight exactly the named varieties.}
#'     \item{`NULL`}{No highlighting; all points drawn in a single colour.}
#'   }
#'   Default is `"default"`.  Ignored for `type = "gmat"`.
#' @param theme       A complete ggplot2 theme object. Default
#'   [ggplot2::theme_bw()].
#' @param return_data Logical. If `TRUE` returns the tidy data frame used to
#'   build the plot rather than the plot itself. Default `FALSE`.
#' @param ...         Additional arguments passed to the background
#'   `geom_point()` call (e.g. `size`, `alpha`, `shape`). Not used for
#'   `"gmat"`.
#'
#' @return A `ggplot` object (when `return_data = FALSE`) or a `data.frame`
#'   (when `return_data = TRUE`).
#'
#' @seealso [randomRegress()], [ggplot2::ggplot()]
#'
#' @examples
#' \dontrun{
#' res <- randomRegress(model, levs = c("N0", "N1", "N2"))
#'
#' # Regression plot with default 6 highlighted varieties
#' plot_randomRegress(res)
#'
#' # Quadrant plot with user-specified highlights
#' plot_randomRegress(res, type = "quadrant",
#'                   highlight = c("Var01", "Var15", "Var22"))
#'
#' # Suppress highlighting
#' plot_randomRegress(res, type = "quadrant", highlight = NULL)
#'
#' # G-matrix correlation heatmap
#' plot_randomRegress(res, type = "gmat")
#'
#' # Retrieve the tidy data frame
#' df <- plot_randomRegress(res, type = "quadrant", return_data = TRUE)
#'
#' # Partial conditioning — default (first conditioning treatment on x)
#' res_part <- randomRegress(model, term = "us(TSite):Variety",
#'                           levs = c("N0","N1","N2"), type = "partial")
#' plot_randomRegress(res_part, type = "regress")              # cond_x = 1L
#'
#' # Show the second conditioning treatment on the x-axis for all panels
#' plot_randomRegress(res_part, type = "regress", cond_x = 2L)
#'
#' # Mixed: second for panels 1 & 3, first for panel 2
#' plot_randomRegress(res_part, type = "regress", cond_x = c(2L, 1L, 2L))
#' }
#'
#' @export
plot_randomRegress <- function(res,
                               type        = c("regress", "quadrant", "gmat"),
                               treatments  = NULL,
                               highlight   = "default",
                               centre      = FALSE,
                               cond_x      = 1L,
                               theme       = ggplot2::theme_bw(),
                               return_data = FALSE,
                               ...) {

  # ---- Validate ----------------------------------------------------------
  type <- match.arg(type)

  required <- c("blups", "beta", "Gmat", "cond_list", "type")
  missing  <- setdiff(required, names(res))
  if (length(missing))
    stop("'res' is missing expected element(s): ",
         paste(missing, collapse = ", "),
         ". Was it produced by randomRegress()?")

  if (!is.null(treatments) && !is.character(treatments))
    stop("'treatments' must be a character vector or NULL.")

  if (!inherits(theme, "theme"))
    stop("'theme' must be a ggplot2 theme object, e.g. ggplot2::theme_bw().")

  if (!is.logical(centre) || length(centre) != 1L)
    stop("'centre' must be a single logical value (TRUE or FALSE).")

  if (!is.numeric(cond_x) || length(cond_x) < 1L ||
      any(is.na(cond_x)) || any(cond_x < 1L) ||
      any(cond_x != as.integer(cond_x)))
    stop("'cond_x' must be a positive integer or integer vector.")

  # ---- Build tidy data ---------------------------------------------------
  df <- switch(type,
    regress  = .rreg_regress_data( res, treatments, centre, cond_x),
    quadrant = .rreg_quadrant_data(res, treatments, centre),
    gmat     = .rreg_gmat_data(    res)
  )

  if (return_data) return(df)

  # ---- Resolve highlights (regress + quadrant only) ----------------------
  hl <- NULL
  if (type %in% c("regress", "quadrant")) {
    if (identical(highlight, "default")) {
      # Always derive highlights from quadrant space
      qdata <- if (type == "quadrant") df else
                 .rreg_quadrant_data(res, treatments, centre)
      hl    <- .rreg_default_highlights(qdata)
    } else if (is.character(highlight) && length(highlight) > 0L) {
      all_vars <- unique(as.character(df$Variety))
      bad      <- setdiff(highlight, all_vars)
      if (length(bad))
        warning("'highlight' varieties not found in data: ",
                paste(bad, collapse = ", "))
      valid <- intersect(highlight, all_vars)
      if (length(valid) > 0L)
        hl <- data.frame(Variety = valid, group = "custom",
                         stringsAsFactors = FALSE)
    }
    # NULL highlight: hl stays NULL
  }

  # ---- Build plot --------------------------------------------------------
  switch(type,
    regress  = .rreg_plot_regress( df, hl, centre, theme, ...),
    quadrant = .rreg_plot_quadrant(df, hl, centre, theme, ...),
    gmat     = .rreg_plot_gmat(    df,             theme, ...)
  )
}
