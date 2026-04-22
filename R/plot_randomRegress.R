# ============================================================
#  plot_randomRegress.R
#  ggplot2 visualisation for randomRegress() output.
# ============================================================

# ---- Data preparation helpers ------------------------------------------

#' @noRd
.rreg_regress_data <- function(res, treatments) {

  blups     <- res$blups
  cond_list <- res$cond_list
  conditioned <- names(Filter(Negate(is.null), cond_list))

  if (!is.null(treatments))
    conditioned <- intersect(conditioned, treatments)
  if (length(conditioned) == 0L)
    stop("No conditioned treatments to plot. Check 'treatments'.")

  rows <- lapply(conditioned, function(lv_j) {
    A_j     <- cond_list[[lv_j]]
    cond_lv <- A_j[1L]                   # first (primary) conditioning treatment
    beta_j  <- res$beta[[lv_j]]          # ns x |A_j| matrix, rownames = sites

    # Map site-specific beta onto every blups row by Site
    site_beta <- setNames(beta_j[, cond_lv], rownames(beta_j))

    data.frame(
      Site       = blups$Site,
      Variety    = blups$Variety,
      x          = blups[[cond_lv]],
      y          = blups[[lv_j]],
      pair_label = paste0(lv_j, " | ", cond_lv),
      beta       = site_beta[as.character(blups$Site)],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' @noRd
.rreg_quadrant_data <- function(res, treatments) {

  blups     <- res$blups
  cond_list <- res$cond_list
  conditioned <- names(Filter(Negate(is.null), cond_list))

  if (!is.null(treatments))
    conditioned <- intersect(conditioned, treatments)
  if (length(conditioned) == 0L)
    stop("No conditioned treatments to plot. Check 'treatments'.")

  # Unconditional (efficiency) treatment — x-axis in every panel
  eff_lv <- names(Filter(is.null, cond_list))[1L]

  rows <- lapply(conditioned, function(lv_j) {
    resp_col <- paste0("resp.", lv_j)
    if (!(resp_col %in% names(blups)))
      stop("Column '", resp_col, "' not found in res$blups.")

    # pair label mirrors the regress convention for easy comparison
    cond_lv <- cond_list[[lv_j]][1L]

    data.frame(
      Site           = blups$Site,
      Variety        = blups$Variety,
      x              = blups[[eff_lv]],
      y              = blups[[resp_col]],
      pair_label     = paste0(lv_j, " | ", cond_lv),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}



#' @noRd
.rreg_gmat_data <- function(res) {

  Gmat <- res$Gmat
  if (is.null(Gmat))
    stop("res$Gmat is NULL.")

  # Convert covariance to correlation
  cormat <- cov2cor(Gmat)
  nms    <- rownames(cormat)
  n      <- length(nms)

  # Full matrix to long (upper + lower + diagonal)
  rows <- vector("list", n * n)
  idx  <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      rows[[idx]] <- data.frame(
        row_var = nms[i],
        col_var = nms[j],
        corr    = cormat[i, j],
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  df <- do.call(rbind, rows)

  # Order factors to preserve matrix layout
  df$row_var <- factor(df$row_var, levels = nms)
  df$col_var <- factor(df$col_var, levels = rev(nms))
  df
}


# ---- Plot builders -------------------------------------------------------

#' @noRd
.rreg_plot_regress <- function(df, theme, ...) {

  # One regression line + beta annotation per panel (pair_label x Site)
  panel_df <- unique(df[, c("pair_label", "Site", "beta")])
  panel_df$beta_label <- paste0("\u03b2 = ", round(panel_df$beta, 2L))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
    ggplot2::geom_point(size = 1.6, alpha = 0.75, colour = "#4E79A7", ...) +
    ggplot2::geom_abline(
      data     = panel_df,
      ggplot2::aes(slope = beta, intercept = 0),
      linetype = "dotted", linewidth = 0.7, colour = "grey20"
    ) +
    ggplot2::geom_text(
      data         = panel_df,
      ggplot2::aes(x = -Inf, y = Inf, label = beta_label),
      hjust        = -0.15,
      vjust        = 1.5,
      size         = 3.2,
      colour       = "grey20",
      fontface     = "italic",
      inherit.aes  = FALSE
    ) +
    ggplot2::facet_grid(pair_label ~ Site, scales = "free") +
    ggplot2::labs(
      x       = "Conditioning treatment BLUP",
      y       = "Conditioned treatment BLUP",
      caption = "Dotted line: site-specific random regression"
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
.rreg_plot_quadrant <- function(df, theme, ...) {

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                        linewidth  = 0.7, colour = "grey30") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted",
                        linewidth  = 0.7, colour = "grey30") +
    ggplot2::geom_point(size = 1.6, alpha = 0.75, colour = "#E15759", ...) +
    ggplot2::facet_grid(pair_label ~ Site, scales = "free") +
    ggplot2::labs(
      x       = "Efficiency (unconditional BLUP)",
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
      low      = "#4575B4",
      mid      = "white",
      high     = "#D73027",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Correlation"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1,
                                            size  = 7),
      axis.text.y  = ggplot2::element_text(size = 7)
    )
  p
}


# ---- Public function -----------------------------------------------------

#' Plot Output from randomRegress()
#'
#' @description
#' Produces one of four ggplot2 visualisations from a [randomRegress()] result
#' object.  All plot types return a `ggplot` object that can be further
#' customised with the standard `+` operator.
#'
#' @details
#' The four `type` options are:
#' \describe{
#'   \item{`"regress"`}{Grid of scatter plots faceted by BLUP pair (rows) and
#'     site (columns).  Each panel plots the raw conditioned-treatment BLUPs
#'     (y) against the conditioning-treatment BLUPs (x) for one site × one
#'     treatment pair.  A dotted random regression line with the site-specific
#'     beta slope passes through the origin.}
#'   \item{`"quadrant"`}{Grid of scatter plots faceted by BLUP pair (rows) and
#'     site (columns).  Each panel plots responsiveness BLUPs (y) against the
#'     conditioning treatment BLUPs (x = efficiency) for one site × one
#'     treatment pair.  Dotted zero reference lines on both axes divide each
#'     panel into four quadrants.}
#'   \item{`"gmat"`}{Heatmap of the G-matrix converted to a correlation
#'     matrix via [stats::cov2cor()].  Fill uses a diverging palette centred
#'     at zero.}
#' }
#'
#' @param res      A list returned by [randomRegress()].
#' @param type     Character string selecting the plot type. One of
#'   `"regress"` (default), `"quadrant"`, or `"gmat"`.
#' @param treatments Character vector restricting which conditioned treatments
#'   are included in the plot.  `NULL` (default) includes all conditioned
#'   treatments.  Ignored for `type = "gmat"`.  When `type = "regress"` the
#'   site-specific \eqn{\hat{\beta}} for the selected treatments is annotated
#'   in the top-left corner of each panel.
#' @param theme    A complete ggplot2 theme object applied to the plot.
#'   Defaults to [ggplot2::theme_bw()].
#' @param return_data Logical.  If `TRUE` the function returns the tidy data
#'   frame used to build the plot instead of the plot itself.  Useful for
#'   bespoke customisation. Default `FALSE`.
#' @param ...      Additional arguments passed to the primary `geom_point()`
#'   call (e.g. `size`, `alpha`, `shape`).  Not used for `"beta"` or
#'   `"gmat"` types.
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
#' # Regression plot (default)
#' plot_randomRegress(res)
#'
#' # Quadrant plot, further customised
#' plot_randomRegress(res, type = "quadrant") +
#'   ggplot2::labs(title = "Efficiency vs Responsiveness")
#'
#' # G-matrix correlation heatmap
#' plot_randomRegress(res, type = "gmat")
#'
#' # Retrieve the tidy data frame instead of the plot
#' df <- plot_randomRegress(res, type = "beta", return_data = TRUE)
#' }
#'
#' @export
plot_randomRegress <- function(res,
                               type        = c("regress", "quadrant", "gmat"),
                               treatments  = NULL,
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

  # ---- Build tidy data ---------------------------------------------------
  df <- switch(type,
    regress  = .rreg_regress_data( res, treatments),
    quadrant = .rreg_quadrant_data(res, treatments),
    gmat     = .rreg_gmat_data(    res)
  )

  if (return_data) return(df)

  # ---- Build plot --------------------------------------------------------
  switch(type,
    regress  = .rreg_plot_regress( df, theme, ...),
    quadrant = .rreg_plot_quadrant(df, theme, ...),
    gmat     = .rreg_plot_gmat(    df, theme, ...)
  )
}
