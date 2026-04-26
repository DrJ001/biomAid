# ============================================================
#  plot_compare.R
#  ggplot2 visualisations for compare() output.
#
#  Three plot types:
#   "dotplot"  — sorted predicted values with shaded criterion band
#   "letters"  — compact letter display (CLD) on a sorted dot plot
#   "heatmap"  — pairwise difference matrix coloured by significance
#
#  Design:
#   - res already carries the criterion as a named column (HSD / LSD /
#     Bonferroni), so no re-computation is needed.
#   - Group detection is automatic: columns that are constant within
#     each block of rows sharing the same criterion value are "by"
#     columns; the remaining factor column is the comparison variable.
#   - Multi-factor by groups produce colon-joined labels and are
#     treated as a single facet variable.
#   - All plot types support return_data = TRUE for bespoke ggplot2
#     customisation.
# ============================================================


# ---- Internal helpers ------------------------------------------------------

## Detect which column is the criterion (HSD / LSD / Bonferroni).
#' @noRd
.pc_crit_col <- function(res) {
  cands <- intersect(names(res), c("HSD", "LSD", "Bonferroni"))
  if (length(cands) == 0L)
    stop("No criterion column (HSD, LSD, or Bonferroni) found in 'res'. ",
         "Was it produced by compare()?")
  cands[1L]
}

## Identify (a) the group columns and (b) the comparison variable.
##
## Strategy:
##   - Remove known non-factor bookkeeping columns.
##   - For each remaining factor column, check whether its values are
##     constant within each unique criterion value.  Columns that ARE
##     constant are "by" (group) columns; the one that is NOT is the
##     comparison variable.
##   - Works for zero, one, or multiple by-columns.
#' @noRd
.pc_detect_cols <- function(res, crit_col) {

  skip <- c("predicted.value", "std.error", "avsed", crit_col)
  candidates <- setdiff(names(res), skip)

  # Split rows by unique criterion value to identify groups
  grp_id <- as.character(res[[crit_col]])

  by_cols   <- character(0L)
  comp_cols <- character(0L)

  for (col in candidates) {
    vals_per_grp <- tapply(as.character(res[[col]]), grp_id,
                           function(x) length(unique(x)))
    if (all(vals_per_grp == 1L))
      by_cols <- c(by_cols, col)
    else
      comp_cols <- c(comp_cols, col)
  }

  # comp_cols should be exactly one column; warn if ambiguous
  if (length(comp_cols) == 0L)
    stop("Could not identify the comparison variable in 'res'. ",
         "All factor columns are constant within groups.")
  if (length(comp_cols) > 1L)
    warning("Multiple candidate comparison columns found: ",
            paste(comp_cols, collapse = ", "),
            ". Using '", comp_cols[1L], "'.",
            call. = FALSE)

  list(by_cols  = by_cols,
       comp_col = comp_cols[1L])
}

## Build a composite group label column ("All" when no by-columns).
#' @noRd
.pc_group_label <- function(res, by_cols) {
  if (length(by_cols) == 0L)
    return(rep("All", nrow(res)))
  if (length(by_cols) == 1L)
    return(as.character(res[[by_cols]]))
  apply(res[, by_cols, drop = FALSE], 1L, paste, collapse = ":")
}

## Compact letter display (CLD) for one group.
##
## Algorithm:
##   1. Sort varieties by predicted value descending.
##   2. Start with all varieties sharing letter "a".
##   3. For each pair (i, j) with i < j (in sorted order):
##      if abs(pred_i - pred_j) > criterion, ensure they do NOT share
##      a letter: if they share all letters, assign a new letter to j
##      and all varieties between i and j that are also not significantly
##      different from j.
##   This is the classical "sweep" absorption algorithm.
#' @noRd
.pc_cld_group <- function(pred_vals, names_vals, criterion) {

  n <- length(pred_vals)
  if (n == 0L) return(character(0L))
  if (n == 1L) return(stats::setNames("a", names_vals))

  ord     <- order(pred_vals, decreasing = TRUE)
  pv_s    <- pred_vals[ord]
  nm_s    <- names_vals[ord]

  # Each variety starts with a set of letters (represented as integer set)
  letter_sets <- vector("list", n)
  for (k in seq_len(n)) letter_sets[[k]] <- 1L   # letter "a" = 1

  next_letter <- 2L

  for (i in seq_len(n - 1L)) {
    for (j in seq(i + 1L, n)) {
      if (abs(pv_s[i] - pv_s[j]) > criterion) {
        # i and j are significantly different — must not share a letter
        shared <- intersect(letter_sets[[i]], letter_sets[[j]])
        if (length(shared) > 0L) {
          # Give j (and its non-sig neighbours above it) a new letter
          new_let <- next_letter
          next_letter <- next_letter + 1L
          # Find all varieties between i+1 and j that are NOT sig diff from j
          for (k in seq(i + 1L, j)) {
            if (abs(pv_s[k] - pv_s[j]) <= criterion)
              letter_sets[[k]] <- union(letter_sets[[k]], new_let)
          }
          # Remove the shared letter from j if i still has it
          for (k in seq(i, j)) {
            if (any(letter_sets[[k]] %in% letter_sets[[i]]) &&
                abs(pv_s[k] - pv_s[j]) > criterion) {
              letter_sets[[k]] <- setdiff(letter_sets[[k]], shared)
              if (length(letter_sets[[k]]) == 0L)
                letter_sets[[k]] <- new_let
            }
          }
        }
      }
    }
  }

  # Convert integer sets to letter strings
  int_to_letters <- function(ints) {
    ints <- sort(unique(ints))
    paste(letters[ints], collapse = "")
  }
  letter_strs <- vapply(letter_sets, int_to_letters, character(1L))
  stats::setNames(letter_strs, nm_s)
}

## Build the tidy data frame for the dotplot type.
#' @noRd
.pc_dotplot_data <- function(res, crit_col, comp_col, by_cols) {

  res$Group    <- .pc_group_label(res, by_cols)
  res$Variety  <- as.character(res[[comp_col]])
  res$pred     <- res$predicted.value
  res$crit     <- res[[crit_col]]

  # Within each group rank by predicted value (1 = highest)
  res$rank <- ave(res$pred, res$Group,
                  FUN = function(x) rank(-x, ties.method = "first"))

  # Top-ranked predicted value per group (for the criterion band)
  top_pred <- tapply(res$pred, res$Group, max, na.rm = TRUE)
  res$top_pred <- as.numeric(top_pred[res$Group])

  # Significant vs not relative to the top-ranked variety
  res$sig <- abs(res$pred - res$top_pred) > res$crit

  res[order(res$Group, res$rank), ]
}

## Build the tidy data frame for the letters type.
#' @noRd
.pc_letters_data <- function(res, crit_col, comp_col, by_cols) {

  base <- .pc_dotplot_data(res, crit_col, comp_col, by_cols)

  # Compute CLD per group
  grps <- unique(base$Group)
  letter_rows <- lapply(grps, function(g) {
    sub  <- base[base$Group == g, ]
    lmap <- .pc_cld_group(sub$pred, sub$Variety, unique(sub$crit))
    data.frame(Group   = g,
               Variety = names(lmap),
               letter  = unname(lmap),
               stringsAsFactors = FALSE)
  })
  letter_df <- do.call(rbind, letter_rows)

  merge(base, letter_df, by = c("Group", "Variety"), all.x = TRUE)
}

## Build the tidy data frame for the heatmap type.
#' @noRd
.pc_heatmap_data <- function(res, crit_col, comp_col, by_cols) {

  res$Group   <- .pc_group_label(res, by_cols)
  res$Variety <- as.character(res[[comp_col]])
  res$pred    <- res$predicted.value
  res$crit    <- res[[crit_col]]

  grps <- unique(res$Group)

  rows <- lapply(grps, function(g) {
    sub  <- res[res$Group == g, ]
    # Sort by predicted value descending for axis ordering
    sub  <- sub[order(sub$pred, decreasing = TRUE), ]
    vars <- sub$Variety
    cr   <- unique(sub$crit)
    pv   <- stats::setNames(sub$pred, vars)

    pairs <- expand.grid(Var_x = vars, Var_y = vars,
                         stringsAsFactors = FALSE)
    pairs$pred_x <- pv[pairs$Var_x]
    pairs$pred_y <- pv[pairs$Var_y]
    pairs$diff   <- abs(pairs$pred_x - pairs$pred_y)
    pairs$sig    <- pairs$diff > cr
    pairs$Group  <- g
    pairs$crit   <- cr

    # Factor levels in descending predicted value order
    pairs$Var_x  <- factor(pairs$Var_x, levels = vars)
    pairs$Var_y  <- factor(pairs$Var_y, levels = rev(vars))
    pairs
  })

  do.call(rbind, rows)
}


# ---- Plot builders ---------------------------------------------------------

#' @noRd
.pc_plot_dotplot <- function(df, crit_col, n_groups, theme, ...) {

  base_col <- "#4E79A7"
  sig_col  <- "#E15759"

  # Shaded criterion band around the top-ranked variety per group
  band_df <- unique(df[, c("Group", "top_pred", "crit")])
  band_df$ymin <- band_df$top_pred - band_df$crit
  band_df$ymax <- band_df$top_pred

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = pred, y = stats::reorder(Variety, pred))) +
    ggplot2::geom_rect(
      data        = band_df,
      ggplot2::aes(xmin = ymin, xmax = ymax,
                   ymin = -Inf, ymax = Inf),
      fill        = "#FFEFC0",
      alpha       = 0.6,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(
      data     = band_df,
      ggplot2::aes(xintercept = top_pred),
      linetype = "dashed",
      linewidth = 0.45,
      colour   = "grey40"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = sig),
      size = 2.4, ...
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = base_col, "TRUE" = sig_col),
      labels = c("FALSE" = paste0("Not sig. from best (", crit_col, ")"),
                 "TRUE"  = paste0("Sig. from best (", crit_col, ")")),
      name   = NULL
    ) +
    ggplot2::labs(
      x       = "Predicted value",
      y       = NULL,
      caption = paste0("Shaded band: ", crit_col,
                       " criterion below the top-ranked variety. ",
                       "Red points fall outside the band.")
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      legend.position  = "bottom",
      panel.spacing    = ggplot2::unit(0.5, "lines"),
      axis.text.y      = ggplot2::element_text(size = 7)
    )

  if (n_groups > 1L)
    p <- p + ggplot2::facet_wrap(~ Group, scales = "free")

  p
}

#' @noRd
.pc_plot_letters <- function(df, crit_col, n_groups, theme, ...) {

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = pred, y = stats::reorder(Variety, pred))) +
    ggplot2::geom_vline(
      data      = unique(df[, c("Group", "top_pred")]),
      ggplot2::aes(xintercept = top_pred),
      linetype  = "dashed",
      linewidth = 0.45,
      colour    = "grey60"
    ) +
    ggplot2::geom_point(colour = "#4E79A7", size = 2.4, ...) +
    ggplot2::geom_text(
      ggplot2::aes(label = letter),
      hjust  = -0.5,
      size   = 3.0,
      colour = "grey20",
      fontface = "bold"
    ) +
    ggplot2::labs(
      x       = "Predicted value",
      y       = NULL,
      caption = paste0("Letters from compact letter display (",
                       crit_col, " criterion). ",
                       "Varieties sharing a letter are not significantly different.")
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.spacing    = ggplot2::unit(0.5, "lines"),
      axis.text.y      = ggplot2::element_text(size = 7)
    )

  if (n_groups > 1L)
    p <- p + ggplot2::facet_wrap(~ Group, scales = "free")

  p
}

#' @noRd
.pc_plot_heatmap <- function(df, crit_col, n_groups, theme, ...) {

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = Var_x, y = Var_y, fill = diff)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_point(
      data        = df[df$sig, ],
      ggplot2::aes(x = Var_x, y = Var_y),
      shape       = 4L,          # cross
      size        = 1.4,
      colour      = "white",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient(
      low  = "#E8F4FC",
      high = "#08306B",
      name = "|Difference|"
    ) +
    ggplot2::labs(
      x       = NULL,
      y       = NULL,
      caption = paste0("Tile colour: absolute pairwise difference. ",
                       "White \u00d7 marks: significant at ",
                       crit_col, " criterion.")
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1,
                                               size  = 7),
      axis.text.y      = ggplot2::element_text(size = 7),
      panel.spacing    = ggplot2::unit(0.5, "lines"),
      legend.position  = "right"
    )

  if (n_groups > 1L)
    p <- p + ggplot2::facet_wrap(~ Group, scales = "free")

  p
}


# ---- Public function -------------------------------------------------------

#' Plot Output from compare()
#'
#' @description
#' Produces one of three ggplot2 visualisations from a [compare()] result
#' object, all correctly handling the `by`-group structure (including
#' multi-factor groups) through automatic faceting.
#'
#' @details
#' The three `type` options are:
#' \describe{
#'   \item{`"dotplot"`}{Varieties sorted by predicted value (highest at top)
#'     within each group.  A shaded band of width equal to the comparison
#'     criterion is drawn below the top-ranked variety — any variety whose
#'     point falls within the band is **not** significantly different from the
#'     best.  Points outside the band are coloured red.  The dashed vertical
#'     line marks the top-ranked predicted value.}
#'   \item{`"letters"`}{Compact letter display (CLD) overlaid on the sorted
#'     dot plot.  Varieties sharing at least one letter are not significantly
#'     different.  Letters are assigned by the standard descending sweep
#'     algorithm using the comparison criterion as the significance threshold.}
#'   \item{`"heatmap"`}{An \eqn{n \times n} tile matrix with varieties sorted
#'     by predicted value on both axes.  Tile colour encodes the absolute
#'     pairwise difference; a white cross (\eqn{\times}) marks pairs that are
#'     significantly different.  Useful for large variety sets where a letter
#'     display becomes unreadable.}
#' }
#'
#' **Group detection** — `plot_compare()` infers which columns are grouping
#' variables and which is the comparison variable directly from `res`.
#' Factor columns that are constant within every group (i.e. within every
#' block of rows sharing the same criterion value) are treated as `by`
#' columns; the remaining factor column is the comparison variable.  When
#' multiple `by` columns are present their values are pasted into a
#' composite facet label.
#'
#' @param res         A data frame returned by [compare()].
#' @param type        Character string selecting the plot type. One of
#'   `"dotplot"` (default), `"letters"`, or `"heatmap"`.
#' @param theme       A complete ggplot2 theme object.
#'   Default [ggplot2::theme_bw()].
#' @param return_data Logical. If `TRUE` returns the tidy data frame used to
#'   build the plot rather than the plot itself. Default `FALSE`.
#' @param ...         Additional arguments passed to the background
#'   `geom_point()` call (e.g. `size`, `alpha`).
#'
#' @return A `ggplot` object (when `return_data = FALSE`) or a `data.frame`
#'   (when `return_data = TRUE`).
#'
#' @seealso [compare()], [ggplot2::ggplot()]
#'
#' @examples
#' \dontrun{
#' res <- compare(model,
#'                term = "Treatment:Site:Variety",
#'                by   = "Site",
#'                type = "HSD")
#'
#' # Sorted dot plot with criterion band
#' plot_compare(res)
#'
#' # Compact letter display
#' plot_compare(res, type = "letters")
#'
#' # Pairwise significance heatmap
#' plot_compare(res, type = "heatmap")
#'
#' # Retrieve the tidy data frame
#' df <- plot_compare(res, return_data = TRUE)
#'
#' # Extend with a custom title
#' plot_compare(res, type = "letters") +
#'   ggplot2::ggtitle("Variety comparison — HSD (Tukey)")
#' }
#'
#' @export
plot_compare <- function(res,
                         type        = c("dotplot", "letters", "heatmap"),
                         theme       = ggplot2::theme_bw(),
                         return_data = FALSE,
                         ...) {

  # ---- Validate -----------------------------------------------------------
  type <- match.arg(type)

  if (!is.data.frame(res))
    stop("'res' must be a data frame returned by compare().")

  required <- c("predicted.value", "std.error", "avsed")
  missing  <- setdiff(required, names(res))
  if (length(missing))
    stop("'res' is missing expected column(s): ",
         paste(missing, collapse = ", "),
         ". Was it produced by compare()?")

  if (!inherits(theme, "theme"))
    stop("'theme' must be a ggplot2 theme object, e.g. ggplot2::theme_bw().")

  if (!is.logical(return_data) || length(return_data) != 1L)
    stop("'return_data' must be a single logical value.")

  # ---- Detect structure ---------------------------------------------------
  crit_col          <- .pc_crit_col(res)
  cols              <- .pc_detect_cols(res, crit_col)
  by_cols           <- cols$by_cols
  comp_col          <- cols$comp_col
  n_groups          <- length(unique(.pc_group_label(res, by_cols)))

  # ---- Build tidy data ----------------------------------------------------
  df <- switch(type,
    dotplot = .pc_dotplot_data(res, crit_col, comp_col, by_cols),
    letters = .pc_letters_data(res, crit_col, comp_col, by_cols),
    heatmap = .pc_heatmap_data(res, crit_col, comp_col, by_cols)
  )

  if (return_data) return(df)

  # ---- Build plot ---------------------------------------------------------
  switch(type,
    dotplot = .pc_plot_dotplot(df, crit_col, n_groups, theme, ...),
    letters = .pc_plot_letters(df, crit_col, n_groups, theme, ...),
    heatmap = .pc_plot_heatmap(df, crit_col, n_groups, theme, ...)
  )
}
