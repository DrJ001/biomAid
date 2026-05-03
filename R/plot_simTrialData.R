# =============================================================================
# plot_simTrialData.R  --  ggplot2 visualisations for simTrialData() output
# =============================================================================
#
# Four plot types:
#
#   "trial"       -- Tiled field layout heatmap: Row x Column grid faceted by
#                    site. fill = any data column (NULL -> Rep); label = any
#                    data column overlaid as text. Numeric fill produces a
#                    continuous viridis scale; factor/character fill produces
#                    a discrete palette. Ideal for inspecting block structure,
#                    treatment strip layout, or the raw yield surface.
#
#   "incidence"   -- Variety x Site presence/absence tile grid. Varieties are
#                    annotated with their site count; sites are annotated with
#                    their variety count. sort = TRUE reorders varieties from
#                    most- to least-connected (top to bottom).
#
#   "correlation" -- Full heatmap of the true genetic correlation matrix
#                    cov2cor(params$G). Diverging RdBu scale centred at 0.
#                    Correlation values printed when the matrix is <= 15 x 15.
#
#   "blup"        -- Heatmap of the true genetic BLUPs (params$g_arr):
#                    Variety x Group (Site or TSite). Diverging RdBu scale
#                    centred at 0 reveals the GEI structure used to generate
#                    the data. sort = TRUE orders varieties by decreasing mean
#                    BLUP (highest at top).
#
# All types support return_data = TRUE.
# =============================================================================

utils::globalVariables(c(
  # trial
  "Row", "Column",
  # incidence
  "Variety", "Site", "Present",
  # correlation
  "Group1", "Group2", "Correlation",
  # blup
  "Group", "BLUP"
))

#' @importFrom stats cov2cor
NULL

# =============================================================================
# MAIN EXPORTED FUNCTION
# =============================================================================

#' Plot simulated trial data from \code{simTrialData()}
#'
#' @description
#' Four complementary plot types for exploring simulated plant breeding trial
#' data returned by \code{\link{simTrialData}()}:
#'
#' \describe{
#'   \item{\code{"trial"}}{A tiled field layout heatmap showing the Row x
#'     Column grid for each site, faceted by site. Tiles can be filled and/or
#'     labelled by any column in \code{res$data}. Use \code{fill = NULL}
#'     (default) to colour by replicate and inspect block structure; use
#'     \code{fill = "yield"} for a continuous yield surface; use
#'     \code{fill = "Variety"} or \code{fill = "Treatment"} to check
#'     randomisation and treatment strip layout.}
#'   \item{\code{"incidence"}}{A Variety x Site presence/absence tile grid
#'     showing which varieties appear at which sites. Each variety label is
#'     annotated with its site count \code{[n]}; each site label is annotated
#'     with its variety count \code{[n]}. Particularly useful for visualising
#'     unbalanced designs generated with \code{incidence = "unbalanced"}.}
#'   \item{\code{"correlation"}}{A full heatmap of the true genetic correlation
#'     matrix \code{cov2cor(params$G)}. A diverging RdBu palette centred at
#'     zero reveals between-environment genetic correlations. Correlation values
#'     are printed when the matrix has 15 or fewer groups.}
#'   \item{\code{"blup"}}{A heatmap of the true genetic BLUPs
#'     (\code{params$g_arr}) with varieties on the y-axis and environments (or
#'     Treatment x Site groups) on the x-axis. The diverging RdBu scale
#'     centred at zero reveals the genotype-by-environment interaction (GEI)
#'     structure embedded in the simulated data.}
#' }
#'
#' @param res         List returned by \code{\link{simTrialData}()}.
#' @param type        Plot type. One of \code{"trial"} (default),
#'   \code{"incidence"}, \code{"correlation"}, or \code{"blup"}.
#' @param fill        \emph{\code{"trial"} type only.} Name of a column in
#'   \code{res$data} to map to the tile fill colour, or \code{NULL} (default)
#'   to colour by replicate (\code{Rep}). Numeric columns use a continuous
#'   viridis scale; character or factor columns use a discrete palette.
#' @param label       \emph{\code{"trial"} type only.} Name of a column in
#'   \code{res$data} to overlay as text on each plot cell, or \code{NULL}
#'   (default) for no labels. Works best with short values such as variety
#'   abbreviations or replicate numbers; \code{"yield"} will show rounded
#'   values.
#' @param sites       Character vector of site names to include, or \code{NULL}
#'   (default) for all sites. Applied to \code{"trial"} and \code{"blup"}
#'   types. For \code{"blup"} in a multi-treatment design, pass the full
#'   \code{TSite} labels (e.g. \code{"T0-Env01"}); partial matching by
#'   substring is used as a fallback.
#' @param sort        Logical. When \code{TRUE} (default):
#'   \itemize{
#'     \item \code{"incidence"}: varieties ordered top-to-bottom by decreasing
#'       site count (most-connected variety at top).
#'     \item \code{"blup"}: varieties ordered top-to-bottom by decreasing mean
#'       BLUP across all groups (highest overall performer at top).
#'   }
#'   No effect for \code{"trial"} or \code{"correlation"}.
#' @param ncol        Integer or \code{NULL}. \emph{\code{"trial"} type only.}
#'   Number of columns passed to \code{\link[ggplot2]{facet_wrap}}. \code{NULL}
#'   (default) lets ggplot2 choose automatically.
#' @param theme       A ggplot2 theme object applied to all plots.
#'   Default \code{\link[ggplot2]{theme_bw}()}.
#' @param return_data Logical. If \code{TRUE}, returns a named list with
#'   elements \code{$plot} (the \code{ggplot} object) and \code{$data} (the
#'   data frame used to build the plot). Default \code{FALSE}.
#' @param ...         Currently unused.
#'
#' @return A \code{ggplot} object, or a named list with elements \code{$plot}
#'   and \code{$data} when \code{return_data = TRUE}.
#'
#' @seealso \code{\link{simTrialData}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' out <- simTrialData(nvar = 20, nsite = 4, seed = 1, verbose = FALSE)
#'
#' # Field layout coloured by replicate (default)
#' plot_simTrialData(out)
#'
#' # Yield heatmap with variety labels
#' plot_simTrialData(out, type = "trial", fill = "yield", label = "Variety")
#'
#' # Rep colouring with treatment strips visible in column layout
#' out_sp <- simTrialData(nvar = 10, nsite = 3,
#'                        treatments = c("T0", "T1", "T2"),
#'                        seed = 5, verbose = FALSE)
#' plot_simTrialData(out_sp, type = "trial", fill = "Treatment")
#'
#' # Incidence for an unbalanced design
#' out2 <- simTrialData(nvar = 30, nsite = 8, incidence = "unbalanced",
#'                      seed = 42, verbose = FALSE)
#' plot_simTrialData(out2, type = "incidence")
#'
#' # True genetic correlation matrix
#' plot_simTrialData(out, type = "correlation")
#'
#' # True BLUP heatmap revealing GEI structure
#' plot_simTrialData(out, type = "blup")
#'
#' # Return data for custom modifications
#' res <- plot_simTrialData(out, type = "blup", return_data = TRUE)
#' res$plot + ggplot2::ggtitle("My custom title")
#' head(res$data)
#' }
plot_simTrialData <- function(res,
                              type        = c("trial", "incidence",
                                              "correlation", "blup"),
                              fill        = NULL,
                              label       = NULL,
                              sites       = NULL,
                              sort        = TRUE,
                              ncol        = NULL,
                              theme       = ggplot2::theme_bw(),
                              return_data = FALSE,
                              ...) {

  type <- match.arg(type)

  # ---- Validate res --------------------------------------------------------
  if (!is.list(res) || !all(c("data", "params") %in% names(res)))
    stop("'res' must be the list returned by simTrialData().")
  if (!is.data.frame(res$data))
    stop("'res$data' must be a data frame.")
  needed_params <- c("G", "incidence", "g_arr")
  missing_p <- setdiff(needed_params, names(res$params))
  if (length(missing_p))
    stop("'res$params' is missing: ", paste(missing_p, collapse = ", "),
         ". Re-run simTrialData() with the current version of biomAid.")

  # ---- Scalar argument checks ----------------------------------------------
  if (!is.logical(sort) || length(sort) != 1L)
    stop("'sort' must be a single logical value (TRUE or FALSE).")
  if (!is.logical(return_data) || length(return_data) != 1L)
    stop("'return_data' must be a single logical value (TRUE or FALSE).")

  # ---- Type-specific argument warnings -------------------------------------
  if (!is.null(fill) && type != "trial")
    warning("'fill' is only used for type = \"trial\" and will be ignored.")
  if (!is.null(label) && type != "trial")
    warning("'label' is only used for type = \"trial\" and will be ignored.")
  if (!is.null(ncol) && type != "trial")
    warning("'ncol' is only used for type = \"trial\" and will be ignored.")
  if (!is.null(sites) && !type %in% c("trial", "blup"))
    warning("'sites' is only used for type = \"trial\" or \"blup\" and will be ignored.")

  # ---- Dispatch ------------------------------------------------------------
  out <- switch(type,
    trial       = .pst_trial(res, fill, label, sites, ncol, theme),
    incidence   = .pst_incidence(res, sort, theme),
    correlation = .pst_correlation(res, theme),
    blup        = .pst_blup(res, sites, sort, theme)
  )

  if (return_data) out else out$plot
}


# =============================================================================
# PRIVATE HELPERS
# =============================================================================

# Discrete colour palette (Tableau-10 inspired)
.pst_pal <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
              "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
              "#9C755F", "#BAB0AC")


# -----------------------------------------------------------------------------
# "trial" — field layout heatmap
# -----------------------------------------------------------------------------

.pst_trial <- function(res, fill, label, sites, ncol, theme) {

  dat <- res$data

  # Subset sites
  if (!is.null(sites)) {
    unknown <- setdiff(sites, levels(dat$Site))
    if (length(unknown))
      warning("Sites not found in data and ignored: ",
              paste(unknown, collapse = ", "))
    dat <- dat[dat$Site %in% sites, ]
    dat$Site <- droplevels(dat$Site)
    if (nrow(dat) == 0L)
      stop("No data remaining after subsetting sites.")
  }

  # Default fill column
  fill_col <- if (is.null(fill)) "Rep" else fill
  if (!fill_col %in% names(dat))
    stop(sprintf("fill column \"%s\" not found in res$data.", fill_col))

  if (!is.null(label) && !label %in% names(dat))
    stop(sprintf("label column \"%s\" not found in res$data.", label))

  is_cont <- is.numeric(dat[[fill_col]])

  # Base plot
  p <- ggplot2::ggplot(
         dat,
         ggplot2::aes(x    = Column,
                      y    = Row,
                      fill = .data[[fill_col]])) +
    ggplot2::geom_tile(colour = "grey85", linewidth = 0.25) +
    ggplot2::scale_y_reverse() +
    ggplot2::facet_wrap(~ Site, scales = "free", ncol = ncol) +
    ggplot2::labs(title = "Field Trial Layout",
                  x     = "Column",
                  y     = "Row",
                  fill  = fill_col) +
    theme +
    ggplot2::theme(
      aspect.ratio     = 1,
      strip.background = ggplot2::element_rect(fill = "#E8EEF4"),
      strip.text       = ggplot2::element_text(face = "bold", size = 9),
      axis.text        = ggplot2::element_text(size = 7),
      panel.grid       = ggplot2::element_blank()
    )

  # Fill scale
  if (is_cont) {
    p <- p + ggplot2::scale_fill_viridis_c(option = "D", name = fill_col)
  } else {
    lvls <- if (is.factor(dat[[fill_col]]))
      levels(dat[[fill_col]])
    else
      unique(as.character(dat[[fill_col]]))
    nL <- length(lvls)
    if (nL <= length(.pst_pal)) {
      p <- p + ggplot2::scale_fill_manual(
        values = setNames(.pst_pal[seq_len(nL)], lvls),
        name   = fill_col
      )
    } else {
      p <- p + ggplot2::scale_fill_discrete(name = fill_col)
    }
  }

  # Optional text overlay
  if (!is.null(label)) {
    txt_colour <- if (is_cont) "white" else "grey10"
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data[[label]]),
      size        = 1.8,
      colour      = txt_colour,
      show.legend = FALSE
    )
  }

  list(plot = p, data = dat)
}


# -----------------------------------------------------------------------------
# "incidence" — variety x site presence/absence grid
# -----------------------------------------------------------------------------

.pst_incidence <- function(res, sort, theme) {

  inc <- res$params$incidence  # nvar x nsite integer matrix

  sites_per_var <- rowSums(inc)
  vars_per_site <- colSums(inc)

  # Variety order: most connected at top (last factor level = top of y-axis)
  var_order <- if (sort)
    rownames(inc)[order(sites_per_var, decreasing = FALSE)]   # low → bottom
  else
    rev(rownames(inc))   # original order, top-to-bottom

  site_order <- colnames(inc)

  # Annotated axis labels
  var_labels  <- setNames(
    paste0(var_order,  " [", sites_per_var[var_order],  "]"),
    var_order
  )
  site_labels <- setNames(
    paste0(site_order, "\n[", vars_per_site[site_order], "]"),
    site_order
  )

  # Long format
  df <- data.frame(
    Variety = rep(rownames(inc), times = ncol(inc)),
    Site    = rep(colnames(inc), each  = nrow(inc)),
    Present = as.integer(as.vector(inc)),
    stringsAsFactors = FALSE
  )
  df$Variety <- factor(df$Variety, levels = var_order)
  df$Site    <- factor(df$Site,    levels = site_order)

  n_obs    <- sum(inc)
  pct_obs  <- 100 * mean(inc)

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = Site, y = Variety,
                      fill = factor(Present))) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(
      values = c("0" = "#E0E0E0", "1" = "#2166AC"),
      labels = c("0" = "Absent",  "1" = "Present"),
      name   = NULL
    ) +
    ggplot2::scale_y_discrete(labels = var_labels) +
    ggplot2::scale_x_discrete(labels = site_labels, position = "top") +
    ggplot2::labs(
      title    = "Variety \u00d7 Site Incidence",
      subtitle = sprintf(
        "%d varieties  \u00d7  %d sites  |  %d observations  |  %.1f%% incidence",
        nrow(inc), ncol(inc), n_obs, pct_obs
      ),
      x = NULL, y = NULL
    ) +
    theme +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(size = 8, hjust = 0.5),
      axis.text.y     = ggplot2::element_text(size = 7),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  list(plot = p, data = df)
}


# -----------------------------------------------------------------------------
# "correlation" — true genetic correlation matrix heatmap
# -----------------------------------------------------------------------------

.pst_correlation <- function(res, theme) {

  G         <- res$params$G
  R         <- cov2cor(G)
  nG        <- nrow(R)
  grp_names <- rownames(R)

  # Full matrix → long format
  # y-axis uses reversed levels so [1,1] element is top-left
  df <- expand.grid(
    Group1 = factor(grp_names, levels = grp_names),
    Group2 = factor(grp_names, levels = rev(grp_names)),
    stringsAsFactors = FALSE
  )
  df$Correlation <- as.vector(R[levels(df$Group1),
                                rev(levels(df$Group2))])

  r_off <- R[upper.tri(R)]

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = Group1, y = Group2, fill = Correlation)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_distiller(
      palette   = "RdBu",
      limits    = c(-1, 1),
      direction = 1,
      name      = "Genetic\ncorrelation"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title    = "True Genetic Correlation Matrix",
      subtitle = sprintf(
        "%d group%s  |  mean r = %.3f  |  range [%.3f, %.3f]",
        nG, if (nG == 1L) "" else "s",
        mean(r_off), min(r_off), max(r_off)
      ),
      x = NULL, y = NULL
    ) +
    theme +
    ggplot2::theme(
      axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1,
                                                size = 8),
      axis.text.y       = ggplot2::element_text(size = 8),
      panel.grid        = ggplot2::element_blank(),
      legend.key.height = ggplot2::unit(1.5, "cm")
    )

  # Overlay correlation values when matrix is small enough to be legible
  if (nG <= 15L) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", Correlation)),
      size   = if (nG <= 6L) 3.0 else 2.2,
      colour = "grey15"
    )
  }

  list(plot = p, data = df)
}


# -----------------------------------------------------------------------------
# "blup" — true genetic BLUP heatmap (variety x group)
# -----------------------------------------------------------------------------

.pst_blup <- function(res, sites, sort, theme) {

  g_arr <- res$params$g_arr  # nvar x ngroup

  # Subset groups (sites or TSites)
  if (!is.null(sites)) {
    keep <- colnames(g_arr) %in% sites
    if (!any(keep)) {
      # Fallback: substring match (handles TSite labels like "T0-Env01")
      keep <- vapply(colnames(g_arr),
                     function(g) any(vapply(sites,
                                            function(s) grepl(s, g, fixed = TRUE),
                                            logical(1L))),
                     logical(1L))
      if (any(keep))
        message("'sites' matched as substrings of group labels.")
    }
    if (!any(keep))
      stop("None of the 'sites' values matched any BLUP group names.\n",
           "  Available groups: ",
           paste(head(colnames(g_arr), 6L), collapse = ", "),
           if (ncol(g_arr) > 6L) " ..." else "")
    g_arr <- g_arr[, keep, drop = FALSE]
  }

  nvar   <- nrow(g_arr)
  ngroup <- ncol(g_arr)

  # Variety order: highest mean BLUP at top (last factor level)
  mean_blup <- rowMeans(g_arr)
  var_order <- if (sort)
    rownames(g_arr)[order(mean_blup, decreasing = FALSE)]   # low → bottom
  else
    rev(rownames(g_arr))

  # Long format
  df <- data.frame(
    Variety = rep(rownames(g_arr), times = ngroup),
    Group   = rep(colnames(g_arr), each  = nvar),
    BLUP    = as.vector(g_arr),
    stringsAsFactors = FALSE
  )
  df$Variety <- factor(df$Variety, levels = var_order)
  df$Group   <- factor(df$Group,   levels = colnames(g_arr))

  # Symmetric limits centred at 0
  max_abs <- max(abs(df$BLUP), na.rm = TRUE)

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = Group, y = Variety, fill = BLUP)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_distiller(
      palette   = "RdBu",
      limits    = c(-max_abs, max_abs),
      direction = 1,
      name      = "True BLUP"
    ) +
    ggplot2::labs(
      title    = "True Genetic BLUPs",
      subtitle = sprintf(
        "%d varieties  \u00d7  %d group%s  |  BLUP range [%.1f, %.1f]",
        nvar, ngroup, if (ngroup == 1L) "" else "s",
        min(df$BLUP), max(df$BLUP)
      ),
      x = NULL, y = NULL
    ) +
    theme +
    ggplot2::theme(
      axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1,
                                                size = 8),
      axis.text.y       = ggplot2::element_text(size = 7),
      panel.grid        = ggplot2::element_blank(),
      legend.key.height = ggplot2::unit(1.5, "cm")
    )

  list(plot = p, data = df)
}
