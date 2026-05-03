#' Simulate Multi-Environment Plant Breeding Trial Data
#'
#' @description
#' Generates a balanced or unbalanced field trial dataset with a realistic
#' genetic covariance structure across environments.
#'
#' **Trial structure:** \code{treatments = NULL} produces a MET-only randomised
#' complete block design. Supply \code{treatments = c("T0", "T1", ...)} for a
#' split-plot design where whole plots are treatment strips and sub-plots are
#' varieties; the genetic structure then operates over all Treatment x Site
#' (\code{TSite}) combinations.
#'
#' **Genetic covariance** is controlled by the \code{G} argument:
#' \itemize{
#'   \item \code{G = "auto"} (default): a random symmetric positive-definite
#'     covariance matrix is generated with pairwise correlations approximately
#'     in the range (\code{g_cor_min}, \code{g_cor_max}), set via
#'     \code{sim.options}. Scales to \code{sigma_genetic}.
#'   \item \code{G = matrix}: a user-supplied J x J symmetric positive-definite
#'     covariance matrix used directly. J = \code{nsite} for MET-only, or
#'     J = \code{nsite * length(treatments)} for multi-treatment designs.
#'     Group order must match \code{outer(treatments, sites)}.
#' }
#' Variety BLUPs are drawn from \eqn{\mathcal{N}(\mathbf{0}, \mathbf{G})}
#' via Cholesky decomposition regardless of which mode is used.
#'
#' **Incidence structure** (\code{incidence} argument) controls which varieties
#' are observed at which sites. All \code{nvar} varieties have a true genetic
#' value at every site but only observed varieties contribute plots.
#'
#' @param nvar        Integer. Maximum number of varieties. Default \code{20}.
#' @param nsite       Integer. Number of environments (sites). Default \code{10}.
#' @param treatments  Character vector of treatment labels (length >= 2), or
#'   \code{NULL} (default) for a MET-only trial with no treatment structure.
#' @param nrep        Integer. Replicates per environment. Default \code{2}.
#' @param G           Controls the true genetic covariance structure. One of:
#'   \describe{
#'     \item{\code{"auto"}}{(Default) Auto-generate a random SPD covariance
#'       matrix with pairwise correlations approximately in the range
#'       (\code{g_cor_min}, \code{g_cor_max}) set via \code{sim.options}.
#'       Scale is controlled by \code{sigma_genetic}.}
#'     \item{Matrix}{A user-supplied J x J symmetric positive-definite
#'       covariance matrix used directly. J = \code{nsite} for MET-only designs
#'       or J = \code{nsite * length(treatments)} for multi-treatment designs.
#'       Group order must match \code{outer(treatments, sites)}.}
#'   }
#' @param incidence   Controls which varieties are observed at which sites.
#'   One of:
#'   \describe{
#'     \item{\code{"balanced"}}{(Default) All \code{nvar} varieties are
#'       observed at every site. Produces identical PEV for all varieties
#'       within a site.}
#'     \item{\code{"unbalanced"}}{Automatically generates a two-tier incidence
#'       structure. Proportions are controlled by five \code{sim.options}
#'       entries: \code{core_pct}, \code{core_min_sites_pct},
#'       \code{reg_min_sites_pct}, \code{reg_max_sites_pct}, and
#'       \code{min_vars_pct} (see below).}
#'     \item{Matrix}{An \code{nvar x nsite} matrix of \code{0}/\code{1} or
#'       \code{TRUE}/\code{FALSE}. Entry \code{[v, s] = 1} means variety
#'       \code{v} is observed at site \code{s}. Every variety must appear in
#'       at least 1 site; every site must have at least 2 varieties.}
#'   }
#' @param seed        Integer or \code{NULL}. Passed to \code{\link{set.seed}}.
#'   Default \code{NULL}.
#' @param verbose     Logical. Print a design summary and suggested ASReml-R
#'   model. Default \code{TRUE}.
#' @param sim.options Named list of optional simulation controls. Unknown names
#'   produce a warning. Recognised elements (with defaults) are:
#'   \describe{
#'     \item{\code{site_mean}}{Numeric. Grand mean yield. Default \code{4500}.}
#'     \item{\code{site_sd}}{Numeric. SD of site mean yields. Default \code{600}.}
#'     \item{\code{treat_effects}}{Numeric vector (length =
#'       \code{length(treatments)}) of fixed treatment effects, or \code{NULL}
#'       for auto-spacing from 0 to \code{site_mean * 0.10}.
#'       Multi-treatment only.}
#'     \item{\code{sigma_genetic}}{Numeric. Target mean genetic SD per group.
#'       Default \code{250}.}
#'     \item{\code{rep_sd}}{Numeric. SD of replicate effects. Default \code{150}.}
#'     \item{\code{row_sd}}{Numeric. SD of row spatial effects. Default \code{80}.}
#'     \item{\code{col_sd}}{Numeric. SD of column spatial effects.
#'       Default \code{60}.}
#'     \item{\code{error_sd}}{Numeric. SD of plot-level residual.
#'       Default \code{350}.}
#'     \item{\code{sep}}{Character. Separator for \code{TSite} labels.
#'       Default \code{"-"}.}
#'     \item{\code{variety_prefix}}{Character. Prefix for variety labels.
#'       Default \code{"Var"}.}
#'     \item{\code{site_prefix}}{Character. Prefix for site labels.
#'       Default \code{"Env"}.}
#'     \item{\code{outfile}}{Character path or \code{NULL}. CSV output.
#'       Default \code{NULL}.}
#'     \item{\code{g_cor_min}}{Numeric in (-1, 1). Minimum pairwise genetic
#'       correlation for the auto-generated G matrix. Used only when
#'       \code{G = "auto"}. Default \code{0.20}.}
#'     \item{\code{g_cor_max}}{Numeric in (-1, 1]. Maximum pairwise genetic
#'       correlation. Must be > \code{g_cor_min}. Used only when
#'       \code{G = "auto"}. Default \code{0.90}.}
#'     \item{\code{core_pct}}{Numeric in (0, 1). Proportion of varieties
#'       designated as "core" entries when \code{incidence = "unbalanced"}.
#'       Default \code{0.20}.}
#'     \item{\code{core_min_sites_pct}}{Numeric in (0, 1]. Minimum proportion
#'       of sites that each core variety must appear in. Default \code{0.75}.}
#'     \item{\code{reg_min_sites_pct}}{Numeric in (0, 1]. Minimum proportion
#'       of sites for regular (non-core) varieties. Default \code{0.40}.}
#'     \item{\code{reg_max_sites_pct}}{Numeric in (0, 1]. Maximum proportion
#'       of sites for regular varieties. Must be >= \code{reg_min_sites_pct}.
#'       Default \code{0.85}.}
#'     \item{\code{min_vars_pct}}{Numeric in (0, 1]. Minimum proportion of
#'       \code{nvar} that every site must contain. Default \code{0.40}.}
#'   }
#'
#' @return A named list:
#' \describe{
#'   \item{\code{data}}{Data frame ordered by \code{Site}, \code{Row},
#'     \code{Column} with columns: \code{Site}, \code{Variety}, \code{Rep},
#'     \code{Row}, \code{Column}, \code{yield} (MET-only), or additionally
#'     \code{Treatment} and \code{TSite} (multi-treatment). Only observed
#'     variety x site combinations are included.}
#'   \item{\code{params}}{Named list of true simulation parameters:
#'     \describe{
#'       \item{\code{G}}{True genetic covariance matrix (J x J).}
#'       \item{\code{site_means}}{True site mean yields (length \code{nsite}).}
#'       \item{\code{incidence}}{Integer \code{nvar x nsite} incidence matrix.}
#'       \item{\code{treat_effects}}{True treatment fixed effects
#'         (multi-treatment only).}
#'       \item{\code{g_arr}}{\code{nvar x ngroup} matrix of true genetic BLUPs.
#'         Always returned. Column names are site names for MET-only designs and
#'         \code{TSite} labels for multi-treatment designs.}
#'     }}
#' }
#'
#' @examples
#' \dontrun{
#' # Simplest call — just run it
#' out <- simTrialData(verbose = FALSE)
#' head(out$data)
#'
#' # MET-only: 30 varieties, 8 sites
#' out <- simTrialData(nvar = 30, nsite = 8, seed = 42)
#' round(cov2cor(out$params$G), 2)   # true genetic correlations
#'
#' # Unbalanced: produces spread in per-variety accuracy
#' out2 <- simTrialData(nvar = 30, nsite = 8,
#'                      incidence = "unbalanced", seed = 42)
#' table(rowSums(out2$params$incidence))   # sites per variety
#'
#' # Custom correlation range
#' out3 <- simTrialData(nvar = 20, nsite = 6, seed = 99,
#'                      sim.options = list(g_cor_min = 0.3, g_cor_max = 0.8))
#' round(cov2cor(out3$params$G), 2)
#'
#' # User-supplied G matrix
#' G_true <- matrix(c(40000, 15000, 15000, 40000), 2, 2)
#' out4 <- simTrialData(nvar = 20, nsite = 2, G = G_true, seed = 10)
#'
#' # User-supplied incidence matrix
#' inc <- matrix(1L, nrow = 20, ncol = 6)
#' inc[1:5, 1:3] <- 0L
#' out5 <- simTrialData(nvar = 20, nsite = 6, incidence = inc, seed = 1)
#'
#' # Multi-treatment split-plot
#' out6 <- simTrialData(nvar = 20, nsite = 6,
#'                      treatments  = c("T0", "T1", "T2"),
#'                      incidence   = "unbalanced",
#'                      seed        = 1,
#'                      sim.options = list(treat_effects = c(0, 150, 350)))
#' }
#'
#' @export
simTrialData <- function(nvar        = 20L,
                         nsite       = 10L,
                         treatments  = NULL,
                         nrep        = 2L,
                         G           = "auto",
                         incidence   = "balanced",
                         seed        = NULL,
                         verbose     = TRUE,
                         sim.options = list()) {

  # ---- Merge sim.options with defaults -------------------------------------
  .defaults <- list(
    site_mean           = 4500,
    site_sd             = 600,
    treat_effects       = NULL,
    sigma_genetic       = 250,
    rep_sd              = 150,
    row_sd              = 80,
    col_sd              = 60,
    error_sd            = 350,
    sep                 = "-",
    variety_prefix      = "Var",
    site_prefix         = "Env",
    outfile             = NULL,
    # auto G controls (used only when G = "auto")
    g_cor_min           = 0.20,
    g_cor_max           = 0.90,
    # unbalanced incidence controls (used only when incidence = "unbalanced")
    core_pct            = 0.20,
    core_min_sites_pct  = 0.75,
    reg_min_sites_pct   = 0.40,
    reg_max_sites_pct   = 0.85,
    min_vars_pct        = 0.40
  )

  unknown <- setdiff(names(sim.options), names(.defaults))
  if (length(unknown))
    warning("Unknown sim.options element(s) ignored: ",
            paste(unknown, collapse = ", "), ".")

  opts <- modifyList(.defaults,
                     sim.options[intersect(names(sim.options),
                                           names(.defaults))])

  # ---- Validation ----------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  has_treatments <- !is.null(treatments) && length(treatments) >= 2L

  if (!is.null(treatments) && length(treatments) < 2L)
    stop("'treatments' must have at least 2 elements, or NULL for MET-only.")

  if (has_treatments && any(grepl(opts$sep, treatments, fixed = TRUE)))
    stop("'sep' (\"", opts$sep,
         "\") must not appear inside any treatment label.")

  ntreat <- if (has_treatments) length(treatments) else 1L
  ngroup <- ntreat * nsite

  if (!is.null(opts$treat_effects) && length(opts$treat_effects) != ntreat)
    stop("'treat_effects' must have length ", ntreat, ".")

  # Validate G argument
  use_G_mode <- if (identical(G, "auto")) "auto" else "user"

  if (use_G_mode == "auto") {
    if (!is.numeric(opts$g_cor_min) || opts$g_cor_min <= -1 || opts$g_cor_min >= 1)
      stop("sim.options$g_cor_min must be in (-1, 1).")
    if (!is.numeric(opts$g_cor_max) || opts$g_cor_max <= -1 || opts$g_cor_max > 1)
      stop("sim.options$g_cor_max must be in (-1, 1].")
    if (opts$g_cor_max <= opts$g_cor_min)
      stop("sim.options$g_cor_max must be > g_cor_min.")
  }

  if (use_G_mode == "user") {
    if (!is.matrix(G) || !is.numeric(G))
      stop("User-supplied 'G' must be a numeric matrix.")
    if (!identical(dim(G), c(ngroup, ngroup)))
      stop(sprintf("'G' must be %d x %d (ngroup x ngroup); got %d x %d.",
                   ngroup, ngroup, nrow(G), ncol(G)))
    if (!isSymmetric(G, tol = 1e-8))
      stop("'G' must be a symmetric matrix.")
    ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev <= 0))
      stop("'G' must be positive definite (all eigenvalues > 0).")
  }

  # Validate unbalanced incidence proportions
  if (identical(incidence, "unbalanced")) {
    if (!is.numeric(opts$core_pct) || opts$core_pct <= 0 || opts$core_pct >= 1)
      stop("sim.options$core_pct must be a number in (0, 1).")
    if (!is.numeric(opts$core_min_sites_pct) ||
        opts$core_min_sites_pct <= 0 || opts$core_min_sites_pct > 1)
      stop("sim.options$core_min_sites_pct must be a number in (0, 1].")
    if (!is.numeric(opts$reg_min_sites_pct) ||
        opts$reg_min_sites_pct <= 0 || opts$reg_min_sites_pct > 1)
      stop("sim.options$reg_min_sites_pct must be a number in (0, 1].")
    if (!is.numeric(opts$reg_max_sites_pct) ||
        opts$reg_max_sites_pct <= 0 || opts$reg_max_sites_pct > 1)
      stop("sim.options$reg_max_sites_pct must be a number in (0, 1].")
    if (opts$reg_max_sites_pct < opts$reg_min_sites_pct)
      stop("sim.options$reg_max_sites_pct must be >= reg_min_sites_pct.")
    if (!is.numeric(opts$min_vars_pct) ||
        opts$min_vars_pct <= 0 || opts$min_vars_pct > 1)
      stop("sim.options$min_vars_pct must be a number in (0, 1].")
  }

  # ---- Labels --------------------------------------------------------------
  ndig_var  <- nchar(as.character(nvar))
  ndig_site <- nchar(as.character(nsite))
  varieties <- paste0(opts$variety_prefix,
                      sprintf(paste0("%0", ndig_var,  "d"), seq_len(nvar)))
  sites     <- paste0(opts$site_prefix,
                      sprintf(paste0("%0", ndig_site, "d"), seq_len(nsite)))

  if (has_treatments) {
    tsite_levs <- as.vector(outer(treatments, sites, paste, sep = opts$sep))
    grp_names  <- tsite_levs
  } else {
    grp_names <- sites
  }

  # ---- Incidence matrix ----------------------------------------------------

  .gen_unbalanced_inc <- function(nv, ns, opts) {
    n_core  <- max(1L, ceiling(nv * opts$core_pct))
    mat     <- matrix(0L, nrow = nv, ncol = ns)
    n_core_sites <- max(2L, ceiling(ns * opts$core_min_sites_pct))
    for (v in seq_len(n_core))
      mat[v, sample.int(ns, min(n_core_sites, ns))] <- 1L
    lo <- max(2L, ceiling(ns * opts$reg_min_sites_pct))
    hi <- max(lo, floor(ns * opts$reg_max_sites_pct))
    for (v in seq(n_core + 1L, nv)) {
      k <- sample.int(hi - lo + 1L, 1L) + lo - 1L
      mat[v, sample.int(ns, k)] <- 1L
    }
    min_per_site <- max(ceiling(nv * opts$min_vars_pct), 2L)
    for (s in seq_len(ns)) {
      deficit <- min_per_site - sum(mat[, s])
      if (deficit > 0L) {
        cands <- which(mat[, s] == 0L)
        if (length(cands) > 0L)
          mat[sample(cands, min(deficit, length(cands))), s] <- 1L
      }
    }
    mat
  }

  .validate_inc <- function(inc, nv, ns) {
    if (!is.matrix(inc))
      stop("User-supplied 'incidence' must be a matrix.")
    if (!identical(dim(inc), c(nv, ns)))
      stop(sprintf(
        "'incidence' must be %d x %d (nvar x nsite); got %d x %d.",
        nv, ns, nrow(inc), ncol(inc)))
    if (is.logical(inc)) storage.mode(inc) <- "integer"
    vals <- unique(as.vector(inc))
    if (!all(vals %in% c(0L, 1L)))
      stop("'incidence' values must be 0/1 or TRUE/FALSE.")
    bad_var  <- which(rowSums(inc) == 0L)
    if (length(bad_var))
      stop("Varieties in row(s) ", paste(bad_var, collapse = ", "),
           " appear in no sites.")
    bad_site <- which(colSums(inc) < 2L)
    if (length(bad_site))
      stop("Sites in column(s) ", paste(bad_site, collapse = ", "),
           " have fewer than 2 varieties.")
    storage.mode(inc) <- "integer"
    inc
  }

  if (identical(incidence, "balanced")) {
    inc_mat <- matrix(1L, nrow = nvar, ncol = nsite)
  } else if (identical(incidence, "unbalanced")) {
    inc_mat <- .gen_unbalanced_inc(nvar, nsite, opts)
  } else if (is.matrix(incidence)) {
    inc_mat <- .validate_inc(incidence, nvar, nsite)
  } else {
    stop("'incidence' must be \"balanced\", \"unbalanced\", or an ",
         nvar, " x ", nsite, " matrix of 0s and 1s.")
  }

  rownames(inc_mat) <- varieties
  colnames(inc_mat) <- sites

  # ---- Layout helper (per-site) --------------------------------------------
  .best_dims <- function(n) {
    sq  <- floor(sqrt(n))
    fac <- which(n %% seq_len(sq) == 0L)
    if (!length(fac)) return(c(rows = 1L, cols = as.integer(n)))
    r <- max(fac)
    c(rows = as.integer(r), cols = as.integer(n %/% r))
  }

  # ---- Genetic covariance matrix -------------------------------------------

  # Auto-generate: random SPD G with correlations in (g_cor_min, g_cor_max)
  .gen_G_auto <- function(d, sigma_genetic, cor_min, cor_max) {
    k     <- max(1L, floor(d / 3L))
    Lraw  <- matrix(runif(d * k, 0.5, 2.0), d, k)
    if (k >= 2L)
      for (f in 2L:k) {
        flip <- sample(d, floor(d * 0.4))
        Lraw[flip, f] <- -Lraw[flip, f]
      }
    G_raw <- Lraw %*% t(Lraw) + diag(rep(mean(diag(Lraw %*% t(Lraw))) * 0.1, d))
    R_raw <- cov2cor(G_raw)
    r_off <- R_raw[upper.tri(R_raw)]
    r_min <- min(r_off);  r_max <- max(r_off)
    r_new <- if (abs(r_max - r_min) > 1e-8)
      (r_off - r_min) / (r_max - r_min) * (cor_max - cor_min) + cor_min
    else
      rep((cor_min + cor_max) / 2, length(r_off))
    R_new <- matrix(0, d, d)
    R_new[upper.tri(R_new)] <- r_new
    R_new <- R_new + t(R_new)
    diag(R_new) <- 1
    ev_min <- min(eigen(R_new, symmetric = TRUE, only.values = TRUE)$values)
    if (ev_min <= 1e-8)
      R_new <- cov2cor(R_new + diag(abs(ev_min) + 1e-6, d))
    sigma_genetic^2 * R_new
  }

  G_sim <- if (use_G_mode == "auto")
    .gen_G_auto(ngroup, opts$sigma_genetic, opts$g_cor_min, opts$g_cor_max)
  else
    G   # user-supplied, already validated

  rownames(G_sim) <- colnames(G_sim) <- grp_names

  # ---- BLUP simulation: rows ~ MVN(0, G_sim) via Cholesky -----------------
  Z <- matrix(rnorm(nvar * ngroup), nrow = nvar, ncol = ngroup)
  U <- Z %*% chol(G_sim)
  rownames(U) <- varieties
  colnames(U) <- grp_names

  # ---- Treatment fixed effects ---------------------------------------------
  if (has_treatments) {
    treat_effects <- opts$treat_effects
    if (is.null(treat_effects))
      treat_effects <- seq(0, opts$site_mean * 0.10, length.out = ntreat)
    names(treat_effects) <- treatments
  }

  # ---- Site means ----------------------------------------------------------
  site_means <- setNames(rnorm(nsite, opts$site_mean, opts$site_sd), sites)

  # ---- Build plot-level observations ---------------------------------------
  plots <- vector("list", nvar * nsite * nrep * ntreat)
  idx   <- 1L

  for (s in seq_len(nsite)) {
    vars_s <- which(inc_mat[, s] == 1L)
    n_s    <- length(vars_s)
    if (n_s == 0L) next

    dims_s <- .best_dims(n_s)

    if (has_treatments) {
      rows_strip_s <- dims_s["rows"]
      cols_strip_s <- dims_s["cols"]
      nrow_e_s     <- nrep   * rows_strip_s
      ncol_e_s     <- ntreat * cols_strip_s
    } else {
      rows_block_s <- dims_s["rows"]
      cols_block_s <- dims_s["cols"]
      nrow_e_s     <- nrep * rows_block_s
      ncol_e_s     <- cols_block_s
    }

    row_eff <- rnorm(nrow_e_s, 0, opts$row_sd)
    col_eff <- rnorm(ncol_e_s, 0, opts$col_sd)
    rep_eff <- rnorm(nrep,     0, opts$rep_sd)

    if (has_treatments) {
      for (r in seq_len(nrep)) {
        treat_order <- sample(ntreat)
        for (t_pos in seq_len(ntreat)) {
          t        <- treat_order[t_pos]
          ts_label <- paste(treatments[t], sites[s], sep = opts$sep)
          u_col    <- match(ts_label, grp_names)
          rows_in_rep   <- ((r     - 1L) * rows_strip_s + 1L) :
                            (r            * rows_strip_s)
          cols_in_strip <- ((t_pos - 1L) * cols_strip_s + 1L) :
                            (t_pos        * cols_strip_s)
          sub_grid <- expand.grid(Row = rows_in_rep, Column = cols_in_strip)
          sub_grid <- sub_grid[order(sub_grid$Row, sub_grid$Column), ]
          var_perm <- sample(vars_s)
          for (p in seq_len(n_s)) {
            v   <- var_perm[p]
            row <- sub_grid$Row[p]
            col <- sub_grid$Column[p]
            yld <- site_means[s]    +
                   treat_effects[t] +
                   U[v, u_col]      +
                   rep_eff[r]       +
                   row_eff[row]     +
                   col_eff[col]     +
                   rnorm(1L, 0, opts$error_sd)
            plots[[idx]] <- data.frame(
              Site      = sites[s],
              Treatment = treatments[t],
              TSite     = ts_label,
              Variety   = varieties[v],
              Rep       = r,
              Row       = row,
              Column    = col,
              yield     = round(yld, 1L),
              stringsAsFactors = FALSE
            )
            idx <- idx + 1L
          }
        }
      }
    } else {
      for (r in seq_len(nrep)) {
        rows_in_block <- ((r - 1L) * rows_block_s + 1L) : (r * rows_block_s)
        sub_grid <- expand.grid(Row    = rows_in_block,
                                Column = seq_len(cols_block_s))
        sub_grid <- sub_grid[order(sub_grid$Row, sub_grid$Column), ]
        var_perm <- sample(vars_s)
        for (p in seq_len(n_s)) {
          v   <- var_perm[p]
          row <- sub_grid$Row[p]
          col <- sub_grid$Column[p]
          yld <- site_means[s] +
                 U[v, s]       +
                 rep_eff[r]    +
                 row_eff[row]  +
                 col_eff[col]  +
                 rnorm(1L, 0, opts$error_sd)
          plots[[idx]] <- data.frame(
            Site    = sites[s],
            Variety = varieties[v],
            Rep     = r,
            Row     = row,
            Column  = col,
            yield   = round(yld, 1L),
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
  }

  dat <- do.call(rbind, plots[seq_len(idx - 1L)])

  # ---- Factor levels and ordering ------------------------------------------
  obs_vars    <- unique(dat$Variety)
  dat$Site    <- factor(dat$Site,    levels = sites)
  dat$Variety <- factor(dat$Variety, levels = varieties[varieties %in% obs_vars])
  dat$Rep     <- factor(dat$Rep,     levels = as.character(seq_len(nrep)))

  if (has_treatments) {
    dat$Treatment <- factor(dat$Treatment, levels = treatments)
    dat$TSite     <- factor(dat$TSite,     levels = tsite_levs)
  }

  dat <- dat[order(dat$Site, dat$Row, dat$Column), ]
  rownames(dat) <- NULL

  # ---- Optional CSV output -------------------------------------------------
  if (!is.null(opts$outfile)) {
    write.csv(dat, opts$outfile, row.names = FALSE)
    if (verbose) cat("CSV written to:", opts$outfile, "\n")
  }

  # ---- Verbose summary -----------------------------------------------------
  if (verbose) {
    R_off       <- cov2cor(G_sim)[upper.tri(G_sim)]
    vars_per_s  <- colSums(inc_mat)
    sites_per_v <- rowSums(inc_mat)
    inc_label   <- if (identical(incidence, "balanced"))   "balanced" else
                   if (identical(incidence, "unbalanced")) "unbalanced (auto)" else
                   "user-supplied matrix"
    G_label     <- if (use_G_mode == "auto")
                     sprintf("auto (cor [%.2f, %.2f])", opts$g_cor_min, opts$g_cor_max)
                   else "user-supplied matrix"
    k_sug       <- min(3L, max(2L, floor(ngroup / 3L)))

    cat("\n=== simTrialData summary ===\n")
    cat(sprintf("  Environments    : %d\n", nsite))
    if (has_treatments)
      cat(sprintf("  Treatments      : %d  ->  %s\n",
                  ntreat, paste(treatments, collapse = ", ")))
    cat(sprintf("  Varieties (max) : %d\n", nvar))
    cat(sprintf("  Replicates      : %d\n", nrep))
    cat(sprintf("  G structure     : %s\n", G_label))
    cat(sprintf("  Groups          : %d  (%s)\n", ngroup,
                if (has_treatments) "Treatment x Site" else "Site"))
    cat(sprintf("  Incidence       : %s\n", inc_label))
    cat(sprintf("  Vars/site       : mean = %.1f,  range [%d, %d]\n",
                mean(vars_per_s), min(vars_per_s), max(vars_per_s)))
    cat(sprintf("  Sites/variety   : mean = %.1f,  range [%d, %d]\n",
                mean(sites_per_v), min(sites_per_v), max(sites_per_v)))
    cat(sprintf("  Total plots     : %d\n", nrow(dat)))
    cat(sprintf("  True G_jj       : mean = %.1f,  range [%.1f, %.1f]\n",
                mean(diag(G_sim)), min(diag(G_sim)), max(diag(G_sim))))
    cat(sprintf("  Genetic corr    : mean = %.3f,  range [%.3f, %.3f]\n",
                mean(R_off), min(R_off), max(R_off)))
    cat("\n  Suggested ASReml-R model:\n")
    if (has_treatments) {
      cat(sprintf("    asreml(yield ~ TSite,\n"))
      cat(sprintf("           random   = ~ fa(TSite, %d):id(Variety) + at(Site):Rep,\n",
                  k_sug))
      cat(sprintf("           data     = dat)\n"))
      cat(sprintf("    randomRegress(model, Env = \"TSite:Variety\",\n"))
      cat(sprintf("                  levs = c(%s), sep = \"%s\")\n",
                  paste(sprintf('"%s"', treatments), collapse = ", "),
                  opts$sep))
    } else {
      cat(sprintf("    asreml(yield ~ Site,\n"))
      cat(sprintf("           random   = ~ fa(Site, %d):id(Variety) + at(Site):Rep,\n",
                  k_sug))
      cat(sprintf("           data     = dat)\n"))
      cat(sprintf("    accuracy(model)\n"))
    }
    cat("============================\n\n")
  }

  # ---- Return --------------------------------------------------------------
  params <- list(G          = G_sim,
                 site_means = site_means,
                 incidence  = inc_mat,
                 g_arr      = U)
  if (has_treatments) {
    params$treat_effects <- treat_effects
  }

  list(data = dat, params = params)
}
