#' Simulate Multi-Environment Plant Breeding Trial Data
#'
#' @description
#' Generates a balanced or unbalanced field trial dataset with a realistic
#' factor-analytic genetic correlation structure across environments.
#'
#' **MET-only** (`treatments = NULL`): a simple randomised complete block
#' design with one observation per observed Variety x Site x Rep combination.
#'
#' **Multi-treatment** (`treatments` is a character vector of length >= 2):
#' a split-plot design where whole plots are treatment strips and sub-plots
#' are varieties, nested within replicate blocks. The genetic structure
#' operates over all Treatment x Site (`TSite`) combinations.
#'
#' **Genetic simulation** uses a factor-analytic (FA) covariance model:
#' \deqn{
#'   \mathbf{u}_v = \boldsymbol{\Lambda}\mathbf{f}_v + \boldsymbol{\delta}_v,
#'   \quad \mathbf{f}_v \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_{k}),
#'   \quad \boldsymbol{\delta}_v \sim \mathcal{N}(\mathbf{0}, \boldsymbol{\Psi})
#' }
#' where \eqn{\boldsymbol{\Lambda}} (\eqn{J \times k}) is the loadings matrix,
#' \eqn{k} = `n_fa`, \eqn{J} is the number of environments (or
#' Treatment x Site combinations), and
#' \eqn{\boldsymbol{\Psi} = \mathrm{diag}(\psi_1,\ldots,\psi_J)} is the
#' specific variance matrix. The true genetic covariance matrix is
#' \eqn{\mathbf{G} = \boldsymbol{\Lambda}\boldsymbol{\Lambda}^\top +
#' \boldsymbol{\Psi}}.
#'
#' **Incidence structure** (`incidence` argument) controls which varieties are
#' observed at which sites. All `nvar` varieties have a true genetic value at
#' every site but only observed varieties contribute plots. Unbalanced designs
#' produce meaningful per-variety variation in prediction error variance —
#' varieties present in fewer sites are predicted less accurately.
#'
#' @param nvar        Integer. Maximum number of varieties. Default `20`.
#' @param nsite       Integer. Number of environments (sites). Default `10`.
#' @param treatments  Character vector of treatment labels (length >= 2), or
#'   `NULL` (default) for a MET-only trial with no treatment structure.
#' @param nrep        Integer. Replicates per environment. Default `2`.
#' @param n_fa        Integer. Number of FA factors used in the data-generating
#'   genetic model. Must be less than the number of groups (environments or
#'   Treatment x Site combinations). Default `2`.
#' @param incidence   Controls which varieties are observed at which sites.
#'   One of:
#'   \describe{
#'     \item{`"balanced"`}{(Default) All `nvar` varieties are observed at every
#'       site. Produces identical PEV for all varieties within a site.}
#'     \item{`"unbalanced"`}{Automatically generates a two-tier incidence
#'       structure. The proportions governing the structure are controlled by
#'       five \code{sim.options} entries: \code{core_pct},
#'       \code{core_min_sites_pct}, \code{reg_min_sites_pct},
#'       \code{reg_max_sites_pct}, and \code{min_vars_pct} (see below).}
#'     \item{Matrix}{An `nvar x nsite` matrix of `0`/`1` or `TRUE`/`FALSE`
#'       supplied by the user. Entry `[v, s] = 1` means variety `v` is
#'       observed at site `s`. Every variety must appear in at least 1 site;
#'       every site must have at least 2 varieties.}
#'   }
#' @param seed        Integer or `NULL`. Passed to [set.seed()]. Default `NULL`.
#' @param verbose     Logical. Print a design summary and suggested ASReml-R
#'   model. Default `TRUE`.
#' @param sim.options Named list of optional simulation controls. Unknown names
#'   produce a warning. Recognised elements (with defaults) are:
#'   \describe{
#'     \item{`site_mean`}{Numeric. Grand mean yield. Default `4500`.}
#'     \item{`site_sd`}{Numeric. SD of site mean yields. Default `600`.}
#'     \item{`treat_effects`}{Numeric vector (length = `length(treatments)`)
#'       of fixed treatment effects, or `NULL` for auto-spacing from 0 to
#'       `site_mean * 0.10`. Multi-treatment only.}
#'     \item{`sigma_genetic`}{Numeric. Target mean genetic SD per group.
#'       Default `250`.}
#'     \item{`loading_min`}{Numeric >= 0. Minimum absolute FA loading before
#'       scaling. Default `0.5`.}
#'     \item{`loading_max`}{Numeric > `loading_min`. Maximum absolute FA
#'       loading before scaling. Default `2.5`.}
#'     \item{`specific_pct`}{Numeric in (0, 1). Specific variances as a
#'       fraction of mean diagonal of Lambda Lambda'. Default `0.20`.}
#'     \item{`rep_sd`}{Numeric. SD of replicate effects. Default `150`.}
#'     \item{`row_sd`}{Numeric. SD of row spatial effects. Default `80`.}
#'     \item{`col_sd`}{Numeric. SD of column spatial effects. Default `60`.}
#'     \item{`error_sd`}{Numeric. SD of plot-level residual. Default `350`.}
#'     \item{`sep`}{Character. Separator for `TSite` labels. Default `"-"`.}
#'     \item{`variety_prefix`}{Character. Prefix for variety labels.
#'       Default `"Var"`.}
#'     \item{`site_prefix`}{Character. Prefix for site labels.
#'       Default `"Env"`.}
  #'     \item{`outfile`}{Character path or `NULL`. CSV output. Default `NULL`.}
#'     \item{`core_pct`}{Numeric in (0, 1). Proportion of varieties designated
#'       as "core" entries when \code{incidence = "unbalanced"}. Core varieties
#'       appear in at least \code{core_min_sites_pct} of sites. Default
#'       \code{0.20}.}
#'     \item{`core_min_sites_pct`}{Numeric in (0, 1]. Minimum proportion of
#'       sites that each core variety must appear in. Default \code{0.75}.}
#'     \item{`reg_min_sites_pct`}{Numeric in (0, 1]. Minimum proportion of
#'       sites for regular (non-core) varieties. Default \code{0.40}.}
#'     \item{`reg_max_sites_pct`}{Numeric in (0, 1]. Maximum proportion of
#'       sites for regular varieties. Must be >= \code{reg_min_sites_pct}.
#'       Default \code{0.85}.}
#'     \item{`min_vars_pct`}{Numeric in (0, 1]. Minimum proportion of
#'       \code{nvar} that every site must contain. Default \code{0.40}.}
#'   }
#'
#' @return A named list:
#' \describe{
#'   \item{`data`}{Data frame ordered by `Site`, `Row`, `Column` with columns:
#'     `Site`, `Variety`, `Rep`, `Row`, `Column`, `yield` (MET-only), or
#'     additionally `Treatment` and `TSite` (multi-treatment). Only observed
#'     variety x site combinations are included.}
#'   \item{`params`}{Named list of true simulation parameters:
#'     \describe{
#'       \item{`G`}{True genetic covariance matrix (J x J).}
#'       \item{`Lambda`}{True loadings matrix (J x k).}
#'       \item{`Psi`}{True specific variances (length J).}
#'       \item{`site_means`}{True site mean yields (length `nsite`).}
#'       \item{`n_fa`}{Number of FA factors used.}
#'       \item{`incidence`}{Integer `nvar x nsite` incidence matrix used.}
#'       \item{`treat_effects`}{True treatment fixed effects
#'         (multi-treatment only).}
#'       \item{`g_arr`}{`nvar x ngroup` matrix of true genetic BLUPs
#'         (multi-treatment only).}
#'     }}
#' }
#'
#' @examples
#' \dontrun{
#' # MET-only balanced: 30 varieties, 8 sites, FA(2) genetic structure
#' out <- simTrialData(nvar = 30, nsite = 8, n_fa = 2, seed = 42)
#' round(cov2cor(out$params$G), 2)   # true genetic correlations
#'
#' # Unbalanced: produces spread in per-variety accuracy
#' out2 <- simTrialData(nvar = 30, nsite = 8, n_fa = 2,
#'                      incidence = "unbalanced", seed = 42)
#' table(colSums(out2$params$incidence))   # varieties per site
#'
#' # User-supplied incidence matrix
#' inc <- matrix(1L, nrow = 20, ncol = 6)
#' inc[1:5, 1:3] <- 0L   # first 5 varieties absent from first 3 sites
#' out3 <- simTrialData(nvar = 20, nsite = 6, n_fa = 2,
#'                      incidence = inc, seed = 1)
#'
#' # Multi-treatment with custom treat_effects
#' out4 <- simTrialData(nvar = 20, nsite = 6,
#'                      treatments  = c("T0", "T1", "T2"),
#'                      n_fa        = 2,
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
                         n_fa        = 2L,
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
    loading_min         = 0.5,
    loading_max         = 2.5,
    specific_pct        = 0.20,
    rep_sd              = 150,
    row_sd              = 80,
    col_sd              = 60,
    error_sd            = 350,
    sep                 = "-",
    variety_prefix      = "Var",
    site_prefix         = "Env",
    outfile             = NULL,
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

  if (n_fa < 1L)
    stop("'n_fa' must be at least 1.")
  if (n_fa >= ngroup)
    stop("'n_fa' (", n_fa, ") must be less than the number of groups (",
         ngroup, ").")

  if (!is.null(opts$treat_effects) && length(opts$treat_effects) != ntreat)
    stop("'treat_effects' must have length ", ntreat, ".")

  # Validate unbalanced incidence proportions (only meaningful when used)
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

  # Auto-generate a two-tier unbalanced incidence matrix.
  # Proportions are controlled by sim.options entries passed via opts.
  .gen_unbalanced_inc <- function(nv, ns, opts) {
    n_core  <- max(1L, ceiling(nv * opts$core_pct))
    mat     <- matrix(0L, nrow = nv, ncol = ns)

    # Core varieties: appear in >= core_min_sites_pct of sites
    n_core_sites <- max(2L, ceiling(ns * opts$core_min_sites_pct))
    for (v in seq_len(n_core))
      mat[v, sample.int(ns, min(n_core_sites, ns))] <- 1L

    # Regular varieties: appear in reg_min_sites_pct to reg_max_sites_pct of sites
    lo <- max(2L, ceiling(ns * opts$reg_min_sites_pct))
    hi <- max(lo, floor(ns * opts$reg_max_sites_pct))
    for (v in seq(n_core + 1L, nv)) {
      k <- sample.int(hi - lo + 1L, 1L) + lo - 1L
      mat[v, sample.int(ns, k)] <- 1L
    }

    # Ensure every site meets minimum variety count (min_vars_pct * nvar)
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

  # Validate a user-supplied incidence matrix.
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

  # Resolve incidence argument
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

  # ---- FA genetic simulation -----------------------------------------------
  Lambda_raw <- matrix(
    runif(ngroup * n_fa, min = opts$loading_min, max = opts$loading_max),
    nrow = ngroup, ncol = n_fa
  )

  if (n_fa >= 2L) {
    for (f in 2L:n_fa) {
      flip <- sample(ngroup, floor(ngroup * 0.4))
      Lambda_raw[flip, f] <- -Lambda_raw[flip, f]
    }
  }

  G_unscaled <- Lambda_raw %*% t(Lambda_raw)
  psi_raw    <- rep(mean(diag(G_unscaled)) * opts$specific_pct, ngroup)
  G_raw      <- G_unscaled + diag(psi_raw)

  scale_f <- opts$sigma_genetic^2 / mean(diag(G_raw))
  Lambda  <- Lambda_raw * sqrt(scale_f)
  Psi     <- psi_raw    * scale_f
  G       <- G_raw      * scale_f

  rownames(G)      <- colnames(G) <- grp_names
  rownames(Lambda) <- grp_names
  names(Psi)       <- grp_names

  # All nvar varieties get true BLUPs at all groups (even unobserved sites)
  F_mat <- matrix(rnorm(n_fa * nvar), nrow = n_fa, ncol = nvar)
  U_fa  <- t(Lambda %*% F_mat)
  U_psi <- matrix(
    rnorm(nvar * ngroup, 0, rep(sqrt(Psi), each = nvar)),
    nrow = nvar, ncol = ngroup
  )
  U <- U_fa + U_psi
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
  # Pre-allocate for the maximum possible number of plots
  plots <- vector("list", nvar * nsite * nrep * ntreat)
  idx   <- 1L

  for (s in seq_len(nsite)) {

    # Varieties present at this site
    vars_s <- which(inc_mat[, s] == 1L)
    n_s    <- length(vars_s)
    if (n_s == 0L) next

    # Per-site layout dimensions based on actual variety count
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
          var_perm <- sample(vars_s)   # only varieties at this site

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
        var_perm <- sample(vars_s)   # only varieties at this site

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
  # Drop varieties that appear in no sites (possible with user matrix)
  obs_vars <- unique(dat$Variety)
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
    R_off      <- cov2cor(G)[upper.tri(G)]
    vars_per_s <- colSums(inc_mat)
    sites_per_v <- rowSums(inc_mat)
    inc_label  <- if (identical(incidence, "balanced"))   "balanced" else
                  if (identical(incidence, "unbalanced")) "unbalanced (auto)" else
                  "user-supplied matrix"

    cat("\n=== simTrialData summary ===\n")
    cat(sprintf("  Environments    : %d\n", nsite))
    if (has_treatments)
      cat(sprintf("  Treatments      : %d  ->  %s\n",
                  ntreat, paste(treatments, collapse = ", ")))
    cat(sprintf("  Varieties (max) : %d\n", nvar))
    cat(sprintf("  Replicates      : %d\n", nrep))
    cat(sprintf("  FA factors (k)  : %d\n", n_fa))
    cat(sprintf("  Groups          : %d  (%s)\n", ngroup,
                if (has_treatments) "Treatment x Site" else "Site"))
    cat(sprintf("  Incidence       : %s\n", inc_label))
    cat(sprintf("  Vars/site       : mean = %.1f,  range [%d, %d]\n",
                mean(vars_per_s), min(vars_per_s), max(vars_per_s)))
    cat(sprintf("  Sites/variety   : mean = %.1f,  range [%d, %d]\n",
                mean(sites_per_v), min(sites_per_v), max(sites_per_v)))
    cat(sprintf("  Total plots     : %d\n", nrow(dat)))
    cat(sprintf("  True G_jj       : mean = %.1f,  range [%.1f, %.1f]\n",
                mean(diag(G)), min(diag(G)), max(diag(G))))
    cat(sprintf("  Genetic corr    : mean = %.3f,  range [%.3f, %.3f]\n",
                mean(R_off), min(R_off), max(R_off)))
    k_sug <- min(3L, n_fa + 1L)
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
  params <- list(G          = G,
                 Lambda     = Lambda,
                 Psi        = Psi,
                 site_means = site_means,
                 n_fa       = n_fa,
                 incidence  = inc_mat)
  if (has_treatments) {
    params$treat_effects <- treat_effects
    params$g_arr         <- U
  }

  list(data = dat, params = params)
}
