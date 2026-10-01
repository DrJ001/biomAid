# ---- Private helper: build the conditioning structure -------------------

#' @noRd
.condList <- function(levs, type, cond) {

  nk <- length(levs)
  cl <- setNames(vector("list", nk), levs)   # all NULL by default

  if (type == "baseline") {
    ## All non-first treatments conditioned on levs[1] only
    for (j in seq(2L, nk))
      cl[[j]] <- levs[1L]

  } else if (type == "sequential") {
    ## Treatment j conditioned on all preceding treatments levs[1:(j-1)]
    ## This is the Gram-Schmidt / LDL' Cholesky decomposition;
    ## TGmat = T G T' is *diagonal* (fully orthogonal components).
    for (j in seq(2L, nk))
      cl[[j]] <- levs[seq_len(j - 1L)]

  } else if (type == "partial") {
    ## Each treatment conditioned on *all other* treatments simultaneously.
    ## Diagonal of TGmat gives the partial genetic variances.
    for (j in seq_len(nk))
      cl[[j]] <- levs[-j]

  } else {    # "custom"
    if (is.null(cond))
      stop("'cond' must be supplied when type = \"custom\".")
    if (!is.list(cond) || is.null(names(cond)))
      stop("'cond' must be a named list with names matching elements of 'levs'.")
    if (length(bad <- setdiff(names(cond), levs)))
      stop("Names in 'cond' not found in 'levs': ", paste(bad, collapse = ", "))
    for (lv in names(cond)) {
      if (!is.null(cond[[lv]])) {
        if (length(bad2 <- setdiff(cond[[lv]], levs)))
          stop("Conditioning levels for '", lv, "' not in 'levs': ",
               paste(bad2, collapse = ", "))
        if (lv %in% cond[[lv]])
          stop("Level '", lv, "' cannot be in its own conditioning set.")
      }
      cl[[lv]] <- cond[[lv]]
    }
  }
  cl
}


# ---- Private helper: parse the random term string -----------------------

#' Parse an ASReml-R random interaction term for randomRegress()
#'
#' Accepts the full left-hand interaction term as the user supplies it, e.g.:
#'   "us(TSite):Variety"
#'   "fa(TSite, 2):Variety"
#'   "corgh(TSite):Variety"
#'   "corh(TSite):Variety"
#'   "diag(TSite):Variety"
#'   "us(TSite):vm(Variety, giv1)"
#'   "us(TSite):ide(Variety)"
#'
#' Returns a list with:
#'   struct      – variance structure keyword (fa / us / corgh / corh / diag)
#'   group_var   – bare group-factor name (e.g. "TSite")
#'   by_var      – bare variety-factor name (e.g. "Variety"), wrappers stripped
#'   by_wrapper  – wrapper on the by-variable if present (e.g. "vm", "ide", or NULL)
#'   by_raw      – full right-hand side as supplied (e.g. "vm(Variety, giv1)")
#'   n_fa        – integer FA order, or NULL for non-FA structures
#'   only_term   – string to pass to predict(only=...)
#'   classify    – string to pass to predict(classify=...)
#'
#' @noRd
.parse_rreg_term <- function(term) {

  # Locate the first ":" at parenthesis depth 0 *before* stripping whitespace,
  # so that we can preserve the original right-hand side for only_term.
  # ASReml-R's predict(only=) must match the term label exactly as stored in
  # the model's random formula — including spaces inside wrappers such as
  # vm(Variety, giv1).  Stripping all whitespace first causes only_term to
  # contain "vm(Variety,giv1)" (no space), which ASReml-R cannot match,
  # making vm() and ide() indistinguishable.
  cp_orig <- .top_colon(term)
  if (is.na(cp_orig))
    stop("'term' must be an interaction of the form 'struct(Group):Variety' or 'struct(Group):wrapper(Variety,...)'.")

  rgt_orig <- substr(term, cp_orig + 1L, nchar(term))   # rhs with original spacing

  term <- gsub("\\s+", "", term)   # strip all whitespace (for parsing only)

  # Locate the colon again in the stripped string (position may shift)
  cp <- .top_colon(term)

  lft <- substr(term, 1L,       cp - 1L)
  rgt <- substr(term, cp + 1L,  nchar(term))

  # ---- Left-hand side: struct(Group [, k]) --------------------------------
  struct <- sub("\\(.*", "", lft)
  if (!struct %in% c("fa", "us", "corgh", "corh", "diag"))
    stop("Unsupported variance structure '", struct, "'. ",
         "Supported: fa, us, corgh, corh, diag.")

  # group variable: first comma-separated token inside lft parens
  group_var <- trimws(strsplit(.inside(lft), ",")[[1L]][1L])

  # FA order (only meaningful for FA terms)
  n_fa <- if (struct == "fa") {
    suppressWarnings(
      as.integer(trimws(strsplit(.inside(lft), ",")[[1L]][2L]))
    )
  } else NULL

  # ---- Right-hand side: Variety | vm(Variety,...) | ide(Variety) ----------
  # by_raw preserves the original spacing so only_term matches ASReml's stored
  # term label exactly (e.g. "vm(Variety, giv1)" not "vm(Variety,giv1)").
  by_raw     <- rgt_orig
  by_wrapper <- NULL

  # Check for a wrapping function: any "word(" pattern (use stripped rgt for
  # parsing — positional logic is simpler on whitespace-free strings)
  if (grepl("^[A-Za-z][A-Za-z0-9_]*\\(", rgt)) {
    by_wrapper <- sub("\\(.*", "", rgt)
    by_var     <- trimws(strsplit(.inside(rgt), ",")[[1L]][1L])
  } else {
    by_var <- rgt
  }

  # ---- Build predict() strings -------------------------------------------
  # classify always uses bare variable names (no wrappers).
  # only uses the bare group name (or full fa() spec) on the left, but keeps
  # any vm() / ide() wrapper on the right — ASReml-R requires it when present.
  # For FA:   only = "fa(Group, k):rhs"  (space after comma required by ASReml)
  # For rest: only = "Group:rhs"
  # IMPORTANT: by_raw retains original spacing so the string matches the term
  # label in the model's random formula exactly.
  only_term <- if (struct == "fa" && !is.null(n_fa))
    sprintf("fa(%s, %d):%s", group_var, n_fa, by_raw)
  else
    paste0(group_var, ":", by_raw)

  classify <- paste0(group_var, ":", by_var)

  list(
    struct     = struct,
    group_var  = group_var,
    by_var     = by_var,
    by_wrapper = by_wrapper,
    by_raw     = by_raw,
    n_fa       = n_fa,
    only_term  = only_term,
    classify   = classify
  )
}


# ---- Private helper: resolve composite group labels into level + section -

#' Resolve G-matrix column labels into decomposed level and section
#'
#' The G-matrix of a term such as `us(TSite):Variety` is indexed by the levels
#' of the grouping factor, and when that factor is a composite the two-way
#' structure survives only in the label strings (e.g. `"N0-Env1"`).  This
#' helper recovers it once, authoritatively, so that neither randomRegress()
#' nor plot_randomRegress() has to split strings ad hoc.
#'
#' Splitting is attempted at the **first** and at the **last** separator, and
#' the candidate retained is the one whose level part contains every member of
#' `levs`.  This is what makes separator-bearing section names such as
#' `"N0-North-West"` safe, in either label order.  If both candidates qualify
#' the labels are genuinely ambiguous and an error is raised rather than a
#' silent choice being made.
#'
#' @param tsnams Character vector of G-matrix column labels.
#' @param levs   Character vector of levels to decompose.
#' @param sep    Separator for composite labels.
#' @param enam   Grouping factor name, used only in error messages.
#' @return List with `level`, `section` (both parallel to `tsnams`),
#'   `composite` (logical) and `level_side` (1 or 2; `NA` when not composite).
#' @noRd
.rreg_labels <- function(tsnams, levs, sep, enam = "the grouping factor") {

  has_sep <- grepl(sep, tsnams, fixed = TRUE)

  # ---- Plain grouping factor: labels are the levels, one section ---------
  if (!any(has_sep)) {
    if (!all(levs %in% tsnams))
      stop("Levels supplied in 'levs' do not exist in ", enam, ": ",
           paste(setdiff(levs, tsnams), collapse = ", "), ".")
    return(list(level      = tsnams,
                section    = rep("Single", length(tsnams)),
                composite  = FALSE,
                level_side = NA_integer_))
  }

  if (!all(has_sep))
    stop("Some levels of ", enam, " contain the separator '", sep,
         "' and some do not, so they cannot be split consistently: ",
         paste(utils::head(tsnams[!has_sep], 3L), collapse = ", "),
         ". Check the 'sep' argument.")

  # ---- Candidate splits: at the first, and at the last, separator -------
  # Done with literal position arithmetic rather than regex, so that a
  # separator containing regex metacharacters (".", "|", "+", ...) needs no
  # escaping and cannot change the meaning of the pattern.
  nsep <- nchar(sep)
  at   <- lapply(gregexpr(sep, tsnams, fixed = TRUE), as.integer)
  first_at <- vapply(at, function(p) p[1L],         integer(1L))
  last_at  <- vapply(at, function(p) p[length(p)],  integer(1L))

  pre1  <- substr(tsnams, 1L, first_at - 1L)
  post1 <- substr(tsnams, first_at + nsep, nchar(tsnams))
  preL  <- substr(tsnams, 1L, last_at - 1L)
  postL <- substr(tsnams, last_at + nsep, nchar(tsnams))

  # A: level on the left  (level itself free of sep)
  # B: level on the right (level itself free of sep)
  cand <- list(
    list(level = pre1,  section = post1, side = 1L),
    list(level = postL, section = preL,  side = 2L)
  )
  ok <- vapply(cand, function(cd) all(levs %in% unique(cd$level)), logical(1L))

  if (!any(ok))
    stop("Levels supplied in 'levs' do not exist in ", enam,
         " on either side of the separator '", sep, "'. Supplied: ",
         paste(levs, collapse = ", "), ".")

  if (all(ok) && !identical(cand[[1L]]$level, cand[[2L]]$level))
    stop("The labels of ", enam, " are ambiguous: every level in 'levs' ",
         "appears on both sides of the separator '", sep, "', so the ",
         "decomposed dimension cannot be identified. Rename the levels of ",
         enam, " so that only one side carries the levels in 'levs'.")

  res <- cand[[which(ok)[1L]]]

  # Each level x section combination must be unique, or downstream lookups
  # silently resolve to the first match and whole sections disappear.
  key <- paste(res$level, res$section, sep = "\r")
  if (anyDuplicated(key)) {
    dup <- unique(key[duplicated(key)])
    stop("The labels of ", enam, " do not form a unique level-by-section ",
         "crossing after splitting on '", sep, "'. Duplicated: ",
         paste(sub("\r", " / ", utils::head(dup, 3L)), collapse = "; "),
         ". Check the 'sep' argument.")
  }

  list(level      = res$level,
       section    = res$section,
       composite  = TRUE,
       level_side = res$side)
}


# ---- Private helper: convert corh/corgh vparameters matrix to G-matrix --

#' Convert a heterogeneous-correlation vparameters matrix to a covariance G-matrix
#'
#' ASReml-R V4 stores `corh` and `corgh` variance parameters in a matrix
#' accessed via `summary(model, vparameters = TRUE)$vparameters[["Group:Variety"]]`
#' (using the bare stripped term, e.g. `"TSite:Variety"`).  This matrix has:
#'   - Genetic variances (sigma^2_j) on the diagonal
#'   - Correlations (r_ij) on the off-diagonal
#'
#' To use it as a proper covariance G-matrix in the regression we need:
#'   G_ij = r_ij * sigma_i * sigma_j    (off-diagonal)
#'   G_jj = sigma^2_j                   (diagonal, unchanged)
#'
#' @param M  Square matrix from vparameters with variances on diagonal and
#'   correlations on off-diagonal.
#' @return   Symmetric covariance matrix of the same dimensions.
#'
#' @noRd
.cor_to_cov_Gmat <- function(M) {
  sds    <- sqrt(diag(M))          # sigma_j for each group level
  R      <- M
  diag(R) <- 1                     # pure correlation matrix (1s on diagonal)
  G      <- diag(sds) %*% R %*% diag(sds)   # G = diag(sigma) R diag(sigma)
  dimnames(G) <- dimnames(M)       # restore row/col names lost by %*%
  G
}


# ---- Main function -------------------------------------------------------

#' Multivariate Random Regression of Variety BLUPs
#'
#' @description
#' Uses the G-matrix from an ASReml-R V4 model to decompose a multivariate set
#' of variety BLUPs into **efficiency** and **responsiveness** components,
#' supporting four conditioning schemes via the `type` argument.
#'
#' The decomposed dimension holds multiple **treatments** or multiple
#' **traits** — levels applied to, or measured on, the same plants within one
#' experiment — nominated in `levs`.  The decomposition is repeated
#' independently within each **section**, ordinarily a site, of which there may
#' be one or many.  Multi-treatment multi-environment data is the most complex
#' case, not the only one; see **Usage regimes** below.
#'
#' For any level \eqn{j} with conditioning set \eqn{A_j}, the multivariate
#' conditional normal distribution gives:
#'
#' \deqn{
#'   \boldsymbol{\beta}_j = \boldsymbol{G}_{A_j A_j}^{-1}\, \boldsymbol{G}_{A_j j}
#'   \qquad
#'   \tilde{u}_j = u_j - \boldsymbol{\beta}_j^\top \boldsymbol{u}_{A_j}
#' }
#'
#' The four built-in conditioning schemes differ only in how \eqn{A_j} is
#' chosen for each level:
#'
#' \describe{
#'   \item{`"baseline"` (default)}{Every non-first level is conditioned on
#'     \code{levs[1]} alone.  Responsiveness BLUPs are orthogonal to the
#'     baseline but may be correlated with each other.  The transformed
#'     G-matrix `TGmat` is block-diagonal.}
#'   \item{`"sequential"`}{Level \eqn{j} is conditioned on all preceding
#'     levels \code{levs[1:(j-1)]}.  This is the Gram-Schmidt
#'     orthogonalisation of the BLUPs, equivalent to the \eqn{LDL^\top}
#'     Cholesky decomposition of the G-matrix.  All components are mutually
#'     orthogonal and `TGmat` is **diagonal**, with the Schur complements on
#'     the diagonal.  The ordering of `levs` matters.}
#'   \item{`"partial"`}{Each level is conditioned on **all other** levels
#'     simultaneously.  The diagonal of `TGmat` gives the partial genetic
#'     variances; off-diagonals are generally non-zero.}
#'   \item{`"custom"`}{The conditioning set for each level is specified
#'     explicitly via the `cond` argument.}
#' }
#'
#' @section Usage regimes:
#' Two dimensions are involved, playing different roles.  The **decomposed**
#' dimension holds the treatments or traits named in `levs`; the **section**
#' dimension is the site or environment the decomposition is repeated within.
#' Which regime applies is determined automatically from the G-matrix column
#' labels.
#'
#' \describe{
#'   \item{Plain grouping factor — one section}{When the labels contain no
#'     `sep` the labels themselves are the decomposed dimension, and a single
#'     section is reported as `"Single"`.  This covers a single-site
#'     multi-treatment model, `us(Treatment):Variety`, and a single-site
#'     multi-trait model, `us(Trait):Variety`.}
#'   \item{Composite grouping factor — several sections}{When the labels contain
#'     `sep` — e.g. a Treatment-by-Site factor `TSite` with levels
#'     `"N0-Env1"` — they are split in two.  The component holding the `levs`
#'     values is the decomposed dimension; the other becomes the section, and
#'     the decomposition is carried out independently within each.  The two
#'     components may appear in either order: both `"N0-Env1"` and `"Env1-N0"`
#'     are recognised.}
#' }
#'
#' Only treatments and traits are decomposable: they are commensurable levels
#' observed on common material within one experiment, so the conditional
#' distribution of one given another is biologically interpretable.
#' **Environments are not.**  Regressing one site's BLUPs on another's would
#' not give an efficiency-responsiveness decomposition, because separate sites
#' are separate experiments; environments belong in the section role.  Where
#' the genetic covariance *between* environments is the question of interest,
#' use [faSummary()] or [fastIC()] instead.
#'
#' Throughout the documentation and output, *level* refers to a level of the
#' decomposed dimension (a member of `levs`) and *section* refers to a level
#' of the second, repeated-over dimension.
#'
#' The word *section* is borrowed from ASReml-R, which uses it for a factor
#' across which a structure is replicated independently — as in
#' `residual = ~ dsum(~ ar1(Row):ar1(Column) | Site)`.  Strictly, ASReml-R's
#' sections partition the **residual** structure whereas these partition the
#' second dimension of the **G** structure; in a typical MET both are `Site`
#' and they coincide, but they are not the same thing by definition.
#' *Stratum* is deliberately avoided, since in experimental design it already
#' denotes an error stratum of the blocking structure.
#'
#' @param model An ASReml-R V4 model object containing a random term that
#'   crosses a grouping factor with the variety factor.  The grouping factor
#'   may index treatments or traits, or be a composite of treatments or traits
#'   with a site factor.
#' @param term Character string giving the **full** random-effect interaction
#'   term exactly as it appears in the model formula, written as
#'   `"<struct>(<Group>):<Variety>"`.  The function parses the structure
#'   keyword, group factor name, and variety factor name automatically, so
#'   the correct strings are used for \code{predict.asreml()}.
#'
#'   Supported variance structures on the left-hand side:
#'   \describe{
#'     \item{`us(TSite)`}{Unstructured G-matrix — the most general form.}
#'     \item{`fa(TSite, k)`}{Factor-analytic of order \eqn{k}.}
#'     \item{`corgh(TSite)`}{Heterogeneous correlation — one correlation
#'       parameter shared across groups with group-specific variances.}
#'     \item{`corh(TSite)`}{Correlation and variance structure for two-group
#'       (two-level) models.}
#'     \item{`diag(TSite)`}{Diagonal — independent genetic variances per group,
#'       zero between-group covariances.}
#'   }
#'
#'   Supported wrappers on the right-hand side (variety factor):
#'   \describe{
#'     \item{`vm(Variety, giv1)`}{Genomic or pedigree relationship matrix via
#'       \code{asreml::vm()}.  Only the bare factor name is used for
#'       \code{predict.asreml()}.}
#'     \item{`ide(Variety)`}{Identity-scaled term via \code{asreml::ide()}.}
#'     \item{`Variety` (no wrapper)}{Plain factor — the default.}
#'   }
#'
#'   Examples:
#'   \preformatted{
#'   term = "us(TSite):Variety"           # treatments within environments
#'   term = "fa(TSite, 2):Variety"
#'   term = "corgh(TSite):Variety"
#'   term = "corh(TSite):Variety"
#'   term = "us(TSite):vm(Variety, giv1)"
#'   term = "us(TSite):ide(Variety)"
#'   term = "us(Treatment):Variety"       # treatments, single site
#'   term = "us(Trait):Variety"           # traits, single site
#'   }
#' @param levs Character vector of length \eqn{\ge 2} naming the levels to
#'   decompose — treatment labels or trait names, depending on what the
#'   grouping factor indexes.  For `type = "baseline"` and
#'   `type = "sequential"` the **first** element is the baseline (efficiency)
#'   level.  For `type = "partial"` the ordering does not affect results.  For
#'   `type = "custom"` the ordering determines which element of `cond` applies
#'   to which level.
#' @param type Character string selecting the conditioning scheme. One of
#'   `"baseline"` (default), `"sequential"`, `"partial"`, or `"custom"`.
#'   See **Description** for full details of each scheme.
#' @param cond Named list required when `type = "custom"`.  Each element name
#'   must be a level from `levs`; each element value is either `NULL` (the
#'   level is unconditional / efficiency) or a character vector of levels from
#'   `levs` that form the conditioning set.  Levels absent from `cond` are
#'   treated as unconditional.  Example for a three-level sequential-style
#'   custom scheme:
#'   \preformatted{
#'   cond = list(T0 = NULL,
#'               T1 = "T0",
#'               T2 = c("T0", "T1"))
#'   }
#' @param sep Character separator used inside the labels of a composite
#'   grouping factor, e.g. `"-"` in `"N0-Env1"`.  Ignored when the grouping
#'   factor's labels contain no separator, in which case a single section
#'   named `"Single"` is reported.  Defaults to `"-"`.
#' @param pev Logical.  If `TRUE` (default) the variance used for HSD
#'   computation is the prediction error variance (PEV) of each responsiveness
#'   BLUP.  If `FALSE` it is the posterior variance
#'   \eqn{\sigma_{j|A_j}^2 - \text{PEV}}.  Ignored for FA models
#'   (HSD is always `NA`).
#' @param ... Additional arguments forwarded to `asreml::predict.asreml()`
#'   (non-FA terms only).
#'
#' @return A named list:
#' \describe{
#'   \item{`blups`}{Data frame with columns: `Site`, `Variety`, one raw BLUP
#'     column per level in `levs`, one `resp.<lev>` column per conditioned
#'     level, and one `HSD.<lev>` column per conditioned level (Tukey's HSD on
#'     the responsiveness scale; `NA` for FA models or absent combinations).
#'     The `Site` column holds the **section** label whatever the stratifying
#'     dimension represents, and is `"Single"` throughout when the grouping
#'     factor is not composite.  The `Variety` column holds the levels of the
#'     variety factor named in `term`, whatever that factor is called in the
#'     model.}
#'   \item{`TGmat`}{Transformed G-matrix \eqn{\boldsymbol{T}\boldsymbol{G}
#'     \boldsymbol{T}^\top}.  Unconditional levels are labelled `eff.<lev>`;
#'     conditioned levels are labelled `resp.<lev>`.  Diagonal for
#'     `type = "sequential"`.}
#'   \item{`Gmat`}{Original G-matrix.}
#'   \item{`beta`}{Named list of length equal to the number of conditioned
#'     levels.  Each element `beta[["<lev>"]]` is an
#'     \eqn{n_s \times |A_j|} matrix of per-section regression coefficients,
#'     with column names equal to the conditioning levels.  `NA` where a
#'     combination is absent from a section.}
#'   \item{`sigmat`}{Numeric matrix of dimensions \eqn{n_s \times n_{\text{cond}}}
#'     containing the scalar conditional genetic variances
#'     \eqn{\sigma_{j|A_j}^2} for each conditioned level in each section.
#'     `NA` where absent.}
#'   \item{`tmat`}{Full transformation matrix \eqn{\boldsymbol{T}}.
#'     Lower-triangular for `type = "sequential"`; sparse (one non-trivial
#'     column per section) for `type = "baseline"`; dense for
#'     `type = "partial"`.}
#'   \item{`cond_list`}{The resolved conditioning structure as a named list,
#'     one element per level in `levs`.}
#'   \item{`type`}{The `type` argument used.}
#'   \item{`sep`}{The `sep` argument used.}
#'   \item{`label_map`}{Data frame with one row per G-matrix column —
#'     `label` (the grouping-factor level as ASReml-R stores it), `level` (the
#'     decomposed level) and `section`.  Composite labels are resolved once
#'     here and this mapping is the authority downstream, so
#'     [plot_randomRegress()] never re-splits label strings.}
#' }
#'
#' @seealso [plot_randomRegress()], [fixedRegress()], [faSummary()],
#'   `asreml::predict.asreml()`
#'
#' @examples
#' \dontrun{
#' ## ---- Treatments within environments (composite grouping factor) ------
#' ## Baseline scheme — unstructured G-matrix
#' res_base <- randomRegress(model, term = "us(TSite):Variety",
#'                           levs = c("N0","N1","N2"))
#'
#' ## Sequential (Cholesky) — fully orthogonal components; TGmat is diagonal
#' res_seq  <- randomRegress(model, term = "fa(TSite, 2):Variety",
#'                           levs = c("N0","N1","N2"), type = "sequential")
#'
#' ## Heterogeneous correlation structure
#' res_cor  <- randomRegress(model, term = "corgh(TSite):Variety",
#'                           levs = c("N0","N1","N2"))
#'
#' ## Genomic relationship matrix (vm wrapper on Variety)
#' res_vm   <- randomRegress(model, term = "us(TSite):vm(Variety, giv1)",
#'                           levs = c("N0","N1","N2"))
#'
#' ## Custom conditioning
#' res_cust <- randomRegress(model, term = "us(TSite):Variety",
#'                           levs = c("N0","N1","N2"),
#'                           type = "custom",
#'                           cond = list(N0 = NULL,
#'                                       N1 = "N0",
#'                                       N2 = c("N0","N1")))
#'
#' ## ---- Single site (plain grouping factor; one section, "Single") -------
#' ## Treatments
#' res_1s   <- randomRegress(model, term = "us(Treatment):Variety",
#'                           levs = c("N0","N1","N2"))
#'
#' ## Traits, each conditioned on all others
#' res_mv   <- randomRegress(model, term = "us(Trait):Variety",
#'                           levs = c("Yield","Protein","Height"),
#'                           type = "partial")
#' }
#'
#' @export
randomRegress <- function(model, term = "us(TSite):Variety", levs = NULL,
                           type = "baseline", cond = NULL,
                           sep = "-", pev = TRUE, ...) {

  # ---- Validate and build conditioning structure -------------------------
  if (is.null(levs) || length(levs) < 2L)
    stop("At least two levels must be supplied in 'levs'.")

  ntreat <- length(levs)
  type   <- match.arg(type, c("baseline", "sequential", "partial", "custom"))

  cond_list   <- .condList(levs, type, cond)
  conditioned <- levs[!vapply(cond_list, is.null, logical(1L))]
  n_cond      <- length(conditioned)
  if (n_cond == 0L)
    stop("No levels have a conditioning set. Check 'type' or 'cond'.")

  # ---- Parse term string -------------------------------------------------
  p     <- .parse_rreg_term(term)
  enam  <- p$group_var    # e.g. "TSite"
  vnam  <- p$by_var       # e.g. "Variety"  (bare, wrappers stripped)
  struct <- p$struct      # e.g. "us", "fa", "corgh", "corh", "diag"

  # Match the term in the model's random formula using the parsed components.
  # We match on struct(Group) first, then verify the RHS wrapper agrees — this
  # prevents silently using the wrong model term when the user supplies, e.g.,
  # "us(TSite):ide(Variety)" against a model fitted with "us(TSite):vm(Variety, giv1)".
  # Without the RHS check both calls pass validation, use the same G-matrix key,
  # and ASReml silently falls back to the same predictions, making vm() and ide()
  # look identical.
  formula_terms <- attr(terms(model$formulae$random), "term.labels")
  rterm <- grep(
    paste0(struct, "\\(", enam),
    formula_terms,
    value = TRUE
  )
  if (length(rterm) == 0L)
    stop("Cannot find a term matching '", term, "' in the model's random formula.")
  rterm <- rterm[1L]

  # Validate that the RHS wrapper in the supplied term matches the model formula.
  # Strip whitespace from both sides for a robust comparison.
  rterm_stripped <- gsub("\\s+", "", rterm)
  rterm_cp       <- .top_colon(rterm_stripped)
  rterm_rhs      <- if (!is.na(rterm_cp))
    substr(rterm_stripped, rterm_cp + 1L, nchar(rterm_stripped))
  else rterm_stripped
  user_rhs <- gsub("\\s+", "", p$by_raw)

  if (!identical(rterm_rhs, user_rhs))
    stop("Term mismatch: the model formula contains '", rterm, "' but the ",
         "supplied 'term' has a different right-hand side ('", p$by_raw, "'). ",
         "Check that the wrapper on the variety factor (e.g. vm(), ide(), or none) ",
         "matches the term as it was specified in the model's random formula.")

  # ---- Extract BLUPs and G-matrix ----------------------------------------
  if (struct == "fa") {
    sumfa <- faSummary(model, term = rterm)
    fag   <- sumfa$gammas[[rterm]]
    # Select by name: faSummary() returns <outer>, <inner>, blup, regblup, pres
    pvals <- sumfa$blups[[rterm]]$blups[, c("blup", fag$outer, fag$inner)]
    names(pvals) <- c("blup", enam, vnam)   # standardise FA column names
    Gmat  <- fag$Gmat
    pred  <- NULL                            # vcov unavailable; HSD will be NA
  } else {
    pred  <- predict(model, classify = p$classify, only = p$only_term,
                     vcov = TRUE, ...)
    # vparameters is always keyed by the bare "Group:Variety" term (p$classify),
    # regardless of the variance structure wrapper.
    raw_vp <- .asreml_vparams(model, p$classify)
    # corh / corgh: diagonal = genetic variances, off-diagonal = correlations.
    # Convert to a proper covariance G-matrix before use.
    Gmat <- if (struct %in% c("corh", "corgh"))
      .cor_to_cov_Gmat(raw_vp)
    else
      raw_vp   # us / diag already return a covariance matrix
    pvals <- pred$pvals
    names(pvals)[names(pvals) == "predicted.value"] <- "blup"
    # Normalise column names to bare variable names (strip any wrappers)
    names(pvals) <- sub(paste0("^", enam, "$"), enam, names(pvals))
    names(pvals) <- sub(paste0("^", vnam,  "$"), vnam,  names(pvals))
  }

  tsnams <- dimnames(Gmat)[[2L]]

  # ---- Resolve level / section from the G-matrix column labels -----------
  lab  <- .rreg_labels(tsnams, levs, sep, enam)
  tnam <- lab$level
  snam <- lab$section

  usnams <- unique(snam)
  ns     <- length(usnams)
  nvar   <- length(unique(as.character(pvals[[vnam]])))
  glev   <- unique(as.character(pvals[[vnam]]))

  resp_nams <- paste0("resp.", conditioned)
  hsd_nams  <- paste0("HSD.",  conditioned)

  # ---- Initialise outputs ------------------------------------------------
  tmat <- diag(nrow(Gmat))

  # beta: named list, one ns x |A_j| matrix per conditioned treatment
  beta <- setNames(vector("list", n_cond), conditioned)
  for (lv in conditioned)
    beta[[lv]] <- matrix(NA_real_, ns, length(cond_list[[lv]]),
                         dimnames = list(usnams, cond_list[[lv]]))

  # sigmat: ns x n_cond matrix of scalar conditional genetic variances
  sigmat <- matrix(NA_real_, ns, n_cond, dimnames = list(usnams, conditioned))

  blist <- vector("list", ns)

  # ---- Main loop: one iteration per site ---------------------------------
  for (i in seq_along(usnams)) {

    inds        <- which(snam == usnams[i])
    names(inds) <- tnam[inds]
    present     <- levs[levs %in% names(inds)]

    raw_df  <- matrix(NA_real_, nvar, ntreat, dimnames = list(NULL, levs))
    resp_df <- matrix(NA_real_, nvar, n_cond, dimnames = list(NULL, resp_nams))
    hsd_df  <- matrix(NA_real_, nvar, n_cond, dimnames = list(NULL, hsd_nams))

    # Fill raw BLUPs for all present treatments
    for (lv in present)
      raw_df[, lv] <- pvals$blup[pvals[[enam]] == tsnams[inds[lv]]]

    # Responsiveness BLUP for each conditioned treatment
    for (ci in seq_len(n_cond)) {

      lv_j <- conditioned[ci]
      A_j  <- cond_list[[lv_j]]

      # Skip if treatment or any member of conditioning set is absent
      if (!(lv_j %in% present) || !all(A_j %in% present)) next

      j_ind  <- inds[lv_j]
      a_inds <- inds[A_j]

      # ---- G sub-blocks --------------------------------------------------
      G_jj <- Gmat[j_ind,  j_ind ]
      G_jA <- Gmat[j_ind,  a_inds, drop = FALSE]   # 1 x |A|
      G_AA <- Gmat[a_inds, a_inds, drop = FALSE]   # |A| x |A|
      G_Aj <- Gmat[a_inds, j_ind,  drop = FALSE]   # |A| x 1

      # ---- Regression coefficients: beta_j = G_AA^{-1} G_Aj -------------
      beta_j <- if (length(A_j) == 1L) {
        drop(G_Aj) / G_AA[1L, 1L]               # scalar shortcut
      } else {
        tryCatch(
          drop(solve(G_AA, G_Aj)),
          error = function(e) {
            warning("G sub-matrix for level '", lv_j, "' in section '",
                    usnams[i], "' is singular. Skipping.")
            NULL
          }
        )
      }
      if (is.null(beta_j)) next

      # ---- Conditional genetic variance: sigma_j = G_jj - G_jA beta_j --
      sig_j <- G_jj - drop(G_jA %*% beta_j)

      # A negative conditional variance signals that the estimated G-matrix
      # is indefinite (not positive definite). This typically arises when a
      # correlation parameter is near its boundary (|r| -> 1), making the
      # G-matrix nearly singular. Schur complements of indefinite matrices
      # can be strongly negative — NOT just floating-point noise.
      # Skip this treatment-site combination and warn the user.
      if (is.na(sig_j) || sig_j <= 0) {
        warning("Non-positive conditional variance (", round(sig_j, 5L),
                ") for level '", lv_j, "' in section '", usnams[i], "'. ",
                "The estimated G-matrix may be indefinite (boundary correlation). ",
                "Check model convergence and varcomp estimates.")
        next
      }

      # Store parameters
      beta[[lv_j]][i, ] <- beta_j
      sigmat[i, ci]     <- sig_j

      # Update transformation matrix
      tmat[j_ind, a_inds] <- -beta_j

      # ---- Responsiveness BLUPs: u_j - u_A %*% beta_j ------------------
      u_A    <- raw_df[, A_j, drop = FALSE]         # nvar x |A|
      resp_j <- raw_df[, lv_j] - drop(u_A %*% beta_j)
      resp_df[, paste0("resp.", lv_j)] <- resp_j

      # ---- PEV via block decomposition (non-FA only) --------------------
      # For FA models pred = NULL; HSD columns remain NA.
      if (!is.null(pred)) {

        j_vcov_i <- which(pvals[[enam]] == tsnams[j_ind])
        a_vcov_i <- lapply(A_j, function(lv)
                       which(pvals[[enam]] == tsnams[inds[lv]]))

        # Stack vcov blocks: [j | A_1 | A_2 | ...]
        vcov_sub <- as.matrix(
          pred$vcov[c(j_vcov_i, unlist(a_vcov_i)),
                    c(j_vcov_i, unlist(a_vcov_i))]
        )

        # Linear combination coefficients: [1, -beta_1, -beta_2, ...]
        # PEV(tilde_u_j) = sum_p sum_q tcoef[p]*tcoef[q] * vcov_sub[B_p, B_q]
        # where B_p = block p of nvar rows/cols
        tcoef <- c(1, -beta_j)
        pev_j <- matrix(0, nvar, nvar)
        for (p in seq_along(tcoef)) {
          pi <- (p - 1L) * nvar + seq_len(nvar)
          for (q in seq_along(tcoef)) {
            qi    <- (q - 1L) * nvar + seq_len(nvar)
            pev_j <- pev_j + tcoef[p] * tcoef[q] * vcov_sub[pi, qi]
          }
        }

        if (!pev) pev_j <- diag(sig_j, nvar) - pev_j

        dv  <- diag(pev_j)
        sed <- outer(dv, dv, "+") - 2 * pev_j
        sed <- sed[lower.tri(sed)]
        sed[sed < 0L] <- NA_real_
        hsd_df[, paste0("HSD.", lv_j)] <-
          (mean(sqrt(sed), na.rm = TRUE) / sqrt(2)) *
          qtukey(0.95, nvar, df = nvar - 2L)
      }
    }

    blist[[i]] <- as.data.frame(cbind(raw_df, resp_df, hsd_df))
  }

  # ---- Transformed G-matrix ----------------------------------------------
  # Labels are rebuilt positionally from the resolved level / section pair.
  # Rewriting the composite string with gsub() instead corrupts any label
  # whose section portion contains a level name, and mangles levels that are
  # substrings of one another (e.g. "N1" inside "N10").
  TGmat  <- tmat %*% Gmat %*% t(tmat)
  uncond <- levs[vapply(cond_list, is.null, logical(1L))]

  prefix <- ifelse(tnam %in% uncond,      "eff.",
            ifelse(tnam %in% conditioned, "resp.", ""))

  tsnams_out <- if (!lab$composite) {
    paste0(prefix, tnam)
  } else if (lab$level_side == 1L) {
    paste0(prefix, tnam, sep, snam)
  } else {
    paste0(snam, sep, prefix, tnam)
  }
  dimnames(TGmat) <- list(tsnams_out, tsnams_out)

  # ---- Assemble blups data frame -----------------------------------------
  blups <- do.call(rbind, blist)
  blups <- cbind(
    data.frame(Site    = rep(usnams, each = nvar),
               Variety = rep(glev,   times = ns),
               stringsAsFactors = FALSE),
    blups
  )

  list(blups     = blups,
       TGmat     = TGmat,
       Gmat      = Gmat,
       beta      = beta,
       sigmat    = sigmat,
       tmat      = tmat,
       cond_list = cond_list,
       type      = type,
       sep       = sep,
       label_map = data.frame(label   = tsnams,
                              level   = tnam,
                              section = snam,
                              stringsAsFactors = FALSE))
}
