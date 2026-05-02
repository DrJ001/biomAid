# =============================================================================
# accuracy.R  —  Mrode (2014) model-based BLUP accuracy for ASReml-R V4
# =============================================================================
# Computes per-group mean accuracy:  r_j = mean[ sqrt(1 - PEV_ij / G_jj) ]
# Optionally also returns Cullis (2006) H² per group.
#
# Supports:  fa(Group, k):id(Var),  diag(Group):id(Var),
#            corgh(Group):id(Var),  corh(Group):id(Var),
#            us(Group):id(Var),     and id(Var) alone.
#
# predict.asreml() only= string used per structure:
#   fa(Group, k):id(Var)  ->  "fa(Group, k):Var"   (space after comma required)
#   diag / corgh / corh / us  ->  "Group:Var"
#   id(Var) alone         ->  "Var"
#
# Usage:
#   accuracy(model)                              # both metrics (default)
#   accuracy(model, metric = "accuracy")         # Mrode accuracy only
#   accuracy(model, metric = "gen.H2")           # Cullis (2006) generalised H2 only
#   accuracy(model, by_variety = TRUE)           # per-variety rows
# =============================================================================


#' Model-Based BLUP Accuracy for ASReml-R V4 Multi-Environment Trials
#'
#' Computes per-group mean Mrode (2014) accuracy and/or Cullis (2006) H\eqn{^2}
#' from an ASReml-R V4 mixed model. Supports \code{fa()}, \code{diag()},
#' \code{corgh()}, \code{corh()}, \code{us()}, and single-environment
#' \code{id()} random structures for the variety term.
#'
#' @param model      A fitted \code{asreml} model object.
#' @param classify   Character. Override the \code{predict()} classify string
#'   (auto-detected from the random formula when \code{NULL}).
#' @param metric     Character vector: \code{"accuracy"}, \code{"gen.H2"},
#'   or both (default). Controls which columns appear in the output.
#' @param pworkspace Character. Passed to \code{predict.asreml()}.
#'   Default \code{"2gb"}.
#' @param by_variety Logical. If \code{TRUE}, return one row per variety x
#'   group instead of group-level summaries. Default \code{FALSE}.
#' @param present    Logical. Controls which variety x group combinations are
#'   included in the output when \code{by_variety = TRUE}.
#'   \itemize{
#'     \item \code{TRUE} (default) — restrict to variety x group combinations
#'       that are \emph{observed} in the data. Unobserved combinations
#'       (varieties absent from a site) are excluded, which can produce holes
#'       in heatmaps for unbalanced trials.
#'     \item \code{FALSE} — include \emph{all} combinations returned by
#'       \code{predict.asreml()}, including unobserved ones. For factor-analytic
#'       models the G-matrix allows prediction of BLUPs for unobserved variety x
#'       site combinations; their PEV is higher than for observed combinations,
#'       reflecting reduced information. For simpler structures (e.g.
#'       \code{diag()}), unobserved combinations may receive the same SE as
#'       observed ones.
#'   }
#'
#' @return A data frame with columns \code{group}, \code{n_vars}, \code{G_jj},
#'   and one or both of \code{mean_acc}/\code{sd_acc} (Mrode accuracy) and
#'   \code{gen.H2} (Cullis H\eqn{^2}). When \code{by_variety = TRUE} the
#'   \code{variety} and \code{pev} columns replace \code{n_vars}/\code{sd_acc}.
#'
#' @examples
#' \dontrun{
#' acc <- accuracy(model_fa)
#' acc_bv <- accuracy(model_fa, by_variety = TRUE)
#' acc_bv_all <- accuracy(model_fa, by_variety = TRUE, present = FALSE)
#' accuracy(model_fa, metric = "gen.H2")
#' }
#'
#' @export
accuracy <- function(model,
                     classify   = NULL,
                     metric     = c("accuracy", "gen.H2"),
                     pworkspace = "2gb",
                     by_variety = FALSE,
                     present    = TRUE) {

  metric <- match.arg(metric, several.ok = TRUE)
  want_H2  <- "gen.H2"   %in% metric
  want_acc <- "accuracy" %in% metric

  p   <- .parse_met_formula(model)
  if (!is.null(classify)) p$classify <- classify

  dat <- .get_model_data(model)
  gl  <- if (!is.null(p$group_var))
           levels(droplevels(dat[[p$group_var]])) else NULL

  # When there is no group structure (single-environment), treat the whole
  # model as one group labelled "all" so it runs through the same loop.
  single_env <- is.null(gl)
  if (single_env) gl <- "all"

  G <- .extract_G_diag(model, if (single_env) NULL else gl, p$group_var)
  if (is.null(G)) stop("[accuracy] Cannot extract G_jj from model.")
  if (single_env) G <- setNames(G, "all")

  # predict(only=) returns pure BLUP SEs — avoids fixed-effect inflation
  # that causes PEV > G_jj and spurious zero accuracies
  need_sed <- want_H2
  out <- predict(model, classify = p$classify, only = p$only_term,
                 sed = need_sed, pworkspace = pworkspace)
  pv  <- out$pvals
  names(pv) <- tolower(names(pv))
  pv$pev    <- pv$std.error^2
  pv$.row   <- seq_len(nrow(pv))   # preserve original row order for SED lookup

  # Mark only the variety-group combinations actually in the data.
  # For FA models, predict() returns rows for all factor-level combinations
  # including unobserved ones; we restrict to observed pairs via data lookup.
  gc <- if (!single_env) tolower(p$group_var) else NULL
  vc <- tolower(p$by_var)

  gd <- if (!single_env) names(dat)[match(tolower(p$group_var), tolower(names(dat)))] else NULL
  vd <- names(dat)[match(tolower(p$by_var), tolower(names(dat)))]

  if (!single_env && !is.null(gd) && !is.null(vd) && present) {
    # Restrict to variety x group combinations observed in the data.
    obs      <- paste(dat[[gd]], dat[[vd]], sep = "\001")
    pv$valid <- paste(pv[[gc]], pv[[vc]], sep = "\001") %in% obs
  } else {
    # Include all combinations with finite, positive SE.
    # For FA models this covers unobserved variety x site combinations whose
    # BLUPs are supported by the G-matrix (present = FALSE use case).
    pv$valid <- is.finite(pv$std.error) & pv$std.error > 0
  }

  # ---- Unified loop (works for single- and multi-environment) ---------------
  rows <- lapply(gl, function(grp) {

    if (single_env) {
      idx <- which(pv$valid)
    } else {
      idx <- which(pv[[gc]] == grp & pv$valid)
    }
    gjj <- unname(G[grp])   # unname() prevents row-name warnings in data.frame

    if (length(idx) < 2L || is.na(gjj) || gjj <= 0) {
      if (by_variety) {
        # Return zero-row frame with the correct columns
        row <- data.frame(group = character(0), variety = character(0),
                          G_jj = numeric(0),   pev     = numeric(0))
        if (want_acc) row$accuracy  <- numeric(0)
        if (want_H2)  row$gen.H2   <- numeric(0)
      } else {
        row <- data.frame(group = grp, n_vars = 0L, G_jj = gjj)
        if (want_acc)  { row$mean_acc <- NA_real_; row$sd_acc <- NA_real_ }
        if (want_H2)     row$gen.H2  <- NA_real_
      }
      return(row)
    }

    # Mrode accuracy
    a <- sqrt(pmax(0, 1 - pv$pev[idx] / gjj))

    # Cullis (2006) generalised H²: 1 - (mean_SED)² / (2 * G_jj)
    h2 <- if (want_H2) {
      if (!is.null(out$sed) && max(idx) <= nrow(out$sed)) {
        ut <- out$sed[idx, idx]
        ut <- ut[upper.tri(ut)]
        ut <- ut[is.finite(ut) & ut > 0]
        if (length(ut)) max(0, min(1, 1 - mean(ut)^2 / (2 * gjj))) else NA_real_
      } else NA_real_
    } else NULL

    if (by_variety) {
      out_bv <- data.frame(group    = grp,
                           variety  = as.character(pv[[vc]][idx]),
                           G_jj     = gjj,
                           pev      = pv$pev[idx],
                           stringsAsFactors = FALSE)
      if (want_acc) out_bv$accuracy <- a
      if (want_H2)  out_bv$gen.H2  <- h2   # group-level statistic (Cullis H²)
      return(out_bv)
    }

    row <- data.frame(group  = grp,
                      n_vars = length(idx),
                      G_jj   = gjj)
    if (want_acc) {
      row$mean_acc <- mean(a, na.rm = TRUE)
      row$sd_acc   <- sd(a,   na.rm = TRUE)
    }
    if (want_H2) row$gen.H2 <- h2
    row
  })

  do.call(rbind, c(rows, make.row.names = FALSE))
}

# =============================================================================
# HELPERS
# =============================================================================

.get_model_data <- function(model) {
  dat <- tryCatch(eval(model$call$data), error = function(e) NULL)
  if (is.data.frame(dat)) return(dat)
  dat <- tryCatch(eval(model$call$data, envir = parent.frame(2L)),
                  error = function(e) NULL)
  if (is.data.frame(dat)) return(dat)
  stop("[accuracy] Cannot retrieve data from model$call$data.")
}

.extract_G_diag <- function(model, group_levels, group_var = NULL) {

  vc  <- summary(model)$varcomp
  nms <- rownames(vc)
  cmp <- vc$component

  if (is.null(group_levels)) {
    ok <- which(!grepl("units|!R$|R!|residual", nms, ignore.case = TRUE))
    if (length(ok) >= 1L) return(setNames(max(0, cmp[ok[1L]]), "all"))
    return(NULL)
  }

  # ---- FA: specific variances (Psi) + factor loadings (Lambda) --------------
  # Psi row names end with "!GroupLevel!var" (FA notation)
  vi_fa <- sapply(group_levels, function(g) {
    i <- which(endsWith(nms, paste0("!", g, "!var")))
    if (length(i) == 1L) i else NA_integer_
  })

  if (!all(is.na(vi_fa))) {
    Psi <- pmax(0, cmp[vi_fa])

    Lc <- list()
    for (f in seq_len(20L)) {
      fi <- sapply(group_levels, function(g) {
        i <- which(endsWith(nms, paste0("!", g, "!fa", f)))
        if (length(i) == 1L) i else NA_integer_
      })
      if (any(is.na(fi))) break
      Lc[[f]] <- cmp[fi]
    }

    G <- if (length(Lc)) rowSums(do.call(cbind, Lc)^2) + Psi else Psi
    return(setNames(G, group_levels))
  }

  # ---- diag / corgh / corh: specific variances -------------------------------
  # Row names end with "!Env_GroupLevel"  (e.g. "Env:Variety!Env_E1")
  vi_diag <- sapply(group_levels, function(g) {
    i <- which(endsWith(nms, paste0("_", g)))
    if (length(i) == 1L) i else NA_integer_
  })

  if (!all(is.na(vi_diag))) {
    G <- pmax(0, cmp[vi_diag])
    return(setNames(G, group_levels))
  }

  # ---- us: unstructured — diagonal entries are "!Env_Gj:Gj" ----------------
  # Off-diagonal entries are "!Env_Gj:Gi" (j > i); we want only j == i.
  vi_us <- sapply(group_levels, function(g) {
    i <- which(endsWith(nms, paste0("_", g, ":", g)))
    if (length(i) == 1L) i else NA_integer_
  })

  if (!all(is.na(vi_us))) {
    G <- pmax(0, cmp[vi_us])
    return(setNames(G, group_levels))
  }

  NULL
}

# Split string on sep at parenthesis depth 0
.split_top <- function(s, sep = "+") {
  res <- character(0L); depth <- 0L; start <- 1L
  for (i in seq_len(nchar(s))) {
    ch <- substr(s, i, i)
    if      (ch == "(")           depth <- depth + 1L
    else if (ch == ")")           depth <- depth - 1L
    else if (ch == sep && !depth) { res <- c(res, substr(s, start, i-1L)); start <- i+1L }
  }
  c(res, substr(s, start, nchar(s)))
}

# Position of first ":" at depth 0
.top_colon <- function(s) {
  depth <- 0L
  for (i in seq_len(nchar(s))) {
    ch <- substr(s, i, i)
    if      (ch == "(") depth <- depth + 1L
    else if (ch == ")") depth <- depth - 1L
    else if (ch == ":" && !depth) return(i)
  }
  NA_integer_
}

# Content inside outermost parens, e.g. "id(Var)" -> "Var"
.inside <- function(s) {
  op <- regexpr("\\(", s)
  if (op < 0L) return(s)
  substr(s, op + 1L, nchar(s) - 1L)
}

.parse_met_formula <- function(model) {

  raw <- tryCatch(paste(deparse(model$call$random), collapse = ""),
                  error = function(e)
                    stop("[accuracy] Cannot read model$call$random."))

  rnd   <- gsub("\\s+", "", sub("^~\\s*", "", raw))
  terms <- .split_top(rnd, "+")

  # Supported multi-environment variance structures
  met <- grep("^(fa|diag|corh|corgh|us)\\(", terms, value = TRUE)

  if (!length(met)) {
    # Single-environment: e.g.  id(Variety)  or  Variety
    by <- .inside(terms[1L])
    return(list(group_var  = NULL,
                by_var     = by,
                n_fa       = NULL,
                only_term  = by,         # only="Variety" works for id(Variety)
                classify   = by))
  }

  cp  <- .top_colon(met[1L])
  lft <- substr(met[1L], 1L, cp - 1L)
  rgt <- substr(met[1L], cp + 1L, nchar(met[1L]))

  type <- sub("\\(.*", "", lft)
  grp  <- trimws(strsplit(.inside(lft), ",")[[1L]][1L])
  by   <- .inside(rgt)   # strip id() wrapper if present

  nfa  <- if (type == "fa")
            suppressWarnings(as.integer(trimws(strsplit(.inside(lft), ",")[[1L]][2L])))
          else NULL

  # Build only= string — must match what ASReml stores internally:
  #   fa:  "fa(Group, k):Var"  (space after comma is required)
  #   all others: "Group:Var"
  only <- if (type == "fa" && !is.null(nfa))
            sprintf("fa(%s, %d):%s", grp, nfa, by)
          else
            paste0(grp, ":", by)

  list(group_var  = grp,
       by_var     = by,
       n_fa       = nfa,
       only_term  = only,
       classify   = paste0(grp, ":", by))
}
