#' @keywords internal
"_PACKAGE"

#' @importFrom stats coef contr.helmert lm model.matrix p.adjust pchisq pf
#'   predict qt qtukey residuals rnorm runif sd setNames terms
#' @importFrom utils combn write.csv
NULL

# ---- Internal wrappers for ASReml-R calls --------------------------------
# These thin wrappers live in the biomAid namespace so that
# testthat::local_mocked_bindings() can substitute them in tests without
# requiring an ASReml-R licence in CI.

#' @noRd
.fa_asreml <- function(model, ...) {
  if (!requireNamespace("ASExtras4", quietly = TRUE))
    stop("Package 'ASExtras4' is required for FA analysis. ",
         "Please install it from the ASReml-R website.", call. = FALSE)
  ASExtras4::fa.asreml(model, ...)
}

#' @noRd
.asreml_vparams <- function(model, term) {
  summary(model, vparameters = TRUE)$vparameters[[term]]
}
