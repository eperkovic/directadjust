#' R1 Entner: certify adjustment sets for a single treatment
#'
#' Implements the first rule of Entner, Hoyer and Spirtes (2013), the
#' paper's Theorem 1, for a single treatment: a candidate set \eqn{Z} is
#' certified as an adjustment set relative to \eqn{(X, Y)} when some witness
#' covariate \eqn{W_a} satisfies \eqn{W_a \not\perp Y \mid Z} and
#' \eqn{W_a \perp Y \mid Z \cup \{X\}}. All certified sets are returned.
#'
#' The covariates must be known not to be descendants of the treatment (the
#' pretreatment assumption); \eqn{X} must precede \eqn{Y}.
#'
#' @inheritParams r1_build
#' @param x Treatment column name.
#' @param max_size Largest size of the candidate set \eqn{Z} (drawn from the
#'   covariates excluding the witness). Defaults to all of them.
#'
#' @return An object of class `r1_sets`; see [r1_build()].
#'
#' @seealso [r1_build()] and [r1_combine()] for multiple treatments, and
#'   [c_equivalent()] for extending the certified sets through
#'   c-equivalence (the paper's Theorem 3).
#'
#' @references Entner, D., Hoyer, P., and Spirtes, P. (2013). Data-driven
#'   covariate selection for nonparametric estimation of causal effects.
#'   *AISTATS*.
#'
#' @examples
#' # C confounds treatment and outcome; W is a witness
#' set.seed(1)
#' n <- 4000
#' C <- rnorm(n); W <- rnorm(n)
#' X <- C + W + rnorm(n)
#' Y <- X + C + rnorm(n)
#' d <- data.frame(C, W, X, Y)
#' r1_entner(d, "X", "Y")   # certifies sets containing C
#'
#' @export
r1_entner <- function(data, x, y,
                      covariates = setdiff(colnames(as.data.frame(data)),
                                           c(x, y)),
                      alpha = 0.05, alpha_indep = alpha,
                      max_size = NULL, tester = NULL) {
  if (is.null(tester)) tester <- make_ci_tester(data, alpha, alpha_indep)
  W <- covariates
  if (length(W) < 1)
    stop("r1_entner needs at least one covariate (a witness).")
  if (is.null(max_size)) max_size <- length(W) - 1L
  found <- list()
  for (wa in W) {
    rest <- setdiff(W, wa)
    for (sz in 0:min(max_size, length(rest))) {
      for (Z in .combn_list(rest, sz)) {
        Z <- as.character(Z)
        if (tester$dep(wa, y, Z) && tester$indep(wa, y, c(Z, x)))
          found[[length(found) + 1L]] <- sort(Z)
      }
    }
  }
  sets <- unique(found)
  structure(list(found = length(sets) > 0, sets = sets,
                 n_tests = tester$count(), x = x, y = y,
                 method = "R1 Entner"),
            class = "r1_sets")
}
