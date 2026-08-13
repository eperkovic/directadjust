#' R1 Build: certify adjustment sets for two treatments in one pass
#'
#' Implements the Build procedure (Theorem 5 of the paper) for two treatments
#' `x1` and `x2` with `x1` preceding `x2` in the causal order. The procedure
#' searches over ordered pairs of witness covariates \eqn{(W_a, W_b)} and
#' candidate sets \eqn{Z} drawn from the remaining covariates, and certifies
#' \eqn{Z} as a valid adjustment set for the joint effect of
#' \eqn{(X_1, X_2)} on \eqn{Y} when
#' \deqn{W_a \not\perp Y \mid Z, \quad W_a \perp Y \mid Z \cup \{X_1\},}
#' \deqn{W_b \not\perp Y \mid Z \cup \{X_1\}, \quad
#'       W_b \perp Y \mid Z \cup \{X_1, X_2\},}
#' with dependence and independence judged by the tester at its level. All
#' passing sets are returned (no early stopping), matching the paper's
#' simulation protocol.
#'
#' The covariates must be known not to be descendants of the treatments
#' (pretreatment covariates), and `x1` must precede `x2` causally; both are
#' assumptions of Theorem 5, not testable by the procedure.
#'
#' @param data A numeric data frame or matrix.
#' @param x1,x2 Column names of the two treatments, `x1` causally first.
#' @param y Column name of the outcome.
#' @param covariates Character vector of candidate (pretreatment) covariate
#'   names. Defaults to all remaining columns.
#' @param alpha Significance level for declaring dependence (default 0.05,
#'   as in the paper's simulations). Ignored when `tester` is supplied.
#' @param alpha_indep Threshold for declaring independence (p-value above
#'   it). Defaults to `alpha`; setting it higher (e.g. 0.1) certifies a set
#'   only on decisive evidence, as in Entner, Hoyer and Spirtes (2013).
#'   Ignored when `tester` is supplied.
#' @param max_size Largest size of the candidate set \eqn{Z} (drawn from the
#'   covariates excluding the two witnesses). Defaults to all of them.
#' @param tester Optional tester from [make_ci_tester()], e.g. to share a
#'   test counter across procedures.
#'
#' @return An object of class `r1_sets`: a list with the certified `sets`
#'   (list of character vectors, possibly empty), `found` (logical),
#'   `n_tests` performed, and the query (`x`, `y`).
#'
#' @note The paper's Theorem 5 covers any number of treatments; this
#'   implementation specializes to two, matching the simulation study. The
#'   general case follows the same pattern, chaining one witness condition
#'   pair per treatment along the causal ordering.
#'
#' @seealso [r1_combine()] for the per-treatment variant, and
#'   [estimate_joint_effect()] to estimate the effect over certified sets.
#'
#' @references LaPlante, S., Triantafillou, S., and Perkovic, E.
#'   Data-driven adjustment for multiple treatments. *Journal of Causal
#'   Inference* (to appear). Entner, D., Hoyer, P., and Spirtes, P. (2013).
#'   Data-driven covariate selection for nonparametric estimation of causal
#'   effects. *AISTATS*.
#'
#' @examples
#' # C confounds both treatments and the outcome; W1, W2 are witnesses
#' set.seed(1)
#' n <- 4000
#' C  <- rnorm(n); W1 <- rnorm(n); W2 <- rnorm(n)
#' X1 <- C + W1 + rnorm(n)
#' X2 <- C + W2 + rnorm(n)
#' Y  <- X1 + X2 + C + rnorm(n)
#' d  <- data.frame(C, W1, W2, X1, X2, Y)
#' r1_build(d, "X1", "X2", "Y")   # certifies sets containing C
#'
#' @export
r1_build <- function(data, x1, x2, y,
                     covariates = setdiff(colnames(as.data.frame(data)),
                                          c(x1, x2, y)),
                     alpha = 0.05, alpha_indep = alpha,
                     max_size = NULL, tester = NULL) {
  if (is.null(tester)) tester <- make_ci_tester(data, alpha, alpha_indep)
  W <- covariates
  if (length(W) < 2)
    stop("r1_build needs at least two covariates (two witnesses).")
  if (is.null(max_size)) max_size <- length(W) - 2L
  found <- list()
  for (wa in W) for (wb in setdiff(W, wa)) {
    rest <- setdiff(W, c(wa, wb))
    for (sz in 0:min(max_size, length(rest))) {
      for (Z in .combn_list(rest, sz)) {
        Z <- as.character(Z)
        ok <- tester$dep(wa, y, Z) &&
              tester$indep(wa, y, c(Z, x1)) &&
              tester$dep(wb, y, c(Z, x1)) &&
              tester$indep(wb, y, c(Z, x1, x2))
        if (ok) found[[length(found) + 1L]] <- sort(Z)
      }
    }
  }
  sets <- unique(found)
  structure(list(found = length(sets) > 0, sets = sets,
                 n_tests = tester$count(), x = c(x1, x2), y = y,
                 method = "R1 Build"),
            class = "r1_sets")
}
