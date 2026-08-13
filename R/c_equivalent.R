#' Test c-equivalence of two candidate adjustment sets
#'
#' Implements the probabilistic criteria for confounding equivalence (the
#' paper's Theorem 3, after Pearl and Paz): sets \eqn{Z} and \eqn{T} are
#' declared c-equivalent relative to \eqn{(X, Y)} when either
#' \describe{
#'   \item{(i)}{\eqn{X \perp Z \mid T} and \eqn{Y \perp T \mid Z \cup X}, or}
#'   \item{(ii)}{\eqn{X \perp T \mid Z} and \eqn{Y \perp Z \mid T \cup X}.}
#' }
#' If \eqn{T} is c-equivalent to a certified adjustment set (for example one
#' found by [r1_entner()] or [r1_combine()]), then \eqn{T} is an adjustment
#' set as well, which extends the reach of the data-driven rules; see the
#' paper's Examples 2, 4, and 9.
#'
#' Joint independence statements involving sets are decided through the
#' chain decomposition: \eqn{A \perp B \mid C} holds when every
#' \eqn{A_i \perp B_j \mid C \cup A_{<i} \cup B_{<j}} is declared
#' independent, each by a Fisher-z test.
#'
#' @inheritParams r1_build
#' @param x Treatment column name(s).
#' @param z,t The two candidate sets (character vectors; either may be
#'   `character(0)` for the empty set).
#'
#' @return A list with `equivalent` (logical), `criterion` (`"i"`, `"ii"`,
#'   or `NA` when neither holds), and `n_tests`.
#'
#' @references Pearl, J., and Paz, A. (2014). Confounding equivalence in
#'   causal inference. *Journal of Causal Inference* 2(1):75--93.
#'
#' @examples
#' # Z is an adjustment set; the empty set is c-equivalent to it when
#' # X and Z are marginally independent (criterion i)
#' set.seed(1)
#' n <- 4000
#' Z <- rnorm(n)
#' X <- rnorm(n)
#' Y <- X + Z + rnorm(n)
#' d <- data.frame(Z, X, Y)
#' c_equivalent(d, "X", "Y", z = "Z", t = character(0))
#'
#' @export
c_equivalent <- function(data, x, y, z, t,
                         alpha = 0.05, alpha_indep = alpha,
                         tester = NULL) {
  if (is.null(tester)) tester <- make_ci_tester(data, alpha, alpha_indep)
  z <- as.character(z); t <- as.character(t); x <- as.character(x)

  # joint independence A _||_ B | C via the chain decomposition; TRUE only
  # when every pairwise test is declared independent
  joint_indep <- function(A, B, C) {
    if (length(A) == 0 || length(B) == 0) return(TRUE)
    for (i in seq_along(A)) for (j in seq_along(B)) {
      cond <- c(C, A[seq_len(i - 1)], B[seq_len(j - 1)])
      if (!tester$indep(A[i], B[j], cond)) return(FALSE)
    }
    TRUE
  }

  crit_i  <- joint_indep(x, z, t) && joint_indep(y, t, c(z, x))
  crit_ii <- !crit_i &&
             joint_indep(x, t, z) && joint_indep(y, z, c(t, x))

  list(equivalent = crit_i || crit_ii,
       criterion  = if (crit_i) "i" else if (crit_ii) "ii" else NA_character_,
       n_tests    = tester$count())
}
