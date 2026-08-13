#' Estimate the joint treatment effect over certified adjustment sets
#'
#' Given the output of [r1_build()] or [r1_combine()], estimates the joint
#' effect of raising both treatments by one unit,
#' \eqn{E[Y \mid do(x_1 + 1, x_2 + 1)] - E[Y \mid do(x_1, x_2)]}, by the sum
#' of the two treatment coefficients in the linear regression of the outcome
#' on both treatments and a certified set, averaged over all certified sets.
#' This matches the estimator of the paper's simulation study; it is
#' appropriate for linear models and serves as a plug-in illustration
#' otherwise.
#'
#' @param data The data used for certification.
#' @param sets An `r1_sets` object, or a list of character vectors.
#' @param x1,x2,y Column names; taken from the `r1_sets` object when omitted.
#'   For the single-treatment output of [r1_entner()], only `x1` applies and
#'   the estimate is that treatment's coefficient.
#'
#' @return A list with the per-set estimates (`by_set`) and their mean
#'   (`estimate`), or `NULL` when no set was certified.
#'
#' @seealso [r1_build()] and [r1_combine()], which produce the `r1_sets`
#'   objects this function consumes.
#'
#' @examples
#' # true joint effect of raising both treatments by one unit is 2
#' set.seed(1)
#' n <- 4000
#' C  <- rnorm(n); W1 <- rnorm(n); W2 <- rnorm(n)
#' X1 <- C + W1 + rnorm(n)
#' X2 <- C + W2 + rnorm(n)
#' Y  <- X1 + X2 + C + rnorm(n)
#' d  <- data.frame(C, W1, W2, X1, X2, Y)
#' sets <- r1_build(d, "X1", "X2", "Y")
#' estimate_joint_effect(d, sets)$estimate
#'
#' @export
estimate_joint_effect <- function(data, sets, x1 = NULL, x2 = NULL, y = NULL) {
  if (inherits(sets, "r1_sets")) {
    xs <- sets$x
    if (is.null(y)) y <- sets$y
    sets <- sets$sets
  } else {
    xs <- c(x1, x2)
  }
  if (!is.null(x1) || !is.null(x2)) xs <- c(x1, x2)
  if (length(sets) == 0) return(NULL)
  data <- as.data.frame(data)
  by_set <- vapply(sets, function(Z) {
    f <- stats::lm(stats::reformulate(c(xs, Z), response = y), data = data)
    sum(stats::coef(f)[xs])
  }, numeric(1))
  list(by_set = by_set, estimate = mean(by_set))
}
