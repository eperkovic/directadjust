#' Conditional independence tester with a test counter
#'
#' Wraps a conditional independence test around a data set, keeping a
#' running count of the tests performed. The returned object is what
#' [r1_entner()], [r1_build()], [r1_combine()], and [c_equivalent()]
#' consume, so a single tester can be shared between procedures to obtain a
#' combined test count.
#'
#' The certification procedures are agnostic to the test: they consume only
#' the dependence/independence decisions, so any conditional independence
#' test appropriate for the data can be supplied through `test`. The
#' Gaussian Fisher-z default (and the linear-regression estimator in
#' [estimate_joint_effect()]) are where Gaussian/linear assumptions enter.
#'
#' @param data A numeric data frame or matrix; columns are variables.
#' @param alpha Significance level for declaring dependence: dependence is
#'   declared when the test's p-value is at most `alpha`. The paper's
#'   simulations use `alpha = 0.05`.
#' @param alpha_indep Threshold for declaring independence: independence is
#'   declared when the p-value exceeds `alpha_indep`. Defaults to `alpha`,
#'   in which case every test is decisive. Setting `alpha_indep` larger
#'   than `alpha` (for example 0.1, as in Entner, Hoyer and Spirtes, 2013)
#'   leaves an undecided region in between, so a candidate set is certified
#'   only when every required condition is decisively in the right
#'   direction.
#' @param test Optional p-value function `function(i, j, S, data)`: two
#'   column names, a character vector of conditioning names, and the data;
#'   must return the p-value of the conditional independence test of `i`
#'   and `j` given `S`. Defaults to the Gaussian Fisher-z test
#'   ([pcalg::gaussCItest()]). Supply, for example, a G-square test for
#'   discrete data or a kernel-based test for nonlinear continuous data.
#'
#' @return A list with components
#'   \describe{
#'     \item{`dep(i, j, S)`}{function of two variable names and a character
#'       vector of conditioning names, returning `TRUE` when dependence is
#'       declared (p-value at most `alpha`);}
#'     \item{`indep(i, j, S)`}{likewise, `TRUE` when independence is
#'       declared (p-value above `alpha_indep`);}
#'     \item{`count()`}{function returning the number of tests performed so
#'       far.}
#'   }
#'
#' @examples
#' d <- data.frame(a = rnorm(200))
#' d$b <- d$a + rnorm(200)
#' ct <- make_ci_tester(d)
#' ct$indep("a", "b", character(0))   # FALSE, they are dependent
#' ct$dep("a", "b", character(0))     # TRUE
#' ct$count()                          # 2
#'
#' @export
make_ci_tester <- function(data, alpha = 0.05, alpha_indep = alpha,
                           test = NULL) {
  data <- as.data.frame(data)
  nm <- colnames(data)
  count <- 0L
  if (is.null(test)) {
    C <- stats::cor(data)
    n <- nrow(data)
    pval <- function(i, j, S)
      pcalg::gaussCItest(match(i, nm), match(j, nm), match(S, nm),
                         list(C = C, n = n))
  } else {
    pval <- function(i, j, S) test(i, j, S, data)
  }
  pv <- function(i, j, S) {
    count <<- count + 1L
    pval(i, j, S)
  }
  list(dep   = function(i, j, S = character(0)) pv(i, j, S) <= alpha,
       indep = function(i, j, S = character(0)) pv(i, j, S) > alpha_indep,
       count = function() count)
}

# all subsets of size m of a character vector, as a list (internal)
.combn_list <- function(v, m) {
  if (m == 0) return(list(character(0)))
  if (length(v) < m) return(list())
  asplit(utils::combn(v, m), 2)
}
