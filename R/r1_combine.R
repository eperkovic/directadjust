#' R1 Combine: assemble adjustment sets from single-treatment certificates
#'
#' Implements the Combine procedure (Theorem 7 of the paper). For each
#' treatment \eqn{X_i} separately, the procedure collects every pair of a
#' witness covariate \eqn{W_a} and a set \eqn{T} with
#' \deqn{W_a \not\perp Y \mid T \quad\text{and}\quad
#'       W_a \perp Y \mid T \cup \{X_i\},}
#' reduces the passing \eqn{T}'s to the inclusion-minimal ones, prunes each
#' to a minimal adjustment set \eqn{Z_i \subseteq T} (see `prune` below),
#' and returns all unions \eqn{Z_1 \cup Z_2} over pairs of pruned
#' per-treatment sets, with treatments removed. The candidate pool for the
#' second treatment includes the first treatment as a possible member of
#' \eqn{T} (never as a witness), reflecting that `x1` is a non-descendant
#' of `x2`.
#'
#' Theorem 7 requires each \eqn{Z_i} to be a *minimal* adjustment set
#' relative to \eqn{(X_i, Y)}; the paper's Example 12 shows that combining
#' non-minimal sets can produce an invalid union. The default
#' `prune = "lemma8"` therefore applies the paper's Lemma 8: \eqn{T} is
#' accepted as minimal when for every \eqn{T_j \in T} both
#' \eqn{X_i \not\perp T_j \mid T \setminus \{T_j\}} and
#' \eqn{Y \not\perp T_j \mid [T \setminus \{T_j\}] \cup \{X_i\}} are
#' declared dependent; otherwise the smallest subset \eqn{Z \subset T}
#' satisfying Lemma 8's conditions (iii)-(v) is used, and the candidate is
#' dropped when no such subset is certified. The alternative
#' `prune = "witness"` prunes \eqn{T} to the smallest subset still passing
#' the two witness conditions; this is cheaper but does not verify
#' minimality, so on structures like the paper's Example 12 it can certify
#' an invalid set.
#'
#' As with [r1_build()], covariates must be pretreatment and `x1` must
#' precede `x2` causally; these are assumptions, not tested.
#'
#' @inheritParams r1_build
#' @param max_size Largest size of the per-treatment candidate set \eqn{T}.
#'   Defaults to all available candidates.
#' @param prune How each certified \eqn{T} is reduced to a minimal
#'   adjustment set: `"lemma8"` (default) tests minimality through the
#'   paper's Lemma 8; `"witness"` uses a faster heuristic (smallest
#'   subset still passing the witness conditions).
#'
#' @return An object of class `r1_sets`; see [r1_build()].
#'
#' @note The paper's Theorem 7 covers any number of treatments; this
#'   implementation specializes to two, matching the simulation study. The
#'   general case pares and unions one certified set per treatment in the
#'   same way.
#'
#' @examples
#' set.seed(1)
#' n <- 4000
#' C  <- rnorm(n); W1 <- rnorm(n); W2 <- rnorm(n)
#' X1 <- C + W1 + rnorm(n)
#' X2 <- C + W2 + rnorm(n)
#' Y  <- X1 + X2 + C + rnorm(n)
#' d  <- data.frame(C, W1, W2, X1, X2, Y)
#' r1_combine(d, "X1", "X2", "Y")
#'
#' @export
r1_combine <- function(data, x1, x2, y,
                       covariates = setdiff(colnames(as.data.frame(data)),
                                            c(x1, x2, y)),
                       alpha = 0.05, alpha_indep = alpha,
                       max_size = NULL, tester = NULL,
                       prune = c("lemma8", "witness")) {
  prune <- match.arg(prune)
  if (is.null(tester)) tester <- make_ci_tester(data, alpha, alpha_indep)
  W <- covariates
  if (length(W) < 1)
    stop("r1_combine needs at least one covariate (a witness).")

  collect <- function(x_i, extra) {
    hits <- list()
    for (wa in W) {
      pool <- c(setdiff(W, wa), extra)
      cap <- if (is.null(max_size)) length(pool) else min(max_size, length(pool))
      for (sz in 0:cap) for (Ti in .combn_list(pool, sz)) {
        Ti <- as.character(Ti)
        if (tester$dep(wa, y, Ti) && tester$indep(wa, y, c(Ti, x_i)))
          hits[[length(hits) + 1L]] <- list(W = wa, T = Ti)
      }
    }
    hits
  }
  # keep one witness per distinct T, then keep inclusion-minimal T's
  dedupe_min <- function(h) {
    h <- h[!duplicated(vapply(h, function(e) paste(sort(e$T), collapse = ","), ""))]
    Ts <- lapply(h, `[[`, "T")
    keep <- vapply(seq_along(Ts), function(i)
      !any(vapply(seq_along(Ts), function(j)
        j != i && length(Ts[[j]]) < length(Ts[[i]]) &&
        all(Ts[[j]] %in% Ts[[i]]), TRUE)), TRUE)
    h[keep]
  }

  # joint independence A _||_ B | C via the chain decomposition
  joint_indep <- function(A, B, C) {
    if (length(A) == 0 || length(B) == 0) return(TRUE)
    for (i in seq_along(A)) for (j in seq_along(B)) {
      cond <- c(C, A[seq_len(i - 1)], B[seq_len(j - 1)])
      if (!tester$indep(A[i], B[j], cond)) return(FALSE)
    }
    TRUE
  }
  # Lemma 8 conditions (i)/(ii) resp. (iii)/(iv): every element of S is
  # tied to both x_i and y given the rest of S
  elementwise_dep <- function(S, x_i) {
    all(vapply(S, function(v) {
      rest <- setdiff(S, v)
      tester$dep(x_i, v, rest) && tester$dep(y, v, c(rest, x_i))
    }, TRUE))
  }
  # smallest minimal adjustment set Z within T, per Lemma 8; NULL when no
  # subset is certified
  minimal_Z_lemma8 <- function(Ti, x_i) {
    if (elementwise_dep(Ti, x_i)) return(Ti)          # (i)-(ii): T minimal
    if (length(Ti) >= 1) {
      for (sz in 0:(length(Ti) - 1L)) for (S in .combn_list(Ti, sz)) {
        S <- as.character(S)
        if (!elementwise_dep(S, x_i)) next            # (iii)-(iv)
        rest <- setdiff(Ti, S)
        if (joint_indep(y, rest, c(S, x_i)) ||        # (v)
            joint_indep(x_i, rest, S)) return(S)
      }
    }
    NULL
  }
  # simulation-study protocol: smallest subset of T still passing the
  # witness conditions; does not verify minimality
  minimal_Z_witness <- function(wa, Ti, x_i) {
    for (sz in 0:length(Ti)) for (S in .combn_list(Ti, sz)) {
      S <- as.character(S)
      if (tester$dep(wa, y, S) && tester$indep(wa, y, c(S, x_i))) return(S)
    }
    NULL
  }
  reduce <- function(wa, Ti, x_i) {
    if (prune == "lemma8") minimal_Z_lemma8(Ti, x_i)
    else minimal_Z_witness(wa, Ti, x_i)
  }

  h1 <- collect(x1, extra = character(0))
  h2 <- if (length(h1) > 0) collect(x2, extra = x1) else list()
  sets <- list()
  if (length(h1) > 0 && length(h2) > 0) {
    m1 <- dedupe_min(h1)
    m2 <- dedupe_min(h2)
    for (a in m1) for (b in m2) {
      Z1 <- reduce(a$W, a$T, x1)
      Z2 <- reduce(b$W, b$T, x2)
      if (is.null(Z1) || is.null(Z2)) next
      sets[[length(sets) + 1L]] <- sort(setdiff(union(Z1, Z2), c(x1, x2)))
    }
    sets <- unique(sets)
  }
  structure(list(found = length(sets) > 0, sets = sets,
                 n_tests = tester$count(), x = c(x1, x2), y = y,
                 method = "R1 Combine"),
            class = "r1_sets")
}

#' Print certified adjustment sets
#'
#' Prints the procedure that produced the object, the treatments and outcome
#' it applies to, the number of conditional independence tests performed, and
#' each certified set. The empty set prints as `{ }`.
#'
#' @param x An `r1_sets` object, as returned by [r1_entner()], [r1_build()],
#'   or [r1_combine()].
#' @param ... Ignored, present for consistency with the generic.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' set.seed(1)
#' n <- 4000
#' C  <- rnorm(n); W1 <- rnorm(n); W2 <- rnorm(n)
#' X1 <- C + W1 + rnorm(n)
#' X2 <- C + W2 + rnorm(n)
#' Y  <- X1 + X2 + C + rnorm(n)
#' print(r1_build(data.frame(C, W1, W2, X1, X2, Y), "X1", "X2", "Y"))
#'
#' @export
print.r1_sets <- function(x, ...) {
  cat(x$method, "for joint effect of {", paste(x$x, collapse = ", "),
      "} on", x$y, "\n")
  if (!x$found) {
    cat("no adjustment set certified;", x$n_tests, "CI tests performed\n")
  } else {
    cat(length(x$sets), "certified set(s);", x$n_tests, "CI tests performed\n")
    for (s in x$sets)
      cat("  {", if (length(s) == 0) "" else paste(s, collapse = ", "), "}\n")
  }
  invisible(x)
}
