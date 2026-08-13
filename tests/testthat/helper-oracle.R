# A d-separation oracle for deterministic tests. `edges` is a two-column
# character matrix of (parent, child) pairs. The returned function plugs
# into make_ci_tester(test = ...): it reports p-value 1 when i and j are
# d-separated given S in the DAG, and 0 otherwise, so test decisions are
# exact and independent of any generated data.
make_dsep_oracle <- function(edges) {
  parents <- function(v) edges[edges[, 2] == v, 1]
  ancestors <- function(nodes) {
    anc <- character(0); stack <- nodes
    while (length(stack)) {
      v <- stack[[1]]; stack <- stack[-1]
      if (v %in% anc) next
      anc <- c(anc, v); stack <- c(stack, parents(v))
    }
    anc
  }
  function(i, j, S, data) {
    if (i %in% S || j %in% S) return(1)
    keep <- ancestors(c(i, j, S))
    e <- edges[edges[, 1] %in% keep & edges[, 2] %in% keep, , drop = FALSE]
    und <- rbind(e, e[, 2:1, drop = FALSE])
    for (v in unique(e[, 2])) {                     # marry parents
      pa <- unique(e[e[, 2] == v, 1])
      if (length(pa) > 1) {
        prs <- t(utils::combn(pa, 2))
        und <- rbind(und, prs, prs[, 2:1, drop = FALSE])
      }
    }
    und <- und[!(und[, 1] %in% S) & !(und[, 2] %in% S), , drop = FALSE]
    seen <- character(0); stack <- i
    while (length(stack)) {
      v <- stack[[1]]; stack <- stack[-1]
      if (v == j) return(0)
      if (v %in% seen) next
      seen <- c(seen, v)
      stack <- c(stack, und[und[, 1] == v, 2])
    }
    1
  }
}

# dummy data frame carrying only the column names the tester needs
oracle_frame <- function(vars) {
  d <- as.data.frame(matrix(stats::rnorm(5 * length(vars)), 5, length(vars)))
  names(d) <- vars
  d
}
