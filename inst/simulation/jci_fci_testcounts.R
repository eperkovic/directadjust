# jci_fci_testcounts.R — fill the FCI + GAC row of the test-count table
# WITHOUT rerunning the full simulation. Uses the same seed and the same
# generation code as jci_simulation_rebuild.R; since the analysis methods
# consume no randomness, the datasets regenerated here are identical to the
# full run's. Only plain FCI with a counting CI test runs on each dataset.
#   source("jci_fci_testcounts.R")     # a few minutes

suppressMessages(library(pcalg))

N_MODELS <- 100
NS       <- c(100, 1000, 10000)
set.seed(20260716)                     # same seed as the full run

p   <- 15
Wv  <- 6:12
X1i <- 13; X2i <- 14; Yi <- 15
obs <- c(Wv, X1i, X2i, Yi)
obs_name <- c(paste0("W", 1:7), "X1", "X2", "Y")

gen_dag <- function(force_effect) {
  A <- matrix(0L, p, p)
  for (i in 2:p) {
    k <- sample(0:min(i - 1, 3), 1)
    if (k > 0) A[i, sample(seq_len(i - 1), k)] <- 1L
  }
  A[Yi, c(X1i, X2i)] <- if (force_effect) 1L else 0L
  A
}
gen_sem <- function(A) {
  B <- matrix(0, p, p)
  nz <- which(A == 1L, arr.ind = TRUE)
  B[nz] <- sample(c(-1, 1), nrow(nz), TRUE) * runif(nrow(nz), 0.1, 1.5)
  list(B = B, s2 = runif(p, 0.5, 1))
}
sim_data <- function(sem, n) {
  V <- matrix(0, n, p)
  for (i in 1:p)
    V[, i] <- V %*% sem$B[i, ] + rnorm(n, 0, sqrt(sem$s2[i]))
  colnames(V) <- paste0("V", 1:p)
  V
}
dag_dsep <- function(A, x, y, S) {
  anc <- unique(c(x, y, S)); repeat {
    more <- unique(c(anc, unlist(lapply(anc, function(v) which(A[v, ] == 1L)))))
    if (length(more) == length(anc)) break
    anc <- more
  }
  M <- A[anc, anc, drop = FALSE]
  U <- (M + t(M)) > 0
  for (i in seq_along(anc)) {
    pa <- which(M[i, ] == 1L)
    if (length(pa) > 1) U[pa, pa] <- TRUE
  }
  diag(U) <- FALSE
  keep <- !(anc %in% S)
  U <- U[keep, keep, drop = FALSE]; ids <- anc[keep]
  xi <- match(x, ids); yi <- match(y, ids)
  if (is.na(xi) || is.na(yi)) return(TRUE)
  reach <- xi; frontier <- xi
  while (length(frontier) > 0) {
    v <- frontier[1]; frontier <- frontier[-1]
    nb <- setdiff(which(U[v, ]), reach)
    if (yi %in% nb) return(FALSE)
    reach <- c(reach, nb); frontier <- c(frontier, nb)
  }
  TRUE
}
valid_set <- function(A, Zfull) {
  Ac <- A; Ac[Yi, c(X1i, X2i)] <- 0L
  dag_dsep(Ac, X1i, Yi, c(Zfull, X2i)) && dag_dsep(Ac, X2i, Yi, c(Zfull, X1i))
}
combn_list <- function(v, m) {
  if (m == 0) return(list(integer(0)))
  if (length(v) < m) return(list())
  asplit(combn(v, m), 2)
}
dag_set_exists <- function(A) {
  for (sz in 0:length(Wv)) for (Z in combn_list(Wv, sz))
    if (valid_set(A, as.integer(Z))) return(TRUE)
  FALSE
}

# ---- regenerate the same models (same RNG sequence as the full run) ---------
models <- list(S1 = list(), S2 = list(), S3 = list())
while (any(vapply(models, length, 1L) < N_MODELS)) {
  need_eff <- length(models$S1) < N_MODELS || length(models$S2) < N_MODELS
  if (need_eff) {
    A <- gen_dag(force_effect = TRUE)
    key <- if (dag_set_exists(A)) "S1" else "S2"
  } else {
    A <- gen_dag(force_effect = FALSE)
    key <- "S3"
  }
  if (length(models[[key]]) < N_MODELS)
    models[[key]] <- c(models[[key]], list(list(A = A, sem = gen_sem(A))))
}
cat("models regenerated\n")

# ---- count FCI's tests on every dataset --------------------------------------
res <- list()
for (setting in names(models)) for (mi in seq_along(models[[setting]])) {
  m <- models[[setting]][[mi]]
  for (n in NS) {
    dat <- sim_data(m$sem, n)[, obs, drop = FALSE]
    colnames(dat) <- obs_name
    C <- cor(dat)
    cnt <- 0L
    counting_test <- function(x, y, S, suffStat) {
      cnt <<- cnt + 1L
      gaussCItest(x, y, S, suffStat)
    }
    fit <- tryCatch(suppressWarnings(
      pcalg::fci(list(C = C, n = n), counting_test, alpha = 0.05,
                 labels = obs_name, verbose = FALSE)),
      error = function(e) NULL)
    res[[length(res) + 1]] <- data.frame(
      setting = setting, model = mi, n = n,
      ntests = if (is.null(fit)) NA else cnt)
  }
  if (mi %% 20 == 0) cat(setting, mi, "done\n")
}
res <- do.call(rbind, res)
saveRDS(res, "jci_fci_testcounts.rds")

cat("\nFCI + GAC test counts (min / median / mean / max):\n")
print(aggregate(ntests ~ setting + n, res,
                function(x) c(min = min(x), median = median(x),
                              mean = round(mean(x)), max = max(x))))
