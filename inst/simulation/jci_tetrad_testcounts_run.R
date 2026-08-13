# jci_tetrad_testcounts_run.R — counts-only pass: exact CI-test counts for
# the Tetrad arms (GFCI with tiers, and Tetrad's FCI with tiers) on the SAME
# 300 models x 3 sample sizes as the main rebuild.
#
# Dataset identity: the main run's datasets are reproduced exactly by
# replaying the RNG stream - same seed, same generation loop (verified
# against jci_models.rds before anything runs), same sim_data order. The
# methods themselves consume no RNG, so the datasets here are bit-identical
# to the ones the main run's counts and outcomes came from.
#
# Output: jci_tetrad_testcounts.rds + aggregate tables formatted like the
# the supplement's test-count table (min / median / mean / max per setting at each n).
# Runtime: roughly 15-25 minutes (1,800 Tetrad searches).
#
#   source("jci_tetrad_testcounts_run.R")

library(reticulate)
use_virtualenv("r-tetrad", required = TRUE)

py_run_string("
import jpype
import jpype.imports
import pytetrad.tools.translate as tr

from edu.cmu.tetrad.data import Knowledge
import edu.cmu.tetrad.search as _search

def _get_class(mod, *names):
    for n in names:
        c = getattr(mod, n, None)
        if c is not None:
            return c
    return None

GFCI_CLS = _get_class(_search, 'GFci', 'Gfci', 'GFCI')
FCI_CLS  = _get_class(_search, 'Fci', 'FCI')

try:
    from edu.cmu.tetrad.search.test import IndTestFisherZ
    IFACE = 'edu.cmu.tetrad.search.test.IndependenceTest'
except ImportError:
    from edu.cmu.tetrad.search import IndTestFisherZ
    IFACE = 'edu.cmu.tetrad.search.IndependenceTest'

try:
    from edu.cmu.tetrad.search.score import SemBicScore
except ImportError:
    from edu.cmu.tetrad.search import SemBicScore

class CountingTest:
    def __init__(self, inner):
        self.inner = inner
        self.count = 0
    def checkIndependence(self, *args):
        self.count += 1
        return self.inner.checkIndependence(*args)
    def __getattr__(self, name):
        return getattr(self.inner, name)

def tetrad_search_with_count(df, alpha, tiers, algo='gfci'):
    data  = tr.pandas_data_to_tetrad(df)
    inner = IndTestFisherZ(data, float(alpha))
    counter = CountingTest(inner)
    proxy = jpype.JProxy(IFACE, inst=counter)
    kn = Knowledge()
    for tier, names in enumerate(tiers):
        if isinstance(names, str):
            names = [names]
        for nm in names:
            kn.addToTier(tier, nm)
    if algo == 'gfci':
        try:
            score = SemBicScore(data, True)
        except TypeError:
            score = SemBicScore(data)
        alg = GFCI_CLS(proxy, score)
    else:
        alg = FCI_CLS(proxy)
    alg.setKnowledge(kn)
    alg.search()
    return counter.count
")

# ---- exact replay of the main run's generation and data stream --------------
N_MODELS <- 100
NS       <- c(100, 1000, 10000)
set.seed(20260716)                       # the main run's seed

p        <- 15
lat      <- 1:5
Wv       <- 6:12
X1i <- 13; X2i <- 14; Yi <- 15
obs      <- c(Wv, X1i, X2i, Yi)
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
  for (i in 1:p) V[, i] <- V %*% sem$B[i, ] + rnorm(n, 0, sqrt(sem$s2[i]))
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

# hard check: the replayed models must equal the saved ones bit for bit
saved <- readRDS("jci_models.rds")
stopifnot(isTRUE(all.equal(models, saved)))
cat("model replay verified against jci_models.rds\n")

# ---- the counting pass -------------------------------------------------------
tiers <- list(as.list(paste0("W", 1:7)), "X1", "X2", "Y")
res <- list()
for (setting in names(models)) for (mi in seq_along(models[[setting]])) {
  m <- models[[setting]][[mi]]
  for (n in NS) {
    dat <- sim_data(m$sem, n)[, obs, drop = FALSE]   # same RNG order as main run
    d2 <- as.data.frame(dat); colnames(d2) <- obs_name
    n_gfci <- tryCatch(py$tetrad_search_with_count(r_to_py(d2), 0.05, tiers, "gfci"),
                       error = function(e) NA)
    n_tfci <- tryCatch(py$tetrad_search_with_count(r_to_py(d2), 0.05, tiers, "fci"),
                       error = function(e) NA)
    res[[length(res) + 1]] <- data.frame(
      setting = setting, model = mi, n = n,
      method = c("GFCI + GAC", "FCItiers (Tetrad)"),
      ntests = c(n_gfci, n_tfci))
  }
  if (mi %% 10 == 0) cat(setting, mi, "of", N_MODELS, "\n")
}
res <- do.call(rbind, res)
saveRDS(res, "jci_tetrad_testcounts.rds")

cat("\ncounts by method, setting, n (min / median / mean / max):\n")
print(aggregate(ntests ~ method + n + setting, res,
                function(x) c(min = min(x), median = median(x),
                              mean = round(mean(x)), max = max(x))))
cat("done; saved to jci_tetrad_testcounts.rds\n")
