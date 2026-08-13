# jci_fig8_covariates.R — replicate the published Figure 8: running time
# versus the NUMBER OF COVARIATES (5, 7, 10) at fixed n = 1,000.
#
# REWRITTEN (17 Jul) so that every arm runs its FULL pipeline, matching the
# main rebuild. Timing is decomposed for the graph-based arms:
#   secs_search = discovery (FCI / tFCI / GFCI search, incl. amat extraction)
#   secs_cert   = decision + GAC set LISTING via dagitty (adjustmentSets,
#                 type = "all"), i.e. the x + y split
#   secs_total  = secs_search + secs_cert   <- what the figure plots
# The R1 arms run the full search (Build: complete enumeration, all sets;
# Combine: enumeration + inclusion-minimal reduction + CI-minimal pruning +
# unions), which IS their certification, so secs_cert is NA and
# secs_total = secs_search.
#   source("jci_fig8_covariates.R")   # self-contained; tens of minutes,
#                                     # dominated by the 10-covariate cells
suppressMessages({ library(pcalg); library(causalDisco); library(ggplot2)
                   library(dagitty) })

N_MODELS <- 50          # models per covariate count (effect arm, unscreened)
N_FIXED  <- 1000        # the published figure's sample size
W_COUNTS <- c(5, 7, 10)
INCLUDE_TETRAD <- TRUE  # adds BOTH Tetrad arms (FCItiers and GFCI, with tiers)
                        # via the counting-proxy path, matching the main
                        # rebuild; needs the r-tetrad virtualenv
INCLUDE_TFCI <- FALSE   # the causalDisco tFCI arm is NOT part of the
                        # manuscript's comparison (FCItiers via Tetrad is);
                        # keep FALSE unless comparing implementations
if (INCLUDE_TETRAD) {
  reticulate::use_virtualenv("r-tetrad", required = TRUE)
  reticulate::py_run_string("
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

def tetrad_run(df, alpha, tiers, algo='gfci'):
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
    g = alg.search()
    return {'count': counter.count, 'graph': str(g)}
")
}
set.seed(20260717)

combn_list <- function(v, m) {
  if (m == 0) return(list(integer(0)))
  if (length(v) < m) return(list())
  asplit(combn(v, m), 2)
}

# PAG amat (pcalg codes) -> dagitty pag { } object; generic in dim/names
pag_to_dagitty <- function(amat) {
  nm <- colnames(amat)
  tok <- function(mi, mj) {
    key <- paste0(mi, mj)
    switch(key,
      "32" = "->", "23" = "<-", "22" = "<->", "33" = "--",
      "12" = "@->", "21" = "<-@", "11" = "@-@",
      stop("unmapped PAG edge marks: ", key))
  }
  p <- ncol(amat); stmts <- nm
  for (i in 1:(p - 1)) for (j in (i + 1):p)
    if (amat[i, j] != 0 || amat[j, i] != 0)
      stmts <- c(stmts, paste(nm[i], tok(amat[j, i], amat[i, j]), nm[j]))
  dagitty::dagitty(paste0("pag { ", paste(stmts, collapse = " ; "), " }"))
}

# Tetrad's printed graph -> pcalg PAG amat; generic in the name vector
tetrad_string_to_amat <- function(txt, obs_name) {
  amat <- matrix(0L, length(obs_name), length(obs_name),
                 dimnames = list(obs_name, obs_name))
  edge_re <- "^\\s*\\d+\\.\\s+(\\S+)\\s+(<->|o->|<-o|o-o|-->|<--|---)\\s+(\\S+)"
  mk <- c("<" = 2L, ">" = 2L, "o" = 1L, "-" = 3L)
  for (ln in unlist(strsplit(txt, "\n"))) {
    m <- regmatches(ln, regexec(edge_re, ln))[[1]]
    if (length(m) != 4) next
    from <- m[2]; e <- m[3]; to <- m[4]
    if (!(from %in% obs_name) || !(to %in% obs_name)) next
    amat[to, from] <- mk[[substr(e, 1, 1)]]
    amat[from, to] <- mk[[substr(e, nchar(e), nchar(e))]]
  }
  amat
}

# decision + GAC set listing (the certification component, timed as "y").
# Mirrors pag_decision in the rebuild, generalized to nW; dagitty certifier.
pag_cert <- function(amat, nW) {
  xs <- c(nW + 1, nW + 2); ys <- nW + 3
  mX <- amat[ys, xs]
  if (any(mX == 2) || all(mX == 0)) return("negative")
  if (any(mX == 3)) {
    g <- tryCatch(pag_to_dagitty(amat), error = function(e) NULL)
    if (is.null(g)) return("unknown")
    sets <- tryCatch(
      dagitty::adjustmentSets(g, exposure = c("X1", "X2"), outcome = "Y",
                              type = "all"),
      error = function(e) NULL)
    return(if (!is.null(sets) && length(sets) > 0) "found" else "unknown")
  }
  "unknown"
}

run_for_W <- function(nW) {
  p   <- 5 + nW + 3                    # 5 latents + W's + X1, X2, Y
  Wv  <- 6:(5 + nW)
  X1i <- p - 2; X2i <- p - 1; Yi <- p
  obs <- c(Wv, X1i, X2i, Yi)
  obs_name <- c(paste0("W", 1:nW), "X1", "X2", "Y")
  Wd <- seq_len(nW); X1d <- nW + 1; X2d <- nW + 2; Yd <- nW + 3

  gen_dag <- function() {
    A <- matrix(0L, p, p)
    for (i in 2:p) {
      k <- sample(0:min(i - 1, 3), 1)
      if (k > 0) A[i, sample(seq_len(i - 1), k)] <- 1L
    }
    A[Yi, c(X1i, X2i)] <- 1L           # effect arm
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
  make_tester <- function(dat) {       # R1 methods' level, as in the code
    C <- cor(dat); n <- nrow(dat)
    function(i, j, S) gaussCItest(i, j, S, list(C = C, n = n)) > 0.01
  }

  # ---- FULL R1 Build: complete enumeration, returns all passing sets ------
  r1_build_full <- function(ci) {
    found <- list()
    for (W1 in Wd) for (W2 in setdiff(Wd, W1)) {
      rest <- setdiff(Wd, c(W1, W2))
      for (sz in 0:min(nW - 2, length(rest))) for (Z in combn_list(rest, sz)) {
        Z <- as.integer(Z)
        if (!ci(W1, Yd, Z) && ci(W1, Yd, c(Z, X1d)) &&
            !ci(W2, Yd, c(Z, X1d)) && ci(W2, Yd, c(Z, X1d, X2d)))
          found[[length(found) + 1]] <- Z
      }
    }
    unique(found)
  }

  # ---- FULL R1 Combine: enumeration + minimal reduction + pruning + unions
  r1_combine_full <- function(ci) {
    collect <- function(Xd, extra) {
      hits <- list()
      for (Wi in Wd) {
        pool <- c(setdiff(Wd, Wi), extra)
        for (sz in 0:min(nW - 1, length(pool)))
          for (Ti in combn_list(pool, sz)) {
            Ti <- as.integer(Ti)
            if (!ci(Wi, Yd, Ti) && ci(Wi, Yd, c(Ti, Xd)))
              hits[[length(hits) + 1]] <- list(W = Wi, T = Ti)
          }
      }
      hits
    }
    minimal_Z <- function(Wi, Ti, Xd) {
      for (mz in 0:length(Ti)) for (S in combn_list(Ti, mz)) {
        S <- as.integer(S)
        if (!ci(Wi, Yd, S) && ci(Wi, Yd, c(S, Xd))) return(S)
      }
      NULL
    }
    h1 <- collect(X1d, extra = integer(0))
    if (length(h1) == 0) return(list())
    h2 <- collect(X2d, extra = X1d)
    if (length(h2) == 0) return(list())
    skey <- function(v) paste(v, collapse = ",")
    dedupe_min <- function(h) {
      h <- h[!duplicated(vapply(h, function(x) skey(x$T), ""))]
      Ts <- lapply(h, `[[`, "T")
      keep <- vapply(seq_along(Ts), function(i)
        !any(vapply(seq_along(Ts), function(j)
          j != i && length(Ts[[j]]) < length(Ts[[i]]) &&
          all(Ts[[j]] %in% Ts[[i]]), TRUE)), TRUE)
      h[keep]
    }
    m1 <- dedupe_min(h1); m2 <- dedupe_min(h2)
    Zs <- list()
    for (a in m1) for (b in m2) {
      Z1 <- minimal_Z(a$W, a$T, X1d); Z2 <- minimal_Z(b$W, b$T, X2d)
      if (is.null(Z1) || is.null(Z2)) next
      Zs[[length(Zs) + 1]] <- sort(setdiff(union(Z1, Z2), c(X1d, X2d)))
    }
    unique(Zs)
  }

  # ---- discovery searches, each returning a PAG amat (or NULL) ------------
  tier_knowledge <- function(d2) {
    if (nW == 5) knowledge(d2, tier(
      covariates ~ c(W1, W2, W3, W4, W5),
      treat1 ~ X1, treat2 ~ X2, outcome ~ Y))
    else if (nW == 7) knowledge(d2, tier(
      covariates ~ c(W1, W2, W3, W4, W5, W6, W7),
      treat1 ~ X1, treat2 ~ X2, outcome ~ Y))
    else knowledge(d2, tier(
      covariates ~ c(W1, W2, W3, W4, W5, W6, W7, W8, W9, W10),
      treat1 ~ X1, treat2 ~ X2, outcome ~ Y))
  }
  fci_search <- function(dat) {
    fit <- tryCatch(suppressWarnings(
      pcalg::fci(list(C = cor(dat), n = nrow(dat)), gaussCItest,
                 alpha = 0.05, labels = obs_name, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) NULL else fit@amat
  }
  tfci_search <- function(dat) {
    d2 <- as.data.frame(dat); colnames(d2) <- obs_name
    fit <- tryCatch(
      disco(d2, tfci(engine = "causalDisco", test = "fisher_z", alpha = 0.05),
            knowledge = tier_knowledge(d2)),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    edges_df <- tryCatch(as.data.frame(fit$caugi@edges), error = function(e) NULL)
    if (is.null(edges_df)) return(NULL)
    mark_left  <- function(ch) ifelse(ch == "-", 3L,
                               ifelse(ch == "<", 2L,
                               ifelse(ch == "o", 1L, 0L)))
    mark_right <- function(ch) ifelse(ch == "-", 3L,
                               ifelse(ch == ">", 2L,
                               ifelse(ch == "o", 1L, 0L)))
    amat <- matrix(0L, length(obs_name), length(obs_name),
                   dimnames = list(obs_name, obs_name))
    for (i in seq_len(nrow(edges_df))) {
      from <- as.character(edges_df$from[i]); to <- as.character(edges_df$to[i])
      e <- as.character(edges_df$edge[i])
      if (!(from %in% obs_name) || !(to %in% obs_name)) next
      amat[to,   from] <- mark_left(substr(e, 1, 1))
      amat[from, to]   <- mark_right(substr(e, nchar(e), nchar(e)))
    }
    amat
  }
  tetrad_search <- function(dat, algo) {
    d2 <- as.data.frame(dat); colnames(d2) <- obs_name
    tiers_nw <- c(list(as.list(paste0("W", seq_len(nW)))),
                  list("X1", "X2", "Y"))
    out <- tryCatch(
      reticulate::py$tetrad_run(reticulate::r_to_py(d2), 0.05, tiers_nw, algo),
      error = function(e) { cat("  [tetrad]", conditionMessage(e), "\n"); NULL })
    if (is.null(out)) return(NULL)
    tryCatch(tetrad_string_to_amat(paste(out$graph, collapse = "\n"), obs_name),
             error = function(e) NULL)
  }
  gfci_search  <- function(dat) tetrad_search(dat, "gfci")
  tfciT_search <- function(dat) tetrad_search(dat, "fci")

  # timed x + y wrapper for the graph arms
  timed_arm <- function(search_fun, dat) {
    t_search <- system.time(amat <- search_fun(dat))["elapsed"]
    t_cert <- if (is.null(amat)) NA_real_ else
      system.time(pag_cert(amat, nW))["elapsed"]
    c(search = unname(t_search), cert = unname(t_cert))
  }

  out <- list()
  if (INCLUDE_TETRAD) invisible(gfci_search(   # warm-up so the JVM boot does
    sim_data(gen_sem(gen_dag()), 200)[, obs])) # not land in model 1's timing
  for (m in seq_len(N_MODELS)) {
    A <- gen_dag(); sem <- gen_sem(A)
    dat <- sim_data(sem, N_FIXED)[, obs, drop = FALSE]
    t_build   <- system.time(r1_build_full(make_tester(dat)))["elapsed"]
    t_combine <- system.time(r1_combine_full(make_tester(dat)))["elapsed"]
    a_fci   <- timed_arm(fci_search,  dat)
    a_tfci  <- if (INCLUDE_TFCI) timed_arm(tfci_search, dat) else NULL
    a_tfciT <- if (INCLUDE_TETRAD) timed_arm(tfciT_search, dat) else NULL
    a_gfci  <- if (INCLUDE_TETRAD) timed_arm(gfci_search,  dat) else NULL
    rows <- data.frame(
      nW = nW, model = m,
      method = c("R1 Build", "R1 Combine", "FCI + GAC",
                 if (INCLUDE_TFCI) "tFCI + GAC",
                 if (INCLUDE_TETRAD) c("FCItiers (Tetrad)", "GFCI + GAC")),
      secs_search = c(t_build, t_combine, a_fci["search"],
                      if (INCLUDE_TFCI) a_tfci["search"],
                      if (INCLUDE_TETRAD) c(a_tfciT["search"], a_gfci["search"])),
      secs_cert   = c(NA, NA, a_fci["cert"],
                      if (INCLUDE_TFCI) a_tfci["cert"],
                      if (INCLUDE_TETRAD) c(a_tfciT["cert"], a_gfci["cert"])))
    rows$secs_total <- rows$secs_search +
      ifelse(is.na(rows$secs_cert), 0, rows$secs_cert)
    out[[m]] <- rows
    if (m %% 10 == 0) cat("W =", nW, ":", m, "of", N_MODELS, "\n")
  }
  do.call(rbind, out)
}
# NOTE on the warm-up draw: it consumes RNG, so this run's graphs differ from
# earlier sweeps'. Fine here - the sweep's models were never persisted, and
# all arms see identical data WITHIN the run, which is the comparison shown.

res8 <- do.call(rbind, lapply(W_COUNTS, run_for_W))
saveRDS(res8, "jci_fig8_times.rds")

summ <- do.call(rbind, by(res8, list(res8$nW, res8$method), function(d)
  data.frame(nW = d$nW[1], method = d$method[1],
             med = median(d$secs_total),
             lo = quantile(d$secs_total, .1), hi = quantile(d$secs_total, .9))))
summ$method <- droplevels(factor(summ$method,
  levels = c("R1 Build", "R1 Combine", "FCI + GAC", "tFCI + GAC",
             "FCItiers (Tetrad)", "GFCI + GAC")))

fig8 <- ggplot(summ, aes(nW, med, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, color = NA) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = W_COUNTS) +
  theme_minimal() +
  labs(x = "number of covariates",
       y = "seconds per dataset, full pipeline (n = 1,000)",
       color = NULL, fill = NULL)
ggsave("fig8_covariates_rebuild.pdf", fig8, width = 7, height = 4.5)

# the x + y decomposition: discovery vs decision + GAC listing
cat("\nTOTAL time (search + certification):\n")
print(aggregate(secs_total ~ method + nW, res8,
                function(x) c(median = round(median(x), 3),
                              mean = round(mean(x), 3))))
cat("\nsearch component (x):\n")
print(aggregate(secs_search ~ method + nW, res8,
                function(x) c(median = round(median(x), 3),
                              mean = round(mean(x), 3))))
cat("\ncertification component (y; graph arms only):\n")
cert_rows <- subset(res8, !is.na(secs_cert))
if (nrow(cert_rows) > 0)
  print(aggregate(secs_cert ~ method + nW, cert_rows,
                  function(x) c(median = round(median(x), 3),
                                mean = round(mean(x), 3))))
cat("done; figure in fig8_covariates_rebuild.pdf\n")
