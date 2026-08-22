# jci_simulation_rebuild.R — the simulation study of LaPlante,
#
# TIERED-KNOWLEDGE VALIDITY GUARD. The tiers imply hard constraints on the
# marks of any graph a knowledge-aware search returns: on an edge between
# nodes of different tiers, the endpoint at the later tier must be an
# arrowhead (the later node cannot be an ancestor of the earlier one). Not
# every released Tetrad build enforces this, and a build that drops such
# marks silently corrupts the Setting 1/2 screen. This script therefore
# verifies every returned graph against these constraints: the ground-truth
# screen refuses to run if any repair is needed, and the data arms repair
# and tally (a correct build needs zero repairs). Verified clean under
# Tetrad 7.6.11; see session_info.txt. Note that py-tetrad's package version
# does not identify the bundled Tetrad engine -- use the jar version and
# SHA-256 recorded in session_info.txt.

# Triantafillou & Perkovic, "Data-Driven Adjustment for Multiple Treatments"
# (Journal of Causal Inference), Section 5.
#
# Design:
#   * random DAGs: 15 nodes; causal order = 5 latents, 7 covariates W, X1, X2, Y;
#     |pa(Vi)| ~ Unif{0..min(#predecessors,3)}, parents uniform among predecessors;
#     then X1->Y and X2->Y edges forced by arm: present in S1/S2, absent in S3
#   * settings: S1 = a fully observed valid adjustment set exists;
#     S2 = none exists; S3 = no-effect arm, no screening.
#     Validity of Z: cut X1->Y, X2->Y in the true DAG, then
#     X1 dsep Y | Z u {X2}  and  X2 dsep Y | Z u {X1}
#   * SEM: b_ij ~ Unif(+-[0.1,1.5]), sigma_i^2 ~ Unif[.5,1];
#     true effect = B[Y,X1] + B[Y,X2]
#   * R1 Build: full enumeration, |Z| <= 5, returns ALL passing sets
#   * R1 Combine: per treatment all passing (witness, T), |T| <= 6,
#     inclusion-minimal T's, minimal Z_i per T, all unions Z1uZ2
#   * FCItiers: tiered background knowledge; negative iff the learned PAG has
#     an arrowhead into X1 or X2 on its edge with Y, or neither X shares an
#     edge with Y; adjustment sets = all inclusion-minimal GAC sets
#   * outcomes: correct = found & AT LEAST ONE set valid; estimation = mean over
#     found sets of joint lm coefficient sum; negative -> 0; unknown -> none
#   * CI-test counters per method and an MSE panel alongside the bias plot
#
# Requires: pcalg, causalDisco, ggplot2, dagitty.
# First run: set TEST_MODE <- TRUE.

suppressMessages({ library(pcalg); library(causalDisco); library(ggplot2)
                   library(dagitty) })   # dagitty certifies adjustment sets (see gac_bug_repro.R)

TEST_MODE <- FALSE # TRUE: 10 models/setting, n in {100, 1000}
N_MODELS  <- if (TEST_MODE) 10 else 100
NS        <- if (TEST_MODE) c(100, 1000) else c(100, 1000, 10000)
# ---------------------------------------------------------------------------
# SETTING CLASSIFICATION - the single most consequential design choice, made
# explicit here because the original implementation contains both variants and they
# produce very different Setting 1 populations:
#
#   SCREENING = "pag"  (the study's primary configuration; Section 5 text:
#       "we note if there is an adjustment set relative to ({X1,X2},Y) that
#        we can learn from observed data ... we consider the PAG that
#        represents in/dependencies from such data as well as knowledge that
#        W < X1 < X2 < Y")
#     Setting 1 <=> the GROUND-TRUTH TIERED PAG (oracle FCI on the true DAG
#     with the latents marked latent + the four tiers) has a GAC-certifiable
#     adjustment set; Setting 2 <=> effect-arm graph where it has none.
#     Identifiability is learnable from observables by construction, which
#     is why the R1 bars reach ~85% in Setting 1.
#
#   SCREENING = "dag"  (sensitivity variant)
#     Setting 1 <=> some observed Z is valid in the TRUE DAG (edge-cut
#     d-separation). Weaker: keeps models where the set exists in the DAG
#     but cannot be certified from data, so Setting 1 success rates are
#     systematically lower (R1 methods plateau near 50%).
#
# Outputs of "dag" runs carry the suffix "_dagscreen" and are kept as the
# sensitivity variant; the "pag" run owns the unsuffixed (primary) files.
# Setting 3 (no-effect arm) is unscreened in both variants.
#
# SCREENING KNOWLEDGE: the ground-truth screening uses oracle FCI on the
# true DAG with the four tiers as background knowledge.
SCREENING <- "pag"
SCREEN_KNOWLEDGE <- "paper"
# NOTE: model generation and data draw the same RNG stream regardless of the
# testing levels, so every run below uses IDENTICAL graphs and datasets and
# results are comparable arm by arm across runs.
# TESTING LEVELS -------------------------------------------------------------
# R1_ALPHA_DEP : R1 methods declare DEPENDENCE when p <= this
# R1_ALPHA_IND : R1 methods declare INDEPENDENCE when p > this (setting it
#                above R1_ALPHA_DEP leaves an undecided gap, the EHS design;
#                setting it equal makes every test decisive)
# COMP_ALPHA   : test level for ALL graph-based arms (FCI, FCItiers, GFCI,
#                and the R-native tFCI)
# RUN_TAG      : appended to every output filename (with "_"); keeps runs
#                from overwriting each other
#
# The manuscript's configuration: every method tests at alpha = 0.05.
# With RUN_TAG = "" a full run writes jci_models.rds and jci_sim_results.rds,
# the files shipped in inst/results/. To explore other levels, change the
# alphas and set a RUN_TAG so the shipped files are not overwritten.
R1_ALPHA_DEP <- 0.05
R1_ALPHA_IND <- 0.05
COMP_ALPHA   <- 0.05
RUN_TAG      <- ""
INCLUDE_GFCI <- TRUE    # TRUE adds a "GFCI + GAC" arm via py-tetrad/reticulate
                         # (needs a Java JDK >= 17 and, once,
                         #  reticulate::py_install("py-tetrad", pip = TRUE)).
                         # Adding the arm consumes no RNG, so the other arms
                         # reproduce exactly. Smoke-test gfci_gac() on one
                         # dataset before a full run.
INCLUDE_TFCI <- FALSE   # the R-native causalDisco tFCI arm is NOT part of
                        # the manuscript's comparison (FCItiers via Tetrad
                        # is); keep FALSE unless comparing implementations
INCLUDE_TETRAD_FCI <- TRUE   # "FCItiers (Tetrad)": Tetrad's FCI with tiered
                             # Knowledge - the same library call the original
                             # study made, with exact test counts via the
                             # counting proxy.
if (INCLUDE_GFCI || INCLUDE_TETRAD_FCI || SCREENING == "pag") {
  # One-time setup (old system Pythons won't do; reticulate needs its own):
  #   reticulate::install_python("3.11:latest")
  #   reticulate::virtualenv_create("r-tetrad", version = "3.11:latest")
  #   reticulate::virtualenv_install("r-tetrad", "py-tetrad")
  # plus a Java JDK >= 17 with JAVA_HOME set (Sys.getenv("JAVA_HOME")).
  reticulate::use_virtualenv("r-tetrad", required = TRUE)
  # Counting proxy: Tetrad accepts any object implementing its
  # IndependenceTest interface; we hand it a jpype.JProxy around Tetrad's own
  # IndTestFisherZ that counts checkIndependence() calls and delegates the
  # rest. Returns the graph string AND the exact test count.
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

# ---- ground-truth tiered PAG (oracle FCI), for the primary screening ----
from edu.cmu.tetrad.graph import GraphNode, EdgeListGraph, NodeType
try:
    from edu.cmu.tetrad.search.test import MsepTest as _MSEP
except ImportError:
    try:
        from edu.cmu.tetrad.search import MsepTest as _MSEP
    except ImportError:
        from edu.cmu.tetrad.search import IndTestDSep as _MSEP

def gt_pag(directed, latents, tiers):
    '''Oracle tiered PAG of a true DAG: nodes in `latents` are marked
    LATENT, FCI runs with the d-separation oracle (MsepTest) and the tier
    knowledge, mirroring Tetrad_wrapper('graphbg','fci',...) in the
    original study. Returns the PAG's edge string.'''
    g = EdgeListGraph()
    nodes = {}
    names = set()
    for a, b in directed:
        names.add(a); names.add(b)
    for nm in sorted(names):
        n = GraphNode(nm)
        if nm in latents:
            n.setNodeType(NodeType.LATENT)
        nodes[nm] = n
        g.addNode(n)
    for a, b in directed:
        g.addDirectedEdge(nodes[a], nodes[b])
    test = _MSEP(g)
    kn = Knowledge()
    for tier, nms in enumerate(tiers):
        if isinstance(nms, str):
            nms = [nms]
        for nm in nms:
            kn.addToTier(tier, nm)
    alg = FCI_CLS(test)
    alg.setKnowledge(kn)
    return str(alg.search())
")
  tiers_py <- list(as.list(paste0("W", 1:7)), "X1", "X2", "Y")

  # report which Tetrad actually answered (goes in session_info / the log)
  invisible(tryCatch(reticulate::py_run_string("
import importlib.metadata as _md
try: print('py-tetrad version:', _md.version('py-tetrad'))
except Exception: print('py-tetrad version: unknown')
try:
    from edu.cmu.tetrad.util import Version as _V
    print('Tetrad jar version:', str(_V.currentViewableVersion()))
except Exception: print('Tetrad jar version: unknown')
"), error = function(e) invisible(NULL)))

# ---- TIER-MARK GUARD --------------------------------------------------------
# The tiers are hard knowledge: on any edge between nodes of different tiers,
# the mark at the LATER-tier node must be an arrowhead (the later node cannot
# be an ancestor of the earlier one in any tier-consistent MAG). A circle
# there means the search failed to apply knowledge it was given; the repair
# (circle -> arrowhead) is sound in every tier-consistent MAG. A tail there
# contradicts the tiers outright. amat coding: amat[i, j] = mark at j.
tier_of_node <- function(nm)
  ifelse(grepl("^W", nm), 1L, ifelse(nm == "X1", 2L, ifelse(nm == "X2", 3L, 4L)))
tier_guard_tally <- new.env()
tier_guard_tally$stamps <- 0L; tier_guard_tally$conflicts <- 0L
check_tier_marks <- function(amat) {
  nm <- rownames(amat); tiers <- tier_of_node(nm)
  stamps <- 0L; conflicts <- 0L
  for (a in 1:(nrow(amat) - 1)) for (b in (a + 1):ncol(amat)) {
    if (amat[a, b] == 0 && amat[b, a] == 0) next
    if (tiers[a] == tiers[b]) next
    later <- if (tiers[a] < tiers[b]) b else a
    other <- if (tiers[a] < tiers[b]) a else b
    m <- amat[other, later]                     # mark at the later-tier node
    if (m == 1) { stamps <- stamps + 1L; amat[other, later] <- 2L }
    else if (m == 3) conflicts <- conflicts + 1L
  }
  list(amat = amat, stamps = stamps, conflicts = conflicts)
}
tier_guard_report <- function()
  cat(sprintf(paste0("tier-mark guard: %d repaired mark(s), %d conflict(s) ",
                     "across all data-arm PAGs (both should be 0 on a ",
                     "correct Tetrad/causalDisco)\n"),
      tier_guard_tally$stamps, tier_guard_tally$conflicts))
  # warm-up (deterministic data, NO RNG consumed) so the JVM boot does not
  # land in the first dataset's recorded runtime
  dummy <- as.data.frame(matrix(sin(1:80), 20, 4))
  colnames(dummy) <- c("W1", "X1", "X2", "Y")
  invisible(tryCatch(
    reticulate::py$tetrad_run(reticulate::r_to_py(dummy), 0.05,
                              list("W1", "X1", "X2", "Y"), "fci"),
    error = function(e) NULL))
}
SUFFIX <- paste0(
  if (SCREENING == "dag") "_dagscreen"
  else if (SCREEN_KNOWLEDGE == "legacy") "_legacyscreen" else "",
  if (nzchar(RUN_TAG)) paste0("_", RUN_TAG) else "")
out_file <- function(name)   # suffix all outputs so runs never overwrite
  sub("\\.(pdf|rds)$", paste0(SUFFIX, ".\\1"), name)
set.seed(20260716)

p        <- 15
lat      <- 1:5
Wv       <- 6:12
X1i <- 13; X2i <- 14; Yi <- 15
obs      <- c(Wv, X1i, X2i, Yi)
obs_name <- c(paste0("W", 1:7), "X1", "X2", "Y")

# ---- generation ---------------------------------------------------------------
# Random DAGs: random parents from predecessors, then the X->Y edges are
# FORCED by the arm (probeffect 1 or 0). X1->X2 can occur; X's never point at W's.
gen_dag <- function(force_effect) {
  A <- matrix(0L, p, p)              # A[i, j] = 1  <=>  j -> i (parent j)
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

true_effect <- function(sem) {       # joint do(): direct coefficients only
  sem$B[Yi, X1i] + sem$B[Yi, X2i]    # (X1->X2->Y is cut by do(X2); no other
}                                    #  mediators possible in this design)

# ---- d-separation in the true DAG (moralization) --------------------------------
# A is parent coding: A[i, j] = 1 <=> j -> i.
dag_dsep <- function(A, x, y, S) {
  # ancestors of {x, y} u S
  anc <- unique(c(x, y, S)); repeat {
    more <- unique(c(anc, unlist(lapply(anc, function(v) which(A[v, ] == 1L)))))
    if (length(more) == length(anc)) break
    anc <- more
  }
  M <- A[anc, anc, drop = FALSE]     # moralize the ancestral subgraph
  U <- (M + t(M)) > 0
  for (i in seq_along(anc)) {        # marry parents
    pa <- which(M[i, ] == 1L)
    if (length(pa) > 1) U[pa, pa] <- TRUE
  }
  diag(U) <- FALSE
  keep <- !(anc %in% S)              # delete S, test connectivity x ~ y
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

# validity of an observed set Z (V-indices): edge-cut d-separation criterion —
# cut the X->Y edges, then X1 dsep Y | Z u {X2} and X2 dsep Y | Z u {X1}
valid_set <- function(A, Zfull) {
  Ac <- A; Ac[Yi, c(X1i, X2i)] <- 0L
  dag_dsep(Ac, X1i, Yi, c(Zfull, X2i)) && dag_dsep(Ac, X2i, Yi, c(Zfull, X1i))
}

# setting classification defaults to the DAG:
# S1 <=> some Z subset of the observed W's is valid in the true DAG
dag_set_exists <- function(A) {
  for (sz in 0:length(Wv)) for (Z in combn_list(Wv, sz))
    if (valid_set(A, as.integer(Z))) return(TRUE)
  FALSE
}

# primary (PAG) screening: does the ground-truth tiered PAG certify an
# adjustment set? Oracle FCI on the true
# DAG (latents marked latent, four tiers as knowledge), then GAC set listing
# on the resulting PAG (via dagitty; see the gac-bug note above).
# Smoke test:  A <- gen_dag(TRUE); gt_pag_has_set(A)
gt_pag_has_set <- function(A) {
  nm_all <- c(paste0("L", 1:5), obs_name)     # V1..V5 latent, then observed
  idx <- which(A == 1L, arr.ind = TRUE)       # A[i, j] = 1  <=>  j -> i
  edges <- lapply(seq_len(nrow(idx)), function(k)
    list(nm_all[idx[k, 2]], nm_all[idx[k, 1]]))
  screen_tiers <- if (SCREEN_KNOWLEDGE == "legacy")
    list(as.list(c(paste0("W", 1:7), "X1")), "X2", "Y")   # original code
  else tiers_py                                           # manuscript text
  txt <- reticulate::py$gt_pag(edges, as.list(paste0("L", 1:5)), screen_tiers)
  txt <- paste(txt, collapse = "\n")
  edge_lines <- grep("^\\s*[0-9]+\\.", strsplit(txt, "\n")[[1]], value = TRUE)
  if (any(grepl("\\bL[1-5]\\b", edge_lines)))
    stop("latent nodes leaked into the ground-truth PAG; check the MsepTest ",
         "import path in the python block")
  amat <- tetrad_string_to_amat(txt)
  chk <- check_tier_marks(amat)
  if (chk$stamps > 0 || chk$conflicts > 0)
    stop("SCREEN GUARD: the Tetrad in use returned an oracle tiered PAG with ",
         chk$stamps, " missing tier-forced arrowhead(s) and ", chk$conflicts,
         " tier-conflicting tail(s), so its handling of tiered background ",
         "knowledge cannot be trusted and the Setting 1/2 screen would be ",
         "corrupted. Upgrade the py-tetrad installation (the bundled Tetrad ",
         "jar must apply tiered knowledge; version 7.6.11 is verified) and ",
         "rerun:\n",
         "  pip install --force-reinstall git+https://github.com/cmu-phil/py-tetrad")
  # cheap necessary condition before the dagitty call: amenability requires a
  # visible edge out of each treatment toward Y, so a circle mark at either
  # treatment on its Y-edge (or a missing X-Y edge with no possibly-directed
  # proper path, the common oracle case: no edge at all) cannot be repaired
  xs <- c(8, 9); ys <- 10
  mX <- amat[ys, xs]
  if (any(mX == 1) || any(mX == 2)) return(FALSE)
  # minimal instead of all: emptiness is equivalent (a set exists iff a
  # minimal one does) and the listing is far cheaper; matches the original
  # pcalgAdjustmentSets(..., 'pag', 'minimal') call
  sets <- tryCatch(
    dagitty::adjustmentSets(pag_to_dagitty(amat),
                            exposure = c("X1", "X2"), outcome = "Y",
                            type = "minimal"),
    error = function(e) NULL)
  !is.null(sets) && length(sets) > 0
}

screen_S1 <- function(A)
  if (SCREENING == "pag") gt_pag_has_set(A) else dag_set_exists(A)

combn_list <- function(v, m) {
  if (m == 0) return(list(integer(0)))
  if (length(v) < m) return(list())
  asplit(combn(v, m), 2)
}

# ---- exact tFCI test counting via trace() ----------------------------------------
# causalDisco's tfci engine takes its test by NAME ("fisher_z"), so unlike
# pcalg::fci there is no test function to wrap. fisher_z runs through the
# internal cor_test (ls(getNamespace("causalDisco"), pattern = "test")), which
# we count by tracing it for the duration of the disco() call. SANITY CHECK on
# first use: the counts should be the same order as the wrapped FCI + GAC
# counter on the same data (same skeleton phase; tiers forbid some edges, so
# tFCI may test somewhat fewer).
tfci_counter <- new.env(parent = emptyenv()); tfci_counter$n <- 0L
tfci_count_hit <- function() tfci_counter$n <- tfci_counter$n + 1L

with_tfci_counter <- function(expr) {
  tfci_counter$n <- 0L
  traced <- tryCatch({
    suppressMessages(trace("cor_test", where = getNamespace("causalDisco"),
                           tracer = quote(tfci_count_hit()), print = FALSE))
    TRUE
  }, error = function(e) FALSE)
  on.exit(if (traced) suppressMessages(
    tryCatch(untrace("cor_test", where = getNamespace("causalDisco")),
             error = function(e) invisible(NULL))))
  val <- expr
  list(value = val, ntests = if (traced && tfci_counter$n > 0)
                               tfci_counter$n else NA)
}

# ---- CI testing with a counter --------------------------------------------------
make_tester <- function(dat) {
  C <- cor(dat); n <- nrow(dat)
  count <- 0L
  pv <- function(i, j, S) {
    count <<- count + 1L
    gaussCItest(i, j, S, list(C = C, n = n))
  }
  dep <- function(i, j, S = integer(0)) pv(i, j, S) <= R1_ALPHA_DEP
  ind <- function(i, j, S = integer(0)) pv(i, j, S) >  R1_ALPHA_IND
  list(dep = dep, ind = ind, count = function() count)
}

# ---- R1 Build (Theorem 5), k = 2 -------------------------------------------------
# indices within the OBSERVED data matrix (columns ordered as obs_name)
Wd <- 1:7; X1d <- 8; X2d <- 9; Yd <- 10

# R1 Build: FULL enumeration (no early exit) over
# ordered witness pairs and Z with |Z| <= nobs-5 = 5; return ALL passing Z's.
r1_build <- function(tester) {
  dep <- tester$dep; ind <- tester$ind
  found <- list()
  for (W1 in Wd) for (W2 in setdiff(Wd, W1)) {
    rest <- setdiff(Wd, c(W1, W2))
    for (sz in 0:min(5, length(rest))) for (Z in combn_list(rest, sz)) {
      Z <- as.integer(Z)
      ok <- dep(W1, Yd, c(Z)) &&              # (i)  i=1: W1 dep Y | Z
            ind(W1, Yd, c(Z, X1d)) &&         # (ii) i=1
            dep(W2, Yd, c(Z, X1d)) &&         # (i)  i=2: given Z u {X1}
            ind(W2, Yd, c(Z, X1d, X2d))       # (ii) i=2
      if (ok) found[[length(found) + 1]] <- Z
    }
  }
  if (length(found) == 0) return(list(found = FALSE))
  list(found = TRUE, Zs = unique(found))
}

# ---- R1 Combine (Theorem 7), k = 2 -----------------------------------------------
# R1 Combine: per treatment collect ALL passing
# (witness, T) with |T| <= 6 (pool for X2 includes X1 as a member, never witness);
# reduce T's to inclusion-minimal; prune each T to the smallest subset still
# passing the witness conditions; return all unions Z1 u Z2.

r1_combine <- function(tester) {
  dep <- tester$dep; ind <- tester$ind
  collect <- function(Xd, extra) {
    hits <- list()
    for (Wi in Wd) {
      pool <- c(setdiff(Wd, Wi), extra)
      for (sz in 0:min(6, length(pool))) for (Ti in combn_list(pool, sz)) {
        Ti <- as.integer(Ti)
        if (dep(Wi, Yd, Ti) && ind(Wi, Yd, c(Ti, Xd)))
          hits[[length(hits) + 1]] <- list(W = Wi, T = Ti)
      }
    }
    hits
  }
  # Lemma 8 minimality pruning (the manuscript's configuration; matches the
  # package's r1_combine(prune = "lemma8")).
  joint_ind <- function(a, B, C) {
    if (length(B) == 0) return(TRUE)
    for (j in seq_along(B))
      if (!ind(a, B[j], c(C, B[seq_len(j - 1)]))) return(FALSE)
    TRUE
  }
  elementwise_dep <- function(S, Xd) {
    all(vapply(S, function(v) {
      rest <- setdiff(S, v)
      dep(Xd, v, rest) && dep(Yd, v, c(rest, Xd))
    }, TRUE))
  }
  minimal_Z <- function(Wi, Ti, Xd) {
    if (elementwise_dep(Ti, Xd)) return(Ti)            # (i)-(ii): T minimal
    if (length(Ti) >= 1) {
      for (sz in 0:(length(Ti) - 1L)) for (S in combn_list(Ti, sz)) {
        S <- as.integer(S)
        if (!elementwise_dep(S, Xd)) next              # (iii)-(iv)
        rest <- setdiff(Ti, S)
        if (joint_ind(Yd, rest, c(S, Xd)) ||           # (v)
            joint_ind(Xd, rest, S)) return(S)
      }
    }
    NULL
  }
  h1 <- collect(X1d, extra = integer(0))
  if (length(h1) == 0) return(list(found = FALSE))
  h2 <- collect(X2d, extra = X1d)
  if (length(h2) == 0) return(list(found = FALSE))
  # dedupe by T (keep witness of first hit), then keep inclusion-minimal T's
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
  m1 <- dedupe_min(h1)
  m2 <- dedupe_min(h2)
  Zs <- list()
  for (a in m1) for (b in m2) {
    Z1 <- minimal_Z(a$W, a$T, X1d); Z2 <- minimal_Z(b$W, b$T, X2d)
    if (is.null(Z1) || is.null(Z2)) next
    Zs[[length(Zs) + 1]] <- sort(setdiff(union(Z1, Z2), c(X1d, X2d)))
  }
  if (length(Zs) == 0) return(list(found = FALSE))
  list(found = TRUE, Zs = unique(Zs))
}

# ---- FCI + GAC --------------------------------------------------------------------
# tiers = TRUE : FCItiers via causalDisco tfci (the paper's comparison method).
# tiers = FALSE: plain pcalg::fci, NO background knowledge — A/B arm to answer
# whether the tails at X (hence GAC amenability) come from the tier knowledge
# or from FCI's own orientation rules (visibility, discriminating paths).
fci_gac <- function(dat, tiers = TRUE) {
  if (!tiers) {
    C <- cor(dat); n <- nrow(dat)
    # NB: causalDisco masks pcalg::fci with its own method constructor,
    # hence the explicit namespace (the unqualified call errors with
    # "'arg' must be NULL or a character vector" from match.arg).
    # A counting wrapper around the CI test gives the exact number of tests
    # FCI performs (skeleton + possible-d-sep phases included).
    cnt <- 0L
    counting_test <- function(x, y, S, suffStat) {
      cnt <<- cnt + 1L
      gaussCItest(x, y, S, suffStat)
    }
    fit <- tryCatch(suppressWarnings(
      pcalg::fci(list(C = C, n = n), counting_test, alpha = COMP_ALPHA,
                 labels = colnames(dat), verbose = FALSE)),
      error = function(e) {
        cat("  [fci error]", conditionMessage(e), "\n"); NULL })
    if (is.null(fit)) return(list(outcome = "unknown", ntests = NA))
    amat <- fit@amat                   # pcalg coding: amat[i, j] = mark at j
    res <- pag_decision(amat)
    res$ntests <- cnt
    return(res)
  }
  # FCItiers via causalDisco CRAN 1.0.1 disco() interface, exactly as in
  # tiers declared via knowledge() + tier(); result edges at fit$caugi@edges.
  #
  # METHODOLOGICAL CAVEAT: pcalg::gac(type = "pag") is designed
  # for proper PAGs; a tier-augmented PAG is not one, so GAC on tfci output
  # carries no established guarantee. This mirrors the PAPER's documented
  # comparison method (FCItiers + GAC), whose theoretical shakiness is part of the
  # paper's own argument: no complete adjustment criteria exist for
  # knowledge-augmented PAGs. Any comparison method misfires are
  # caught honestly downstream, since found sets are validated against the
  # TRUE DAG in is_correct_set().
  d2 <- as.data.frame(dat)                       # columns W1..W7, X1, X2, Y
  kn <- knowledge(d2, tier(
    covariates ~ c(W1, W2, W3, W4, W5, W6, W7),
    treat1     ~ X1,
    treat2     ~ X2,
    outcome    ~ Y))
  # All graph-based arms test at COMP_ALPHA.
  wrapped <- with_tfci_counter(tryCatch(
    disco(d2, tfci(engine = "causalDisco", test = "fisher_z",
                   alpha = COMP_ALPHA),
          knowledge = kn),
    error = function(e) {
      cat("  [tfci error]", conditionMessage(e), "\n"); NULL }))
  fit <- wrapped$value
  if (is.null(fit)) return(list(outcome = "unknown", ntests = wrapped$ntests))
  edges_df <- tryCatch(as.data.frame(fit$caugi@edges), error = function(e) NULL)
  if (is.null(edges_df)) return(list(outcome = "unknown", ntests = wrapped$ntests))
  # caugi edge strings -> pcalg PAG mark codes (0 none, 1 circle, 2 arrow, 3 tail)
  mark_left  <- function(ch) ifelse(ch == "-", 3L,
                             ifelse(ch == "<", 2L,
                             ifelse(ch == "o", 1L, 0L)))
  mark_right <- function(ch) ifelse(ch == "-", 3L,
                             ifelse(ch == ">", 2L,
                             ifelse(ch == "o", 1L, 0L)))
  amat <- matrix(0L, length(obs_name), length(obs_name),
                 dimnames = list(obs_name, obs_name))
  for (i in seq_len(nrow(edges_df))) {
    from <- as.character(edges_df$from[i])
    to   <- as.character(edges_df$to[i])
    e    <- as.character(edges_df$edge[i])
    if (!(from %in% obs_name) || !(to %in% obs_name)) next
    L <- substr(e, 1, 1); R <- substr(e, nchar(e), nchar(e))
    # pcalg convention: amat[i, j] = mark at node j on the edge between i and j
    amat[to,   from] <- mark_left(L)
    amat[from, to]   <- mark_right(R)
  }
  chk <- check_tier_marks(amat)       # tier-mark guard: repair + tally
  tier_guard_tally$stamps <- tier_guard_tally$stamps + chk$stamps
  tier_guard_tally$conflicts <- tier_guard_tally$conflicts + chk$conflicts
  amat <- chk$amat
  res <- pag_decision(amat)
  res$ntests <- wrapped$ntests       # exact via the cor_test trace
  res
}

# ---- GFCI + GAC (optional arm) ----------------------------------------------
# GFCI with tiered background knowledge via Tetrad, reached through py-tetrad
# and reticulate. There is no native R GFCI. Tetrad's Knowledge object accepts
# tiers, but NO completeness theory exists for GFCI-with-tiers, and the
# certification caveat is identical to the tFCI arm's: the output is not
# guaranteed to be a PAG, so dagitty/GAC on it is a best-effort construction.
# The arm tests empirically what the paper argues theoretically.
# Smoke test before a full run:
#   m <- readRDS("jci_models.rds")$S1[[1]]
#   d <- sim_data(m$sem, 1000)[, obs]; colnames(d) <- obs_name
#   gfci_gac(d)
tetrad_gac <- function(dat, algo = c("gfci", "fci")) {
  algo <- match.arg(algo)
  d2 <- as.data.frame(dat)
  out <- tryCatch(
    reticulate::py$tetrad_run(reticulate::r_to_py(d2), COMP_ALPHA, tiers_py,
                              algo),
    error = function(e) {
      cat("  [tetrad/", algo, "]", conditionMessage(e), "\n"); NULL })
  if (is.null(out)) return(list(outcome = "unknown", ntests = NA))
  amat <- tryCatch(
    tetrad_string_to_amat(paste(out$graph, collapse = "\n")),
    error = function(e) NULL)
  if (is.null(amat)) return(list(outcome = "unknown", ntests = out$count))
  chk <- check_tier_marks(amat)       # tier-mark guard: repair + tally
  tier_guard_tally$stamps <- tier_guard_tally$stamps + chk$stamps
  tier_guard_tally$conflicts <- tier_guard_tally$conflicts + chk$conflicts
  amat <- chk$amat
  res <- pag_decision(amat)           # same decision rule + dagitty certifier
  res$ntests <- out$count             # exact, via the counting proxy
  res
}
gfci_gac <- function(dat) tetrad_gac(dat, "gfci")

# Parse Tetrad's printed graph ("Graph Edges: 1. W1 o-> X1 ...") into a pcalg
# PAG amat (pcalg codes: 1 circle, 2 arrow, 3 tail; amat[i, j] = mark at j).
tetrad_string_to_amat <- function(txt) {
  amat <- matrix(0L, length(obs_name), length(obs_name),
                 dimnames = list(obs_name, obs_name))
  edge_re <- "^\\s*\\d+\\.\\s+(\\S+)\\s+(<->|o->|<-o|o-o|-->|<--|---)\\s+(\\S+)"
  mk <- c("<" = 2L, ">" = 2L, "o" = 1L, "-" = 3L)
  for (ln in unlist(strsplit(txt, "\n"))) {
    m <- regmatches(ln, regexec(edge_re, ln))[[1]]
    if (length(m) != 4) next
    from <- m[2]; e <- m[3]; to <- m[4]
    if (!(from %in% obs_name) || !(to %in% obs_name)) next
    amat[to, from] <- mk[[substr(e, 1, 1)]]              # mark at `from`
    amat[from, to] <- mk[[substr(e, nchar(e), nchar(e))]] # mark at `to`
  }
  amat
}

# Decision flow mirrors the original study's FCI arm (decision logic
# inlined there). mX = marks AT X1/X2 on their edges with Y (amat[i, j] = mark
# at j, so mark at X is amat[Y, X]):
#   any mark == 2 (arrowhead into an X)  -> negative, estimate 0
#   else any mark == 3 (tail at an X)     -> found iff GAC sets exist,
#        otherwise UNKNOWN with estimate 0 (their effects(...) = 0 branch)
#   else both marks == 0 (no X-Y edges)   -> negative, estimate 0
#   else (circles)                        -> unknown, NO estimate
# Note: the tail pre-check is implied by GAC amenability, so
# this flow and "GAC sets exist" coincide on the found call; the flow only
# changes the estimate-0 bookkeeping. The learned amat travels with every
# return so anomalies can be saved and dissected.
#
# ADJUSTMENT-SET CERTIFIER: pcalg::gac(type = "pag") has a bug —
# for joint X it passes PAGs whose OTHER treatment has a circle-start edge to
# Y, violating amenability (isolated in gac_bug_repro.R). Sets are therefore
# certified with dagitty::isAdjustmentSet / adjustmentSets, an independent
# implementation of the same criterion (Perković, Textor, Kalisch & Maathuis).
# pcalg's verdict is still computed and disagreements are tallied (gacdis).

pag_to_dagitty <- function(amat) {
  nm <- colnames(amat); if (is.null(nm)) nm <- c(paste0("W", 1:7), "X1", "X2", "Y")
  # dagitty edge token from the two end marks (pcalg codes: 1 circle, 2 arrow,
  # 3 tail); mark at i = amat[j, i], mark at j = amat[i, j]
  tok <- function(mi, mj) {
    key <- paste0(mi, mj)
    switch(key,
      "32" = "->", "23" = "<-", "22" = "<->", "33" = "--",
      "12" = "@->", "21" = "<-@", "11" = "@-@",
      stop("unmapped PAG edge marks: ", key))
  }
  p <- ncol(amat); stmts <- nm                     # declare all nodes
  for (i in 1:(p - 1)) for (j in (i + 1):p)
    if (amat[i, j] != 0 || amat[j, i] != 0)
      stmts <- c(stmts, paste(nm[i], tok(amat[j, i], amat[i, j]), nm[j]))
  dagitty::dagitty(paste0("pag { ", paste(stmts, collapse = " ; "), " }"))
}

pag_decision <- function(amat) {
  xs <- c(8, 9); ys <- 10
  arrowX <- any(amat[ys, xs] == 2 & amat[xs, ys] != 0)
  mX <- amat[ys, xs]
  markX <- paste(mX, collapse = "")
  base <- list(arrowX = arrowX, markX = markX, amat = amat)
  if (any(mX == 2))
    return(c(list(outcome = "negative", est0 = TRUE), base))
  W_nm <- paste0("W", 1:7)
  dagitty_sets <- function() {           # ALL valid sets, via dagitty
    g <- tryCatch(pag_to_dagitty(amat), error = function(e) {
      cat("  [pag_to_dagitty]", conditionMessage(e), "\n"); NULL })
    if (is.null(g)) return(list())
    sets <- tryCatch(
      dagitty::adjustmentSets(g, exposure = c("X1", "X2"), outcome = "Y",
                              type = "all"),
      error = function(e) NULL)
    if (is.null(sets)) {                 # fallback: enumerate subsets
      sets <- list()
      for (sz in 0:7) for (Z in combn_list(1:7, sz))
        if (isTRUE(tryCatch(
              dagitty::isAdjustmentSet(g, W_nm[as.integer(Z)],
                                       exposure = c("X1", "X2"), outcome = "Y"),
              error = function(e) FALSE)))
          sets[[length(sets) + 1]] <- W_nm[as.integer(Z)]
    }
    lapply(sets, function(s) match(as.character(s), W_nm))
  }
  pcalg_gac_any <- function() {          # buggy reference, for the tally only
    for (sz in 0:7) for (Z in combn_list(1:7, sz))
      if (isTRUE(tryCatch(gac(amat, xs, ys, as.integer(Z), type = "pag")$gac,
                          error = function(e) FALSE))) return(TRUE)
    FALSE
  }
  if (any(mX == 3)) {
    found <- dagitty_sets()
    gacdis <- (length(found) > 0) != pcalg_gac_any()
    if (length(found) > 0)
      return(c(list(outcome = "found", Zs = found, gacdis = gacdis), base))
    return(c(list(outcome = "unknown", est0 = TRUE, gacdis = gacdis), base))
  }
  if (all(mX == 0))
    return(c(list(outcome = "negative", est0 = TRUE), base))
  c(list(outcome = "unknown"), base)   # circles: no estimate
}

# ---- outcome classification vs the truth --------------------------------------------
# Scoring: a found set is correct iff valid in the TRUE DAG (edge-cut
# d-separation criterion); the run is "correct set" iff AT LEAST ONE found set
# is valid.
any_correct <- function(A, Zs) {
  any(vapply(Zs, function(Z) valid_set(A, obs[Z]), TRUE))
}

# ---- driver ---------------------------------------------------------------------------
run_all <- function() {
  models <- list(S1 = list(), S2 = list(), S3 = list())
  tried <- 0L
  # S1/S2 from the effect arm (X->Y edges forced present), screened by the
  # true-DAG criterion; S3 = no-effect arm, unscreened.
  while (any(vapply(models, length, 1L) < N_MODELS)) {
    tried <- tried + 1L
    need_eff <- length(models$S1) < N_MODELS || length(models$S2) < N_MODELS
    if (need_eff) {
      A <- gen_dag(force_effect = TRUE)
      key <- if (screen_S1(A)) "S1" else "S2"
    } else {
      A <- gen_dag(force_effect = FALSE)
      key <- "S3"
    }
    if (length(models[[key]]) < N_MODELS) {
      models[[key]] <- c(models[[key]], list(list(A = A, sem = gen_sem(A))))
      cat(sprintf("generation: S1 %d/%d  S2 %d/%d  S3 %d/%d  (candidates tried: %d)\n",
                  length(models$S1), N_MODELS, length(models$S2), N_MODELS,
                  length(models$S3), N_MODELS, tried))
    }
  }
  res <- list()
  for (setting in names(models)) for (mi in seq_along(models[[setting]])) {
    m <- models[[setting]][[mi]]
    for (n in NS) {
      dat <- sim_data(m$sem, n)[, obs, drop = FALSE]
      colnames(dat) <- obs_name
      te <- true_effect(m$sem)
      methods_run <- c("R1 Build", "R1 Combine",
                       if (INCLUDE_TFCI) "tFCI + GAC",
                       "FCI + GAC",
                       if (INCLUDE_TETRAD_FCI) "FCItiers (Tetrad)",
                       if (INCLUDE_GFCI) "GFCI + GAC")
      for (method in methods_run) {
        tester <- make_tester(dat)
        t0 <- Sys.time()
        out <- switch(method,
          "R1 Build"   = { r <- r1_build(tester)
                           list(outcome = if (r$found) "found" else "unknown",
                                Zs = r$Zs, ntests = tester$count()) },
          "R1 Combine" = { r <- r1_combine(tester)
                           list(outcome = if (r$found) "found" else "unknown",
                                Zs = r$Zs, ntests = tester$count()) },
          "tFCI + GAC" = { r <- fci_gac(dat, tiers = TRUE)
                           r },              # ntests via the cor_test trace
          "FCI + GAC"  = { r <- fci_gac(dat, tiers = FALSE)
                           r },              # ntests counted inside
          "FCItiers (Tetrad)" = tetrad_gac(dat, "fci"),
          "GFCI + GAC"        = tetrad_gac(dat, "gfci"))
                                             # both count exactly via the proxy
        secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        # runtime covers the full pipeline per method: search + (for the FCI
        # arms) discovery, decision, and set certification — comparable to
        # the per-method times the original study recorded.
        # (gac_bug_repro.R holds the pcalg::gac reproducible example;
        #  the gacdis column tallies certifier disagreements)
        outcome <- out$outcome
        est <- NA_real_
        nsets <- 0L; nvalid <- 0L        # per-set TP/FP, as the original study kept
        if (outcome == "found") {
          valids <- vapply(out$Zs, function(Z) valid_set(m$A, obs[Z]), TRUE)
          nsets  <- length(valids)
          nvalid <- sum(valids)
          outcome <- if (nvalid > 0) "correct set" else "none correct"
          # Effect estimate: MEAN over all found sets of the joint
          # regression coefficient sum
          df <- as.data.frame(dat)
          ests <- vapply(out$Zs, function(Z) {
            f <- lm(reformulate(c("X1", "X2", obs_name[Z]), "Y"), data = df)
            sum(coef(f)[c("X1", "X2")])
          }, numeric(1))
          est <- mean(ests)
        } else if (outcome == "negative" || isTRUE(out$est0)) est <- 0
        res[[length(res) + 1]] <- data.frame(
          setting = setting, model = mi, n = n, method = method,
          outcome = outcome, est = est, true = te, ntests = out$ntests,
          arrowX = if (!is.null(out$arrowX)) out$arrowX else NA,
          markX  = if (!is.null(out$markX))  out$markX  else NA_character_,
          gacdis = if (!is.null(out$gacdis)) out$gacdis else NA,
          secs   = secs, nsets = nsets, nvalid = nvalid)
      }
    }
    cat(setting, mi, "done\n")
  }
  # persist the generated models (DAG + SEM) so any case can be replicated
  saveRDS(models, out_file("jci_models.rds"))
  cat(sprintf("saved %d models to %s\n", sum(lengths(models)),
              out_file("jci_models.rds")))
  do.call(rbind, res)
}

res <- run_all()
saveRDS(res, out_file("jci_sim_results.rds"))

# ---- outputs: Fig 6/7 analogues + NEW test-count table + NEW MSE panel ---------------
res$outcome <- factor(res$outcome,
  levels = c("correct set", "none correct", "negative", "unknown"))
res$method <- factor(res$method,
  levels = c("R1 Build", "R1 Combine", "FCI + GAC", "tFCI + GAC",
             "FCItiers (Tetrad)", "GFCI + GAC"))
res$method <- droplevels(res$method)     # drops unused levels when arms are off
print(with(res, table(setting, method, outcome, n)))

fig6 <- ggplot(res, aes(factor(n), fill = outcome)) +
  geom_bar(position = "stack") +
  facet_grid(setting ~ method) + theme_minimal() +
  labs(x = "sample size", y = "count")
ggsave(out_file("fig6_rebuild.pdf"), fig6, width = 9, height = 6)

est <- subset(res, !is.na(est))
est$diff <- est$est - est$true
fig7 <- ggplot(est, aes(factor(n), diff)) +
  geom_boxplot() + facet_grid(setting ~ method) + theme_minimal() +
  labs(x = "sample size", y = "estimated - true effect")
ggsave(out_file("fig7_rebuild.pdf"), fig7, width = 9, height = 6)

# TPR panel: among runs that returned sets,
# the share of returned sets that are valid in the true DAG
tpr_rows <- subset(res, nsets > 0)
tpr_rows$tpr <- tpr_rows$nvalid / tpr_rows$nsets
fig_tpr <- ggplot(tpr_rows, aes(factor(n), tpr)) + geom_boxplot() +
  facet_grid(setting ~ method) + theme_minimal() +
  labs(x = "sample size", y = "share of returned sets that are valid")
ggsave(out_file("fig_tpr_rebuild.pdf"), fig_tpr, width = 9, height = 6)

# timing panel: per-dataset runtime
fig_time <- ggplot(res, aes(factor(n), secs)) + geom_boxplot() +
  scale_y_log10() + facet_grid(setting ~ method) + theme_minimal() +
  labs(x = "sample size", y = "seconds per dataset (log scale)")
ggsave(out_file("fig_times_rebuild.pdf"), fig_time, width = 9, height = 6)

# paper-matching layouts: Figures 6/7 show Settings 1-2 only,
# with Setting 3 in Supplement D as Figures 10/11, and Figure 8 = run times
s12 <- subset(res, setting != "S3")
ggsave(out_file("fig8_rebuild.pdf"),
  ggplot(s12, aes(factor(n), secs)) + geom_boxplot() + scale_y_log10() +
    facet_grid(setting ~ method) + theme_minimal() +
    labs(x = "sample size", y = "seconds per dataset (log scale)"),
  width = 9, height = 5)
s3 <- subset(res, setting == "S3")
ggsave(out_file("fig10_rebuild.pdf"),
  ggplot(s3, aes(factor(n), fill = outcome)) + geom_bar(position = "stack") +
    facet_grid(. ~ method) + theme_minimal() +
    labs(x = "sample size", y = "count"),
  width = 9, height = 3.2)
est3 <- subset(s3, !is.na(est)); est3$diff <- est3$est - est3$true
ggsave(out_file("fig11_rebuild.pdf"),
  ggplot(est3, aes(factor(n), diff)) + geom_boxplot() +
    facet_grid(. ~ method) + theme_minimal() +
    labs(x = "sample size", y = "estimated - true effect"),
  width = 9, height = 3.2)

# CI-test counts, full distribution, and MSE
print(aggregate(ntests ~ method + n + setting, res,
                function(x) c(min = min(x), median = median(x),
                              mean = round(mean(x)), max = max(x))))
# NEW: per-method runtimes in seconds (R1-4, the brute-force comment)
print(aggregate(secs ~ method + n + setting, res,
                function(x) c(min = round(min(x), 3), median = round(median(x), 3),
                              mean = round(mean(x), 3), max = round(max(x), 3))))
# TPR among found runs: fraction of returned sets that are valid
found_rows <- subset(res, nsets > 0)
if (nrow(found_rows) > 0)
  print(aggregate(cbind(tpr = nvalid / nsets) ~ method + n + setting,
                  found_rows, function(x) round(mean(x), 3)))
print(aggregate(I(diff^2) ~ method + n + setting, est, mean))
# DIAGNOSTIC: how often does the learned PAG carry an arrowhead into a
# treatment on an X-Y edge (the Tetrad-divergence hypothesis)? Reported for
# BOTH FCI arms so the A/B answers where the marks come from.
fci_rows <- subset(res, method %in% c("tFCI + GAC", "FCI + GAC",
                                      "FCItiers (Tetrad)", "GFCI + GAC") &
                        !is.na(arrowX))
if (nrow(fci_rows) > 0)
  print(aggregate(arrowX ~ method + n + setting, fci_rows, mean))
# AMENABILITY / KNOWLEDGE A/B: mark pattern at (X1, X2) on their Y-edges
# (0 none, 1 circle, 2 arrow, 3 tail). "33" = both tails (certifiable);
# any 1 = circle (amenability fails regardless of the rest of the PAG).
# Tier knowledge places arrowheads AT X from earlier-tier W's, and FCI's rules
# propagate those into tails on X->Y (visibility); if "FCI plain" shows circles
# where "FCI+GAC" shows tails, the tails are knowledge-induced.
if (nrow(fci_rows) > 0)
  for (mm in unique(fci_rows$method)) {
    cat("\nmarkX for", mm, ":\n")
    print(with(subset(fci_rows, method == mm & !is.na(markX)),
               table(setting, markX, n)))
  }
# pcalg::gac BUG IMPACT: fraction of tail-bearing FCI runs where pcalg's gac
# and dagitty disagree on whether any adjustment set exists
dis_rows <- subset(fci_rows, !is.na(gacdis))
if (nrow(dis_rows) > 0) {
  cat("\npcalg-gac vs dagitty disagreement rate (the gac joint-X bug):\n")
  print(aggregate(gacdis ~ method + n + setting, dis_rows, mean))
}
cat("done; results in", out_file("jci_sim_results.rds"),
    "- figures/models suffixed", if (nzchar(SUFFIX)) SUFFIX else "(none)", "\n")
cat(sprintf(
  "config: R1 dep <= %.2f, R1 indep > %.2f, comparison methods at %.2f, screening %s/%s\n",
  R1_ALPHA_DEP, R1_ALPHA_IND, COMP_ALPHA, SCREENING, SCREEN_KNOWLEDGE))
cat("NOTE: all runs share identical graphs and datasets (same RNG stream),\n",
    "so any arm with unchanged settings must reproduce exactly across runs.\n")
tier_guard_report()
