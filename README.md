# directadjust
[![DOI](https://zenodo.org/badge/1333114488.svg)](https://doi.org/10.5281/zenodo.21919198)

Data-driven covariate adjustment, certified directly from conditional
independence tests, without learning a causal graph as an intermediate
object. The package implements all three routes of

> LaPlante, S., Triantafillou, S., and Perković, E.
> *Data-Driven Adjustment for Multiple Treatments.*
> Journal of Causal Inference.

- `r1_entner()` — single-treatment certification (Theorem 1, after Entner,
  Hoyer and Spirtes, 2013);
- `c_equivalent()` — extends certified sets through confounding equivalence
  (Theorem 3);
- `r1_build()` and `r1_combine()` — certification for two treatments
  (Theorems 5 and 7; the theorems cover any number of treatments, the
  implementation specializes to two, matching the paper's simulations).

`estimate_joint_effect()` turns certified sets into effect estimates, and
`make_ci_tester()` provides the shared testing engine.

## Installation

```r
# install.packages("remotes")
remotes::install_github("eperkovic/directadjust", build_vignettes = TRUE)
```

Requires `pcalg`. The simulation scripts additionally use `dagitty`,
`ggplot2`, and (for the Tetrad comparison method arms) `reticulate` with
`py-tetrad`.

## Usage

```r
library(directadjust)

# C confounds both treatments and the outcome; W1, W2 act as witnesses
set.seed(1)
n  <- 4000
C  <- rnorm(n); W1 <- rnorm(n); W2 <- rnorm(n)
X1 <- C + W1 + rnorm(n)
X2 <- C + W2 + rnorm(n)
Y  <- X1 + X2 + C + rnorm(n)
d  <- data.frame(C, W1, W2, X1, X2, Y)

r1_entner(d, "X1", "Y", covariates = c("C", "W1", "W2"))  # one treatment
b <- r1_build(d, "X1", "X2", "Y")      # two treatments, one-pass
k <- r1_combine(d, "X1", "X2", "Y")    # two treatments, per-treatment
estimate_joint_effect(d, k)            # averages over certified sets
```

Both multi-treatment procedures assume the covariates are pretreatment (no
descendants of the treatments) and that `x1` precedes `x2` in the causal
order. These are assumptions of the theorems, not testable by the
procedures; see the paper.

`r1_combine()` reduces each per-treatment set to a minimal adjustment set
through the paper's Lemma 8 by default, as Theorem 7 requires;
`prune = "witness"` uses a faster heuristic instead (see `?r1_combine`
for the distinction).

The default testing level is `alpha = 0.05`, the level of the paper's
simulations. Dependence and independence can also be declared at different
thresholds (dependence when p is at most `alpha`, independence when p
exceeds `alpha_indep`), leaving an undecided region in between so sets are
certified only on decisive evidence:

```r
r1_build(d, "X1", "X2", "Y", alpha = 0.01, alpha_indep = 0.1)
```

The procedures consume only independence decisions, so any conditional
independence test can be plugged in; Gaussian assumptions enter only
through the default Fisher-z test and the linear estimator:

```r
my_test <- function(i, j, S, data) { ... }   # returns a p-value
ct <- make_ci_tester(d, alpha = 0.01, test = my_test)
r1_build(d, "X1", "X2", "Y", tester = ct)
ct$count()                                   # tests accumulate across calls
```

## Vignette

- **Section 5, reproduced** — regenerates the paper's simulation figures and
  tables from the shipped results:
  `vignette("simulation-results", package = "directadjust")`

## Repository layout

- `R/`, `man/`, `tests/` — the package itself.
- `inst/simulation/` — the complete Section 5 pipeline: graph generation,
  ground-truth screening, all method arms including the Tetrad comparison
  methods, decision outcomes, effect estimates, exact CI-test counts, and
  runtimes. The configuration header of `jci_simulation_rebuild.R`
  documents every design choice, including the provenance of the original
  study's settings.
- `inst/results/` — the saved simulation results (`.rds`). Figures are not
  committed; they are regenerated from these files on demand, either by
  knitting the vignette or with `inst/simulation/jci_replot_manuscript_arms.R`.
  `jci_sim_results.rds` and `jci_models.rds` are the main study,
  `jci_fig8_times.rds` the covariate sweep of Figure 8.
- `vignettes/` — regenerates the paper's figures and tables from
  `inst/results/`: `vignette("simulation-results", package = "directadjust")`.

## Reproducing the simulation study

Everything needed to rerun is in this repository; no other folder from the
project is required.

**1. R side.** R (>= 4.0; the study ran on R 4.x). Install the package and
the simulation dependencies:

```r
install.packages(c("pcalg", "dagitty", "ggplot2", "causalDisco",
                   "reticulate", "knitr", "rmarkdown", "testthat"))
remotes::install_github("eperkovic/directadjust", build_vignettes = TRUE)
```

`pcalg` needs Bioconductor pieces on a fresh machine:
`BiocManager::install(c("graph", "RBGL"))` first if it fails to install.
Note: adjustment-set validity is certified with `dagitty`, not
`pcalg::gac`, deliberately — `gac(type = "pag")` has a bug for joint
treatments (see `inst/simulation/gac_bug_repro.R`).

**2. Tetrad side** (only for the comparison arms `FCItiers` and
`GFCI + GAC`; the R1 arms run without it). A Java JDK >= 17 on the PATH
with `JAVA_HOME` set, then a one-time Python setup managed by reticulate —
system Pythons will not do:

```r
reticulate::install_python("3.11:latest")
reticulate::virtualenv_create("r-tetrad", version = "3.11:latest")
reticulate::virtualenv_install("r-tetrad", "py-tetrad")
```

`py-tetrad` bundles the Tetrad jar it runs, so its pip version pins the
Tetrad version.

**3. Record the exact versions.** After a successful run, snapshot the
environment next to the results:

```r
source(system.file("simulation", "session_info.R", package = "directadjust"))
```

This writes `session_info.txt` (R, package, Python, py-tetrad, and Java
versions) so future reruns can reconstruct the environment exactly.

**4. Run.** `inst/simulation/jci_simulation_rebuild.R` top to bottom. Set
`TEST_MODE <- TRUE` first for a fast small-scale pass; the full run takes
hours. The script ships in the manuscript's configuration (every method at
alpha = 0.05, Lemma 8 pruning) and writes `jci_models.rds` and
`jci_sim_results.rds`, the files under `inst/results/`. Figures come from
`jci_replot_manuscript_arms.R` or the vignette.

## License

MIT. See `LICENSE`.
