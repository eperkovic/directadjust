# session_info.R — record the software environment of a simulation run into
# session_info.txt. Run from the repo root (or anywhere) in a FRESH R session:
#   source(system.file("simulation", "session_info.R", package = "directadjust"))
#
# NOTE: py-tetrad's declared package version does not identify the bundled
# Tetrad engine (it has read "0.1.2" across many different jars). The engine
# is therefore identified below by the jar's reported version, byte size,
# and SHA-256.

out <- c(sprintf("directadjust simulation environment, captured %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
         strrep("-", 70), "", "## R", "", R.version.string,
         paste0("Platform: ", R.version$platform), "", "## R packages", "")

pkgs <- c("directadjust", "pcalg", "dagitty", "ggplot2", "causalDisco",
          "reticulate", "graph", "RBGL", "knitr", "rmarkdown", "testthat")
for (p in pkgs) {
  v <- tryCatch(as.character(utils::packageVersion(p)),
                error = function(e) "NOT INSTALLED")
  out <- c(out, sprintf("%-14s %s", p, v))
}

## ---- Python / py-tetrad / Tetrad jar --------------------------------------
reticulate::use_virtualenv("r-tetrad", required = TRUE)
reticulate::py_run_string("
import sys, os
import importlib.metadata as _md
import pytetrad
_jar = os.path.join(os.path.dirname(pytetrad.__file__),
                    'resources', 'tetrad-current.jar')
_info = {'pyver': sys.version.split()[0],
         'declared': _md.version('py-tetrad'),
         'jar': _jar,
         'jar_bytes': str(os.path.getsize(_jar))}
try:
    import pytetrad.tools.translate      # boots the JVM
    from edu.cmu.tetrad.util import Version
    _info['tetrad'] = str(Version.currentViewableVersion())
except Exception as _e:
    _info['tetrad'] = 'unknown (' + str(_e) + ')'
")
info <- reticulate::py$`_info`

sha <- tryCatch({
  cmd <- if (nzchar(Sys.which("shasum"))) c("shasum", "-a", "256")
         else c("sha256sum")
  strsplit(system2(cmd[1], c(cmd[-1], info$jar), stdout = TRUE), " ")[[1]][1]
}, error = function(e) "unavailable")

out <- c(out, "", "## Python / py-tetrad / Tetrad", "",
  "python:            r-tetrad virtualenv",
  paste0("python version:    ", info$pyver),
  paste0("py-tetrad package: ", info$declared,
         "  (NOTE: this version string never changes upstream;"),
  "                   identify the engine by the jar below)",
  paste0("Tetrad version:    ", info$tetrad),
  "Tetrad jar:        pytetrad/resources/tetrad-current.jar (in the py-tetrad install)",
  paste0("jar size (bytes):  ", info$jar_bytes),
  paste0("jar SHA-256:       ", sha))

## ---- Java -----------------------------------------------------------------
jh <- Sys.getenv("JAVA_HOME")
jbin <- if (nzchar(jh)) file.path(jh, "bin", "java") else "java"
jver <- tryCatch(system2(jbin, "-version", stdout = TRUE, stderr = TRUE),
                 error = function(e) "unavailable")
out <- c(out, "", "## Java", "", jver)

## ---- run configuration ----------------------------------------------------
out <- c(out, "", "## Simulation configuration (inst/simulation/jci_simulation_rebuild.R)", "",
  "seed 20260716; N_MODELS 100/setting; n in {100, 1000, 10000}",
  "R1 dep <= 0.05, R1 indep > 0.05, comparison alpha 0.05",
  "screening: pag/paper (oracle tiered PAG via Tetrad MsepTest FCI + dagitty)",
  "tier-mark guard: verifies every knowledge-forced endpoint mark on every",
  "  returned graph; a correct Tetrad build reports 0 repairs, 0 conflicts")

writeLines(out, "session_info.txt")
cat("written to session_info.txt\n")
