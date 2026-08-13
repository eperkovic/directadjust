# session_info.R — snapshot the exact software environment of a simulation
# run, so a future rerun can reconstruct it. Writes session_info.txt to the
# working directory (copy it into inst/results/ next to the .rds it belongs
# to, or commit it at the repo root).
#
# Covers: R version and platform, loaded/installed R package versions,
# the reticulate Python (path and version), py-tetrad's pip version, and
# the Java that py-tetrad will use.

out <- "session_info.txt"
con <- file(out, "w")
w <- function(...) writeLines(paste0(...), con)

w("directadjust simulation environment, captured ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"))
w(strrep("-", 70))

# --- R ---
w("\n## R\n")
w(R.version.string)
w("Platform: ", R.version$platform)

pkgs <- c("directadjust", "pcalg", "dagitty", "ggplot2", "causalDisco",
          "reticulate", "graph", "RBGL", "knitr", "rmarkdown", "testthat")
w("\n## R packages\n")
for (p in pkgs) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "NOT INSTALLED")
  w(sprintf("%-14s %s", p, v))
}

# --- Python / py-tetrad (only if the r-tetrad virtualenv exists) ---
w("\n## Python / py-tetrad\n")
ok <- tryCatch({
  reticulate::use_virtualenv("r-tetrad", required = TRUE)
  cfg <- reticulate::py_config()
  w("python:    ", cfg$python)
  w("version:   ", as.character(cfg$version))
  pt <- tryCatch(reticulate::py_run_string(
    "import importlib.metadata as md; v = md.version('py-tetrad')")$v,
    error = function(e) "NOT FOUND")
  w("py-tetrad: ", pt)
  TRUE
}, error = function(e) { w("r-tetrad virtualenv not available: ", conditionMessage(e)); FALSE })

# --- Java ---
w("\n## Java\n")
jv <- tryCatch(system2("java", "-version", stdout = TRUE, stderr = TRUE),
               error = function(e) "java not on PATH")
w(paste(jv, collapse = "\n"))
w("JAVA_HOME: ", Sys.getenv("JAVA_HOME", unset = "(unset)"))

close(con)
message("Wrote ", normalizePath(out))
