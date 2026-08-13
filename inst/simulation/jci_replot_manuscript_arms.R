# jci_replot_manuscript_arms.R — regenerate the manuscript figures
# from the SAVED results, restricted to the manuscript's five arms:
#   R1 Build, R1 Combine, FCI + GAC, FCItiers + GAC, GFCI + GAC
# (the causalDisco tFCI arm is dropped from display; its rows remain in the
# rds files). Also renames "FCItiers (Tetrad)" -> "FCItiers + GAC".
# No simulation is rerun - this only replots. Takes seconds.
#   source("jci_replot_manuscript_arms.R")

suppressMessages(library(ggplot2))

ARM_LEVELS <- c("R1 Build", "R1 Combine", "FCI + GAC",
                "FCItiers + GAC", "GFCI + GAC")

# palette of the original figures + bottom-up stacking, used for every figure that
# replaces a published panel, so the visual diff against the paper is minimal
ORIG_FILL <- scale_fill_manual(values = c(
  "correct set"  = "#0072BD",   # blue
  "none correct" = "#7E7E7E",   # gray
  "negative"     = "#D95319",   # orange-red
  "unknown"      = "#EDB120"))  # yellow
STACK <- geom_bar(position = position_stack(reverse = TRUE))

prep <- function(res) {
  res$method <- as.character(res$method)
  res <- subset(res, method != "tFCI + GAC")
  res$method[res$method == "FCItiers (Tetrad)"] <- "FCItiers + GAC"
  res$method <- factor(res$method, levels = ARM_LEVELS)
  res$outcome <- factor(res$outcome,
    levels = c("correct set", "none correct", "negative", "unknown"))
  res
}

replot_main <- function(rds, suffix = "") {
  res <- prep(readRDS(rds))
  fn <- function(name) sub("\\.pdf$", paste0(suffix, ".pdf"), name)

  ggsave(fn("fig6_rebuild.pdf"),
    ggplot(res, aes(factor(n), fill = outcome)) +
      STACK + ORIG_FILL +
      facet_grid(setting ~ method) + theme_minimal() +
      labs(x = "sample size", y = "count", fill = NULL),
    width = 9, height = 6)

  est <- subset(res, !is.na(est)); est$diff <- est$est - est$true
  ggsave(fn("fig7_rebuild.pdf"),
    ggplot(est, aes(factor(n), diff)) +
      geom_boxplot() + facet_grid(setting ~ method) + theme_minimal() +
      labs(x = "sample size", y = "estimated - true effect"),
    width = 9, height = 6)

  tpr <- subset(res, nsets > 0); tpr$tpr <- tpr$nvalid / tpr$nsets
  ggsave(fn("fig_tpr_rebuild.pdf"),
    ggplot(tpr, aes(factor(n), tpr)) + geom_boxplot() +
      facet_grid(setting ~ method) + theme_minimal() +
      labs(x = "sample size", y = "share of returned sets that are valid"),
    width = 9, height = 6)

  ggsave(fn("fig_times_rebuild.pdf"),
    ggplot(res, aes(factor(n), secs)) + geom_boxplot() + scale_y_log10() +
      facet_grid(setting ~ method) + theme_minimal() +
      labs(x = "sample size", y = "seconds per dataset (log scale)"),
    width = 9, height = 6)

  s12 <- subset(res, setting != "S3")
  ggsave(fn("fig8_rebuild.pdf"),
    ggplot(s12, aes(factor(n), secs)) + geom_boxplot() + scale_y_log10() +
      facet_grid(setting ~ method) + theme_minimal() +
      labs(x = "sample size", y = "seconds per dataset (log scale)"),
    width = 9, height = 5)

  s3 <- subset(res, setting == "S3")
  ggsave(fn("fig10_rebuild.pdf"),
    ggplot(s3, aes(factor(n), fill = outcome)) +
      STACK + ORIG_FILL + facet_grid(. ~ method) +
      theme_minimal() + labs(x = "sample size", y = "count", fill = NULL),
    width = 9, height = 3.2)
  est3 <- subset(s3, !is.na(est)); est3$diff <- est3$est - est3$true
  ggsave(fn("fig11_rebuild.pdf"),
    ggplot(est3, aes(factor(n), diff)) + geom_boxplot() +
      facet_grid(. ~ method) + theme_minimal() +
      labs(x = "sample size", y = "estimated - true effect"),
    width = 9, height = 3.2)
  invisible(NULL)
}

# Side-by-side companion: the reproduction of the original
# Figure 6 layout (Settings 1-2, the three original methods, the original
# default colors), for direct comparison against the original figure.
replot_fig6_original_layout <- function(rds, out = "fig6_original_layout.pdf") {
  res <- prep(readRDS(rds))
  # the published figure's third panel is labeled "FCI+GAC" but the method
  # behind it is FCItiers (Tetrad's FCI with the tier knowledge), so the
  # side-by-side uses the FCItiers + GAC arm
  res <- subset(res, setting != "S3" &
                     method %in% c("R1 Build", "R1 Combine", "FCItiers + GAC"))
  res$method <- droplevels(res$method)
  ggsave(out,
    ggplot(res, aes(factor(n), fill = outcome)) +
      STACK +
      scale_fill_manual(values = c(
          "correct set"  = "#0072BD", "none correct" = "#7E7E7E",
          "negative"     = "#D95319", "unknown"      = "#EDB120"),
        labels = c("correct set"  = "Correct Set",
                   "none correct" = "None Correct",
                   "negative"     = "Negative",
                   "unknown"      = "Unknown")) +
      facet_grid(setting ~ method,
                 labeller = labeller(setting = c(S1 = "Setting 1",
                                                 S2 = "Setting 2"))) +
      theme_minimal() + theme(legend.position = "bottom") +
      labs(x = "sample size", y = "count", fill = NULL),
    width = 9, height = 4.5)
  invisible(NULL)
}

replot_main("jci_sim_results.rds")
replot_fig6_original_layout("jci_sim_results.rds")
