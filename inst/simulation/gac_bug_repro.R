# gac_bug_repro.R — minimal distillation of the anomaly (case 2 and friends):
# a PAG where X1 -> Y is visible but X2 o-> Y is a circle edge. For the joint
# X = {X1, X2} this PAG is NOT adjustment amenable (the one-edge path
# X2 o-> Y is proper, possibly directed, and does not start with a visible
# edge), so gac must return FALSE for every Z. Run and see.
#
#   source("gac_bug_repro.R")
#
# amat coding (pcalg amat.pag): amat[i, j] = mark AT j on edge i-j
# (0 none, 1 circle, 2 arrowhead, 3 tail)

suppressMessages(library(pcalg))

nodes <- c("W", "X1", "X2", "Y")
m <- matrix(0L, 4, 4, dimnames = list(nodes, nodes))

# W -> X1   (tail at W, arrowhead at X1) — makes X1 -> Y visible,
#            since W is not adjacent to Y
m["W", "X1"] <- 2L; m["X1", "W"] <- 3L
# X1 -> Y   (tail at X1, arrowhead at Y)
m["X1", "Y"] <- 2L; m["Y", "X1"] <- 3L
# X2 o-> Y  (circle at X2, arrowhead at Y)  <- the amenability killer
m["X2", "Y"] <- 2L; m["Y", "X2"] <- 1L

cat("PAG:  W -> X1 -> Y,  X2 o-> Y\n")
print(m)

run <- function(x, z, label, expect) {
  g <- pcalg:::gac(m, x = x, y = 4, z = z, type = "pag")
  cat(sprintf("%-52s gac = %-5s (should be %s)\n", label, g$gac, expect))
  invisible(g)
}
pcalg:::isAmenable(m, x = c(2, 3), y = 4, type = "pag")

library(dagitty)
?isAdjustmentSet

cat("\n")
g1 <- run(c(2, 3), integer(0), "joint X = {X1, X2}, Z = {}", "FALSE")
g2 <- run(c(2, 3), 1L,         "joint X = {X1, X2}, Z = {W}", "FALSE")
g3 <- run(3,       integer(0), "X2 alone,           Z = {}", "FALSE")
g4 <- run(2,       integer(0), "X1 alone,           Z = {}", "TRUE")

cat("\nfull gac() output for the joint call (which condition passes?):\n")
print(g1)

# If g1/g2 come back TRUE while g3 is FALSE, the amenability check is being
# applied per-treatment only where a directed edge leaves the treatment, and
# circle-start paths out of the OTHER treatment are skipped in the joint case.
# The real learned PAG from the simulation is one line away for comparison:
#   A <- readRDS("jci_anomaly_cases.rds")[[2]]$amat
#   gac(A, x = c(8, 9), y = 10, z = integer(0), type = "pag")


library(dagitty)
g <- dagitty("pag { W ; X1 ; X2 ; Y ; W -> X1 ; X1 -> Y ; X2 @-> Y }")
isAdjustmentSet(g, c(), exposure = c("X1","X2"), outcome = "Y")        # expect FALSE
adjustmentSets(g, exposure = c("X1","X2"), outcome = "Y", type = "all") # expect empty
