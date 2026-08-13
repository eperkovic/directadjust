# Both procedures should certify the confounder in a model where C confounds
# both treatments and the outcome, and W1, W2 serve as witnesses.
make_toy <- function(n = 4000, seed = 1) {
  set.seed(seed)
  C  <- rnorm(n); W1 <- rnorm(n); W2 <- rnorm(n)
  X1 <- C + W1 + rnorm(n)
  X2 <- C + W2 + rnorm(n)
  Y  <- X1 + X2 + C + rnorm(n)
  data.frame(C, W1, W2, X1, X2, Y)
}

test_that("r1_build certifies a set containing the confounder", {
  d <- make_toy()
  r <- r1_build(d, "X1", "X2", "Y")
  expect_true(r$found)
  expect_true(any(vapply(r$sets, function(s) "C" %in% s, TRUE)))
})

test_that("r1_combine certifies a set containing the confounder", {
  d <- make_toy()
  r <- r1_combine(d, "X1", "X2", "Y")
  expect_true(r$found)
  expect_true(any(vapply(r$sets, function(s) "C" %in% s, TRUE)))
})

test_that("estimates over certified sets are near the truth (joint effect 2)", {
  d <- make_toy(n = 8000)
  r <- r1_combine(d, "X1", "X2", "Y")
  est <- estimate_joint_effect(d, r)
  expect_true(abs(est$estimate - 2) < 0.1)
})

test_that("a shared tester accumulates counts across procedures", {
  d <- make_toy(n = 1000)
  ct <- make_ci_tester(d, alpha = 0.01)
  r1 <- r1_build(d, "X1", "X2", "Y", tester = ct)
  n1 <- ct$count()
  r2 <- r1_combine(d, "X1", "X2", "Y", tester = ct)
  expect_gt(ct$count(), n1)
})

test_that("r1_entner certifies a set containing the confounder", {
  d <- make_toy()
  r <- r1_entner(d, "X1", "Y", covariates = c("C", "W1", "W2"))
  expect_true(r$found)
  expect_true(any(vapply(r$sets, function(s) "C" %in% s, TRUE)))
})

test_that("c_equivalent finds the empty set equivalent under marginal independence", {
  set.seed(2)
  n <- 4000
  Z <- rnorm(n); X <- rnorm(n); Y <- X + Z + rnorm(n)
  d <- data.frame(Z, X, Y)
  r <- c_equivalent(d, "X", "Y", z = "Z", t = character(0))
  expect_true(r$equivalent)
  expect_identical(r$criterion, "i")
})

test_that("a stricter independence threshold certifies no more sets", {
  d <- make_toy(n = 1000)
  r_common <- r1_build(d, "X1", "X2", "Y", alpha = 0.01)
  r_asym   <- r1_build(d, "X1", "X2", "Y", alpha = 0.01, alpha_indep = 0.1)
  expect_true(all(vapply(r_asym$sets, function(s)
    any(vapply(r_common$sets, function(t) identical(s, t), TRUE)), TRUE)))
})

# ---- oracle regression tests against the paper's own examples --------------

# Example 12 (Figure 6a): the example constructed to show why Theorem 7
# needs minimal adjustment sets. The naive union {Z1, Z2} is invalid; the
# correct certified set is {Z2}.
ex12_edges <- rbind(
  c("X1", "Y"),  c("X2", "Y"),  c("Z1", "W1"), c("W2", "Z2"),
  c("W2", "X2"), c("Z2", "Y"),  c("Z2", "X2"), c("U1", "X1"),
  c("U1", "W1"), c("U2", "Z1"), c("U2", "Z2"), c("U3", "Z1"),
  c("U3", "Y"))

test_that("Lemma 8 pruning certifies only the valid set in Example 12", {
  d  <- oracle_frame(c("W1", "W2", "Z1", "Z2", "X1", "X2", "Y"))
  ct <- make_ci_tester(d, test = make_dsep_oracle(ex12_edges))
  r  <- r1_combine(d, "X1", "X2", "Y",
                   covariates = c("W1", "W2", "Z1", "Z2"), tester = ct)
  expect_true(r$found)
  expect_true(all(vapply(r$sets, function(s) identical(s, "Z2"), TRUE)))
})

test_that("the simulation-protocol pruning reproduces the naive union", {
  d  <- oracle_frame(c("W1", "W2", "Z1", "Z2", "X1", "X2", "Y"))
  ct <- make_ci_tester(d, test = make_dsep_oracle(ex12_edges))
  r  <- r1_combine(d, "X1", "X2", "Y",
                   covariates = c("W1", "W2", "Z1", "Z2"),
                   tester = ct, prune = "witness")
  expect_true(any(vapply(r$sets,
                         function(s) identical(s, c("Z1", "Z2")), TRUE)))
})

# Examples 6 and 9 (Figure 5a): Build fails, Combine finds {Z} using X1 in
# the pool for X2.
ex9_edges <- rbind(
  c("X1", "X2"), c("X2", "Y"), c("X1", "Y"),
  c("U", "W"),   c("U", "X1"), c("Z", "X1"), c("Z", "X2"))

test_that("Build fails and Combine finds {Z} on the Figure 5a model", {
  d  <- oracle_frame(c("W", "Z", "X1", "X2", "Y"))
  ct <- make_ci_tester(d, test = make_dsep_oracle(ex9_edges))
  rb <- r1_build(d, "X1", "X2", "Y", covariates = c("W", "Z"), tester = ct)
  expect_false(rb$found)
  ct2 <- make_ci_tester(d, test = make_dsep_oracle(ex9_edges))
  rc  <- r1_combine(d, "X1", "X2", "Y", covariates = c("W", "Z"),
                    tester = ct2)
  expect_true(rc$found)
  expect_true(all(vapply(rc$sets, function(s) identical(s, "Z"), TRUE)))
})
