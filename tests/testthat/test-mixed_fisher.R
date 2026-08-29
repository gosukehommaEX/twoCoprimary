# Tests for the Fisher option of the mixed continuous and binary combination.
#
# This option is simulation based, so it lives in its own file: the tolerances
# are Monte Carlo tolerances and the runtime is dominated by the number of
# replications, unlike the four asymptotic methods.
#
# Two defects were repaired in version 1.1.0. The asymptotic branch was not
# inside an else, so it overwrote the simulated result and a request for
# "Fisher" returned the ASc value; and the latent binary variable was centred at
# qnorm(1 - p_j) while both groups were dichotomised at qnorm(1 - p2), which
# pinned the control response probability at 0.5 and reversed the effect.

test_that("Fisher is a distinct method from the asymptotic ones", {
  args <- list(n1 = 60, n2 = 60, delta = 0.5, sd = 1,
               p1 = 0.60, p2 = 0.40, rho = 0.5, alpha = 0.025)

  asc <- do.call(power2MixedContinuousBinary, c(args, list(Test = "ASc")))
  set.seed(101)
  fisher <- do.call(power2MixedContinuousBinary,
                    c(args, list(Test = "Fisher", nMC = 5000)))

  expect_false(isTRUE(all.equal(fisher$powerCoprimary, asc$powerCoprimary)))
  expect_equal(fisher$nMC, 5000)
  expect_true(is.na(asc$nMC))
})

test_that("the Fisher result depends on the number of replications", {
  args <- list(n1 = 60, n2 = 60, delta = 0.5, sd = 1,
               p1 = 0.60, p2 = 0.40, rho = 0.5, alpha = 0.025, Test = "Fisher")

  set.seed(202)
  small <- do.call(power2MixedContinuousBinary, c(args, list(nMC = 1000)))
  set.seed(202)
  large <- do.call(power2MixedContinuousBinary, c(args, list(nMC = 20000)))

  expect_false(isTRUE(all.equal(small$powerCoprimary, large$powerCoprimary)))
  expect_true(abs(small$powerCoprimary - large$powerCoprimary) < 0.05)
})

test_that("the Fisher result is reproducible given a seed", {
  args <- list(n1 = 50, n2 = 50, delta = 0.5, sd = 1,
               p1 = 0.60, p2 = 0.40, rho = 0.3, alpha = 0.025,
               Test = "Fisher", nMC = 2000)

  set.seed(303)
  a <- do.call(power2MixedContinuousBinary, args)$powerCoprimary
  set.seed(303)
  b <- do.call(power2MixedContinuousBinary, args)$powerCoprimary

  expect_identical(a, b)
})

test_that("the binary power recovers the marginal Fisher power", {
  # With a negligible effect on the continuous endpoint the co-primary power is
  # driven by the binary endpoint alone, and the marginal binary power must not
  # exceed it
  set.seed(404)
  res <- power2MixedContinuousBinary(n1 = 60, n2 = 60, delta = 0.5, sd = 1,
                                     p1 = 0.70, p2 = 0.40, rho = 0.2,
                                     alpha = 0.025, Test = "Fisher", nMC = 5000)

  expect_true(res$powerCoprimary <= res$powerCont + 1e-12)
  expect_true(res$powerCoprimary <= res$powerBin + 1e-12)
  expect_true(res$powerBin > 0.5)
})

test_that("the treatment effect points in the right direction", {
  # The binary power must rise as p1 moves away from p2
  set.seed(505)
  low <- power2MixedContinuousBinary(60, 60, 0.5, 1, p1 = 0.45, p2 = 0.40,
                                     rho = 0.2, alpha = 0.025,
                                     Test = "Fisher", nMC = 5000)$powerBin
  set.seed(505)
  high <- power2MixedContinuousBinary(60, 60, 0.5, 1, p1 = 0.75, p2 = 0.40,
                                      rho = 0.2, alpha = 0.025,
                                      Test = "Fisher", nMC = 5000)$powerBin

  expect_true(high > low)
})

test_that("Fisher reproduces Table S5 of Sozu et al. (2012)", {
  # Supporting Information Section C, Table 5: sample size per group and the
  # empirical overall power achieved there, alpha = 0.025, target 0.8, sd = 1.
  # The two cheapest cells are used; dev/check_mixed_fisher_fix.R covers all 24.
  cells <- data.frame(
    delta = c(0.521, 0.521),
    p1 = c(0.99, 0.99),
    p2 = c(0.85, 0.85),
    rho = c(0.5, 0.8),
    n = c(75, 74),
    power_published = c(0.802, 0.801)
  )

  nMC <- 20000
  se <- sqrt(0.8 * 0.2 / nMC)

  for (k in seq_len(nrow(cells))) {
    set.seed(600 + k)
    got <- power2MixedContinuousBinary(
      n1 = cells$n[k], n2 = cells$n[k],
      delta = cells$delta[k], sd = 1,
      p1 = cells$p1[k], p2 = cells$p2[k], rho = cells$rho[k],
      alpha = 0.025, Test = "Fisher", nMC = nMC
    )$powerCoprimary

    expect_lt(abs(got - cells$power_published[k]), 4 * se)
  }
})
