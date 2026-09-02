# Regression tests for the plot method and for the argument validation added
# in 1.1.1. Each test corresponds to a call that failed in 1.1.0.
#
# Every plot call fixes n_points, because each point of a curve runs a complete
# power evaluation or a complete sample size search and the default of fifty
# would spend most of the check time here for no extra coverage. The exact
# binary design is chosen to need about a dozen subjects per group for the same
# reason.

ss_cont <- ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
                         rho = 0.5, r = 1, alpha = 0.025, beta = 0.2,
                         known_var = TRUE)
pw_cont <- power2Continuous(n1 = 100, n2 = 100, delta1 = 0.5, delta2 = 0.5,
                            sd1 = 1, sd2 = 1, rho = 0.5, alpha = 0.025,
                            known_var = TRUE)
ss_bin_ap <- ss2BinaryApprox(p11 = 0.40, p12 = 0.35, p21 = 0.25, p22 = 0.20,
                             rho1 = 0.5, rho2 = 0.5, r = 1,
                             alpha = 0.025, beta = 0.2, Test = "AN")

test_that("power_curve works for sample size objects of every endpoint type", {
  pdf(NULL); on.exit(dev.off())
  expect_s3_class(plot(ss_cont, type = "power_curve", n_points = 5),
                  "data.frame")
  expect_s3_class(plot(ss_bin_ap, type = "power_curve", n_points = 5),
                  "data.frame")
  ss_mcb <- ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.60, p2 = 0.40,
                                     rho = 0.5, r = 1, alpha = 0.025,
                                     beta = 0.2, Test = "AN")
  expect_s3_class(plot(ss_mcb, type = "power_curve", n_points = 5),
                  "data.frame")
  ss_mcc <- ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                                    mu1 = -50, mu2 = 0, sd = 250,
                                    rho1 = 0.5, rho2 = 0.5, r = 1,
                                    alpha = 0.025, beta = 0.2)
  expect_s3_class(plot(ss_mcc, type = "power_curve", n_points = 5),
                  "data.frame")
})

test_that("sample_size_rho works for power objects, which carry no r or beta", {
  pdf(NULL); on.exit(dev.off())
  out <- plot(pw_cont, type = "sample_size_rho", n_points = 5)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("rho", "n2") %in% names(out)))
  expect_true(all(is.finite(out$n2)))
})

test_that("plot dispatches exact binary objects to the exact functions", {
  pdf(NULL); on.exit(dev.off())
  ss_ex <- ss2BinaryExact(p11 = 0.80, p12 = 0.70, p21 = 0.30, p22 = 0.20,
                          rho1 = 0.3, rho2 = 0.3, r = 1,
                          alpha = 0.025, beta = 0.2, Test = "Fisher")
  out_ss <- plot(ss_ex, type = "power_curve", n_points = 4)
  expect_s3_class(out_ss, "data.frame")
  pw_ex <- power2BinaryExact(n1 = 20, n2 = 20, p11 = 0.80, p12 = 0.70,
                             p21 = 0.30, p22 = 0.20, rho1 = 0.3, rho2 = 0.3,
                             alpha = 0.025, Test = "Fisher")
  out_pw <- plot(pw_ex, type = "power_curve", n_points = 4)
  expect_s3_class(out_pw, "data.frame")
  # The exact function, not the approximate one, must have produced the curve
  direct <- vapply(out_pw$n2, function(n) power2BinaryExact(
    n1 = n, n2 = n, p11 = 0.80, p12 = 0.70, p21 = 0.30, p22 = 0.20,
    rho1 = 0.3, rho2 = 0.3, alpha = 0.025, Test = "Fisher")$powerCoprimary,
    numeric(1))
  expect_equal(out_pw$power, direct)
})

test_that("plot rejects single endpoint objects with an explanatory message", {
  pdf(NULL); on.exit(dev.off())
  expect_error(
    plot(ss1Continuous(delta = 0.5, sd = 1, r = 1, alpha = 0.025, beta = 0.2)),
    "single endpoint"
  )
})

test_that("power2BinaryApprox validates Test", {
  expect_error(
    power2BinaryApprox(n1 = 100, n2 = 100, p11 = 0.60, p12 = 0.50,
                       p21 = 0.40, p22 = 0.30, rho1 = 0.3, rho2 = 0.3,
                       alpha = 0.025, Test = "Fisher"),
    "Test must be one of"
  )
})

test_that("ss2MixedCountContinuous rejects a reversed effect direction", {
  expect_error(
    ss2MixedCountContinuous(r1 = 1.25, r2 = 1.0, nu = 0.8, t = 1,
                            mu1 = -50, mu2 = 0, sd = 250,
                            rho1 = 0.5, rho2 = 0.5, r = 1,
                            alpha = 0.025, beta = 0.2),
    "r1 must be less than r2"
  )
  expect_error(
    ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                            mu1 = 0, mu2 = -50, sd = 250,
                            rho1 = 0.5, rho2 = 0.5, r = 1,
                            alpha = 0.025, beta = 0.2),
    "mu1 must be less than mu2"
  )
})

test_that("a large standardized effect with unknown variance is handled", {
  res <- ss2Continuous(delta1 = 4, delta2 = 4, sd1 = 1, sd2 = 1, rho = 0.5,
                       r = 1, alpha = 0.025, beta = 0.2, known_var = FALSE)
  expect_true(is.finite(res$n2))
  expect_true(res$n2 >= 1)
  zero <- power2Continuous(n1 = 1, n2 = 1, delta1 = 4, delta2 = 4,
                           sd1 = 1, sd2 = 1, rho = 0.5, alpha = 0.025,
                           known_var = FALSE)
  expect_equal(zero$powerCoprimary, 0)
})

test_that("the arcsine sample size attains the target power when r is not 1", {
  arcsine_power <- function(p1, p2, n1, n2, alpha) {
    se <- 0.5 * sqrt(1 / n1 + 1 / n2)
    pnorm((asin(sqrt(p1)) - asin(sqrt(p2))) / se - qnorm(1 - alpha))
  }
  for (rv in c(1, 2, 0.5)) {
    res <- ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = rv, alpha = 0.025,
                           beta = 0.1, Test = "AS")
    got <- arcsine_power(0.6, 0.4, res$n1, res$n2, 0.025)
    expect_gte(got, 0.90)
    expect_lt(got, 0.93)
  }
})

test_that("the arcsine continuity correction does not reduce the sample size", {
  for (rv in c(1, 2)) {
    as_n <- ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = rv, alpha = 0.025,
                            beta = 0.1, Test = "AS")$n2
    asc_n <- ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = rv, alpha = 0.025,
                             beta = 0.1, Test = "ASc")$n2
    expect_gte(asc_n, as_n)
  }
})

test_that("the two endpoint search still returns the minimum sample size", {
  for (tst in c("AN", "ANc", "AS", "ASc")) {
    for (rv in c(1, 2)) {
      res <- ss2BinaryApprox(p11 = 0.60, p12 = 0.50, p21 = 0.40, p22 = 0.30,
                             rho1 = 0.3, rho2 = 0.3, r = rv,
                             alpha = 0.025, beta = 0.2, Test = tst)
      at_n <- power2BinaryApprox(n1 = res$n1, n2 = res$n2,
                                 p11 = 0.60, p12 = 0.50, p21 = 0.40,
                                 p22 = 0.30, rho1 = 0.3, rho2 = 0.3,
                                 alpha = 0.025, Test = tst)$powerCoprimary
      below <- power2BinaryApprox(n1 = ceiling(rv * (res$n2 - 1)),
                                  n2 = res$n2 - 1,
                                  p11 = 0.60, p12 = 0.50, p21 = 0.40,
                                  p22 = 0.30, rho1 = 0.3, rho2 = 0.3,
                                  alpha = 0.025, Test = tst)$powerCoprimary
      expect_gte(at_n, 0.8)
      expect_lt(below, 0.8)
    }
  }
})

test_that("design_table forwards nMC to the mixed continuous-binary path", {
  # Checked by seeding rather than by timing: an elapsed time is not a stable
  # assertion on a loaded machine, and the same seed with the same nMC must
  # give the same answer as the direct call only if nMC actually arrives
  grid_mcb <- data.frame(delta = 0.5, sd = 1, p1 = 0.60, p2 = 0.40)
  set.seed(4321)
  from_table <- design_table(
    param_grid = grid_mcb, rho_values = 0.3, r = 1, alpha = 0.025, beta = 0.2,
    endpoint_type = "mixed_cont_binary", Test = "Fisher",
    nMC = 200)[["rho_0.3"]]
  set.seed(4321)
  direct <- ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.60, p2 = 0.40,
                                     rho = 0.3, r = 1, alpha = 0.025,
                                     beta = 0.2, Test = "Fisher",
                                     nMC = 200)$N
  expect_equal(from_table, direct)
})
