# Invariance tests that apply across all five endpoint type combinations.
#
# These check structural properties that any correct implementation must have,
# so they need no reference values and they catch the class of defect that the
# unknown variance branch of power2Continuous() had, where a standard deviation
# failed to cancel.
#
# Only the four asymptotic methods of the mixed continuous and binary
# combination are exercised here. Its Fisher option is simulation based, so it
# satisfies these identities only up to Monte Carlo error and is tested
# separately in test-mixed_fisher.R.

# ==============================================================================
# Independence at zero correlation
# ==============================================================================
# Under the intersection union test the two marginal test statistics are
# independent when the endpoints are uncorrelated, so the co-primary power is
# the product of the two marginal powers.

test_that("two continuous endpoints factorise at zero correlation", {
  res <- power2Continuous(60, 60, 0.5, 0.5, 2, 3, rho = 0, alpha = 0.025,
                          known_var = TRUE)
  expect_equal(res$powerCoprimary, res$power1 * res$power2, tolerance = 1e-10)
})

test_that("two binary endpoints factorise at zero correlation", {
  for (tst in c("AN", "ANc", "AS", "ASc")) {
    res <- power2BinaryApprox(80, 80, 0.55, 0.60, 0.30, 0.35,
                              rho1 = 0, rho2 = 0, alpha = 0.025, Test = tst)
    expect_equal(res$powerCoprimary, res$power1 * res$power2,
                 tolerance = 1e-10, info = tst)
  }

  for (tst in c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")) {
    res <- power2BinaryExact(25, 25, 0.55, 0.60, 0.30, 0.35,
                             rho1 = 0, rho2 = 0, alpha = 0.025, Test = tst)
    expect_equal(res$powerCoprimary, res$power1 * res$power2,
                 tolerance = 1e-10, info = tst)
  }
})

test_that("mixed continuous and binary endpoints factorise at zero correlation", {
  for (tst in c("AN", "ANc", "AS", "ASc")) {
    res <- power2MixedContinuousBinary(80, 80, delta = 0.5, sd = 1,
                                       p1 = 0.60, p2 = 0.40, rho = 0,
                                       alpha = 0.025, Test = tst)
    expect_equal(res$powerCoprimary, res$powerCont * res$powerBin,
                 tolerance = 1e-10, info = tst)
  }
})

test_that("mixed count and continuous endpoints factorise at zero correlation", {
  res <- power2MixedCountContinuous(
    n1 = 150, n2 = 150, r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
    mu1 = -50, mu2 = 0, sd = 250, rho1 = 0, rho2 = 0, alpha = 0.025
  )
  expect_equal(res$powerCoprimary, res$powerCount * res$powerCont,
               tolerance = 1e-10)
})

# ==============================================================================
# Scale invariance of the continuous component
# ==============================================================================
# Power depends on the continuous endpoint only through the standardized effect,
# so multiplying the effect and the standard deviation by the same constant must
# leave every reported power unchanged.

test_that("mixed continuous and binary power is invariant to the scale", {
  for (tst in c("AN", "ANc", "AS", "ASc")) {
    base <- power2MixedContinuousBinary(80, 80, delta = 0.5, sd = 1,
                                        p1 = 0.60, p2 = 0.40, rho = 0.4,
                                        alpha = 0.025, Test = tst)
    for (cc in c(7, 250)) {
      scaled <- power2MixedContinuousBinary(80, 80, delta = cc * 0.5, sd = cc,
                                            p1 = 0.60, p2 = 0.40, rho = 0.4,
                                            alpha = 0.025, Test = tst)
      expect_equal(scaled$powerCont, base$powerCont, tolerance = 1e-12,
                   info = paste(tst, cc))
      expect_equal(scaled$powerBin, base$powerBin, tolerance = 1e-12,
                   info = paste(tst, cc))
      expect_equal(scaled$powerCoprimary, base$powerCoprimary,
                   tolerance = 1e-12, info = paste(tst, cc))
    }
  }
})

test_that("mixed count and continuous power is invariant to the scale", {
  base <- power2MixedCountContinuous(
    n1 = 150, n2 = 150, r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
    mu1 = -50, mu2 = 0, sd = 250, rho1 = 0.5, rho2 = 0.5, alpha = 0.025
  )
  for (cc in c(0.004, 20)) {
    scaled <- power2MixedCountContinuous(
      n1 = 150, n2 = 150, r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
      mu1 = cc * -50, mu2 = 0, sd = cc * 250,
      rho1 = 0.5, rho2 = 0.5, alpha = 0.025
    )
    expect_equal(scaled$powerCount, base$powerCount, tolerance = 1e-12)
    expect_equal(scaled$powerCont, base$powerCont, tolerance = 1e-12)
    expect_equal(scaled$powerCoprimary, base$powerCoprimary, tolerance = 1e-12)
  }
})

test_that("mixed count and continuous power depends only on the mean difference", {
  base <- power2MixedCountContinuous(
    n1 = 150, n2 = 150, r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
    mu1 = -50, mu2 = 0, sd = 250, rho1 = 0.5, rho2 = 0.5, alpha = 0.025
  )
  shifted <- power2MixedCountContinuous(
    n1 = 150, n2 = 150, r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
    mu1 = 950, mu2 = 1000, sd = 250, rho1 = 0.5, rho2 = 0.5, alpha = 0.025
  )
  expect_equal(shifted$powerCoprimary, base$powerCoprimary, tolerance = 1e-12)
})

# ==============================================================================
# The co-primary power never exceeds either marginal power
# ==============================================================================

test_that("co-primary power is bounded by the marginal powers", {
  res <- power2Continuous(50, 50, 0.5, 0.4, 1, 1, 0.3, 0.025, known_var = TRUE)
  expect_lte(res$powerCoprimary, min(res$power1, res$power2) + 1e-12)

  res <- power2BinaryApprox(80, 80, 0.55, 0.60, 0.30, 0.35, 0.3, 0.3,
                            0.025, "AN")
  expect_lte(res$powerCoprimary, min(res$power1, res$power2) + 1e-12)

  res <- power2MixedContinuousBinary(80, 80, 0.5, 1, 0.60, 0.40, 0.4,
                                     0.025, "AN")
  expect_lte(res$powerCoprimary, min(res$powerCont, res$powerBin) + 1e-12)

  res <- power2MixedCountContinuous(150, 150, 1.0, 1.25, 0.8, 1,
                                    -50, 0, 250, 0.5, 0.5, 0.025)
  expect_lte(res$powerCoprimary, min(res$powerCount, res$powerCont) + 1e-12)
})

# ==============================================================================
# Monotonicity in the sample size
# ==============================================================================

test_that("power increases with the sample size for every endpoint type", {
  ns <- c(40, 80, 160)

  p_cont <- vapply(ns, function(n) {
    power2Continuous(n, n, 0.4, 0.4, 1, 1, 0.3, 0.025,
                     known_var = TRUE)$powerCoprimary
  }, numeric(1))
  expect_true(all(diff(p_cont) > 0))

  p_bin <- vapply(ns, function(n) {
    power2BinaryApprox(n, n, 0.55, 0.60, 0.30, 0.35, 0.3, 0.3,
                       0.025, "AN")$powerCoprimary
  }, numeric(1))
  expect_true(all(diff(p_bin) > 0))

  p_mcb <- vapply(ns, function(n) {
    power2MixedContinuousBinary(n, n, 0.4, 1, 0.60, 0.40, 0.4,
                                0.025, "AN")$powerCoprimary
  }, numeric(1))
  expect_true(all(diff(p_mcb) > 0))

  p_mcc <- vapply(c(100, 150, 250), function(n) {
    power2MixedCountContinuous(n, n, 1.0, 1.25, 0.8, 1,
                               -50, 0, 250, 0.5, 0.5, 0.025)$powerCoprimary
  }, numeric(1))
  expect_true(all(diff(p_mcc) > 0))
})

# ==============================================================================
# The sample size functions agree with the power functions
# ==============================================================================
# The returned sample size must be the smallest one reaching the target, so the
# power at that size clears the target and the power one step below does not.

test_that("the sample size is the smallest one reaching the target", {
  target <- 0.8

  ss <- ss2MixedContinuousBinary(
    delta = 0.5, sd = 1, p1 = 0.60, p2 = 0.40, rho = 0.4,
    r = 1, alpha = 0.025, beta = 0.2, Test = "AN"
  )
  at <- power2MixedContinuousBinary(ss$n1, ss$n2, 0.5, 1, 0.60, 0.40, 0.4,
                                    0.025, "AN")$powerCoprimary
  below <- power2MixedContinuousBinary(ss$n1 - 1, ss$n2 - 1, 0.5, 1,
                                       0.60, 0.40, 0.4, 0.025, "AN")$powerCoprimary
  expect_gte(at, target)
  expect_lt(below, target)

  ss <- ss2MixedCountContinuous(
    r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
    mu1 = -50, mu2 = 0, sd = 250,
    rho1 = 0.5, rho2 = 0.5, r = 1, alpha = 0.025, beta = 0.2
  )
  at <- power2MixedCountContinuous(ss$n1, ss$n2, 1.0, 1.25, 0.8, 1,
                                   -50, 0, 250, 0.5, 0.5, 0.025)$powerCoprimary
  below <- power2MixedCountContinuous(ss$n1 - 1, ss$n2 - 1, 1.0, 1.25, 0.8, 1,
                                      -50, 0, 250, 0.5, 0.5, 0.025)$powerCoprimary
  expect_gte(at, target)
  expect_lt(below, target)
})

