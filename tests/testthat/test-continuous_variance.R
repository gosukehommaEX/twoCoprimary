# Regression tests for the unknown variance branch of power2Continuous().
#
# The test statistic is already divided by sd1 and sd2, so the Wishart matrix
# entering the critical value must be that of the standardized endpoints.
# Drawing it from the variance-covariance matrix left a factor of sd_k in
# sqrt(W_kk / nu), which cancelled only when sd1 = sd2 = 1.

test_that("power with unknown variance is invariant to the scale", {
  designs <- list(
    list(n1 = 50, n2 = 50, delta1 = 0.5, delta2 = 0.5, rho = 0.0),
    list(n1 = 50, n2 = 100, delta1 = 0.4, delta2 = 0.6, rho = 0.3),
    list(n1 = 40, n2 = 40, delta1 = 0.6, delta2 = 0.4, rho = 0.8)
  )

  for (d in designs) {
    values <- vapply(c(1, 5, 250), function(cc) {
      set.seed(1234)
      power2Continuous(
        n1 = d$n1, n2 = d$n2,
        delta1 = cc * d$delta1, delta2 = cc * d$delta2,
        sd1 = cc, sd2 = cc, rho = d$rho, alpha = 0.025,
        known_var = FALSE, nMC = 5000
      )$powerCoprimary
    }, numeric(1))

    expect_equal(values[2], values[1], tolerance = 1e-10)
    expect_equal(values[3], values[1], tolerance = 1e-10)
  }
})

test_that("power with unknown variance depends only on the standardized effect", {
  set.seed(1234)
  a <- power2Continuous(60, 60, delta1 = 1.5, delta2 = 1.5,
                        sd1 = 3, sd2 = 3, rho = 0.3, alpha = 0.025,
                        known_var = FALSE, nMC = 5000)$powerCoprimary
  set.seed(1234)
  b <- power2Continuous(60, 60, delta1 = 0.5, delta2 = 0.5,
                        sd1 = 1, sd2 = 1, rho = 0.3, alpha = 0.025,
                        known_var = FALSE, nMC = 5000)$powerCoprimary

  expect_equal(a, b, tolerance = 1e-10)
})

test_that("unequal standard deviations behave like unequal standardized effects", {
  set.seed(9876)
  a <- power2Continuous(60, 60, delta1 = 1.0, delta2 = 1.5,
                        sd1 = 2, sd2 = 3, rho = 0.4, alpha = 0.025,
                        known_var = FALSE, nMC = 5000)$powerCoprimary
  set.seed(9876)
  b <- power2Continuous(60, 60, delta1 = 0.5, delta2 = 0.5,
                        sd1 = 1, sd2 = 1, rho = 0.4, alpha = 0.025,
                        known_var = FALSE, nMC = 5000)$powerCoprimary

  expect_equal(a, b, tolerance = 1e-10)
})

test_that("power with unknown variance is below the known variance value", {
  for (rho in c(0, 0.5)) {
    known <- power2Continuous(40, 40, 0.5, 0.5, 4, 4, rho, 0.025,
                              known_var = TRUE)$powerCoprimary
    set.seed(2468)
    unknown <- power2Continuous(40, 40, 0.5, 0.5, 4, 4, rho, 0.025,
                                known_var = FALSE, nMC = 8000)$powerCoprimary

    expect_true(unknown < known)
    expect_true(unknown > known - 0.1)
  }
})

test_that("the known variance branch is unaffected by the scale", {
  a <- power2Continuous(50, 50, 12.5, 12.5, 25, 25, 0.3, 0.025,
                        known_var = TRUE)$powerCoprimary
  b <- power2Continuous(50, 50, 0.5, 0.5, 1, 1, 0.3, 0.025,
                        known_var = TRUE)$powerCoprimary

  expect_equal(a, b, tolerance = 1e-12)
})

test_that("sample size with unknown variance exceeds the known variance value", {
  set.seed(1357)
  known <- ss2Continuous(0.5, 0.5, 2, 2, rho = 0.3, r = 1,
                         alpha = 0.025, beta = 0.2,
                         known_var = TRUE)[["N"]]
  set.seed(1357)
  unknown <- ss2Continuous(0.5, 0.5, 2, 2, rho = 0.3, r = 1,
                           alpha = 0.025, beta = 0.2,
                           known_var = FALSE, nMC = 2000)[["N"]]

  expect_true(unknown >= known)
})
