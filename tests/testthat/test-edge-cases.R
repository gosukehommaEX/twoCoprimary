# Behaviour at the edges of the parameter space, and the agreement between a
# sample size function and the power function of the same family.
#
# Each block states one property and checks it on the smallest set of cases
# that can distinguish it, using the exact arguments that failed before the
# repair rather than a sweep. The sweeps themselves, over thousands of argument
# combinations, are in dev/audit_all_functions.R and dev/fuzz_all_functions.R,
# which are not shipped; what is here is the part that must stay inside the
# time a CRAN check allows.

# ==============================================================================
# The continuity correction at the smallest sample sizes
# ==============================================================================

test_that("the corrected arcsine power is a probability where the correction leaves (0, 1)", {
  # p12 + c1 is exactly zero here, which made the corrected variance vanish and
  # the reported power one half
  a <- power2BinaryApprox(n1 = 1, n2 = 1, p11 = 0.6, p12 = 0.5, p21 = 0.4,
                          p22 = 0.3, rho1 = 0, rho2 = 0, alpha = 0.025,
                          Test = "ASc")
  expect_equal(a$powerCoprimary, 0)

  # p12 + c1 is negative here, which gave NaN
  b <- power2BinaryApprox(n1 = 1, n2 = 1, p11 = 0.6, p12 = 0.4, p21 = 0.4,
                          p22 = 0.2, rho1 = 0, rho2 = 0, alpha = 0.025,
                          Test = "ASc")
  expect_equal(b$powerCoprimary, 0)

  # p21 + c2 reaches one at n = 5 with a control probability of 0.9
  d <- power2BinaryApprox(n1 = 5, n2 = 5, p11 = 0.99, p12 = 0.95, p21 = 0.90,
                          p22 = 0.85, rho1 = 0, rho2 = 0, alpha = 0.025,
                          Test = "ASc")
  expect_equal(d$powerCoprimary, 0)

  # One subject further on the correction is admissible again and the power is
  # a genuine, small number
  e <- power2BinaryApprox(n1 = 6, n2 = 6, p11 = 0.99, p12 = 0.95, p21 = 0.90,
                          p22 = 0.85, rho1 = 0, rho2 = 0, alpha = 0.025,
                          Test = "ASc")
  expect_true(is.finite(e$powerCoprimary))
  expect_gt(e$powerCoprimary, 0)
  expect_lt(e$powerCoprimary, 0.5)
})

test_that("the mixed corrected arcsine power is a probability at one subject per group", {
  res <- power2MixedContinuousBinary(n1 = 1, n2 = 1, delta = 0.5, sd = 1,
                                     p1 = 0.3, p2 = 0.1, rho = 0.5,
                                     alpha = 0.025, Test = "ASc")
  expect_equal(res$powerBin, 0)
  expect_equal(res$powerCoprimary, 0)
})

test_that("the power is zero, not undefined, where the variance cannot be estimated", {
  mixed <- power2MixedContinuousBinary(n1 = 1, n2 = 1, delta = 0.5, sd = 1,
                                       p1 = 0.6, p2 = 0.4, rho = 0.5,
                                       alpha = 0.025, Test = "Fisher",
                                       nMC = 50)
  expect_equal(mixed$powerCoprimary, 0)

  cont <- power2Continuous(n1 = 1, n2 = 1, delta1 = 4, delta2 = 4, sd1 = 1,
                           sd2 = 1, rho = 0.5, alpha = 0.025,
                           known_var = FALSE, nMC = 50)
  expect_equal(cont$powerCoprimary, 0)
})

# ==============================================================================
# The rejection region when a group has a single subject
# ==============================================================================

test_that("the unconditional tests build a rejection region when a group has one subject", {
  # Only one distinct responder count is then ordered, and the grid of null
  # probabilities used to be simplified from a matrix to a vector
  for (test in c("Z-pool", "Boschloo")) {
    for (nn in list(c(1, 1), c(1, 4), c(4, 1))) {
      RR <- rr1Binary(nn[1], nn[2], 0.025, Test = test, n_grid = 20)
      expect_true(is.logical(RR))
      expect_false(any(is.na(RR)))
      expect_equal(nrow(RR), nn[1] + 1)
      expect_equal(ncol(RR), nn[2] + 1)
    }
  }
})

# ==============================================================================
# The correlation at the Frechet-Hoeffding bound
# ==============================================================================

test_that("the bivariate binomial is still a distribution at the upper Prentice bound", {
  N <- 8
  p <- 0.5
  P <- outer(0:N, 0:N, function(a, b) dbibinom(N, a, b, p, p, 1))
  expect_true(all(is.finite(P)))
  expect_equal(sum(P), 1)
  expect_equal(rowSums(P), dbinom(0:N, N, p))
  expect_equal(colSums(P), dbinom(0:N, N, p))
  y <- 0:N
  covariance <- sum(P * outer(y, y)) - sum(rowSums(P) * y) * sum(colSums(P) * y)
  expect_equal(covariance / (N * p * (1 - p)), 1)
})

test_that("the co-primary power is the smaller marginal when the test statistics coincide", {
  for (test in c("AN", "ANc", "AS", "ASc")) {
    res <- power2BinaryApprox(n1 = 120, n2 = 120, p11 = 0.6, p12 = 0.6,
                              p21 = 0.4, p22 = 0.4, rho1 = 1, rho2 = 1,
                              alpha = 0.025, Test = test)
    expect_true(is.finite(res$powerCoprimary))
    expect_equal(res$powerCoprimary, min(res$power1, res$power2))
  }
  ex <- power2BinaryExact(n1 = 15, n2 = 15, p11 = 0.6, p12 = 0.6, p21 = 0.4,
                          p22 = 0.4, rho1 = 1, rho2 = 1, alpha = 0.025,
                          Test = "Fisher")
  expect_equal(ex$powerCoprimary, ex$power1)
})

test_that("the co-primary power is the Bonferroni bound at a correlation of minus one", {
  res <- power2BinaryApprox(n1 = 120, n2 = 120, p11 = 0.6, p12 = 0.4,
                            p21 = 0.4, p22 = 0.6, rho1 = -1, rho2 = -1,
                            alpha = 0.025, Test = "AN")
  expect_true(is.finite(res$powerCoprimary))
  expect_equal(res$powerCoprimary, max(0, res$power1 + res$power2 - 1))
})

# ==============================================================================
# The size and the power of the binary family solve the same inequality
# ==============================================================================

test_that("the published single endpoint sizes of Sozu et al. Table S5 are reproduced", {
  published <- data.frame(
    p1 = c(0.99, 0.99, 0.99, 0.95, 0.95, 0.95),
    p2 = c(0.95, 0.90, 0.85, 0.90, 0.85, 0.80),
    ANc = c(333, 121, 72, 474, 160, 88),
    ASc = c(300, 103, 60, 464, 153, 83)
  )
  for (i in seq_len(nrow(published))) {
    for (test in c("ANc", "ASc")) {
      got <- ss1BinaryApprox(p1 = published$p1[i], p2 = published$p2[i], r = 1,
                             alpha = 0.025, beta = 0.2, Test = test)$n2
      expect_equal(got, published[[test]][i])
    }
  }
})

test_that("the single endpoint size is the smallest reaching the target power", {
  power1_at <- function(p1, p2, r, m, test) {
    power2BinaryApprox(n1 = ceiling(r * m), n2 = m, p11 = p1, p12 = p1,
                       p21 = p2, p22 = p2, rho1 = 0, rho2 = 0, alpha = 0.025,
                       Test = test)[["power1"]]
  }
  for (test in c("AN", "ANc", "AS", "ASc")) {
    for (r in c(1, 1.5)) {
      res <- ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = r, alpha = 0.025,
                             beta = 0.2, Test = test)
      expect_gte(power1_at(0.6, 0.4, r, res$n2, test), 0.8)
      expect_lt(power1_at(0.6, 0.4, r, res$n2 - 1, test), 0.8)
    }
  }
})

test_that("the continuity corrected size returns for a design that made the iteration cycle", {
  res <- ss1BinaryApprox(p1 = 0.95, p2 = 0.05, r = 5, alpha = 0.025,
                         beta = 0.2, Test = "ANc")
  expect_true(is.finite(res$n2))
  expect_gte(res$n2, 1)
})

test_that("the Fisher size is not left above a smaller one that also reaches the target", {
  exact_power <- function(m) {
    RR <- rr1Binary(m, m, 0.025, Test = "Fisher")
    sum(dbinom(0:m, m, 0.9) * pbinom(rowSums(RR) - 1, m, 0.3))
  }
  res <- ss1BinaryApprox(p1 = 0.9, p2 = 0.3, r = 1, alpha = 0.025, beta = 0.2,
                         Test = "Fisher")
  expect_gte(exact_power(res$n2), 0.8)
  expect_lt(exact_power(res$n2 - 1), 0.8)
})

# ==============================================================================
# The joint probability against the bounds its marginals allow
# ==============================================================================

test_that("the co-primary power is never a small negative number", {
  # The lower Prentice bound with nearly equal marginals sent pbivnorm a few
  # times 1e-19 below zero
  b1 <- corrbound2Binary(0.55, 0.50)
  b2 <- corrbound2Binary(0.45, 0.40)
  rho <- max(b1[["L_bound"]], b2[["L_bound"]])
  res <- power2BinaryApprox(n1 = 1, n2 = 1, p11 = 0.55, p12 = 0.50, p21 = 0.45,
                            p22 = 0.40, rho1 = rho, rho2 = rho, alpha = 0.025,
                            Test = "ANc")
  expect_gte(res$powerCoprimary, 0)

  mix <- power2MixedContinuousBinary(n1 = 1, n2 = 300, delta = 0.5, sd = 1,
                                     p1 = 0.99, p2 = 0.85, rho = -0.99,
                                     alpha = 0.025, Test = "ASc")
  expect_gte(mix$powerCoprimary, 0)
})

test_that("the co-primary power is a number at a very large standardized effect", {
  # The critical values are then some 2500 standard deviations out, where
  # pbivnorm returns NaN
  res <- power2Continuous(n1 = 50, n2 = 50, delta1 = 5, delta2 = 5, sd1 = 0.01,
                          sd2 = 0.01, rho = -0.999, alpha = 0.025,
                          known_var = TRUE)
  expect_true(is.finite(res$powerCoprimary))
  expect_equal(res$powerCoprimary, 1)
})

test_that("the simulated co-primary power does not exceed an exact marginal", {
  set.seed(20260902)
  res <- power2Continuous(n1 = 1, n2 = 3, delta1 = 5, delta2 = 5, sd1 = 1,
                          sd2 = 1.2, rho = 0.999, alpha = 0.025,
                          known_var = FALSE, nMC = 200)
  expect_lte(res$powerCoprimary, min(res$power1, res$power2))
  expect_gte(res$powerCoprimary, max(0, res$power1 + res$power2 - 1))
})

test_that("the unknown variance power exists at one degree of freedom", {
  # rWishart refuses a single degree of freedom for a two by two scale matrix
  set.seed(20260902)
  res <- power2Continuous(n1 = 1, n2 = 2, delta1 = 1, delta2 = 1, sd1 = 1,
                          sd2 = 1, rho = 0.5, alpha = 0.025, known_var = FALSE,
                          nMC = 200)
  expect_true(is.finite(res$powerCoprimary))
  expect_lte(res$powerCoprimary, min(res$power1, res$power2))
})

# ==============================================================================
# The correlation bounds for a count and a continuous endpoint
# ==============================================================================

test_that("the count and continuous bounds do not depend on the location or the scale", {
  # A correlation cannot depend on either, but the quadrature was taken in the
  # untransformed variable and returned zero for a mean far from the origin
  base <- corrbound2MixedCountContinuous(lambda = 1.25, nu = 0.8, mu = 0,
                                         sd = 1)
  expect_lt(base[["L_bound"]], 0)
  expect_gt(base[["U_bound"]], 0)
  expect_equal(corrbound2MixedCountContinuous(1.25, 0.8, -50, 0.5), base)
  expect_equal(corrbound2MixedCountContinuous(1.25, 0.8, 100, 250), base)
})

test_that("a correlation of zero is accepted at a small continuous scale", {
  res <- power2MixedCountContinuous(n1 = 3, n2 = 2, r1 = 1, r2 = 1.25, nu = 5,
                                    t = 1, mu1 = -50, mu2 = 0, sd = 0.5,
                                    rho1 = 0, rho2 = 0, alpha = 0.025)
  expect_true(is.finite(res$powerCoprimary))
})

# ==============================================================================
# The smallest sample size any single endpoint function may return
# ==============================================================================

test_that("a single endpoint sample size is at least one subject per group", {
  # The closed forms return zero when the target power does not exceed the size
  # of the test
  expect_equal(ss1Continuous(delta = 0.5, sd = 1, r = 1, alpha = 0.1,
                             beta = 0.9)$n2, 1)
  expect_equal(ss1Count(r1 = 1, r2 = 1.25, nu = 0.8, t = 1, r = 1, alpha = 0.1,
                        beta = 0.9)$n2, 1)
  for (test in c("AN", "ANc", "AS", "ASc")) {
    res <- ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = 0.25, alpha = 0.1,
                           beta = 0.9, Test = test)
    expect_gte(res$n2, 1)
    expect_gte(res$n1, 1)
  }
})

# ==============================================================================
# Argument validation added in 1.1.1
# ==============================================================================

test_that("power2Continuous validates its arguments like the other power functions", {
  expect_error(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 2, 0.025), "rho")
  expect_error(power2Continuous(100, 100, 0.5, 0.5, -1, 1, 0.5, 0.025), "sd")
  expect_error(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 0.5, 2), "alpha")
  expect_error(power2Continuous(100, 0, 0.5, 0.5, 1, 1, 0.5, 0.025), "n2")
  expect_error(power2Continuous(100, 2.5, 0.5, 0.5, 1, 1, 0.5, 0.025), "n2")
  expect_error(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 0.5, 0.025,
                                known_var = "yes"), "known_var")
})

test_that("power2BinaryApprox validates its arguments like power2BinaryExact", {
  expect_error(power2BinaryApprox(0, 100, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025,
                                  "AN"), "n1")
  expect_error(power2BinaryApprox(100, 2.5, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025,
                                  "AN"), "n2")
  expect_error(power2BinaryApprox(100, 100, 1, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025,
                                  "AN"), "probabilities")
  expect_error(power2BinaryApprox(100, 100, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 1,
                                  "AN"), "alpha")
})

test_that("ss1Count requires the treatment rate to be the lower one", {
  expect_error(ss1Count(r1 = 1.25, r2 = 1.0, nu = 0.8, t = 1, r = 1,
                        alpha = 0.025, beta = 0.2), "less than")
  expect_error(ss1Count(r1 = 1.0, r2 = 1.0, nu = 0.8, t = 1, r = 1,
                        alpha = 0.025, beta = 0.2), "less than")
})

# ==============================================================================
# The plot method at the edges
# ==============================================================================

test_that("the effect contour axes are standardized effect sizes", {
  pdf(NULL); on.exit(dev.off())
  obj <- power2Continuous(n1 = 60, n2 = 60, delta1 = 1.0, delta2 = 1.2,
                          sd1 = 2, sd2 = 3, rho = 0.4, alpha = 0.025,
                          known_var = TRUE)
  dat <- plot(obj, type = "effect_contour", n_points = 4)
  expect_equal(range(dat$delta1), c(0.2, 1.0))
  expect_equal(range(dat$delta2), c(0.2, 1.0))
  expect_true(all(is.finite(dat$power)))
})

test_that("the power curve window does not invert for a design of a few per group", {
  pdf(NULL); on.exit(dev.off())
  obj <- power2Continuous(n1 = 5, n2 = 5, delta1 = 1.5, delta2 = 1.5, sd1 = 1,
                          sd2 = 1, rho = 0.5, alpha = 0.025, known_var = TRUE)
  dat <- plot(obj, type = "power_curve", n_points = 5)
  expect_true(all(diff(dat$n2) > 0))
  expect_gte(min(dat$n2), 2)
})

# ==============================================================================
# Structural contract, at an allocation ratio that is not an integer
# ==============================================================================

test_that("every sample size result satisfies N = n1 + n2 and n1 = ceiling(r n2)", {
  results <- list(
    ss1Continuous(delta = 0.5, sd = 1, r = 1.5, alpha = 0.025, beta = 0.2),
    ss1Count(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1, r = 1.5, alpha = 0.025,
             beta = 0.2),
    ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = 1.5, alpha = 0.025, beta = 0.2,
                    Test = "ANc"),
    ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1, rho = 0.5,
                  r = 1.5, alpha = 0.025, beta = 0.2, known_var = TRUE),
    ss2BinaryApprox(p11 = 0.6, p12 = 0.5, p21 = 0.4, p22 = 0.3, rho1 = 0.3,
                    rho2 = 0.3, r = 1.5, alpha = 0.025, beta = 0.2,
                    Test = "ASc"),
    ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.6, p2 = 0.4,
                             rho = 0.5, r = 1.5, alpha = 0.025, beta = 0.2,
                             Test = "AN")
  )
  for (res in results) {
    expect_s3_class(res, "twoCoprimary")
    expect_equal(nrow(res), 1L)
    expect_equal(res$N, res$n1 + res$n2)
    expect_equal(res$n1, ceiling(res$r * res$n2))
  }
})
