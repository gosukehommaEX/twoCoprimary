# Test sample size calculation functions

# ss1Continuous tests
test_that("ss1Continuous returns valid sample size", {
  result <- ss1Continuous(delta = 0.5, sd = 1, r = 1, alpha = 0.025, beta = 0.1)

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("n1", "n2", "N") %in% names(result)))
  expect_true(result$n1 > 0)
  expect_true(result$n2 > 0)
  expect_equal(result$N, result$n1 + result$n2)
})

test_that("ss1Continuous sample size increases with smaller effect", {
  n1 <- ss1Continuous(delta = 0.5, sd = 1, r = 1, alpha = 0.025, beta = 0.1)$N
  n2 <- ss1Continuous(delta = 0.3, sd = 1, r = 1, alpha = 0.025, beta = 0.1)$N

  expect_true(n2 > n1)
})

# ss1BinaryApprox tests
test_that("ss1BinaryApprox returns valid sample size", {
  result <- ss1BinaryApprox(p1 = 0.5, p2 = 0.3, r = 1, alpha = 0.025, beta = 0.1, Test = "AN")

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("n1", "n2", "N") %in% names(result)))
  expect_true(result$n1 > 0)
  expect_true(result$n2 > 0)
})

test_that("ss1BinaryApprox works with different test methods", {
  tests <- c("AN", "ANc", "AS", "ASc")

  for (test in tests) {
    result <- ss1BinaryApprox(p1 = 0.5, p2 = 0.3, r = 1, alpha = 0.025, beta = 0.1, Test = test)
    expect_s3_class(result, "twoCoprimary")
    expect_true(result$N > 0)
  }
})

# ss1Count tests
test_that("ss1Count returns valid sample size", {
  result <- ss1Count(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1, r = 1, alpha = 0.025, beta = 0.1)

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("n1", "n2", "N") %in% names(result)))
  expect_true(result$n1 > 0)
  expect_true(result$n2 > 0)
})

# ss2Continuous tests
test_that("ss2Continuous returns valid sample size with known variance", {
  result <- ss2Continuous(
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, r = 1,
    alpha = 0.025, beta = 0.1,
    known_var = TRUE
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("n1", "n2", "N") %in% names(result)))
  expect_true(result$n1 > 0)
  expect_true(result$n2 > 0)
})

test_that("ss2Continuous sample size decreases with correlation", {
  # Higher correlation should reduce required sample size (with equal individual power)
  n_rho0 <- ss2Continuous(
    delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
    rho = 0, r = 1, alpha = 0.025, beta = 0.1, known_var = TRUE
  )$N

  n_rho05 <- ss2Continuous(
    delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
    rho = 0.5, r = 1, alpha = 0.025, beta = 0.1, known_var = TRUE
  )$N

  expect_true(n_rho05 < n_rho0)
})

# ss2BinaryApprox tests
test_that("ss2BinaryApprox returns valid sample size", {
  result <- ss2BinaryApprox(
    p11 = 0.5, p12 = 0.4,
    p21 = 0.3, p22 = 0.2,
    rho1 = 0.5, rho2 = 0.5,
    r = 1, alpha = 0.025, beta = 0.1,
    Test = "AN"
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("n1", "n2", "N") %in% names(result)))
  expect_true(result$n1 > 0)
  expect_true(result$n2 > 0)
})

# ss2BinaryExact tests (use larger tolerances for target power)
test_that("ss2BinaryExact returns valid sample size", {
  # Use simpler parameters for efficiency
  result <- ss2BinaryExact(
    p11 = 0.6, p12 = 0.5,
    p21 = 0.4, p22 = 0.3,
    rho1 = 0.3, rho2 = 0.3,
    r = 1, alpha = 0.025, beta = 0.2,  # Use beta=0.2 for smaller sample
    Test = "Chisq"
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("n1", "n2", "N") %in% names(result)))
  expect_true(result$n1 > 0)
  expect_true(result$n2 > 0)
})

# ss2MixedContinuousBinary tests
test_that("ss2MixedContinuousBinary returns valid sample size", {
  result <- ss2MixedContinuousBinary(
    delta = 0.5, sd = 1,
    p1 = 0.6, p2 = 0.4,
    rho = 0.5, r = 1,
    alpha = 0.025, beta = 0.1,
    Test = "AN"
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("n1", "n2", "N") %in% names(result)))
  expect_true(result$n1 > 0)
  expect_true(result$n2 > 0)
})

# ss2MixedCountContinuous tests
test_that("ss2MixedCountContinuous returns valid sample size", {
  result <- ss2MixedCountContinuous(
    r1 = 1.0, r2 = 1.25,
    nu = 0.8, t = 1,
    mu1 = -50, mu2 = 0, sd = 250,
    rho1 = 0.5, rho2 = 0.5,
    r = 1, alpha = 0.025, beta = 0.1
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("n1", "n2", "N") %in% names(result)))
  expect_true(result$n1 > 0)
  expect_true(result$n2 > 0)
})

# Validation tests
test_that("sample size functions validate input parameters", {
  # Invalid alpha
  expect_error(ss1Continuous(delta = 0.5, sd = 1, r = 1, alpha = 1.5, beta = 0.1))

  # Invalid beta
  expect_error(ss1Continuous(delta = 0.5, sd = 1, r = 1, alpha = 0.025, beta = 1.5))

  # Invalid probabilities
  expect_error(ss1BinaryApprox(p1 = 1.5, p2 = 0.3, r = 1, alpha = 0.025, beta = 0.1, Test = "AN"))

  # Invalid allocation ratio
  expect_error(ss2Continuous(
    delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
    rho = 0.3, r = -1, alpha = 0.025, beta = 0.1, known_var = TRUE
  ))
})

test_that("calculated sample size achieves target power", {
  # Calculate sample size
  ss_result <- ss2Continuous(
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, r = 1,
    alpha = 0.025, beta = 0.1,
    known_var = TRUE
  )

  # Calculate power with the obtained sample size
  power_result <- power2Continuous(
    n1 = ss_result$n1, n2 = ss_result$n2,
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, alpha = 0.025,
    known_var = TRUE
  )

  # Power should be at least the target (1 - beta)
  expect_gte(power_result$powerCoprimary, 0.9)
})
