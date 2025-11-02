# Test power calculation functions

# power2Continuous tests
test_that("power2Continuous returns valid power with known variance", {
  result <- power2Continuous(
    n1 = 50, n2 = 50,
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, alpha = 0.025,
    known_var = TRUE
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("power1", "power2", "powerCoprimary") %in% names(result)))
  expect_true(result$power1 >= 0 && result$power1 <= 1)
  expect_true(result$power2 >= 0 && result$power2 <= 1)
  expect_true(result$powerCoprimary >= 0 && result$powerCoprimary <= 1)
  expect_true(result$powerCoprimary <= min(result$power1, result$power2))
})

test_that("power2Continuous increases with sample size", {
  powers <- numeric(3)
  for (i in 1:3) {
    n <- 30 * i
    result <- power2Continuous(
      n1 = n, n2 = n,
      delta1 = 0.5, delta2 = 0.5,
      sd1 = 1, sd2 = 1,
      rho = 0.3, alpha = 0.025,
      known_var = TRUE
    )
    powers[i] <- result$powerCoprimary
  }

  # Power should increase with sample size
  expect_true(all(diff(powers) > 0))
})

# power2BinaryApprox tests
test_that("power2BinaryApprox returns valid power", {
  result <- power2BinaryApprox(
    n1 = 100, n2 = 50,
    p11 = 0.5, p12 = 0.4,
    p21 = 0.3, p22 = 0.2,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025, Test = "AN"
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("power1", "power2", "powerCoprimary") %in% names(result)))
  expect_true(result$powerCoprimary >= 0 && result$powerCoprimary <= 1)
  expect_true(result$powerCoprimary <= min(result$power1, result$power2))
})

test_that("power2BinaryApprox works with different test methods", {
  tests <- c("AN", "ANc", "AS", "ASc")

  for (test in tests) {
    result <- power2BinaryApprox(
      n1 = 50, n2 = 50,
      p11 = 0.6, p12 = 0.5,
      p21 = 0.4, p22 = 0.3,
      rho1 = 0.5, rho2 = 0.5,
      alpha = 0.025, Test = test
    )

    expect_s3_class(result, "twoCoprimary")
    expect_true(result$powerCoprimary >= 0 && result$powerCoprimary <= 1)
  }
})

# power2BinaryExact tests (use small samples for efficiency)
test_that("power2BinaryExact returns valid power", {
  result <- power2BinaryExact(
    n1 = 20, n2 = 20,  # Small sample for efficiency
    p11 = 0.5, p12 = 0.4,
    p21 = 0.3, p22 = 0.2,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025, Test = "Chisq"
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("power1", "power2", "powerCoprimary") %in% names(result)))
  expect_true(result$powerCoprimary >= 0 && result$powerCoprimary <= 1)
  expect_true(result$powerCoprimary <= min(result$power1, result$power2))
})

# power2MixedContinuousBinary tests
test_that("power2MixedContinuousBinary returns valid power for AN", {
  result <- power2MixedContinuousBinary(
    n1 = 50, n2 = 50,
    delta = 0.5, sd = 1,
    p1 = 0.6, p2 = 0.4,
    rho = 0.5, alpha = 0.025,
    Test = "AN"
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("powerCont", "powerBin", "powerCoprimary") %in% names(result)))
  expect_true(result$powerCoprimary >= 0 && result$powerCoprimary <= 1)
  expect_true(result$powerCoprimary <= min(result$powerCont, result$powerBin))
})

# power2MixedCountContinuous tests
test_that("power2MixedCountContinuous returns valid power", {
  result <- power2MixedCountContinuous(
    n1 = 100, n2 = 100,
    r1 = 1.0, r2 = 1.25,
    nu = 0.8, t = 1,
    mu1 = -50, mu2 = 0, sd = 250,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025
  )

  expect_s3_class(result, "twoCoprimary")
  expect_true(all(c("powerCount", "powerCont", "powerCoprimary") %in% names(result)))
  expect_true(result$powerCoprimary >= 0 && result$powerCoprimary <= 1)
  expect_true(result$powerCoprimary <= min(result$powerCount, result$powerCont))
})

# General validation tests
test_that("power functions validate correlation bounds", {
  # Binary endpoints
  expect_error(power2BinaryApprox(
    n1 = 50, n2 = 50,
    p11 = 0.5, p12 = 0.4,
    p21 = 0.3, p22 = 0.2,
    rho1 = 1.5, rho2 = 0.5,  # Invalid correlation
    alpha = 0.025, Test = "AN"
  ))

  # Count-continuous endpoints
  expect_error(power2MixedCountContinuous(
    n1 = 100, n2 = 100,
    r1 = 1.0, r2 = 1.25,
    nu = 0.8, t = 1,
    mu1 = -50, mu2 = 0, sd = 250,
    rho1 = 1.5, rho2 = 0.5,  # Invalid correlation
    alpha = 0.025
  ))
})

test_that("power functions validate input parameters", {
  # Invalid probability
  expect_error(power2BinaryApprox(
    n1 = 50, n2 = 50,
    p11 = 1.5, p12 = 0.4,  # Invalid probability
    p21 = 0.3, p22 = 0.2,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025, Test = "AN"
  ))

  # Negative sample size
  expect_error(power2BinaryExact(
    n1 = -10, n2 = 20,  # Negative n1
    p11 = 0.5, p12 = 0.4,
    p21 = 0.3, p22 = 0.2,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025, Test = "Chisq"
  ))

  # Invalid test method
  expect_error(power2BinaryApprox(
    n1 = 50, n2 = 50,
    p11 = 0.5, p12 = 0.4,
    p21 = 0.3, p22 = 0.2,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025, Test = "InvalidTest"
  ))
})
