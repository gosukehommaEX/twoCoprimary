# Test rejection region function

test_that("rr1Binary returns matrix of correct dimensions", {
  # Use small sample sizes for efficiency
  n1 <- 5
  n2 <- 5
  alpha <- 0.025

  result <- rr1Binary(n1, n2, alpha, Test = "Chisq")

  expect_type(result, "logical")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n1 + 1, n2 + 1))
})

test_that("rr1Binary works with different test methods", {
  n1 <- 5
  n2 <- 5
  alpha <- 0.025

  # Test all methods
  tests <- c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")

  for (test in tests) {
    result <- rr1Binary(n1, n2, alpha, Test = test)
    expect_true(is.matrix(result))
    expect_equal(dim(result), c(n1 + 1, n2 + 1))
    expect_type(result, "logical")
  }
})

test_that("rr1Binary rejection region is monotone", {
  # Small sample for efficiency
  n1 <- 8
  n2 <- 8
  alpha <- 0.025

  RR <- rr1Binary(n1, n2, alpha, Test = "Chisq")

  # Check that rejection region is in the upper-right corner
  # If (i,j) is in rejection region and i' > i, j' < j, then (i',j') should also be in rejection
  # This is a simplified monotonicity check
  expect_true(all(!RR[1, ]))  # Bottom row should not reject
  expect_true(all(!RR[, ncol(RR)]))  # Rightmost column should not reject
})

test_that("rr1Binary validates input", {
  expect_error(rr1Binary(n1 = 0, n2 = 5, alpha = 0.025, Test = "Chisq"))
  expect_error(rr1Binary(n1 = 5, n2 = -1, alpha = 0.025, Test = "Chisq"))
  expect_error(rr1Binary(n1 = 5, n2 = 5, alpha = 1.5, Test = "Chisq"))
  expect_error(rr1Binary(n1 = 5, n2 = 5, alpha = 0.025, Test = "InvalidTest"))
})

test_that("rr1Binary type I error is controlled", {
  # Check that under H0, rejection rate is approximately alpha
  # Using small sample for efficiency
  n1 <- 10
  n2 <- 10
  alpha <- 0.05  # Use larger alpha for better precision with small samples
  p0 <- 0.5  # Common probability under H0

  RR <- rr1Binary(n1, n2, alpha, Test = "Fisher")

  # Calculate rejection probability under H0
  reject_prob <- 0
  for (i in 0:n1) {
    for (j in 0:n2) {
      if (RR[i + 1, j + 1]) {
        reject_prob <- reject_prob + dbinom(i, n1, p0) * dbinom(j, n2, p0)
      }
    }
  }

  # Should be approximately equal to alpha (allowing some tolerance)
  expect_lt(reject_prob, alpha * 1.5)  # Not too liberal
})
