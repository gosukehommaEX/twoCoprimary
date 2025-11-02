# Test bivariate binomial distribution function

test_that("dbibinom returns valid probabilities", {
  # Basic test with small N
  result <- dbibinom(N = 10, y1 = 5, y2 = 5, p1 = 0.5, p2 = 0.5, rho = 0.5)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("dbibinom probabilities sum to 1", {
  # Small N for efficiency
  N <- 10
  p1 <- 0.3
  p2 <- 0.5
  rho <- 0.3

  total_prob <- sum(outer(0:N, 0:N, function(x, y) {
    dbibinom(N, x, y, p1, p2, rho)
  }))

  expect_equal(total_prob, 1, tolerance = 1e-6)
})

test_that("dbibinom handles independence (rho = 0)", {
  # When rho = 0, should match product of independent binomials
  N <- 10
  y1 <- 5
  y2 <- 6
  p1 <- 0.4
  p2 <- 0.5

  result <- dbibinom(N, y1, y2, p1, p2, rho = 0)
  expected <- dbinom(y1, N, p1) * dbinom(y2, N, p2)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("dbibinom validates input", {
  expect_error(dbibinom(N = 10, y1 = 11, y2 = 5, p1 = 0.5, p2 = 0.5, rho = 0.5))
  expect_error(dbibinom(N = 10, y1 = -1, y2 = 5, p1 = 0.5, p2 = 0.5, rho = 0.5))
  expect_error(dbibinom(N = 10, y1 = 5, y2 = 5, p1 = 0, p2 = 0.5, rho = 0.5))
  expect_error(dbibinom(N = 10, y1 = 5, y2 = 5, p1 = 1, p2 = 0.5, rho = 0.5))
})

test_that("dbibinom checks correlation bounds", {
  # Correlation outside valid bounds should error
  expect_error(dbibinom(N = 10, y1 = 5, y2 = 5, p1 = 0.3, p2 = 0.5, rho = 1.5))
  expect_error(dbibinom(N = 10, y1 = 5, y2 = 5, p1 = 0.3, p2 = 0.5, rho = -1.5))
})
