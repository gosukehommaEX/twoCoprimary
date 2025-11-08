# Test correlation bound functions

# ==============================================================================
# corrbound2Binary tests
# ==============================================================================

test_that("corrbound2Binary returns valid bounds", {
  # Basic test
  result <- corrbound2Binary(p1 = 0.3, p2 = 0.5)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_named(result, c("L_bound", "U_bound"))
  expect_true(result[1] <= result[2])
  expect_true(result[1] >= -1 && result[1] <= 1)
  expect_true(result[2] >= -1 && result[2] <= 1)
})

test_that("corrbound2Binary handles equal probabilities", {
  # When p1 = p2, upper bound should be 1
  result <- corrbound2Binary(p1 = 0.4, p2 = 0.4)
  expect_equal(unname(result[2]), 1, tolerance = 1e-10)
})

test_that("corrbound2Binary handles complementary probabilities", {
  # When p1 + p2 = 1, lower bound should be -1
  result <- corrbound2Binary(p1 = 0.3, p2 = 0.7)
  expect_equal(unname(result[1]), -1, tolerance = 1e-10)
})

test_that("corrbound2Binary validates input", {
  expect_error(corrbound2Binary(p1 = 0, p2 = 0.5))
  expect_error(corrbound2Binary(p1 = 1, p2 = 0.5))
  expect_error(corrbound2Binary(p1 = 0.5, p2 = -0.1))
  expect_error(corrbound2Binary(p1 = 0.5, p2 = 1.5))
})

test_that("corrbound2Binary returns symmetric bounds for symmetric inputs", {
  # For p1 = p2, bounds should be symmetric around 0
  result <- corrbound2Binary(p1 = 0.5, p2 = 0.5)
  expect_equal(unname(result[1]), unname(-result[2]), tolerance = 1e-10)
})

# ==============================================================================
# corrbound2MixedCountContinuous tests
# ==============================================================================

test_that("corrbound2MixedCountContinuous returns valid bounds", {
  # Basic test
  result <- corrbound2MixedCountContinuous(lambda = 1.25, nu = 0.8, mu = 0, sd = 250)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_named(result, c("L_bound", "U_bound"))
  expect_true(result[1] <= result[2])
  expect_true(result[1] >= -1 && result[1] <= 1)
  expect_true(result[2] >= -1 && result[2] <= 1)
})

test_that("corrbound2MixedCountContinuous bounds are symmetric around 0", {
  # For symmetric parameters, bounds should be symmetric
  result <- corrbound2MixedCountContinuous(lambda = 1.5, nu = 1.0, mu = 0, sd = 100)
  expect_equal(abs(unname(result[1])), abs(unname(result[2])), tolerance = 1e-10)
})

test_that("corrbound2MixedCountContinuous validates input", {
  expect_error(corrbound2MixedCountContinuous(lambda = -1, nu = 0.8, mu = 0, sd = 250))
  expect_error(corrbound2MixedCountContinuous(lambda = 1.25, nu = 0, mu = 0, sd = 250))
  expect_error(corrbound2MixedCountContinuous(lambda = 1.25, nu = 0.8, mu = 0, sd = -1))
})

test_that("corrbound2MixedCountContinuous handles Poisson limit", {
  # As nu approaches infinity, should converge to Poisson case
  result_large_nu <- corrbound2MixedCountContinuous(
    lambda = 2.0, nu = 1000, mu = 0, sd = 100
  )

  # Bounds should be close to symmetric for large nu
  expect_equal(abs(unname(result_large_nu[1])), abs(unname(result_large_nu[2])), tolerance = 0.01)
})
