# Test correlation bound functions

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

test_that("corrbound2MixedCountContinuous validates input", {
  expect_error(corrbound2MixedCountContinuous(lambda = -1, nu = 0.8, mu = 0, sd = 250))
  expect_error(corrbound2MixedCountContinuous(lambda = 1.25, nu = 0, mu = 0, sd = 250))
  expect_error(corrbound2MixedCountContinuous(lambda = 1.25, nu = 0.8, mu = 0, sd = -1))
})
