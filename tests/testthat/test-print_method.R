# Test print method for twoCoprimary objects

# ==============================================================================
# Print method tests for different endpoint types
# ==============================================================================

test_that("print.twoCoprimary works for continuous endpoints", {
  result <- power2Continuous(
    n1 = 50, n2 = 50,
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, alpha = 0.025,
    known_var = TRUE
  )

  # Should print without error
  expect_output(print(result), "n1")
  expect_output(print(result), "n2")
  expect_output(print(result), "powerCoprimary")
})

test_that("print.twoCoprimary works for binary endpoints", {
  result <- power2BinaryApprox(
    n1 = 100, n2 = 50,
    p11 = 0.5, p12 = 0.4,
    p21 = 0.3, p22 = 0.2,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025, Test = "AN"
  )

  expect_output(print(result), "n1")
  expect_output(print(result), "p \\(group 1\\)")
  expect_output(print(result), "Test")
})

test_that("print.twoCoprimary works for sample size results", {
  result <- ss2Continuous(
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, r = 1,
    alpha = 0.025, beta = 0.1,
    known_var = TRUE
  )

  expect_output(print(result), "n1")
  expect_output(print(result), "n2")
  expect_output(print(result), "N")
})

test_that("print.twoCoprimary works for mixed continuous-binary endpoints", {
  result <- power2MixedContinuousBinary(
    n1 = 50, n2 = 50,
    delta = 0.5, sd = 1,
    p1 = 0.6, p2 = 0.4,
    rho = 0.5, alpha = 0.025,
    Test = "AN"
  )

  expect_output(print(result), "n1")
  expect_output(print(result), "powerCont")
  expect_output(print(result), "powerBin")
})

test_that("print.twoCoprimary works for count-continuous endpoints", {
  result <- power2MixedCountContinuous(
    n1 = 100, n2 = 100,
    r1 = 1.0, r2 = 1.25,
    nu = 0.8, t = 1,
    mu1 = -50, mu2 = 0, sd = 250,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025
  )

  expect_output(print(result), "n1")
  expect_output(print(result), "powerCount")
  expect_output(print(result), "powerCont")
})

# ==============================================================================
# Print method behavior tests
# ==============================================================================

test_that("print.twoCoprimary returns invisibly", {
  result <- power2Continuous(
    n1 = 50, n2 = 50,
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, alpha = 0.025,
    known_var = TRUE
  )

  # Print should return the object invisibly
  expect_invisible(print(result))
  returned <- withVisible(print(result))
  expect_false(returned$visible)
  expect_identical(returned$value, result)
})

test_that("print.twoCoprimary formats numeric values appropriately", {
  result <- ss2Continuous(
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, r = 1,
    alpha = 0.025, beta = 0.1,
    known_var = TRUE
  )

  # Capture output
  output <- capture.output(print(result))
  output_text <- paste(output, collapse = "\n")

  # Check that sample sizes are printed as integers
  expect_match(output_text, "n1.*=.*\\d+")
  expect_match(output_text, "n2.*=.*\\d+")
})

# ==============================================================================
# Print method tests for unified interface functions
# ==============================================================================

test_that("print.twoCoprimary works for twoCoprimary2Continuous power calculation", {
  result <- twoCoprimary2Continuous(
    n1 = 100, n2 = 100,
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, alpha = 0.025,
    known_var = TRUE
  )

  expect_output(print(result), "n1")
  expect_output(print(result), "powerCoprimary")
})

test_that("print.twoCoprimary works for twoCoprimary2Continuous sample size calculation", {
  result <- twoCoprimary2Continuous(
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, power = 0.8, r = 1,
    alpha = 0.025, known_var = TRUE
  )

  expect_output(print(result), "n1")
  expect_output(print(result), "n2")
  expect_output(print(result), "N")
})

test_that("print.twoCoprimary works for twoCoprimary2BinaryApprox", {
  result <- twoCoprimary2BinaryApprox(
    n1 = 100, n2 = 100,
    p11 = 0.5, p12 = 0.4,
    p21 = 0.3, p22 = 0.2,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025, Test = "AN"
  )

  expect_output(print(result), "n1")
  expect_output(print(result), "Test")
})

test_that("print.twoCoprimary works for twoCoprimary2MixedContinuousBinary", {
  result <- twoCoprimary2MixedContinuousBinary(
    n1 = 80, n2 = 80,
    delta = 0.5, sd = 1,
    p1 = 0.6, p2 = 0.4,
    rho = 0.5, alpha = 0.025,
    Test = "AN"
  )

  expect_output(print(result), "n1")
  expect_output(print(result), "powerCont")
  expect_output(print(result), "powerBin")
})

test_that("print.twoCoprimary works for twoCoprimary2MixedCountContinuous", {
  result <- twoCoprimary2MixedCountContinuous(
    n1 = 200, n2 = 200,
    r1 = 1.0, r2 = 1.25,
    nu = 0.8, t = 1,
    mu1 = -50, mu2 = 0, sd = 250,
    rho1 = 0.5, rho2 = 0.5,
    alpha = 0.025
  )

  expect_output(print(result), "n1")
  expect_output(print(result), "powerCount")
  expect_output(print(result), "powerCont")
})

# ==============================================================================
# Print method tests for different parameter combinations
# ==============================================================================

test_that("print.twoCoprimary handles single endpoint calculations", {
  result <- ss1Continuous(delta = 0.5, sd = 1, r = 1, alpha = 0.025, beta = 0.1)

  expect_output(print(result), "n1")
  expect_output(print(result), "n2")
  expect_output(print(result), "N")
})

test_that("print.twoCoprimary handles unbalanced designs", {
  result <- power2Continuous(
    n1 = 100, n2 = 50,  # Unbalanced
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = 0.3, alpha = 0.025,
    known_var = TRUE
  )

  expect_output(print(result), "n1.*=.*100")
  expect_output(print(result), "n2.*=.*50")
})
