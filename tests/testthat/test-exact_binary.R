# Regression tests for the exact binary machinery: tie handling in the
# unconditional tests, the n_grid argument, the compiled bivariate binomial
# kernel, and the matrix product form of the co-primary power.

# ==============================================================================
# .tie_last: positions of the last member of each tie group
# ==============================================================================

test_that(".tie_last identifies tie groups in a sorted vector", {
  tl <- twoCoprimary:::.tie_last

  expect_equal(tl(numeric(0)), integer(0))
  expect_equal(tl(1), 1L)
  expect_equal(tl(c(3, 2, 1)), c(1L, 2L, 3L))
  expect_equal(tl(c(3, 3, 1)), c(2L, 2L, 3L))
  expect_equal(tl(c(3, 3, 3)), c(3L, 3L, 3L))
  expect_equal(tl(c(5, 4, 4, 4, 2, 1, 1)), c(1L, 4L, 4L, 4L, 5L, 7L, 7L))
})

test_that(".tie_last uses a relative tolerance", {
  tl <- twoCoprimary:::.tie_last

  # Values far below the absolute tolerance of fpCompare must stay distinct
  expect_equal(tl(c(1e-14, 2e-14)), c(1L, 2L))

  # Values agreeing to relative 1e-10 are one group
  expect_equal(tl(c(1, 1 + 1e-14)), c(2L, 2L))
})

# ==============================================================================
# rr1Binary: tie handling and the n_grid argument
# ==============================================================================

test_that("rr1Binary rejection region matches the Exact package", {
  skip_if_not_installed("Exact")
  message("test-exact_binary.R: comparing against Exact ",
          as.character(packageVersion("Exact")))

  n1 <- 10
  n2 <- 10
  methods <- c("Z-pool" = "z-pooled", "Boschloo" = "boschloo")

  for (tst in names(methods)) {
    p_ex <- matrix(1, n1 + 1, n2 + 1)
    for (x1 in 0:n1) {
      for (x2 in 0:n2) {
        tab <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2, byrow = TRUE)
        p_ex[x1 + 1, x2 + 1] <- Exact::exact.test(
          tab, alternative = "greater", method = methods[[tst]],
          np.interval = FALSE, npNumbers = 100,
          ref.pvalue = FALSE, to.plot = FALSE
        )$p.value
      }
    }

    for (a in c(0.01, 0.025, 0.05)) {
      # Outcomes whose p-value sits within the difference between the two
      # nuisance grids are not comparable, so they are excluded
      comparable <- abs(p_ex - a) > 2e-4
      expect_equal(
        rr1Binary(n1, n2, a, tst)[comparable],
        (p_ex < a)[comparable],
        info = paste(tst, "at alpha =", a)
      )
    }
  }
})

test_that("rr1Binary rejection region is stable in n_grid", {
  for (tst in c("Z-pool", "Boschloo")) {
    reference <- rr1Binary(20, 20, 0.025, tst, n_grid = 400)
    for (g in c(25, 100, 200)) {
      expect_equal(rr1Binary(20, 20, 0.025, tst, n_grid = g), reference,
                   info = paste(tst, "at n_grid =", g))
    }
  }
})

test_that("rr1Binary rejects invalid n_grid", {
  expect_error(rr1Binary(10, 10, 0.025, "Z-pool", n_grid = 5), "at least 10")
  expect_error(rr1Binary(10, 10, 0.025, "Z-pool", n_grid = 0), "at least 10")
  expect_error(rr1Binary(10, 10, 0.025, "Z-pool", n_grid = -100), "at least 10")
  expect_error(rr1Binary(10, 10, 0.025, "Z-pool", n_grid = 100.5), "single integer")
  expect_error(rr1Binary(10, 10, 0.025, "Z-pool", n_grid = c(50, 100)), "single integer")
  expect_error(rr1Binary(10, 10, 0.025, "Z-pool", n_grid = NA_real_), "single integer")
  expect_error(rr1Binary(10, 10, 0.025, "Z-pool", n_grid = Inf), "single integer")
})

test_that("the n_grid default reproduces the documented sample sizes", {
  for (tst in c("Z-pool", "Boschloo")) {
    default <- ss2BinaryExact(0.54, 0.54, 0.25, 0.25, 0.3, 0.3,
                              1, 0.025, 0.1, tst)[["N"]]
    explicit <- ss2BinaryExact(0.54, 0.54, 0.25, 0.25, 0.3, 0.3,
                               1, 0.025, 0.1, tst, n_grid = 100)[["N"]]
    expect_equal(default, 142)
    expect_equal(explicit, 142)
  }
})

test_that("n_grid does not affect the tests that never order the outcomes", {
  for (tst in c("Chisq", "Fisher", "Fisher-midP")) {
    expect_equal(rr1Binary(15, 15, 0.025, tst, n_grid = 25),
                 rr1Binary(15, 15, 0.025, tst, n_grid = 400),
                 info = tst)
  }
})

# ==============================================================================
# dbibinom: the compiled kernel
# ==============================================================================

test_that("the compiled kernel agrees with the R reference implementation", {
  g_r <- twoCoprimary:::.dbibinom_g_r
  g_cpp <- twoCoprimary:::dbibinom_g

  for (N in c(1, 2, 5, 20, 40)) {
    for (p1 in c(0.25, 0.6)) {
      for (p2 in c(0.4, 0.75)) {
        bounds <- corrbound2Binary(p1, p2)
        for (fr in c(0, 0.5, 1)) {
          rho <- bounds[1] + fr * (bounds[2] - bounds[1])
          z <- rho * sqrt(p2 * (1 - p2) / (p1 * (1 - p1)))
          gam <- z / (1 - z)
          xi <- p2 + gam * (p2 - p1)

          y1 <- rep(0:N, each = N + 1)
          y2 <- rep(0:N, times = N + 1)

          expect_equal(
            g_cpp(as.integer(N), as.integer(y1), as.integer(y2), xi, gam),
            g_r(N, y1, y2, xi, gam),
            tolerance = 1e-12,
            info = paste("N =", N, "p1 =", p1, "p2 =", p2, "rho =", round(rho, 3))
          )
        }
      }
    }
  }
})

test_that("dbibinom is a probability mass function", {
  for (N in c(5, 20, 40)) {
    for (p1 in c(0.25, 0.54)) {
      for (p2 in c(0.4, 0.7)) {
        bounds <- corrbound2Binary(p1, p2)
        for (fr in c(0, 0.5, 0.95)) {
          rho <- bounds[1] + fr * (bounds[2] - bounds[1])
          pm <- outer(0:N, 0:N, function(x, y) dbibinom(N, x, y, p1, p2, rho))

          expect_true(all(pm >= -1e-12))
          expect_equal(sum(pm), 1, tolerance = 1e-10)
          expect_equal(rowSums(pm), dbinom(0:N, N, p1), tolerance = 1e-10)
          expect_equal(colSums(pm), dbinom(0:N, N, p2), tolerance = 1e-10)
        }
      }
    }
  }
})

test_that("dbibinom handles scalar and boundary arguments", {
  expect_length(dbibinom(10, 0, 0, 0.4, 0.5, 0.2), 1)
  expect_length(dbibinom(10, rep(0, 5), 0:4, 0.4, 0.5, 0.2), 5)
  expect_length(dbibinom(1, c(0, 1), c(1, 0), 0.4, 0.5, 0.2), 2)
  expect_true(all(dbibinom(10, 0:10, 0:10, 0.4, 0.5, 0.2) >= 0))
})

# ==============================================================================
# power2BinaryExact: the matrix product form of the co-primary power
# ==============================================================================

test_that("the co-primary power equals the direct double sum", {
  n1 <- 12
  n2 <- 12
  p11 <- 0.55; p12 <- 0.65; p21 <- 0.30; p22 <- 0.40
  rho1 <- 0.2
  rho2 <- 0.2

  for (tst in c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")) {
    RR <- rr1Binary(n1, n2, 0.025, tst)
    pmass1 <- outer(0:n1, 0:n1, function(x, y) dbibinom(n1, x, y, p11, p12, rho1))
    pmass2 <- outer(0:n2, 0:n2, function(x, y) dbibinom(n2, x, y, p21, p22, rho2))

    # Direct evaluation of the double sum over the rejection region
    A <- row(RR)[RR]
    C <- col(RR)[RR]
    direct <- sum(t(pmass1[A, ])[A, ] * t(pmass2[C, ])[C, ])

    result <- power2BinaryExact(n1, n2, p11, p12, p21, p22,
                                rho1, rho2, 0.025, tst)

    expect_equal(result$powerCoprimary, direct, tolerance = 1e-12, info = tst)
    expect_true(result$powerCoprimary <= min(result$power1, result$power2) + 1e-12)
  }
})
