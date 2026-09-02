# Exhaustive sweep of every exported function over its whole documented domain,
# with the grids made dense at the edges.
#
# Why this file exists. dev/audit_all_functions.R checks named properties on
# chosen arguments, and the reproduce scripts check the package against the
# published tables. Both passed while six defects remained, every one of them
# at an edge of the parameter space that no chosen argument happened to reach:
# a continuity correction carrying a probability out of the unit interval, a
# corrected variance vanishing, a t statistic with no degrees of freedom, a
# dependence parameter diverging at the Frechet-Hoeffding bound, and a matrix
# of null probabilities simplified to a vector when a group has one subject.
# Finding those by reading the sources does not scale and did not work. This
# file replaces reading with enumeration.
#
# What it asserts. For every call: it must not raise an unexpected error, it
# must not warn, and the value must satisfy the contract of its type. For a
# result object that means the class and one row, every numeric column finite,
# every power in the unit interval, the co-primary power between the Frechet
# bounds of its marginals, and, for a sample size, N equal to n1 + n2 with n1
# equal to ceiling(r n2). For a rejection region it means a logical matrix of
# the right size with no NA, monotone in both arguments, and every row a
# prefix, which is the shape the marginal power formulas assume. Calls that
# should be refused are listed separately and must be refused.
#
# Every call runs under an elapsed time limit, so a loop that does not
# terminate is reported as a finding instead of freezing the session.
#
# Run with:  source("dev/fuzz_all_functions.R")
#
# Results go to dev/out/fuzz_all_functions.log and .csv. The grids are the
# named vectors just below; shrink them for a quicker pass.

library(twoCoprimary)

TIME_LIMIT <- 20        # seconds allowed for a single call
N_GRID_FUZZ <- 20       # nuisance grid for the unconditional tests
NMC_FUZZ <- 200         # replications for the Monte Carlo paths

t_start_all <- Sys.time()
out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(out_dir, "fuzz_all_functions.log")
csv_path <- file.path(out_dir, "fuzz_all_functions.csv")
if (file.exists(log_path)) file.remove(log_path)

say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_path, append = TRUE)
}

rows <- list()
n_calls <- 0L
record <- function(fn, case, status, detail = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    fn = fn, case = case, status = status, detail = detail,
    stringsAsFactors = FALSE
  )
}

# Run one call under a time limit, capturing errors, warnings and messages, and
# validate the value. Anything other than a silent, contract-satisfying return
# is recorded.
probe <- function(fn, case, expr, validate = NULL, must_fail = FALSE,
                  limit = TIME_LIMIT) {
  n_calls <<- n_calls + 1L
  warns <- character(0)
  msgs <- character(0)
  setTimeLimit(elapsed = limit)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
  val <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) structure(list(msg = conditionMessage(e)),
                                  class = "fuzz_error")
  )
  setTimeLimit(cpu = Inf, elapsed = Inf)

  if (inherits(val, "fuzz_error")) {
    if (must_fail) {
      record(fn, case, "REFUSED", val$msg)
    } else {
      record(fn, case, "ERROR", gsub("\n", " ", val$msg))
    }
    return(invisible(NULL))
  }
  if (must_fail) {
    record(fn, case, "NOT REFUSED", "an invalid argument returned a value")
    return(invisible(val))
  }
  if (length(warns)) {
    record(fn, case, "WARNING", gsub("\n", " ", warns[1]))
    return(invisible(val))
  }
  bad <- if (is.null(validate)) "" else validate(val)
  if (nzchar(bad)) record(fn, case, "VALUE", bad) else record(fn, case, "OK", "")
  invisible(val)
}

# --------------------------------------------------------------------------
# Validators
# --------------------------------------------------------------------------

# Every result object the package returns
v_obj <- function(x) {
  if (!identical(class(x), c("twoCoprimary", "data.frame"))) {
    return(paste("class is", paste(class(x), collapse = "/")))
  }
  if (nrow(x) != 1L) return(sprintf("nrow is %d", nrow(x)))
  num <- x[, vapply(x, is.numeric, logical(1)), drop = FALSE]
  num <- num[, setdiff(names(num), "nMC"), drop = FALSE]
  nf <- names(num)[!vapply(num, function(v) is.finite(v), logical(1))]
  if (length(nf)) return(paste("non-finite:", paste(nf, collapse = ", ")))
  pw <- names(num)[grepl("^power", names(num))]
  if (length(pw)) {
    v <- unlist(num[, pw, drop = FALSE])
    names(v) <- pw
    if (any(v < 0 | v > 1)) {
      return(paste("power outside the unit interval:",
                   paste(sprintf("%s=%.8g", pw, v), collapse = ", ")))
    }
    if ("powerCoprimary" %in% pw) {
      marg <- v[setdiff(pw, "powerCoprimary")]
      j <- v[["powerCoprimary"]]
      if (length(marg) && any(j > marg + 1e-7)) {
        return(sprintf("joint %.8f exceeds a marginal (%s)", j,
                       paste(sprintf("%.8f", marg), collapse = ", ")))
      }
      if (length(marg) && j < sum(marg) - 1 - 1e-7) {
        return(sprintf("joint %.8f below the Bonferroni bound %.8f", j,
                       sum(marg) - 1))
      }
    }
  }
  if (all(c("n1", "n2", "N") %in% names(x))) {
    if (x$n2 < 1 || x$n2 != round(x$n2)) return("n2 is not a positive integer")
    if (x$n1 < 1 || x$n1 != round(x$n1)) return("n1 is not a positive integer")
    if (x$N != x$n1 + x$n2) return("N is not n1 + n2")
    if ("r" %in% names(x) && x$n1 != ceiling(x$r * x$n2)) {
      return("n1 is not ceiling(r n2)")
    }
  }
  ""
}

# A rejection region
v_rr <- function(n1, n2) {
  force(n1); force(n2)
  function(RR) {
    if (!is.matrix(RR) || !is.logical(RR)) return("not a logical matrix")
    if (!identical(dim(RR), c(as.integer(n1) + 1L, as.integer(n2) + 1L))) {
      return(paste("dimensions", paste(dim(RR), collapse = " by ")))
    }
    if (any(is.na(RR))) return("contains NA")
    if (nrow(RR) > 1L &&
        any(RR[-nrow(RR), , drop = FALSE] & !RR[-1, , drop = FALSE])) {
      return("not monotone in the group 1 count")
    }
    if (ncol(RR) > 1L &&
        any(RR[, -1, drop = FALSE] & !RR[, -ncol(RR), drop = FALSE])) {
      return("not monotone in the group 2 count")
    }
    k <- rowSums(RR)
    for (i in seq_len(nrow(RR))) {
      if (any(RR[i, ] != c(rep(TRUE, k[i]), rep(FALSE, ncol(RR) - k[i])))) {
        return("a row is not an initial run, which the marginal power assumes")
      }
    }
    ""
  }
}

# A pair of correlation bounds
v_bound <- function(b) {
  if (!is.numeric(b) || length(b) != 2L) return("not a numeric pair")
  if (!identical(names(b), c("L_bound", "U_bound"))) {
    return(paste("names are", paste(names(b), collapse = "/")))
  }
  if (!all(is.finite(b))) return("non-finite")
  if (b[[1]] > b[[2]]) return("the lower bound exceeds the upper bound")
  if (b[[1]] < -1 - 1e-12 || b[[2]] > 1 + 1e-12) return("outside [-1, 1]")
  if (b[[1]] > 1e-12 || b[[2]] < -1e-12) return("the bounds exclude zero")
  ""
}

section <- function(title) {
  say("")
  say(strrep("-", 78))
  say(title)
  say(strrep("-", 78))
  invisible(Sys.time())
}
section_done <- function(t0, label) {
  n_now <- length(rows)
  say(sprintf("  %s: %d calls in %.1f s", label, n_now,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

say("twoCoprimary version : ", as.character(utils::packageVersion("twoCoprimary")))
say("R version            : ", R.version.string)
say("run at               : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("time limit per call  : ", TIME_LIMIT, " s")

# ==========================================================================
t0 <- section("A. corrbound2Binary over the whole open unit square")
# ==========================================================================
P_FINE <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999)
for (p1 in P_FINE) {
  for (p2 in P_FINE) {
    probe("corrbound2Binary", sprintf("p = (%g, %g)", p1, p2),
          corrbound2Binary(p1, p2), v_bound)
  }
}
section_done(t0, "corrbound2Binary")

# ==========================================================================
t0 <- section("B. corrbound2MixedCountContinuous over the parameter space")
# ==========================================================================
# The quadrature runs once per point of the negative binomial support, so a
# very small dispersion with a large mean is expensive. The support is measured
# first and the case is recorded rather than run when it is large.
for (lambda in c(0.05, 0.5, 1.25, 3, 10)) {
  for (nu in c(0.05, 0.5, 3, 50)) {
    for (mu in c(-100, 0, 100)) {
      for (sd in c(0.1, 1, 250)) {
        y_max <- stats::qnbinom(0.9999, mu = lambda, size = nu)
        case <- sprintf("lambda = %g, nu = %g, mu = %g, sd = %g",
                        lambda, nu, mu, sd)
        if (y_max > 400) {
          record("corrbound2MixedCountContinuous", case, "SKIPPED",
                 sprintf("support reaches %d, too many quadratures", y_max))
          next
        }
        probe("corrbound2MixedCountContinuous", case,
              corrbound2MixedCountContinuous(lambda, nu, mu, sd), v_bound)
      }
    }
  }
}
section_done(t0, "corrbound2MixedCountContinuous")

# ==========================================================================
t0 <- section("C. dbibinom at the bounds of its dependence parameter")
# ==========================================================================
v_pmf <- function(N, p1, p2, rho) {
  force(N); force(p1); force(p2); force(rho)
  function(P) {
    if (!all(is.finite(P))) return("contains a non-finite probability")
    if (any(P < -1e-13)) return(sprintf("smallest entry %.3e", min(P)))
    if (abs(sum(P) - 1) > 1e-10) return(sprintf("total mass %.14f", sum(P)))
    y <- 0:N
    if (max(abs(rowSums(P) - stats::dbinom(y, N, p1))) > 1e-10) {
      return("the first marginal is not binomial")
    }
    if (max(abs(colSums(P) - stats::dbinom(y, N, p2))) > 1e-10) {
      return("the second marginal is not binomial")
    }
    v <- sqrt(N * p1 * (1 - p1) * N * p2 * (1 - p2))
    got <- (sum(P * outer(y, y)) - sum(rowSums(P) * y) * sum(colSums(P) * y)) / v
    if (abs(got - rho) > 1e-8) {
      return(sprintf("the recovered correlation is %.10f, not %.10f", got, rho))
    }
    ""
  }
}
for (N in c(1, 2, 3, 5, 8, 12)) {
  for (pp in list(c(0.01, 0.01), c(0.01, 0.99), c(0.1, 0.5), c(0.5, 0.5),
                  c(0.5, 0.9), c(0.9, 0.9), c(0.99, 0.5), c(0.3, 0.7))) {
    b <- corrbound2Binary(pp[1], pp[2])
    for (rho in c(b[[1]], b[[1]] / 2, 0, b[[2]] / 2, b[[2]])) {
      probe("dbibinom",
            sprintf("N = %d, p = (%g, %g), rho = %.6f", N, pp[1], pp[2], rho),
            outer(0:N, 0:N, function(a, c2) dbibinom(N, a, c2, pp[1], pp[2], rho)),
            v_pmf(N, pp[1], pp[2], rho))
    }
  }
}
section_done(t0, "dbibinom")

# ==========================================================================
t0 <- section("D. rr1Binary over every test, size and pair of group sizes")
# ==========================================================================
NN_RR <- c(1, 2, 3, 4, 5, 6, 8, 10, 15)
ALPHA_RR <- c(0.001, 0.01, 0.025, 0.05, 0.1, 0.5, 0.9)
RR_TESTS <- c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")
for (n1 in NN_RR) {
  for (n2 in NN_RR) {
    for (a in ALPHA_RR) {
      for (tst in RR_TESTS) {
        probe("rr1Binary",
              sprintf("%-11s n = (%d, %d), alpha = %g", tst, n1, n2, a),
              rr1Binary(n1, n2, a, Test = tst, n_grid = N_GRID_FUZZ),
              v_rr(n1, n2))
      }
    }
  }
}
section_done(t0, "rr1Binary")

# The two exact unconditional tests must hold their size, and the conditional
# test must too. This is quadratic in the grid, so it runs on a subset.
t0 <- section("D2. Null rejection probability against the nominal level")
for (n1 in c(2, 5, 8, 12)) {
  for (n2 in c(2, 5, 8, 12)) {
    for (a in c(0.01, 0.025, 0.05)) {
      for (tst in c("Fisher", "Z-pool", "Boschloo")) {
        RR <- rr1Binary(n1, n2, a, Test = tst, n_grid = 100)
        th <- seq(0.005, 0.995, by = 0.005)
        sz <- vapply(th, function(p) {
          sum(outer(stats::dbinom(0:n1, n1, p), stats::dbinom(0:n2, n2, p)) * RR)
        }, numeric(1))
        case <- sprintf("%-9s n = (%2d, %2d), alpha = %g", tst, n1, n2, a)
        if (tst == "Fisher") {
          # Conditioning makes the level exact with no grid to discretise
          if (max(sz) > a + 1e-10) {
            record("rr1Binary size", case, "VALUE",
                   sprintf("largest size %.6f exceeds alpha", max(sz)))
          } else {
            record("rr1Binary size", case, "OK", "")
          }
        } else {
          # The unconditional tests maximise over the nuisance parameter on a
          # finite grid, so a size marginally above alpha is a statement about
          # the grid rather than about the test
          if (max(sz) > a * 1.02) {
            record("rr1Binary size", case, "VALUE",
                   sprintf("largest size %.6f exceeds alpha by more than 2 per cent",
                           max(sz)))
          } else {
            record("rr1Binary size", case, "OK", "")
          }
        }
        n_calls <- n_calls + 1L
      }
    }
  }
}
section_done(t0, "rr1Binary size")

# ==========================================================================
t0 <- section("E. power2BinaryApprox over every test, size and correlation")
# ==========================================================================
# The probability sets reach both edges of the unit interval and both the
# equal-marginal case, where the upper Prentice bound is one.
PROB_SETS <- list(
  c(0.60, 0.50, 0.40, 0.30), c(0.60, 0.40, 0.40, 0.20),
  c(0.30, 0.25, 0.10, 0.05), c(0.99, 0.95, 0.90, 0.85),
  c(0.55, 0.50, 0.45, 0.40), c(0.90, 0.90, 0.70, 0.70),
  c(0.01, 0.05, 0.005, 0.02), c(0.999, 0.99, 0.99, 0.95)
)
NN_PW <- c(1, 2, 3, 4, 5, 6, 10, 50, 500)
for (pp in PROB_SETS) {
  b1 <- corrbound2Binary(pp[1], pp[2])
  b2 <- corrbound2Binary(pp[3], pp[4])
  lo <- max(b1[[1]], b2[[1]])
  hi <- min(b1[[2]], b2[[2]])
  for (rho in c(lo, 0, hi)) {
    for (tst in c("AN", "ANc", "AS", "ASc")) {
      for (n1 in NN_PW) {
        for (n2 in NN_PW) {
          probe("power2BinaryApprox",
                sprintf("%-4s p = (%g, %g, %g, %g), rho = %.4f, n = (%d, %d)",
                        tst, pp[1], pp[2], pp[3], pp[4], rho, n1, n2),
                power2BinaryApprox(n1 = n1, n2 = n2, p11 = pp[1], p12 = pp[2],
                                   p21 = pp[3], p22 = pp[4], rho1 = rho,
                                   rho2 = rho, alpha = 0.025, Test = tst),
                v_obj)
        }
      }
    }
  }
}
for (a in c(0.001, 0.01, 0.1, 0.5, 0.9)) {
  probe("power2BinaryApprox", sprintf("alpha = %g", a),
        power2BinaryApprox(n1 = 40, n2 = 40, p11 = 0.6, p12 = 0.5, p21 = 0.4,
                           p22 = 0.3, rho1 = 0.3, rho2 = 0.3, alpha = a,
                           Test = "ANc"), v_obj)
}
section_done(t0, "power2BinaryApprox")

# ==========================================================================
t0 <- section("F. power2BinaryExact over every test and small group sizes")
# ==========================================================================
for (pp in list(c(0.60, 0.50, 0.40, 0.30), c(0.90, 0.90, 0.70, 0.70),
                c(0.99, 0.95, 0.90, 0.85), c(0.05, 0.10, 0.01, 0.02))) {
  b1 <- corrbound2Binary(pp[1], pp[2])
  b2 <- corrbound2Binary(pp[3], pp[4])
  hi <- min(b1[[2]], b2[[2]])
  for (rho in c(0, hi / 2, hi)) {
    for (tst in RR_TESTS) {
      for (n1 in c(1, 2, 3, 4, 6, 10)) {
        for (n2 in c(1, 2, 3, 4, 6, 10)) {
          probe("power2BinaryExact",
                sprintf("%-11s p = (%g, %g, %g, %g), rho = %.4f, n = (%d, %d)",
                        tst, pp[1], pp[2], pp[3], pp[4], rho, n1, n2),
                power2BinaryExact(n1 = n1, n2 = n2, p11 = pp[1], p12 = pp[2],
                                  p21 = pp[3], p22 = pp[4], rho1 = rho,
                                  rho2 = rho, alpha = 0.025, Test = tst,
                                  n_grid = N_GRID_FUZZ),
                v_obj)
        }
      }
    }
  }
}
section_done(t0, "power2BinaryExact")

# ==========================================================================
t0 <- section("G. power2Continuous, both variance assumptions")
# ==========================================================================
for (kv in c(TRUE, FALSE)) {
  for (n1 in c(1, 2, 3, 5, 50, 500)) {
    for (n2 in c(1, 2, 3, 5, 50, 500)) {
      for (dd in list(c(0.01, 0.01), c(0.5, 0.4), c(5, 5))) {
        for (ss in list(c(0.01, 0.01), c(1, 1.2), c(100, 100))) {
          for (rho in c(-0.999, -0.5, 0, 0.5, 0.999)) {
            probe("power2Continuous",
                  sprintf("known_var = %s, n = (%d, %d), delta = (%g, %g), sd = (%g, %g), rho = %g",
                          kv, n1, n2, dd[1], dd[2], ss[1], ss[2], rho),
                  power2Continuous(n1 = n1, n2 = n2, delta1 = dd[1],
                                   delta2 = dd[2], sd1 = ss[1], sd2 = ss[2],
                                   rho = rho, alpha = 0.025, known_var = kv,
                                   nMC = NMC_FUZZ),
                  v_obj)
          }
        }
      }
    }
  }
}
section_done(t0, "power2Continuous")

# ==========================================================================
t0 <- section("H. power2MixedContinuousBinary over every test")
# ==========================================================================
for (tst in c("AN", "ANc", "AS", "ASc", "Fisher")) {
  nn <- if (tst == "Fisher") c(1, 2, 3, 30) else c(1, 2, 3, 4, 5, 6, 30, 300)
  for (n1 in nn) {
    for (n2 in nn) {
      for (pp in list(c(0.60, 0.40), c(0.30, 0.10), c(0.99, 0.85),
                      c(0.95, 0.80), c(0.02, 0.01))) {
        for (rho in c(-0.99, 0, 0.99)) {
          probe("power2MixedContinuousBinary",
                sprintf("%-7s n = (%d, %d), p = (%g, %g), rho = %g",
                        tst, n1, n2, pp[1], pp[2], rho),
                power2MixedContinuousBinary(n1 = n1, n2 = n2, delta = 0.5,
                                            sd = 1, p1 = pp[1], p2 = pp[2],
                                            rho = rho, alpha = 0.025,
                                            Test = tst, nMC = NMC_FUZZ),
                v_obj)
        }
      }
    }
  }
}
for (ds in list(c(0.001, 1), c(1, 0.001), c(1000, 1), c(1, 1000))) {
  probe("power2MixedContinuousBinary",
        sprintf("delta = %g, sd = %g", ds[1], ds[2]),
        power2MixedContinuousBinary(n1 = 40, n2 = 40, delta = ds[1],
                                    sd = ds[2], p1 = 0.6, p2 = 0.4, rho = 0.5,
                                    alpha = 0.025, Test = "ANc"), v_obj)
}
section_done(t0, "power2MixedContinuousBinary")

# ==========================================================================
t0 <- section("I. power2MixedCountContinuous over the parameter space")
# ==========================================================================
# Each call validates both correlation bounds, and each bound is a quadrature
# over the support of the negative binomial distribution, so the support is
# measured first and an expensive combination is recorded rather than run.
for (rr in list(c(1.0, 1.25), c(1.0, 2.0), c(0.05, 0.06))) {
  for (nu in c(0.2, 0.8, 5)) {
    for (mus in list(c(-50, 0), c(0, 0), c(-0.001, 0))) {
      for (sd in c(0.5, 250)) {
        case0 <- sprintf("rates = (%g, %g), nu = %g, mu = (%g, %g), sd = %g",
                         rr[1], rr[2], nu, mus[1], mus[2], sd)
        y_max <- max(stats::qnbinom(0.9999, mu = rr, size = nu))
        if (y_max > 400) {
          record("power2MixedCountContinuous", case0, "SKIPPED",
                 sprintf("support reaches %d, too many quadratures", y_max))
          next
        }
        b1 <- corrbound2MixedCountContinuous(rr[1], nu, mus[1], sd)
        b2 <- corrbound2MixedCountContinuous(rr[2], nu, mus[2], sd)
        lo <- max(b1[[1]], b2[[1]])
        hi <- min(b1[[2]], b2[[2]])
        for (rho in c(lo, 0, hi)) {
          for (nn in list(c(1, 1), c(3, 2), c(300, 300))) {
            probe("power2MixedCountContinuous",
                  sprintf("%s, rho = %.4f, n = (%d, %d)", case0, rho, nn[1], nn[2]),
                  power2MixedCountContinuous(n1 = nn[1], n2 = nn[2], r1 = rr[1],
                                             r2 = rr[2], nu = nu, t = 1,
                                             mu1 = mus[1], mu2 = mus[2],
                                             sd = sd, rho1 = rho, rho2 = rho,
                                             alpha = 0.025),
                  v_obj)
          }
        }
      }
    }
  }
}
# The follow-up time enters only through lambda = r t, so it is swept separately
for (tt in c(0.1, 1, 5)) {
  probe("power2MixedCountContinuous", sprintf("t = %g", tt),
        power2MixedCountContinuous(n1 = 300, n2 = 300, r1 = 1.0, r2 = 1.25,
                                   nu = 0.8, t = tt, mu1 = -50, mu2 = 0,
                                   sd = 250, rho1 = 0, rho2 = 0, alpha = 0.025),
        v_obj)
}
section_done(t0, "power2MixedCountContinuous")

# ==========================================================================
t0 <- section("J. The single endpoint sample size functions")
# ==========================================================================
R_GRID <- c(0.25, 0.5, 1, 1.5, 2, 3, 5)
BETA_GRID <- c(0.01, 0.1, 0.2, 0.5, 0.9)
for (r in R_GRID) {
  for (beta in BETA_GRID) {
    for (a in c(0.001, 0.025, 0.1, 0.4)) {
      probe("ss1Continuous", sprintf("r = %g, alpha = %g, beta = %g", r, a, beta),
            ss1Continuous(delta = 0.5, sd = 1, r = r, alpha = a, beta = beta),
            v_obj)
      probe("ss1Count", sprintf("r = %g, alpha = %g, beta = %g", r, a, beta),
            ss1Count(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1, r = r, alpha = a,
                     beta = beta), v_obj)
      for (tst in c("AN", "ANc", "AS", "ASc")) {
        for (pp in list(c(0.6, 0.4), c(0.99, 0.95), c(0.95, 0.05),
                        c(0.05, 0.01), c(0.51, 0.50))) {
          probe("ss1BinaryApprox",
                sprintf("%-4s p = (%g, %g), r = %g, alpha = %g, beta = %g",
                        tst, pp[1], pp[2], r, a, beta),
                ss1BinaryApprox(p1 = pp[1], p2 = pp[2], r = r, alpha = a,
                                beta = beta, Test = tst), v_obj)
        }
      }
    }
  }
}
for (r in c(0.5, 1, 2)) {
  for (pp in list(c(0.9, 0.3), c(0.75, 0.35))) {
    probe("ss1BinaryApprox",
          sprintf("Fisher p = (%g, %g), r = %g", pp[1], pp[2], r),
          ss1BinaryApprox(p1 = pp[1], p2 = pp[2], r = r, alpha = 0.025,
                          beta = 0.2, Test = "Fisher"), v_obj)
  }
}
section_done(t0, "single endpoint sample sizes")

# ==========================================================================
t0 <- section("K. The two endpoint sample size functions")
# ==========================================================================
for (r in c(0.5, 1, 1.5, 2, 3)) {
  for (beta in c(0.05, 0.2, 0.5)) {
    probe("ss2Continuous", sprintf("r = %g, beta = %g, known variance", r, beta),
          ss2Continuous(delta1 = 0.5, delta2 = 0.45, sd1 = 1, sd2 = 1.1,
                        rho = 0.5, r = r, alpha = 0.025, beta = beta,
                        known_var = TRUE), v_obj)
    probe("ss2Continuous", sprintf("r = %g, beta = %g, unknown variance", r, beta),
          ss2Continuous(delta1 = 0.5, delta2 = 0.45, sd1 = 1, sd2 = 1.1,
                        rho = 0.5, r = r, alpha = 0.025, beta = beta,
                        known_var = FALSE, nMC = NMC_FUZZ), v_obj)
    for (tst in c("AN", "ANc", "AS", "ASc")) {
      for (pp in list(c(0.60, 0.50, 0.40, 0.30), c(0.99, 0.95, 0.90, 0.85),
                      c(0.30, 0.25, 0.10, 0.05))) {
        probe("ss2BinaryApprox",
              sprintf("%-4s p = (%g, %g, %g, %g), r = %g, beta = %g",
                      tst, pp[1], pp[2], pp[3], pp[4], r, beta),
              ss2BinaryApprox(p11 = pp[1], p12 = pp[2], p21 = pp[3],
                              p22 = pp[4], rho1 = 0.3, rho2 = 0.3, r = r,
                              alpha = 0.025, beta = beta, Test = tst), v_obj)
      }
      probe("ss2MixedContinuousBinary",
            sprintf("%-4s r = %g, beta = %g", tst, r, beta),
            ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.6, p2 = 0.4,
                                     rho = 0.5, r = r, alpha = 0.025,
                                     beta = beta, Test = tst), v_obj)
    }
    probe("ss2MixedCountContinuous", sprintf("r = %g, beta = %g", r, beta),
          ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                                  mu1 = -50, mu2 = 0, sd = 250, rho1 = 0.5,
                                  rho2 = 0.5, r = r, alpha = 0.025,
                                  beta = beta), v_obj)
    for (tst in RR_TESTS) {
      probe("ss2BinaryExact", sprintf("%-11s r = %g, beta = %g", tst, r, beta),
            ss2BinaryExact(p11 = 0.80, p12 = 0.70, p21 = 0.30, p22 = 0.20,
                           rho1 = 0.3, rho2 = 0.3, r = r, alpha = 0.025,
                           beta = beta, Test = tst, n_grid = N_GRID_FUZZ),
            v_obj)
    }
  }
}
probe("ss2MixedContinuousBinary", "Fisher, r = 1, beta = 0.2",
      ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.6, p2 = 0.4,
                               rho = 0.5, r = 1, alpha = 0.025, beta = 0.2,
                               Test = "Fisher", nMC = NMC_FUZZ), v_obj)
section_done(t0, "two endpoint sample sizes")

# ==========================================================================
t0 <- section("L. The unified interfaces, design_table, plot and print")
# ==========================================================================
shapes <- list(
  ss_cont = ss2Continuous(0.5, 0.5, 1, 1, 0.5, 1, 0.025, 0.2, TRUE),
  pw_cont = power2Continuous(60, 60, 0.5, 0.5, 1, 1, 0.5, 0.025, TRUE),
  ss_binap = ss2BinaryApprox(0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 1, 0.025, 0.2, "AN"),
  pw_binap = power2BinaryApprox(80, 80, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "AN"),
  ss_binex = ss2BinaryExact(0.8, 0.7, 0.3, 0.2, 0.3, 0.3, 1, 0.025, 0.2, "Chisq"),
  pw_binex = power2BinaryExact(20, 20, 0.8, 0.7, 0.3, 0.2, 0.3, 0.3, 0.025, "Fisher"),
  ss_mcb = ss2MixedContinuousBinary(0.5, 1, 0.6, 0.4, 0.5, 1, 0.025, 0.2, "AN"),
  pw_mcb = power2MixedContinuousBinary(80, 80, 0.5, 1, 0.6, 0.4, 0.5, 0.025, "AN"),
  ss_mcc = ss2MixedCountContinuous(1, 1.25, 0.8, 1, -50, 0, 250, 1, 0.5, 0.5, 0.025, 0.2),
  pw_mcc = power2MixedCountContinuous(300, 300, 1, 1.25, 0.8, 1, -50, 0, 250, 0.5, 0.5, 0.025),
  ss1_cont = ss1Continuous(0.5, 1, 1, 0.025, 0.2),
  ss1_count = ss1Count(1, 1.25, 0.8, 1, 1, 0.025, 0.2),
  ss1_bin = ss1BinaryApprox(0.6, 0.4, 1, 0.025, 0.2, "AN")
)

is_single <- function(o) {
  !(all(c("delta1", "delta2", "sd1", "sd2") %in% names(o)) ||
      all(c("p11", "p12", "p21", "p22") %in% names(o)) ||
      all(c("delta", "sd", "p1", "p2") %in% names(o)) ||
      all(c("r1", "r2", "nu", "mu1", "mu2") %in% names(o)))
}
is_cont <- function(o) all(c("delta1", "delta2", "sd1", "sd2") %in% names(o))

v_plot <- function(type) {
  force(type)
  want <- switch(type, power_curve = c("n1", "n2", "power"),
                 sample_size_rho = c("rho", "n2"),
                 effect_contour = c("delta1", "delta2", "power"))
  function(d) {
    if (!is.data.frame(d)) return("not a data frame")
    if (!identical(names(d), want)) {
      return(paste("columns", paste(names(d), collapse = ",")))
    }
    if (nrow(d) < 1L) return("no rows")
    if (!all(vapply(d, function(v) all(is.finite(v)), logical(1)))) {
      return("a returned column is not finite")
    }
    ""
  }
}

grDevices::pdf(NULL)
for (nm in names(shapes)) {
  o <- shapes[[nm]]
  for (type in c("power_curve", "sample_size_rho", "effect_contour")) {
    for (np in c(2, 3, 5)) {
      supported <- !is_single(o) && (type != "effect_contour" || is_cont(o))
      probe("plot", sprintf("%-9s %-16s n_points = %d", nm, type, np),
            plot(o, type = type, n_points = np), v_plot(type),
            must_fail = !supported)
    }
  }
  probe("print", nm, {
    txt <- utils::capture.output(back <- print(o))
    if (!identical(back, o)) stop("print did not return its argument unchanged")
    if (length(txt) < 4L) stop("the printout is shorter than four lines")
    o
  }, function(x) "")
}
grDevices::dev.off()

# A single plotted point is an admissible request and must not fail silently
grDevices::pdf(NULL)
probe("plot", "one point on a power curve",
      plot(shapes$ss_cont, type = "power_curve", n_points = 1),
      v_plot("power_curve"))
probe("plot", "one point on a sample size curve",
      plot(shapes$ss_cont, type = "sample_size_rho", n_points = 1),
      v_plot("sample_size_rho"))
grDevices::dev.off()

# The unified interfaces must return exactly what they dispatch to
unified <- list(
  list(id = "continuous power",
       a = quote(twoCoprimary2Continuous(n1 = 60, n2 = 60, delta1 = 0.5,
                                         delta2 = 0.5, sd1 = 1, sd2 = 1,
                                         rho = 0.5, alpha = 0.025)),
       b = quote(power2Continuous(60, 60, 0.5, 0.5, 1, 1, 0.5, 0.025))),
  list(id = "continuous size",
       a = quote(twoCoprimary2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1,
                                         sd2 = 1, rho = 0.5, power = 0.8,
                                         r = 2, alpha = 0.025)),
       b = quote(ss2Continuous(0.5, 0.5, 1, 1, 0.5, 2, 0.025, 0.2))),
  list(id = "binary approx power",
       a = quote(twoCoprimary2BinaryApprox(n1 = 80, n2 = 80, p11 = 0.6,
                                           p12 = 0.5, p21 = 0.4, p22 = 0.3,
                                           rho1 = 0.3, rho2 = 0.3,
                                           alpha = 0.025, Test = "ASc")),
       b = quote(power2BinaryApprox(80, 80, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3,
                                    0.025, "ASc"))),
  list(id = "binary approx size",
       a = quote(twoCoprimary2BinaryApprox(p11 = 0.6, p12 = 0.5, p21 = 0.4,
                                           p22 = 0.3, rho1 = 0.3, rho2 = 0.3,
                                           power = 0.8, r = 1, alpha = 0.025,
                                           Test = "AS")),
       b = quote(ss2BinaryApprox(0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 1, 0.025, 0.2,
                                 "AS"))),
  list(id = "binary exact power",
       a = quote(twoCoprimary2BinaryExact(n1 = 20, n2 = 20, p11 = 0.8,
                                          p12 = 0.7, p21 = 0.3, p22 = 0.2,
                                          rho1 = 0.3, rho2 = 0.3,
                                          alpha = 0.025, Test = "Fisher")),
       b = quote(power2BinaryExact(20, 20, 0.8, 0.7, 0.3, 0.2, 0.3, 0.3,
                                   0.025, "Fisher"))),
  list(id = "binary exact size",
       a = quote(twoCoprimary2BinaryExact(p11 = 0.8, p12 = 0.7, p21 = 0.3,
                                          p22 = 0.2, rho1 = 0.3, rho2 = 0.3,
                                          power = 0.8, r = 1, alpha = 0.025,
                                          Test = "Chisq")),
       b = quote(ss2BinaryExact(0.8, 0.7, 0.3, 0.2, 0.3, 0.3, 1, 0.025, 0.2,
                                "Chisq"))),
  list(id = "mixed continuous binary power",
       a = quote(twoCoprimary2MixedContinuousBinary(n1 = 80, n2 = 80,
                                                    delta = 0.5, sd = 1,
                                                    p1 = 0.6, p2 = 0.4,
                                                    rho = 0.5, alpha = 0.025,
                                                    Test = "ANc")),
       b = quote(power2MixedContinuousBinary(80, 80, 0.5, 1, 0.6, 0.4, 0.5,
                                             0.025, "ANc"))),
  list(id = "mixed continuous binary size",
       a = quote(twoCoprimary2MixedContinuousBinary(delta = 0.5, sd = 1,
                                                    p1 = 0.6, p2 = 0.4,
                                                    rho = 0.5, power = 0.9,
                                                    r = 1, alpha = 0.025,
                                                    Test = "AN")),
       b = quote(ss2MixedContinuousBinary(0.5, 1, 0.6, 0.4, 0.5, 1, 0.025,
                                          0.1, "AN"))),
  list(id = "mixed count continuous power",
       a = quote(twoCoprimary2MixedCountContinuous(n1 = 300, n2 = 300, r1 = 1,
                                                   r2 = 1.25, nu = 0.8, t = 1,
                                                   mu1 = -50, mu2 = 0,
                                                   sd = 250, rho1 = 0.5,
                                                   rho2 = 0.5, alpha = 0.025)),
       b = quote(power2MixedCountContinuous(300, 300, 1, 1.25, 0.8, 1, -50, 0,
                                            250, 0.5, 0.5, 0.025))),
  list(id = "mixed count continuous size",
       a = quote(twoCoprimary2MixedCountContinuous(r1 = 1, r2 = 1.25, nu = 0.8,
                                                   t = 1, mu1 = -50, mu2 = 0,
                                                   sd = 250, rho1 = 0.5,
                                                   rho2 = 0.5, power = 0.8,
                                                   r = 1, alpha = 0.025)),
       b = quote(ss2MixedCountContinuous(1, 1.25, 0.8, 1, -50, 0, 250, 1, 0.5,
                                         0.5, 0.025, 0.2)))
)
for (u in unified) {
  probe("twoCoprimary2*", u$id, {
    a <- eval(u$a, envir = globalenv())
    b <- eval(u$b, envir = globalenv())
    if (!isTRUE(all.equal(as.data.frame(a), as.data.frame(b),
                          check.attributes = FALSE))) {
      stop("the unified interface does not return what it dispatches to")
    }
    a
  }, v_obj)
}

# design_table over all four endpoint types, both modes, and a correlation
# outside the bounds, which must give NA rather than an error
dt <- list(
  list(id = "continuous size",
       e = quote(design_table(data.frame(delta1 = c(0.4, 0.5), delta2 = 0.4,
                                         sd1 = 1, sd2 = 1),
                              rho_values = c(0, 0.5), r = 1, alpha = 0.025,
                              beta = 0.2, endpoint_type = "continuous"))),
  list(id = "continuous power",
       e = quote(design_table(data.frame(n1 = 80, n2 = 80, delta1 = 0.5,
                                         delta2 = 0.4, sd1 = 1, sd2 = 1),
                              rho_values = c(0, 0.5), alpha = 0.025,
                              endpoint_type = "continuous"))),
  list(id = "binary approx size",
       e = quote(design_table(data.frame(p11 = 0.6, p12 = 0.5, p21 = 0.4,
                                         p22 = 0.3),
                              rho_values = c(0, 0.5), r = 1, alpha = 0.025,
                              beta = 0.2, endpoint_type = "binary",
                              Test = "ASc"))),
  list(id = "binary exact power",
       e = quote(design_table(data.frame(n1 = 20, n2 = 20, p11 = 0.8,
                                         p12 = 0.7, p21 = 0.3, p22 = 0.2),
                              rho_values = c(0, 0.5), alpha = 0.025,
                              endpoint_type = "binary", Test = "Fisher"))),
  list(id = "mixed continuous binary size",
       e = quote(design_table(data.frame(delta = 0.5, sd = 1, p1 = 0.6,
                                         p2 = 0.4),
                              rho_values = c(0, 0.5), r = 1, alpha = 0.025,
                              beta = 0.2, endpoint_type = "mixed_cont_binary",
                              Test = "AN"))),
  list(id = "mixed count continuous size",
       e = quote(design_table(data.frame(r1 = 1, r2 = 1.25, nu = 0.8, t = 1,
                                         mu1 = -50, mu2 = 0, sd = 250),
                              rho_values = c(0, 0.5), r = 1, alpha = 0.025,
                              beta = 0.2, endpoint_type = "mixed_count_cont"))),
  list(id = "correlation outside the bounds",
       e = quote(design_table(data.frame(p11 = 0.9, p12 = 0.2, p21 = 0.5,
                                         p22 = 0.1),
                              rho_values = c(0.1, 0.9), r = 1, alpha = 0.025,
                              beta = 0.2, endpoint_type = "binary",
                              Test = "AN"))),
  list(id = "a grid of one row and one correlation",
       e = quote(design_table(data.frame(delta1 = 0.5, delta2 = 0.5, sd1 = 1,
                                         sd2 = 1),
                              rho_values = 0, r = 1, alpha = 0.025, beta = 0.2,
                              endpoint_type = "continuous")))
)
v_dt <- function(x) {
  if (!identical(class(x), c("twoCoprimary_table", "data.frame"))) {
    return(paste("class is", paste(class(x), collapse = "/")))
  }
  if (!any(grepl("^rho_", names(x)))) return("no correlation column")
  ""
}
for (d in dt) {
  probe("design_table", d$id, eval(d$e, envir = globalenv()), v_dt)
}
section_done(t0, "interfaces, tables and methods")

# ==========================================================================
t0 <- section("M. Arguments that must be refused")
# ==========================================================================
bad <- list(
  quote(ss1Continuous(-0.5, 1, 1, 0.025, 0.2)),
  quote(ss1Continuous(0.5, 0, 1, 0.025, 0.2)),
  quote(ss1Continuous(0.5, 1, 0, 0.025, 0.2)),
  quote(ss1Continuous(0.5, 1, 1, 1, 0.2)),
  quote(ss1Continuous(0.5, 1, 1, 0.025, 0)),
  quote(ss1Continuous(c(0.4, 0.5), 1, 1, 0.025, 0.2)),
  quote(ss1Count(0, 1.25, 0.8, 1, 1, 0.025, 0.2)),
  quote(ss1Count(1, 1.25, 0, 1, 1, 0.025, 0.2)),
  quote(ss1Count(1, 1.25, 0.8, 0, 1, 0.025, 0.2)),
  quote(ss1Count(1.25, 1.0, 0.8, 1, 1, 0.025, 0.2)),
  quote(ss1Count(1.0, 1.0, 0.8, 1, 1, 0.025, 0.2)),
  quote(ss1BinaryApprox(1, 0.4, 1, 0.025, 0.2, "AN")),
  quote(ss1BinaryApprox(0.4, 0.6, 1, 0.025, 0.2, "AN")),
  quote(ss1BinaryApprox(0.6, 0.4, 1, 0.025, 0.2, "Boschloo")),
  quote(rr1Binary(0, 10, 0.025, "Fisher")),
  quote(rr1Binary(10.5, 10, 0.025, "Fisher")),
  quote(rr1Binary(10, 10, 0, "Fisher")),
  quote(rr1Binary(10, 10, 0.025, "AN")),
  quote(rr1Binary(10, 10, 0.025, "Boschloo", n_grid = 5)),
  quote(rr1Binary(10, 10, 0.025, "Boschloo", n_grid = 12.5)),
  quote(rr1Binary(10, 10, 0.025, "Boschloo", n_grid = NA)),
  quote(rr1Binary(10, 10, 0.025, "Boschloo", n_grid = Inf)),
  quote(corrbound2Binary(0, 0.5)),
  quote(corrbound2Binary(0.5, 1)),
  quote(corrbound2MixedCountContinuous(0, 0.8, 0, 250)),
  quote(corrbound2MixedCountContinuous(1.25, 0, 0, 250)),
  quote(corrbound2MixedCountContinuous(1.25, 0.8, 0, 0)),
  quote(dbibinom(0, 0, 0, 0.3, 0.5, 0.2)),
  quote(dbibinom(10, c(1, 2), 3, 0.3, 0.5, 0.2)),
  quote(dbibinom(10, 11, 3, 0.3, 0.5, 0.2)),
  quote(dbibinom(10, 3, 5, 0.3, 0.5, 0.99)),
  quote(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 2, 0.025)),
  quote(power2Continuous(100, 100, 0.5, 0.5, -1, 1, 0.5, 0.025)),
  quote(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 0.5, 2)),
  quote(power2Continuous(100, 0, 0.5, 0.5, 1, 1, 0.5, 0.025)),
  quote(power2Continuous(100, 2.5, 0.5, 0.5, 1, 1, 0.5, 0.025)),
  quote(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 0.5, 0.025, known_var = "yes")),
  quote(power2BinaryApprox(0, 100, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "AN")),
  quote(power2BinaryApprox(100, 2.5, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "AN")),
  quote(power2BinaryApprox(100, 100, 1, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "AN")),
  quote(power2BinaryApprox(100, 100, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 1, "AN")),
  quote(power2BinaryApprox(100, 100, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "Boschloo")),
  quote(power2BinaryApprox(100, 100, 0.6, 0.5, 0.4, 0.3, 0.99, 0.3, 0.025, "AN")),
  quote(power2BinaryExact(20, 20, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "AN")),
  quote(power2BinaryExact(0, 20, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "Fisher")),
  quote(power2MixedContinuousBinary(100, 100, 0.5, 1, 0.6, 0.4, 0.5, 0.025, "Boschloo")),
  quote(power2MixedContinuousBinary(100, 100, 0.5, 1, 0.6, 0.4, 1, 0.025, "AN")),
  quote(power2MixedContinuousBinary(100, 100, 0, 1, 0.6, 0.4, 0.5, 0.025, "AN")),
  quote(power2MixedCountContinuous(300, 300, 1, 1.25, 0, 1, -50, 0, 250, 0.5, 0.5, 0.025)),
  quote(power2MixedCountContinuous(300, 300, 1, 1.25, 0.8, 1, -50, 0, 250, 0.95, 0.5, 0.025)),
  quote(ss2Continuous(0.5, 0.5, 1, 1, 1, 1, 0.025, 0.2)),
  quote(ss2Continuous(0.5, 0.5, 1, 1, 0.5, 1, 0.025, 0.2, known_var = "yes")),
  quote(ss2BinaryApprox(0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 1, 0.025, 0.2, "Fisher")),
  quote(ss2BinaryExact(0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 1, 0.025, 0.2, "AN")),
  quote(ss2MixedContinuousBinary(0, 1, 0.6, 0.4, 0.5, 1, 0.025, 0.2, "AN")),
  quote(ss2MixedCountContinuous(1.25, 1.0, 0.8, 1, -50, 0, 250, 1, 0.5, 0.5, 0.025, 0.2)),
  quote(ss2MixedCountContinuous(1.0, 1.25, 0.8, 1, 50, 0, 250, 1, 0.5, 0.5, 0.025, 0.2)),
  quote(design_table(list(delta1 = 0.5), endpoint_type = "continuous")),
  quote(design_table(data.frame(delta1 = 0.5), endpoint_type = "continuous")),
  quote(design_table(data.frame(delta1 = 0.5, delta2 = 0.4, sd1 = 1, sd2 = 1),
                     endpoint_type = "ordinal")),
  quote(twoCoprimary2Continuous(n1 = 100, n2 = 100, delta1 = 0.5, delta2 = 0.4,
                                sd1 = 1, sd2 = 1, rho = 0.3, power = 0.8, r = 1)),
  quote(twoCoprimary2Continuous(n1 = 100, delta1 = 0.5, delta2 = 0.4, sd1 = 1,
                                sd2 = 1, rho = 0.3)),
  quote(twoCoprimary2Continuous(delta1 = 0.5, delta2 = 0.4, sd1 = 1, sd2 = 1,
                                rho = 0.3, power = 0.8)),
  quote(twoCoprimary2Continuous(delta1 = 0.5, delta2 = 0.4, sd1 = 1, sd2 = 1,
                                rho = 0.3))
)
for (e in bad) {
  probe("refused", paste(deparse(e), collapse = " "),
        eval(e, envir = globalenv()), must_fail = TRUE)
}
section_done(t0, "arguments that must be refused")

# ==========================================================================
section("Summary")
# ==========================================================================
res <- do.call(rbind, rows)
write.csv(res, csv_path, row.names = FALSE)

tab <- table(res$status)
say("")
say(sprintf("calls made : %d", nrow(res)))
for (nm in names(tab)) say(sprintf("  %-12s %d", nm, tab[[nm]]))
say("")

problems <- res[!res$status %in% c("OK", "REFUSED", "SKIPPED"), , drop = FALSE]
if (nrow(problems)) {
  say("--- every call that did not behave ---")
  for (i in seq_len(nrow(problems))) {
    say(sprintf("  [%s] %s | %s", problems$status[i], problems$fn[i],
                problems$case[i]))
    say("      ", problems$detail[i])
  }
} else {
  say("No call raised an unexpected error, warned, or returned a value")
  say("outside its contract, and every argument listed as invalid was refused.")
}

skipped <- res[res$status == "SKIPPED", , drop = FALSE]
if (nrow(skipped)) {
  say("")
  say("--- not attempted ---")
  for (i in seq_len(nrow(skipped))) {
    say(sprintf("  %s | %s : %s", skipped$fn[i], skipped$case[i],
                skipped$detail[i]))
  }
}

say("")
say(sprintf("total elapsed : %.1f s",
            as.numeric(difftime(Sys.time(), t_start_all, units = "secs"))))
cat("\nDone. See", log_path, "and", csv_path, "\n")
