# Independent verification of the two continuous endpoint family.
#
# Run with:  source("dev/verify_two_continuous_family.R")
#
# The article validates this endpoint type against Sozu et al. (2011) Table 1,
# which uses a balanced design throughout. Unequal allocation is therefore not
# covered by any published table, and unequal standard deviations are not
# covered either. This script checks both, in the known variance case and the
# unknown variance case, against two independent references: a closed form
# recomputed here with mvtnorm rather than pbivnorm, and a direct patient level
# simulation of the trial.
#
# Everything is written to dev/out/. Console output need not be pasted back.

library(twoCoprimary)
library(mvtnorm)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(out_dir, "continuous.log")
if (file.exists(log_path)) file.remove(log_path)

say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_path, append = TRUE)
}

pass <- 0L; fail <- 0L
check <- function(id, ok, detail = "") {
  if (isTRUE(ok)) pass <<- pass + 1L else fail <<- fail + 1L
  say(sprintf("[%s] %-34s %s", if (isTRUE(ok)) "PASS" else "FAIL", id, detail))
}

say("twoCoprimary version : ", as.character(utils::packageVersion("twoCoprimary")))
say("run at               : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

set.seed(20260830)

# ---------------------------------------------------------------------------
# The design grid. Deliberately asymmetric: unequal standardized effects,
# unequal standard deviations, and allocation ratios on both sides of 1.
# ---------------------------------------------------------------------------
grid <- expand.grid(
  r = c(0.5, 1, 2, 3),
  rho = c(-0.3, 0, 0.5, 0.8),
  case = 1:3,
  KEEP.OUT.ATTRS = FALSE
)
cases <- list(
  list(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,   n2 = 60),  # symmetric
  list(delta1 = 4.0, delta2 = 1.5, sd1 = 8, sd2 = 5,   n2 = 60),  # unequal effect and sd
  list(delta1 = 12,  delta2 = 30,  sd1 = 20, sd2 = 75, n2 = 40)   # different scales
)
grid$delta1 <- vapply(grid$case, function(i) cases[[i]]$delta1, 0)
grid$delta2 <- vapply(grid$case, function(i) cases[[i]]$delta2, 0)
grid$sd1    <- vapply(grid$case, function(i) cases[[i]]$sd1, 0)
grid$sd2    <- vapply(grid$case, function(i) cases[[i]]$sd2, 0)
grid$n2     <- vapply(grid$case, function(i) cases[[i]]$n2, 0)
grid$n1     <- ceiling(grid$r * grid$n2)

# ---------------------------------------------------------------------------
# Part A: known variance against a closed form recomputed from first
# principles with mvtnorm::pmvnorm. This is an exact check with no Monte Carlo
# error, and uses a different bivariate normal routine from the package.
# ---------------------------------------------------------------------------
say("--- Part A: known variance against an independent closed form ---")
ref_known <- function(n1, n2, d1, d2, s1, s2, rho, alpha) {
  se <- sqrt(1 / n1 + 1 / n2)
  w <- c(d1 / (s1 * se), d2 / (s2 * se))          # non-centrality
  z <- stats::qnorm(1 - alpha)
  R <- matrix(c(1, rho, rho, 1), 2)
  as.numeric(mvtnorm::pmvnorm(lower = c(z, z) - w, upper = c(Inf, Inf), corr = R))
}
grid$pkg_known <- NA_real_; grid$ref_known <- NA_real_
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  grid$pkg_known[i] <- power2Continuous(
    n1 = g$n1, n2 = g$n2, delta1 = g$delta1, delta2 = g$delta2,
    sd1 = g$sd1, sd2 = g$sd2, rho = g$rho, alpha = 0.025,
    known_var = TRUE)$powerCoprimary
  grid$ref_known[i] <- ref_known(g$n1, g$n2, g$delta1, g$delta2,
                                 g$sd1, g$sd2, g$rho, 0.025)
}
grid$abs_diff_known <- abs(grid$pkg_known - grid$ref_known)
say(sprintf("  designs compared            : %d", nrow(grid)))
say(sprintf("  max |package - closed form| : %.3e", max(grid$abs_diff_known)))
worst <- grid[which.max(grid$abs_diff_known), ]
say(sprintf("  worst design                : r = %.1f, rho = %.1f, delta = (%g, %g), sd = (%g, %g)",
            worst$r, worst$rho, worst$delta1, worst$delta2, worst$sd1, worst$sd2))
check("A.known.closed.form", max(grid$abs_diff_known) < 1e-6,
      sprintf("max deviation %.3e", max(grid$abs_diff_known)))
say("")

# ---------------------------------------------------------------------------
# Part B: known variance against a direct patient level simulation.
# ---------------------------------------------------------------------------
say("--- Part B: known variance against direct simulation ---")
sim_known <- function(n1, n2, d1, d2, s1, s2, rho, alpha, nsim) {
  S <- matrix(c(s1^2, rho * s1 * s2, rho * s1 * s2, s2^2), 2)
  L <- chol(S)
  draw_mean <- function(n) {
    z <- matrix(stats::rnorm(2 * nsim), ncol = 2) %*% L / sqrt(n)
    z
  }
  m1 <- draw_mean(n1); m1[, 1] <- m1[, 1] + d1; m1[, 2] <- m1[, 2] + d2
  m2 <- draw_mean(n2)
  se <- sqrt(1 / n1 + 1 / n2)
  Z1 <- (m1[, 1] - m2[, 1]) / (s1 * se)
  Z2 <- (m1[, 2] - m2[, 2]) / (s2 * se)
  z <- stats::qnorm(1 - alpha)
  c(power = mean(Z1 > z & Z2 > z), cor = stats::cor(Z1, Z2))
}
nsim <- 200000
sub <- grid[grid$rho %in% c(0, 0.5, 0.8) & grid$r %in% c(0.5, 1, 3), ]
rows <- list()
for (i in seq_len(nrow(sub))) {
  g <- sub[i, ]
  s <- sim_known(g$n1, g$n2, g$delta1, g$delta2, g$sd1, g$sd2, g$rho, 0.025, nsim)
  se3 <- 3 * sqrt(s[["power"]] * (1 - s[["power"]]) / nsim)
  rows[[length(rows) + 1]] <- data.frame(
    r = g$r, rho = g$rho, case = g$case,
    package = g$pkg_known, simulated = s[["power"]],
    within_3SE = abs(g$pkg_known - s[["power"]]) <= se3,
    cor_Z = s[["cor"]])
}
simtab <- do.call(rbind, rows)
write.csv(simtab, file.path(out_dir, "continuous_sim_known.csv"), row.names = FALSE)
say(sprintf("  designs simulated                : %d (nsim = %d)", nrow(simtab), nsim))
say(sprintf("  agreeing with simulation (3 SE)  : %d of %d",
            sum(simtab$within_3SE), nrow(simtab)))
check("B.known.simulation", all(simtab$within_3SE),
      sprintf("max |pkg - sim| = %.5f", max(abs(simtab$package - simtab$simulated))))

# the test statistic correlation should equal rho whatever the allocation
dev_cor <- max(abs(simtab$cor_Z - simtab$rho))
say(sprintf("  max |cor(Z1, Z2) - rho|          : %.5f  (theory says gamma = rho)", dev_cor))
check("B.gamma.equals.rho", dev_cor < 0.01, sprintf("max deviation %.5f", dev_cor))
say("")

# ---------------------------------------------------------------------------
# Part C: unknown variance against a direct patient level simulation with the
# t-test. Individual patients are generated, so the sample covariance matrix
# is not drawn from the Wishart distribution the package assumes.
# ---------------------------------------------------------------------------
say("--- Part C: unknown variance against direct simulation ---")
sim_unknown <- function(n1, n2, d1, d2, s1, s2, rho, alpha, nsim) {
  S <- matrix(c(s1^2, rho * s1 * s2, rho * s1 * s2, s2^2), 2)
  L <- chol(S)
  nu <- n1 + n2 - 2
  tcrit <- stats::qt(1 - alpha, nu)
  rej <- logical(nsim)
  for (b in seq_len(nsim)) {
    x1 <- matrix(stats::rnorm(2 * n1), ncol = 2) %*% L
    x1[, 1] <- x1[, 1] + d1; x1[, 2] <- x1[, 2] + d2
    x2 <- matrix(stats::rnorm(2 * n2), ncol = 2) %*% L
    m1 <- colMeans(x1); m2 <- colMeans(x2)
    v1 <- apply(x1, 2, stats::var); v2 <- apply(x2, 2, stats::var)
    sp <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / nu)
    tt <- (m1 - m2) / (sp * sqrt(1 / n1 + 1 / n2))
    rej[b] <- all(tt > tcrit)
  }
  mean(rej)
}
nsim_t <- 20000
sub_t <- grid[grid$rho %in% c(0, 0.5) & grid$r %in% c(1, 2), ]
rows <- list()
for (i in seq_len(nrow(sub_t))) {
  g <- sub_t[i, ]
  pk <- power2Continuous(n1 = g$n1, n2 = g$n2, delta1 = g$delta1,
                         delta2 = g$delta2, sd1 = g$sd1, sd2 = g$sd2,
                         rho = g$rho, alpha = 0.025,
                         known_var = FALSE, nMC = 50000)$powerCoprimary
  sm <- sim_unknown(g$n1, g$n2, g$delta1, g$delta2, g$sd1, g$sd2,
                    g$rho, 0.025, nsim_t)
  se3 <- 3 * sqrt(sm * (1 - sm) / nsim_t)
  say(sprintf("  r = %.1f rho = %.1f case %d : package = %.4f  simulated = %.4f  (3 SE = %.4f)",
              g$r, g$rho, g$case, pk, sm, se3))
  rows[[length(rows) + 1]] <- data.frame(r = g$r, rho = g$rho, case = g$case,
                                         package = pk, simulated = sm,
                                         within_3SE = abs(pk - sm) <= se3)
}
ttab <- do.call(rbind, rows)
write.csv(ttab, file.path(out_dir, "continuous_sim_unknown.csv"), row.names = FALSE)
check("C.unknown.simulation", all(ttab$within_3SE),
      sprintf("%d of %d within 3 SE", sum(ttab$within_3SE), nrow(ttab)))
say("")

# ---------------------------------------------------------------------------
# Part D: scale invariance. Multiplying both the effect and the standard
# deviation by a constant must leave the power unchanged, in both branches.
# ---------------------------------------------------------------------------
say("--- Part D: scale invariance ---")
scales <- c(1, 7, 250)
dev_known <- 0; dev_unknown <- 0
for (i in seq_len(nrow(sub_t))) {
  g <- sub_t[i, ]
  pk <- vapply(scales, function(k) power2Continuous(
    n1 = g$n1, n2 = g$n2, delta1 = g$delta1 * k, delta2 = g$delta2 * k,
    sd1 = g$sd1 * k, sd2 = g$sd2 * k, rho = g$rho, alpha = 0.025,
    known_var = TRUE)$powerCoprimary, 0)
  dev_known <- max(dev_known, diff(range(pk)))
  set.seed(99)
  pu <- vapply(scales, function(k) { set.seed(99); power2Continuous(
    n1 = g$n1, n2 = g$n2, delta1 = g$delta1 * k, delta2 = g$delta2 * k,
    sd1 = g$sd1 * k, sd2 = g$sd2 * k, rho = g$rho, alpha = 0.025,
    known_var = FALSE, nMC = 20000)$powerCoprimary }, 0)
  dev_unknown <- max(dev_unknown, diff(range(pu)))
}
say(sprintf("  max spread across scales, known variance   : %.3e", dev_known))
say(sprintf("  max spread across scales, unknown variance : %.3e", dev_unknown))
check("D.scale.invariance", dev_known < 1e-12 && dev_unknown < 1e-12,
      sprintf("known %.3e, unknown %.3e", dev_known, dev_unknown))
say("")

# ---------------------------------------------------------------------------
# Part E: does ss2Continuous return the minimum sample size when r is not 1?
# ---------------------------------------------------------------------------
say("--- Part E: minimality of ss2Continuous across allocation ratios ---")
rows <- list()
for (rv in c(0.5, 1, 1.5, 2, 3)) {
  for (ci in 1:3) {
    cs <- cases[[ci]]
    res <- ss2Continuous(delta1 = cs$delta1, delta2 = cs$delta2,
                         sd1 = cs$sd1, sd2 = cs$sd2, rho = 0.5, r = rv,
                         alpha = 0.025, beta = 0.2, known_var = TRUE)
    at_n <- power2Continuous(n1 = res$n1, n2 = res$n2,
                             delta1 = cs$delta1, delta2 = cs$delta2,
                             sd1 = cs$sd1, sd2 = cs$sd2, rho = 0.5,
                             alpha = 0.025, known_var = TRUE)$powerCoprimary
    below <- power2Continuous(n1 = ceiling(rv * (res$n2 - 1)), n2 = res$n2 - 1,
                              delta1 = cs$delta1, delta2 = cs$delta2,
                              sd1 = cs$sd1, sd2 = cs$sd2, rho = 0.5,
                              alpha = 0.025, known_var = TRUE)$powerCoprimary
    ok <- at_n >= 0.8 && below < 0.8
    say(sprintf("  r = %.1f case %d : n1 = %4s n2 = %4s  power(n2) = %.4f  power(n2-1) = %.4f  minimal = %s",
                rv, ci, format(res$n1), format(res$n2), at_n, below,
                if (ok) "YES" else "NO"))
    rows[[length(rows) + 1]] <- data.frame(r = rv, case = ci, n1 = res$n1,
                                           n2 = res$n2, N = res$N,
                                           power_at_n = at_n,
                                           power_below = below, minimal = ok)
  }
}
mtab <- do.call(rbind, rows)
write.csv(mtab, file.path(out_dir, "continuous_minimality.csv"), row.names = FALSE)
check("E.minimality", all(mtab$minimal),
      sprintf("%d of %d designs minimal", sum(mtab$minimal), nrow(mtab)))
say("")

# ---------------------------------------------------------------------------
# Part F: the article's own numbers, recomputed.
# ---------------------------------------------------------------------------
say("--- Part F: the numbers the article reports for this endpoint type ---")
for (rv in c(0, 0.3, 0.5, 0.8)) {
  res <- ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
                       rho = rv, r = 1, alpha = 0.025, beta = 0.2,
                       known_var = TRUE)
  say(sprintf("  rho = %.1f : N = %s", rv, format(res$N)))
}
say("  (the article reports 166, 162, 158, 148 and an 11%% reduction at rho = 0.8)")
say("")

# ---------------------------------------------------------------------------
# Part G: two claims in the article that no existing evidence file covers.
# ---------------------------------------------------------------------------
say("--- Part G: two unverified statements in the article ---")

# G1. "known_var = FALSE ... adds one subject per group relative to the
#      known-variance result at the same design parameters"
kv <- ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1, rho = 0.5,
                    r = 1, alpha = 0.025, beta = 0.2, known_var = TRUE)
set.seed(20260829)
uv <- ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1, rho = 0.5,
                    r = 1, alpha = 0.025, beta = 0.2, known_var = FALSE)
say(sprintf("  known variance   : n per group = %s", format(kv$n2)))
say(sprintf("  unknown variance : n per group = %s", format(uv$n2)))
check("G1.unknown.adds.one", uv$n2 - kv$n2 == 1,
      sprintf("difference is %s subject(s) per group; the article says one",
              format(uv$n2 - kv$n2)))

# G2. "Four of the five endpoint type combinations ... return in well under a
#      second". Timed on the designs the article itself uses.
tm <- c(
  continuous = system.time(ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1,
    sd2 = 1, rho = 0.5, r = 1, alpha = 0.025, beta = 0.2,
    known_var = TRUE))[["elapsed"]],
  binary_approx = system.time(ss2BinaryApprox(p11 = 0.40, p12 = 0.35,
    p21 = 0.25, p22 = 0.20, rho1 = 0.5, rho2 = 0.5, r = 1, alpha = 0.025,
    beta = 0.2, Test = "AN"))[["elapsed"]],
  mixed_cont_binary = system.time(ss2MixedContinuousBinary(delta = 0.5, sd = 1,
    p1 = 0.60, p2 = 0.40, rho = 0.5, r = 1, alpha = 0.025, beta = 0.2,
    Test = "AN"))[["elapsed"]],
  mixed_count_cont = system.time(ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25,
    nu = 0.8, t = 1, mu1 = -50, mu2 = 0, sd = 250, rho1 = 0.5, rho2 = 0.5,
    r = 1, alpha = 0.025, beta = 0.2))[["elapsed"]]
)
for (nm in names(tm)) say(sprintf("  %-20s %.3f s", nm, tm[[nm]]))
write.csv(data.frame(combination = names(tm), seconds = as.numeric(tm)),
          file.path(out_dir, "continuous_timing.csv"), row.names = FALSE)
check("G2.four.combinations.fast", all(tm < 1),
      sprintf("slowest is %.3f s", max(tm)))
say("")

say("=== summary ===")
say("  PASS : ", pass)
say("  FAIL : ", fail)
cat("\nDone. See", log_path, "\n")
