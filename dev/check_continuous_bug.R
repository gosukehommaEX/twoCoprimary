# Verification of the known_var = FALSE fix in power2Continuous().
#
# The test statistic Z_k = delta_k / (sd_k sqrt(1 / n1 + 1 / n2)) is already
# divided by sd_k, so the Wishart matrix entering the critical value must be
# that of the standardized endpoints. Drawing it from the variance-covariance
# matrix leaves a factor of sd_k inside sqrt(W_kk / nu), which cancels only when
# sd1 = sd2 = 1.
#
# Three checks, each of which the previous implementation fails and the fixed
# one should pass.
#
# Part A  Scale invariance. Power depends on the standardized effects
#         delta_k / sd_k, so multiplying delta and sd by the same constant must
#         leave it unchanged. This is exact and needs no reference value.
# Part B  Agreement with a direct simulation of the trial: generate the data,
#         run the two pooled t tests, and count how often both reject.
# Part C  Results at sd1 = sd2 = 1 are identical to the previous implementation,
#         confirming that nothing in the article or the vignettes moves.
#
# Run from the package root:
#   source("dev/check_continuous_bug.R")
#
# Writes dev/out/continuous_bug_scale.csv, dev/out/continuous_bug_sim.csv,
# dev/out/continuous_bug_sd1.csv and dev/out/continuous_bug.log

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "continuous_bug.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_con)
}

say("twoCoprimary version: ", as.character(utils::packageVersion("twoCoprimary")))
say("R version: ", R.version.string)
say("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

probe <- tryCatch({
  power2Continuous(10, 10, 0.5, 0.5, 1, 1, 0.3, 0.025, known_var = TRUE)
  TRUE
}, error = function(e) {
  message("Cannot call power2Continuous(): ", conditionMessage(e))
  message("Reinstall the package, restart the R session and run this again.")
  FALSE
})
if (!probe) {
  close(log_con)
  stop("stale package state; reinstall, restart R and re-run")
}

nMC <- 2e5

# ---------------------------------------------------------------------------
# Part A: scale invariance
# ---------------------------------------------------------------------------

say("--- Part A: scale invariance (nMC = ", nMC, ") ---")

base <- data.frame(
  n1 = c(50, 50, 80, 120, 30),
  n2 = c(50, 100, 80, 60, 30),
  delta1 = c(0.5, 0.4, 0.6, 0.5, 0.8),
  delta2 = c(0.5, 0.6, 0.4, 0.5, 0.7),
  rho = c(0.0, 0.3, 0.5, 0.8, 0.3)
)
scales <- c(1, 5, 20, 250)

scale_res <- data.frame()
for (k in seq_len(nrow(base))) {
  b <- base[k, ]
  vals <- numeric(length(scales))
  for (m in seq_along(scales)) {
    cc <- scales[m]
    set.seed(20260828)
    vals[m] <- power2Continuous(
      n1 = b$n1, n2 = b$n2,
      delta1 = cc * b$delta1, delta2 = cc * b$delta2,
      sd1 = cc, sd2 = cc, rho = b$rho, alpha = 0.025,
      known_var = FALSE, nMC = nMC
    )$powerCoprimary
  }
  scale_res <- rbind(scale_res, data.frame(
    n1 = b$n1, n2 = b$n2, delta1 = b$delta1, delta2 = b$delta2, rho = b$rho,
    scale_1 = vals[1], scale_5 = vals[2], scale_20 = vals[3], scale_250 = vals[4],
    max_abs_deviation = max(abs(vals - vals[1]))
  ))
}
write.csv(scale_res, file.path(out_dir, "continuous_bug_scale.csv"), row.names = FALSE)
print(scale_res)
capture.output(print(scale_res), file = log_con, append = TRUE)
say("")
say("max deviation across all scales : ",
    format(max(scale_res$max_abs_deviation), digits = 4))
say("(with a common seed this should be 0 up to floating point)")
say("")

# ---------------------------------------------------------------------------
# Part B: agreement with a direct simulation of the trial
# ---------------------------------------------------------------------------

say("--- Part B: direct simulation of the trial ---")

simulate_power <- function(n1, n2, delta1, delta2, sd1, sd2, rho, alpha, nsim) {
  Sigma <- matrix(c(sd1 ^ 2, rho * sd1 * sd2, rho * sd1 * sd2, sd2 ^ 2), nrow = 2)
  nu <- n1 + n2 - 2
  crit <- qt(1 - alpha, df = nu)
  both <- 0L
  block <- 2000L
  done <- 0L
  while (done < nsim) {
    b <- min(block, nsim - done)
    rej <- vapply(seq_len(b), function(i) {
      x1 <- mvtnorm::rmvnorm(n1, mean = c(delta1, delta2), sigma = Sigma)
      x2 <- mvtnorm::rmvnorm(n2, mean = c(0, 0), sigma = Sigma)
      ok <- TRUE
      for (k in 1:2) {
        m1 <- mean(x1[, k]); m2 <- mean(x2[, k])
        sp2 <- ((n1 - 1) * var(x1[, k]) + (n2 - 1) * var(x2[, k])) / nu
        tstat <- (m1 - m2) / sqrt(sp2 * (1 / n1 + 1 / n2))
        if (tstat < crit) { ok <- FALSE; break }
      }
      ok
    }, logical(1))
    both <- both + sum(rej)
    done <- done + b
  }
  both / nsim
}

sim_grid <- data.frame(
  n1 = c(50, 50, 40, 60),
  n2 = c(50, 50, 80, 60),
  delta1 = c(2.5, 12.5, 7.5, 100),
  delta2 = c(2.5, 15.0, 6.0, 90),
  sd1 = c(5, 25, 15, 250),
  sd2 = c(5, 25, 15, 250),
  rho = c(0.0, 0.3, 0.5, 0.8)
)

nsim <- 20000
sim_res <- data.frame()
for (k in seq_len(nrow(sim_grid))) {
  g <- sim_grid[k, ]
  set.seed(20260828 + k)
  p_pkg <- power2Continuous(g$n1, g$n2, g$delta1, g$delta2, g$sd1, g$sd2,
                            g$rho, 0.025, known_var = FALSE, nMC = nMC)$powerCoprimary
  set.seed(90000 + k)
  p_sim <- simulate_power(g$n1, g$n2, g$delta1, g$delta2, g$sd1, g$sd2,
                          g$rho, 0.025, nsim)
  se <- sqrt(p_sim * (1 - p_sim) / nsim)
  sim_res <- rbind(sim_res, cbind(g, data.frame(
    power_package = round(p_pkg, 4),
    power_simulated = round(p_sim, 4),
    difference = round(p_pkg - p_sim, 4),
    sim_se = round(se, 4),
    within_3_se = abs(p_pkg - p_sim) < 3 * se
  )))
  say("  ", k, " / ", nrow(sim_grid), ": package=", round(p_pkg, 4),
      " simulated=", round(p_sim, 4), " (3 SE = ", round(3 * se, 4), ")",
      if (abs(p_pkg - p_sim) < 3 * se) "" else "   *** OUTSIDE ***")
}
write.csv(sim_res, file.path(out_dir, "continuous_bug_sim.csv"), row.names = FALSE)
say("")
say("designs agreeing with simulation within 3 SE : ",
    sum(sim_res$within_3_se), " of ", nrow(sim_res))
say("")

# ---------------------------------------------------------------------------
# Part C: nothing moves at sd1 = sd2 = 1
# ---------------------------------------------------------------------------

say("--- Part C: results at sd1 = sd2 = 1 against the previous implementation ---")

power_old <- function(n1, n2, delta1, delta2, sd1, sd2, rho, alpha, nMC) {
  nu <- n1 + n2 - 2
  Z <- c(delta1, delta2) / (c(sd1, sd2) * sqrt(1 / n1 + 1 / n2))
  Sigma <- matrix(c(sd1 ^ 2, rho * sd1 * sd2, rho * sd1 * sd2, sd2 ^ 2), nrow = 2)
  Ws <- rWishart(nMC, df = nu, Sigma = Sigma)
  t_alpha <- qt(1 - alpha, df = nu)
  s <- sqrt(1 / nu)
  mean(pbivnorm::pbivnorm(
    x = -t_alpha * sqrt(Ws[1, 1, ]) * s + Z[1],
    y = -t_alpha * sqrt(Ws[2, 2, ]) * s + Z[2],
    rho = rho
  ))
}

old_res <- data.frame()
for (k in seq_len(nrow(base))) {
  b <- base[k, ]
  set.seed(20260828)
  p_new <- power2Continuous(b$n1, b$n2, b$delta1, b$delta2, 1, 1, b$rho,
                            0.025, known_var = FALSE, nMC = nMC)$powerCoprimary
  set.seed(20260828)
  p_old <- power_old(b$n1, b$n2, b$delta1, b$delta2, 1, 1, b$rho, 0.025, nMC)
  old_res <- rbind(old_res, data.frame(
    n1 = b$n1, n2 = b$n2, delta1 = b$delta1, delta2 = b$delta2, rho = b$rho,
    power_new = p_new, power_old = p_old, abs_diff = abs(p_new - p_old)
  ))
}
write.csv(old_res, file.path(out_dir, "continuous_bug_sd1.csv"), row.names = FALSE)
print(old_res)
capture.output(print(old_res), file = log_con, append = TRUE)
say("")
say("max |new - old| at sd = 1 : ", format(max(old_res$abs_diff), digits = 4))
say("(should be 0 up to floating point, confirming the article is unaffected)")

close(log_con)
