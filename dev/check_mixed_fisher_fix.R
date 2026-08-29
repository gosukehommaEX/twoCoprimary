# Verification of the corrected Fisher option for mixed continuous and binary
# endpoints.
#
# Two defects were repaired. The asymptotic branch is now inside an else, so the
# Monte Carlo result is no longer overwritten by the ASc result. The latent
# binary variable is now centred in both groups and dichotomised at a group
# specific threshold g_j = qnorm(1 - p_j), following Sozu et al. (2012), so the
# simulated response probabilities are p1 and p2. The earlier code centred the
# latent variables at qnorm(1 - p_j) and used a single threshold qnorm(1 - p2),
# which forced the control probability to 0.5 and reversed the effect.
#
# Part A  The simulated response probabilities now equal p1 and p2.
# Part B  Fisher is no longer identical to ASc and now depends on nMC.
# Part C  Validation against Table S5 of the Supporting Information of Sozu
#         et al. (2012), which reports, for each design, the sample size per
#         group and the empirical overall power that Fisher's exact test
#         achieves there. The power at the published sample size is a more
#         stable check than re-running the search, because a Monte Carlo power
#         makes the sequential search jitter by a subject or two.
# Part D  The sequential search reproduces the published sample size on the
#         designs that are cheap enough to search.
#
# Run from the package root:
#   source("dev/check_mixed_fisher_fix.R")
#
# Writes dev/out/mixed_fisher_fix_S5.csv, dev/out/mixed_fisher_fix_search.csv
# and dev/out/mixed_fisher_fix.log

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "mixed_fisher_fix.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_con)
}

say("twoCoprimary version: ", as.character(utils::packageVersion("twoCoprimary")))
say("R version: ", R.version.string)
say("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# ---------------------------------------------------------------------------
# Part A: the simulated response probabilities
# ---------------------------------------------------------------------------

say("--- Part A: response probabilities implied by the latent construction ---")

for (pp in list(c(0.60, 0.40), c(0.99, 0.85), c(0.95, 0.80))) {
  p1 <- pp[1]; p2 <- pp[2]
  gT <- qnorm(1 - p1)
  gC <- qnorm(1 - p2)
  say("  target p1 = ", p1, ", p2 = ", p2,
      "  implied ", format(pnorm(0 - gT), digits = 6),
      " and ", format(pnorm(0 - gC), digits = 6))
}
say("")

# ---------------------------------------------------------------------------
# Part B: Fisher is a distinct, simulation based method again
# ---------------------------------------------------------------------------

say("--- Part B: Fisher against the asymptotic methods ---")

args_common <- list(n1 = 60, n2 = 60, delta = 0.5, sd = 1,
                    p1 = 0.60, p2 = 0.40, rho = 0.5, alpha = 0.025)

for (tst in c("AN", "ANc", "AS", "ASc")) {
  v <- do.call(power2MixedContinuousBinary, c(args_common, list(Test = tst)))
  say("  ", formatC(tst, width = 7), " powerCoprimary = ",
      format(v$powerCoprimary, digits = 6))
}
for (m in c(2000, 20000, 100000)) {
  set.seed(20260828)
  v <- do.call(power2MixedContinuousBinary,
               c(args_common, list(Test = "Fisher", nMC = m)))
  say("  Fisher  nMC = ", formatC(m, width = 6), " powerCoprimary = ",
      format(v$powerCoprimary, digits = 6), "  nMC reported = ", v$nMC)
}
say("")
say("Fisher should now differ from ASc and should settle as nMC grows.")
say("")

# ---------------------------------------------------------------------------
# Part C: Table S5 of Sozu et al. (2012), Supporting Information Section C
# ---------------------------------------------------------------------------
# Sample size per group and empirical overall power for Fisher's exact test,
# alpha = 0.025, overall power 0.8, balanced allocation, sd = 1 so that delta is
# the standardized effect delta1*.

S5 <- data.frame(
  delta = rep(c(0.235, 0.397, 0.521, 0.190, 0.335, 0.457), each = 4),
  p1    = rep(c(0.99, 0.99, 0.99, 0.95, 0.95, 0.95), each = 4),
  p2    = rep(c(0.95, 0.90, 0.85, 0.90, 0.85, 0.80), each = 4),
  rho   = rep(c(0.0, 0.3, 0.5, 0.8), times = 6),
  n_published = c(375, 372, 370, 366,
                  132, 131, 130, 128,
                   76,  75,  75,  74,
                  584, 577, 571, 563,
                  190, 188, 186, 183,
                  102, 101, 100,  98),
  power_published = c(0.801, 0.801, 0.800, 0.800,
                      0.801, 0.802, 0.802, 0.799,
                      0.800, 0.798, 0.802, 0.801,
                      0.800, 0.799, 0.799, 0.801,
                      0.800, 0.801, 0.800, 0.802,
                      0.803, 0.803, 0.803, 0.803)
)

nMC_C <- 20000
say("--- Part C: power at the published sample size (nMC = ", nMC_C, ") ---")

S5$power_package <- NA_real_
for (k in seq_len(nrow(S5))) {
  set.seed(1000 + k)
  S5$power_package[k] <- power2MixedContinuousBinary(
    n1 = S5$n_published[k], n2 = S5$n_published[k],
    delta = S5$delta[k], sd = 1,
    p1 = S5$p1[k], p2 = S5$p2[k], rho = S5$rho[k],
    alpha = 0.025, Test = "Fisher", nMC = nMC_C
  )$powerCoprimary
  say("  ", k, " / ", nrow(S5), ": delta = ", S5$delta[k],
      " (", S5$p1[k], ", ", S5$p2[k], ") rho = ", S5$rho[k],
      "  n = ", S5$n_published[k],
      "  published ", format(S5$power_published[k], nsmall = 3),
      "  package ", format(round(S5$power_package[k], 3), nsmall = 3))
}
S5$difference <- round(S5$power_package - S5$power_published, 4)
write.csv(S5, file.path(out_dir, "mixed_fisher_fix_S5.csv"), row.names = FALSE)

se <- sqrt(0.8 * 0.2 / nMC_C)
say("")
say("Monte Carlo standard error at nMC = ", nMC_C, " : ", format(se, digits = 3))
say("max |package - published| : ", format(max(abs(S5$difference)), digits = 3))
say("cells within 3 standard errors : ",
    sum(abs(S5$difference) < 3 * se), " of ", nrow(S5))
say("")

# ---------------------------------------------------------------------------
# Part D: the sequential search on the cheaper designs
# ---------------------------------------------------------------------------

say("--- Part D: sample size search on the designs with small n ---")

cheap <- S5[S5$n_published <= 135, ]
cheap$n_search <- NA_integer_
for (k in seq_len(nrow(cheap))) {
  set.seed(5000 + k)
  ss <- ss2MixedContinuousBinary(
    delta = cheap$delta[k], sd = 1,
    p1 = cheap$p1[k], p2 = cheap$p2[k], rho = cheap$rho[k],
    r = 1, alpha = 0.025, beta = 0.2, Test = "Fisher", nMC = 20000
  )
  cheap$n_search[k] <- ss$n1
  say("  delta = ", cheap$delta[k], " (", cheap$p1[k], ", ", cheap$p2[k],
      ") rho = ", cheap$rho[k],
      "  published ", cheap$n_published[k], "  search ", ss$n1,
      "  difference ", ss$n1 - cheap$n_published[k])
}
cheap$difference_n <- cheap$n_search - cheap$n_published
write.csv(cheap, file.path(out_dir, "mixed_fisher_fix_search.csv"), row.names = FALSE)
say("")
say("max |search - published| : ", max(abs(cheap$difference_n)), " subjects per group")
say("A Monte Carlo power makes the search jitter, so agreement within a couple")
say("of subjects is the most that can be expected.")

close(log_con)
