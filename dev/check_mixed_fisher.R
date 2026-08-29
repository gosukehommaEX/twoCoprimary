# Two defects in power2MixedContinuousBinary(Test = "Fisher").
#
# Defect 1, structural. The block that computes the asymptotic result is not
# inside an else branch, so it runs for every value of Test and overwrites
# powerCont, powerBin and powerCoprimary that the Fisher branch has just
# computed. Because the asymptotic block then falls through
# if (Test == "AN" || Test == "ANc") ... else ... and if (Test == "AS") ...
# else ..., a request for "Fisher" returns the ASc result. The Monte Carlo
# simulation still runs first, drawing nMC * n1 + nMC * n2 multivariate normal
# vectors whose result is discarded. ss2MixedContinuousBinary() has the
# corresponding guard (if (Test != "Fisher") nMC <- NA), which is what the
# missing else in the power function should have looked like.
#
# Defect 2, numerical. Inside the Fisher branch the latent binary variable of
# group j has mean qnorm(1 - p_j) and both groups are dichotomised at the same
# threshold g = qnorm(1 - p2). The implied response probability is then
# Phi(qnorm(1 - p_j) - qnorm(1 - p2)), which equals 0.5 for the control group
# whatever p2 is, and is not p1 for the treatment group. For p1 > p2 it is also
# smaller than the control probability, so the effect points the wrong way.
#
# Neither defect touches any published number: the article uses Test = "AN" for
# this endpoint type, the vignette uses AN, ANc, AS and ASc, and the existing
# test files are labelled "excluding Fisher test".
#
# Run from the package root:
#   source("dev/check_mixed_fisher.R")
#
# Writes dev/out/mixed_fisher.log

library(twoCoprimary)
library(mvtnorm)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "mixed_fisher.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_con)
}

say("twoCoprimary version: ", as.character(utils::packageVersion("twoCoprimary")))
say("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# ---------------------------------------------------------------------------
# Defect 1: the Fisher result is discarded
# ---------------------------------------------------------------------------

say("--- Defect 1: does Test = 'Fisher' return the ASc result? ---")

args_common <- list(n1 = 60, n2 = 60, delta = 0.5, sd = 1,
                    p1 = 0.60, p2 = 0.40, rho = 0.5, alpha = 0.025)

res <- list()
for (tst in c("AN", "ANc", "AS", "ASc", "Fisher")) {
  set.seed(11)
  res[[tst]] <- do.call(power2MixedContinuousBinary,
                        c(args_common, list(Test = tst, nMC = 2000)))
}

tab <- do.call(rbind, lapply(names(res), function(tst) {
  data.frame(Test = tst,
             powerCont = res[[tst]]$powerCont,
             powerBin = res[[tst]]$powerBin,
             powerCoprimary = res[[tst]]$powerCoprimary,
             nMC = res[[tst]]$nMC)
}))
print(tab)
capture.output(print(tab), file = log_con, append = TRUE)
say("")
say("Fisher identical to ASc : ",
    isTRUE(all.equal(res[["Fisher"]]$powerCoprimary, res[["ASc"]]$powerCoprimary)))
say("nMC reported for Fisher : ", res[["Fisher"]]$nMC,
    "  (a Monte Carlo method should report the number of replications)")
say("")
say("The result is also independent of nMC, which it would not be if the")
say("simulation were being used:")
for (m in c(500, 2000, 20000)) {
  set.seed(11)
  v <- do.call(power2MixedContinuousBinary,
               c(args_common, list(Test = "Fisher", nMC = m)))$powerCoprimary
  say("  nMC = ", m, " -> powerCoprimary = ", format(v, digits = 10))
}
say("")

# ---------------------------------------------------------------------------
# Defect 2: the latent construction does not reproduce p1 and p2
# ---------------------------------------------------------------------------

say("--- Defect 2: response probabilities implied by the latent construction ---")

p1 <- 0.60
p2 <- 0.40
g <- qnorm(1 - p2)
mu_treat <- qnorm(1 - p1)
mu_ctrl <- qnorm(1 - p2)

say("target probabilities       : p1 = ", p1, ", p2 = ", p2)
say("implied for the treatment  : ", format(pnorm(mu_treat - g), digits = 6))
say("implied for the control    : ", format(pnorm(mu_ctrl - g), digits = 6))
say("")
say("The control value is 0.5 for any p2, and the treatment value is below the")
say("control value even though p1 > p2, so the simulated effect has the wrong")
say("sign. A correct construction sets the group mean to g + qnorm(p_j).")
say("")

set.seed(99)
n_draw <- 2e5
sim_treat <- mean(rnorm(n_draw, mu_treat, 1) >= g)
sim_ctrl <- mean(rnorm(n_draw, mu_ctrl, 1) >= g)
say("simulated response rate, treatment : ", format(sim_treat, digits = 4),
    "  (target ", p1, ")")
say("simulated response rate, control   : ", format(sim_ctrl, digits = 4),
    "  (target ", p2, ")")

close(log_con)
