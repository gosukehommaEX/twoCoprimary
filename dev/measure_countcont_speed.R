# Locate the cost in ss2MixedCountContinuous(), which took 4.7 seconds in
# dev/out/continuous.log while the other three closed-form combinations took
# under a tenth of a second.
#
# The hypothesis to test is that almost all of the time is spent in
# corrbound2MixedCountContinuous(), whose Frechet-Hoeffding bounds have no
# closed form for a negative binomial and a normal marginal and are obtained
# by numerical integration over the support of the count distribution. The
# bounds do not depend on n1 or n2, yet power2MixedCountContinuous() calls
# them once per group on every evaluation, and the sequential search calls
# the power function many times.
#
# Run with:  source("dev/measure_countcont_speed.R")

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(out_dir, "countcont_cost.log")
if (file.exists(log_path)) file.remove(log_path)

say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_path, append = TRUE)
}

say("twoCoprimary version : ", as.character(utils::packageVersion("twoCoprimary")))
say("run at               : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# The design used in the article's usage example
args_design <- list(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                    mu1 = -50, mu2 = 0, sd = 250,
                    rho1 = 0.5, rho2 = 0.5, r = 1,
                    alpha = 0.025, beta = 0.2)
lambda1 <- args_design$r1 * args_design$t
lambda2 <- args_design$r2 * args_design$t

# --- 1. the whole sample size search -----------------------------------------
t_ss <- system.time(res <- do.call(ss2MixedCountContinuous, args_design))[["elapsed"]]
say(sprintf("total ss2MixedCountContinuous()      : %6.3f s   (n2 = %s)",
            t_ss, format(res$n2)))

# --- 2. one call to the correlation bounds, per group ------------------------
t_b1 <- system.time(b1 <- corrbound2MixedCountContinuous(
  lambda1, args_design$nu, args_design$mu1, args_design$sd))[["elapsed"]]
t_b2 <- system.time(b2 <- corrbound2MixedCountContinuous(
  lambda2, args_design$nu, args_design$mu2, args_design$sd))[["elapsed"]]
say(sprintf("one corrbound call, group 1          : %6.3f s   bounds [%.4f, %.4f]",
            t_b1, b1[1], b1[2]))
say(sprintf("one corrbound call, group 2          : %6.3f s   bounds [%.4f, %.4f]",
            t_b2, b2[1], b2[2]))

# --- 3. one power evaluation, which calls the bounds twice --------------------
t_pw <- system.time(power2MixedCountContinuous(
  n1 = res$n1, n2 = res$n2, r1 = args_design$r1, r2 = args_design$r2,
  nu = args_design$nu, t = args_design$t, mu1 = args_design$mu1,
  mu2 = args_design$mu2, sd = args_design$sd,
  rho1 = args_design$rho1, rho2 = args_design$rho2,
  alpha = args_design$alpha))[["elapsed"]]
say(sprintf("one power2MixedCountContinuous()     : %6.3f s", t_pw))
say("")

# --- 4. how many power evaluations does the search make? ---------------------
# Counted by wrapping the power function and letting the search call the wrapper.
n_calls <- 0L
counting_power <- function(n1, n2, ...) {
  n_calls <<- n_calls + 1L
  power2MixedCountContinuous(n1 = n1, n2 = n2, ...)
}
search_fun <- get(".ss_sequential_search", envir = asNamespace("twoCoprimary"))
init <- ss1Count(r1 = args_design$r1, r2 = args_design$r2, nu = args_design$nu,
                 t = args_design$t, r = args_design$r,
                 alpha = args_design$alpha, beta = args_design$beta)[["n2"]]
t_cnt <- system.time(
  search_fun(initial_n2 = init, r = args_design$r,
             target_power = 1 - args_design$beta,
             power_fun = counting_power,
             r1 = args_design$r1, r2 = args_design$r2, nu = args_design$nu,
             t = args_design$t, mu1 = args_design$mu1, mu2 = args_design$mu2,
             sd = args_design$sd, rho1 = args_design$rho1,
             rho2 = args_design$rho2, alpha = args_design$alpha)
)[["elapsed"]]
say(sprintf("power evaluations in one search      : %d", n_calls))
say(sprintf("bounds calls implied (2 per power)   : %d", 2L * n_calls))
say(sprintf("predicted bounds time                : %6.3f s",
            n_calls * (t_b1 + t_b2)))
say(sprintf("observed search time                 : %6.3f s", t_cnt))
share <- 100 * n_calls * (t_b1 + t_b2) / max(t_cnt, 1e-9)
say(sprintf("share of the search spent in bounds  : %5.1f %%", share))
say("")

# --- 5. are the bounds invariant to the sample size? -------------------------
say("The bounds depend only on lambda, nu, mu and sd, so they cannot change")
say("with n1 or n2. Confirming that the same two values are recomputed:")
same <- TRUE
for (n in c(20, 60, 200)) {
  bb1 <- corrbound2MixedCountContinuous(lambda1, args_design$nu,
                                        args_design$mu1, args_design$sd)
  bb2 <- corrbound2MixedCountContinuous(lambda2, args_design$nu,
                                        args_design$mu2, args_design$sd)
  same <- same && isTRUE(all.equal(as.numeric(bb1), as.numeric(b1))) &&
    isTRUE(all.equal(as.numeric(bb2), as.numeric(b2)))
}
say(sprintf("  identical across sample sizes      : %s", same))
say("")

write.csv(data.frame(
  quantity = c("total_search", "one_bounds_g1", "one_bounds_g2", "one_power",
               "power_calls", "predicted_bounds_time", "observed_search_time"),
  value = c(t_ss, t_b1, t_b2, t_pw, n_calls, n_calls * (t_b1 + t_b2), t_cnt)
), file.path(out_dir, "countcont_cost.csv"), row.names = FALSE)

say("=== conclusion ===")
if (share > 80) {
  say("The correlation bounds account for most of the run time, and they are")
  say("invariant to the sample size, so the cost is avoidable.")
} else {
  say("The bounds do NOT account for most of the run time. The diagnosis is")
  say("wrong and the cost is elsewhere; do not write the explanation into the")
  say("article until it has been located.")
}

cat("\nDone. See", log_path, "\n")
