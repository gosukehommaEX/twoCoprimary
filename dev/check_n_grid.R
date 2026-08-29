# Verification of the n_grid argument.
#
# n_grid exposes the number of points at which the null tail probability is
# maximized over the nuisance parameter in the Z-pooled and Boschloo exact
# unconditional tests. It was previously hard coded as seq(0, 1, l = 100).
#
# Part A  The default reproduces the sample sizes obtained before the argument
#         existed, so nothing in the article or the vignettes moves.
# Part B  Convergence of the p-values in n_grid, taking a very fine grid as the
#         reference. A coarse grid can only understate the maximum, so the
#         p-values must increase towards the reference.
# Part C  Convergence of the required sample size in n_grid, which is the
#         material for the performance and limitations discussion.
# Part D  The argument rejects values that would silently give a wrong answer.
# Part E  Computation time against n_grid.
#
# Run devtools::document() and reinstall before running this, since the new
# argument needs its help files regenerated.
#
# Run from the package root:
#   source("dev/check_n_grid.R")
#
# Writes dev/out/n_grid_default.csv, dev/out/n_grid_pvalue.csv,
# dev/out/n_grid_samplesize.csv, dev/out/n_grid_timing.csv and
# dev/out/n_grid.log

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "n_grid.log"), open = "wt")
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
  rr1Binary(5, 5, 0.025, "Z-pool", n_grid = 100)
  TRUE
}, error = function(e) {
  message("Cannot call rr1Binary() with n_grid: ", conditionMessage(e))
  message("Run devtools::document(), reinstall, restart R and try again.")
  FALSE
})
if (!probe) {
  close(log_con)
  stop("n_grid not available; document, reinstall, restart R and re-run")
}

# ---------------------------------------------------------------------------
# Part A: the default changes nothing
# ---------------------------------------------------------------------------
# Validation designs of the article: p11 = p12 = 0.54, p21 = p22 = 0.25,
# alpha = 0.025, 1 - beta = 0.9, r = 1.

say("--- Part A: the default reproduces the known sample sizes ---")

known <- data.frame(
  Test = c(rep("Z-pool", 3), rep("Boschloo", 3)),
  rho = rep(c(0, 0.3, 0.5), 2),
  N_expected = c(144, 142, 140, 144, 142, 140),
  stringsAsFactors = FALSE
)
known$N_default <- NA_integer_
known$N_explicit_100 <- NA_integer_
for (k in seq_len(nrow(known))) {
  known$N_default[k] <- ss2BinaryExact(
    0.54, 0.54, 0.25, 0.25, known$rho[k], known$rho[k],
    1, 0.025, 0.1, known$Test[k])[["N"]]
  known$N_explicit_100[k] <- ss2BinaryExact(
    0.54, 0.54, 0.25, 0.25, known$rho[k], known$rho[k],
    1, 0.025, 0.1, known$Test[k], n_grid = 100)[["N"]]
  say("  ", known$Test[k], " rho=", known$rho[k],
      "  expected=", known$N_expected[k],
      " default=", known$N_default[k],
      " explicit=", known$N_explicit_100[k],
      if (known$N_default[k] == known$N_expected[k] &&
          known$N_explicit_100[k] == known$N_expected[k]) "" else "   *** MISMATCH ***")
}
write.csv(known, file.path(out_dir, "n_grid_default.csv"), row.names = FALSE)
say("")
say("designs reproducing the expected sample size : ",
    sum(known$N_default == known$N_expected &
        known$N_explicit_100 == known$N_expected), " of ", nrow(known))
say("")

# ---------------------------------------------------------------------------
# Part B: convergence of the p-values
# ---------------------------------------------------------------------------
# The rejection region is p < alpha, so comparing regions at a sequence of grid
# sizes against a very fine reference shows how many outcomes are classified
# differently by a coarse grid.

say("--- Part B: convergence of the rejection region in n_grid ---")

grids <- c(25, 50, 100, 200, 400, 800)
reference <- 4000

pv_res <- data.frame()
for (nn in c(20, 40)) {
  for (tst in c("Z-pool", "Boschloo")) {
    RR_ref <- rr1Binary(nn, nn, 0.025, tst, n_grid = reference)
    for (g in grids) {
      RR_g <- rr1Binary(nn, nn, 0.025, tst, n_grid = g)
      pv_res <- rbind(pv_res, data.frame(
        n1 = nn, n2 = nn, Test = tst, n_grid = g,
        n_reject = sum(RR_g),
        n_reject_reference = sum(RR_ref),
        n_cells_differ = sum(RR_g != RR_ref),
        extra_rejections = sum(RR_g & !RR_ref)
      ))
    }
    say("  n = ", nn, ", ", tst, ": cells differing from the n_grid = ",
        reference, " reference, by grid size ",
        paste(pv_res$n_cells_differ[pv_res$n1 == nn & pv_res$Test == tst],
              collapse = ", "))
  }
}
write.csv(pv_res, file.path(out_dir, "n_grid_pvalue.csv"), row.names = FALSE)
say("")
say("A coarse grid can only understate the maximum, so every disagreement")
say("should be an extra rejection. extra_rejections equal to n_cells_differ : ",
    all(pv_res$extra_rejections == pv_res$n_cells_differ))
say("")

# ---------------------------------------------------------------------------
# Part C: convergence of the required sample size
# ---------------------------------------------------------------------------

say("--- Part C: required sample size against n_grid ---")

ss_grids <- c(25, 50, 100, 200, 400)
ss_res <- data.frame()
for (tst in c("Z-pool", "Boschloo")) {
  for (g in ss_grids) {
    tt <- system.time(
      N <- ss2BinaryExact(0.54, 0.54, 0.25, 0.25, 0.3, 0.3,
                          1, 0.025, 0.1, tst, n_grid = g)[["N"]]
    )[["elapsed"]]
    ss_res <- rbind(ss_res, data.frame(
      Test = tst, n_grid = g, N = N, seconds = round(tt, 2)
    ))
    say("  ", tst, " n_grid=", g, "  N=", N, "  (", round(tt, 2), " s)")
  }
}
write.csv(ss_res, file.path(out_dir, "n_grid_samplesize.csv"), row.names = FALSE)
say("")

# ---------------------------------------------------------------------------
# Part D: validation
# ---------------------------------------------------------------------------

say("--- Part D: rejected values of n_grid ---")

bad_values <- list(5, 0, -100, 100.5, c(50, 100), NA_real_, Inf)
labels <- c("5", "0", "-100", "100.5", "c(50, 100)", "NA", "Inf")
for (m in seq_along(bad_values)) {
  ok <- tryCatch({
    rr1Binary(10, 10, 0.025, "Z-pool", n_grid = bad_values[[m]])
    FALSE
  }, error = function(e) TRUE)
  say("  n_grid = ", labels[m], " : ", if (ok) "rejected" else "*** ACCEPTED ***")
}
say("")

# ---------------------------------------------------------------------------
# Part E: computation time
# ---------------------------------------------------------------------------

say("--- Part E: time for one rejection region against n_grid ---")

tim_res <- data.frame()
for (nn in c(50, 100)) {
  for (tst in c("Z-pool", "Boschloo")) {
    for (g in c(50, 100, 200, 400)) {
      tt <- system.time(rr1Binary(nn, nn, 0.025, tst, n_grid = g))[["elapsed"]]
      tim_res <- rbind(tim_res, data.frame(
        n1 = nn, n2 = nn, Test = tst, n_grid = g, seconds = round(tt, 3)
      ))
    }
  }
}
write.csv(tim_res, file.path(out_dir, "n_grid_timing.csv"), row.names = FALSE)
print(tim_res)
capture.output(print(tim_res), file = log_con, append = TRUE)

close(log_con)
