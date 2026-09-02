# Verification of the tie handling fix in rr1Binary().
#
# Part A checks the rejection region of the installed package against the one
# implied by the Exact package, which computes the tail probability of the
# ordering statistic directly and therefore includes whole tie groups by
# construction. Cells whose p-value sits within the grid tolerance of alpha are
# reported separately, because rr1Binary() maximizes over seq(0, 1, l = 100)
# while Exact uses its own grid of the same size, and the two can differ by
# order 1e-4.
#
# Part B re-runs the sample sizes that the manuscript reports and compares them
# with the values recorded before the fix.
#
# Requires the package to be reinstalled and the R session restarted first.
#
# Run from the package root:
#   source("dev/verify_exact_binary_ties.R")
#
# Writes dev/out/verify_tie_fix_rr.csv, dev/out/verify_tie_fix_ss.csv and
# dev/out/verify_tie_fix.log

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "verify_tie_fix.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_con)
}

say("twoCoprimary version: ", as.character(utils::packageVersion("twoCoprimary")))
say("R version: ", R.version.string)
say("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# Fail early if the session still holds a stale copy of the package
probe <- tryCatch({
  rr1Binary(5, 5, 0.025, "Z-pool")
  TRUE
}, error = function(e) {
  message("Cannot call rr1Binary(): ", conditionMessage(e))
  message("Reinstall the package, restart the R session ",
          "(Session > Restart R, or Ctrl+Shift+F10) and run this script again.")
  FALSE
})
if (!probe) {
  close(log_con)
  stop("stale package state; reinstall, restart R and re-run")
}

# ---------------------------------------------------------------------------
# Part A: rejection region against the Exact package
# ---------------------------------------------------------------------------

if (!requireNamespace("Exact", quietly = TRUE)) {
  say("--- Part A skipped: package 'Exact' is not installed ---")
} else {
  say("--- Part A: rejection region against Exact ",
      as.character(utils::packageVersion("Exact")), " ---")

  methods <- c("Z-pool" = "z-pooled", "Boschloo" = "boschloo")
  n1 <- 15; n2 <- 15
  alphas <- c(0.01, 0.025, 0.05)

  rr_res <- data.frame()
  for (tst in names(methods)) {

    # p-values from Exact over the whole outcome grid
    p_ex <- matrix(1, n1 + 1, n2 + 1)
    for (x1 in 0:n1) {
      for (x2 in 0:n2) {
        tab <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2, byrow = TRUE)
        p_ex[x1 + 1, x2 + 1] <- tryCatch(
          Exact::exact.test(tab, alternative = "greater",
                            method = methods[[tst]],
                            np.interval = FALSE, npNumbers = 100,
                            ref.pvalue = FALSE, to.plot = FALSE)$p.value,
          error = function(e) NA_real_
        )
      }
    }

    for (a in alphas) {
      RR_pkg <- rr1Binary(n1, n2, a, tst)
      RR_ex <- p_ex < a
      diff_cells <- which(RR_pkg != RR_ex & !is.na(p_ex), arr.ind = TRUE)

      # A disagreement is attributable to the nuisance grid when the Exact
      # p-value lies within the observed grid discrepancy of alpha
      near <- if (nrow(diff_cells) > 0) {
        abs(p_ex[diff_cells] - a) < 2e-4
      } else logical(0)

      rr_res <- rbind(rr_res, data.frame(
        Test = tst, alpha = a,
        n_cells = (n1 + 1) * (n2 + 1),
        n_reject_package = sum(RR_pkg),
        n_reject_Exact = sum(RR_ex, na.rm = TRUE),
        n_disagree = nrow(diff_cells),
        n_disagree_near_alpha = sum(near),
        n_disagree_clear = sum(!near)
      ))
      say("  ", tst, " alpha=", a,
          "  package rejects ", sum(RR_pkg),
          ", Exact rejects ", sum(RR_ex, na.rm = TRUE),
          ", disagreements ", nrow(diff_cells),
          " (", sum(near), " within the grid tolerance of alpha)")
    }
  }

  write.csv(rr_res, file.path(out_dir, "verify_tie_fix_rr.csv"), row.names = FALSE)
  say("")
  say("disagreements not explained by the nuisance grid : ",
      sum(rr_res$n_disagree_clear))
  say("")
}

# ---------------------------------------------------------------------------
# Part B: the sample sizes the manuscript reports
# ---------------------------------------------------------------------------
# Expected values, recorded from verify_tie_impact_on_manuscript.R before the fix was applied.

expected <- data.frame(
  source = c(rep("validation", 16), rep("usage", 2)),
  Test = c(rep("Z-pool", 8), rep("Boschloo", 8), "Z-pool", "Boschloo"),
  rho = c(rep(c(0, 0.3, 0.5, 0.8), 4), 0.5, 0.5),
  r = c(rep(c(1, 1, 1, 1, 2, 2, 2, 2), 2), 1, 1),
  N_before = c(144, 142, 140, 134, 180, 180, 177, 168,
               144, 142, 140, 134, 162, 159, 156, 150,
               368, 368),
  stringsAsFactors = FALSE
)

say("--- Part B: sample sizes reported in the manuscript ---")

pars <- function(src) {
  if (src == "validation") {
    list(p11 = 0.54, p12 = 0.54, p21 = 0.25, p22 = 0.25, alpha = 0.025, beta = 0.1)
  } else {
    list(p11 = 0.40, p12 = 0.35, p21 = 0.25, p22 = 0.20, alpha = 0.025, beta = 0.2)
  }
}

expected$N_after <- NA_integer_
for (k in seq_len(nrow(expected))) {
  q <- pars(expected$source[k])
  expected$N_after[k] <- ss2BinaryExact(
    q$p11, q$p12, q$p21, q$p22,
    expected$rho[k], expected$rho[k],
    expected$r[k], q$alpha, q$beta, expected$Test[k]
  )[["N"]]
  say("  ", k, " / ", nrow(expected), ": ", expected$source[k], " ",
      expected$Test[k], " rho=", expected$rho[k], " r=", expected$r[k],
      "  before=", expected$N_before[k], " after=", expected$N_after[k],
      if (expected$N_after[k] == expected$N_before[k]) "" else "   CHANGED")
}

expected$changed <- expected$N_after != expected$N_before
write.csv(expected, file.path(out_dir, "verify_tie_fix_ss.csv"), row.names = FALSE)

say("")
say("designs whose sample size changed : ", sum(expected$changed),
    " of ", nrow(expected))
if (any(expected$changed)) {
  print(expected[expected$changed, ])
  capture.output(print(expected[expected$changed, ]), file = log_con, append = TRUE)
}

close(log_con)
