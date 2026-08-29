# Verification of the compiled dbibinom() kernel.
#
# The profiling in dev/bench_exact_power.R showed that once the co-primary power
# was rewritten as a matrix product, essentially all of the remaining time in
# ss2BinaryExact() went into building the two bivariate binomial probability
# mass matrices, not into rr1Binary(). The sum over the set M is now formed in
# compiled code.
#
# Part A  The compiled kernel agrees with the R reference implementation, which
#         is retained as .dbibinom_g_r() for exactly this purpose.
# Part B  dbibinom() still sums to one over the whole outcome grid, and its
#         margins are still binomial. These hold for any correct implementation
#         and need no reference.
# Part C  Nothing downstream moves: the sample sizes of the article are
#         unchanged.
# Part D  Timing, before and after.
#
# Run Rcpp::compileAttributes(), devtools::document() and reinstall before
# running this.
#
# Run from the package root:
#   source("dev/check_rcpp_dbibinom.R")
#
# Writes dev/out/rcpp_agreement.csv, dev/out/rcpp_identities.csv,
# dev/out/rcpp_samplesize.csv, dev/out/rcpp_timing.csv and dev/out/rcpp.log

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "rcpp.log"), open = "wt")
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
  dbibinom(5, 0:2, 0:2, 0.4, 0.5, 0.1)
  TRUE
}, error = function(e) {
  message("Cannot call dbibinom(): ", conditionMessage(e))
  message("Run Rcpp::compileAttributes(), devtools::document(), reinstall, ",
          "restart R and try again.")
  FALSE
})
if (!probe) {
  close(log_con)
  stop("compiled code not available; rebuild and re-run")
}

# ---------------------------------------------------------------------------
# Part 0: where the package is actually being loaded from
# ---------------------------------------------------------------------------
# Two warnings appeared on the previous run: no package.rds in Meta/, and S3
# methods declared in NAMESPACE but not found. Both are what R reports when it
# treats a directory that has a NAMESPACE but no installed metadata as an
# installed package. The source tree is one such directory. This section shows
# every copy of the package that is visible on the library path.

say("--- Part 0: installation diagnostic ---")
say("working directory : ", getwd())
say("library paths:")
for (lp in .libPaths()) say("  ", lp)

cands <- character(0)
for (lp in .libPaths()) {
  d <- file.path(lp, "twoCoprimary")
  if (dir.exists(d)) cands <- c(cands, d)
}
say("copies of twoCoprimary on the library path : ", length(cands))
for (d in cands) {
  say("  ", d)
  say("      Meta/package.rds : ", file.exists(file.path(d, "Meta", "package.rds")))
  say("      R/twoCoprimary.rdb : ",
      file.exists(file.path(d, "R", "twoCoprimary.rdb")))
  say("      libs/ : ", dir.exists(file.path(d, "libs")))
  say("      NAMESPACE : ", file.exists(file.path(d, "NAMESPACE")))
}
say("attached from : ", paste(find.package("twoCoprimary", quiet = TRUE),
                              collapse = " | "))

ns <- asNamespace("twoCoprimary")
say("objects present in the loaded namespace:")
for (f in c("print.twoCoprimary", "plot.twoCoprimary", "print.twoCoprimary_table",
            "dbibinom_g", ".dbibinom_g_r", ".tie_last")) {
  say("  ", f, " : ", exists(f, envir = ns, inherits = FALSE))
}
say("S3 methods registered:")
reg <- tryCatch(ls(get(".__S3MethodsTable__.", envir = ns)), error = function(e) character(0))
say("  ", if (length(reg)) paste(reg, collapse = ", ") else "(none)")

say("printing a result object to confirm the S3 method dispatches:")
res_probe <- power2Continuous(50, 50, 0.5, 0.5, 1, 1, 0.3, 0.025, known_var = TRUE)
out_probe <- capture.output(print(res_probe))
say("  print() produced ", length(out_probe), " lines, first line: ",
    if (length(out_probe)) out_probe[1] else "(none)")
say("")

g_r <- twoCoprimary:::.dbibinom_g_r
g_cpp <- twoCoprimary:::dbibinom_g

to_gamma <- function(p1, p2, rho) {
  z <- rho * sqrt(p2 * (1 - p2) / (p1 * (1 - p1)))
  z / (1 - z)
}

# ---------------------------------------------------------------------------
# Part A: compiled kernel against the R reference
# ---------------------------------------------------------------------------

say("--- Part A: compiled kernel against the R reference ---")

settings <- expand.grid(
  N = c(1, 2, 5, 20, 60, 150),
  p1 = c(0.2, 0.5, 0.8),
  p2 = c(0.3, 0.54),
  rho_frac = c(0.0, 0.4, 0.9),
  stringsAsFactors = FALSE
)

agr <- data.frame()
for (k in seq_len(nrow(settings))) {
  st <- settings[k, ]
  b <- corrbound2Binary(st$p1, st$p2)
  rho <- b[1] + st$rho_frac * (b[2] - b[1])
  gam <- to_gamma(st$p1, st$p2, rho)
  xi <- st$p2 + gam * (st$p2 - st$p1)

  y1 <- rep(0:st$N, each = st$N + 1)
  y2 <- rep(0:st$N, times = st$N + 1)

  v_r <- g_r(st$N, y1, y2, xi, gam)
  v_c <- g_cpp(as.integer(st$N), as.integer(y1), as.integer(y2), xi, gam)

  # A relative difference is only meaningful on entries that carry weight. The
  # summand alternates in sign whenever xi + gamma or 1 - xi + gamma is
  # negative, which happens for negative correlations, and entries many orders
  # of magnitude below the largest one then lose relative accuracy in both
  # implementations alike. Absolute agreement is what matters for a probability.
  d <- abs(v_c - v_r)
  big <- max(abs(v_r))
  material <- abs(v_r) > 1e-12 * big

  agr <- rbind(agr, data.frame(
    N = st$N, p1 = st$p1, p2 = st$p2, rho = round(rho, 4),
    n_pairs = length(y1),
    max_abs_diff = max(d),
    max_rel_diff_material = if (any(material)) max(d[material] / abs(v_r)[material]) else 0,
    n_material = sum(material),
    smallest_entry = min(abs(v_r)),
    alternating_signs = (xi + gam < 0) || (1 - xi + gam < 0)
  ))
}
write.csv(agr, file.path(out_dir, "rcpp_agreement.csv"), row.names = FALSE)
say("settings compared : ", nrow(agr))
say("max absolute difference over all settings          : ",
    format(max(agr$max_abs_diff), digits = 4))
say("max relative difference on entries carrying weight : ",
    format(max(agr$max_rel_diff_material), digits = 4))
say("settings whose summand alternates in sign          : ",
    sum(agr$alternating_signs), " of ", nrow(agr))
say("smallest entry seen anywhere                       : ",
    format(min(agr$smallest_entry), digits = 4))
say("")

# Scalar and edge inputs, where the R reference used to collapse to a vector
say("scalar and edge inputs:")
edges <- list(
  list(N = 10, y1 = 0L, y2 = 0L),
  list(N = 10, y1 = 10L, y2 = 10L),
  list(N = 10, y1 = 0L, y2 = 10L),
  list(N = 10, y1 = rep(0L, 5), y2 = 0:4),
  list(N = 1,  y1 = c(0L, 1L), y2 = c(1L, 0L))
)
gam <- to_gamma(0.4, 0.5, 0.2); xi <- 0.5 + gam * (0.5 - 0.4)
for (e in edges) {
  v_c <- g_cpp(as.integer(e$N), as.integer(e$y1), as.integer(e$y2), xi, gam)
  say("  N=", e$N, " y1=(", paste(e$y1, collapse = ","), ") y2=(",
      paste(e$y2, collapse = ","), ") -> ",
      paste(format(v_c, digits = 6), collapse = ", "))
}
say("")

# ---------------------------------------------------------------------------
# Part B: identities that any correct implementation must satisfy
# ---------------------------------------------------------------------------

say("--- Part B: the distribution sums to one and has binomial margins ---")

ident <- data.frame()
for (N in c(5, 20, 60)) {
  for (p1 in c(0.25, 0.54)) {
    for (p2 in c(0.4, 0.7)) {
      b <- corrbound2Binary(p1, p2)
      for (fr in c(0, 0.5, 0.95)) {
        rho <- b[1] + fr * (b[2] - b[1])
        pm <- outer(0:N, 0:N, function(x, y) dbibinom(N, x, y, p1, p2, rho))
        ident <- rbind(ident, data.frame(
          N = N, p1 = p1, p2 = p2, rho = round(rho, 4),
          total_minus_1 = sum(pm) - 1,
          max_margin1_error = max(abs(rowSums(pm) - dbinom(0:N, N, p1))),
          max_margin2_error = max(abs(colSums(pm) - dbinom(0:N, N, p2)))
        ))
      }
    }
  }
}
write.csv(ident, file.path(out_dir, "rcpp_identities.csv"), row.names = FALSE)
say("max |sum - 1|          : ", format(max(abs(ident$total_minus_1)), digits = 4))
say("max margin 1 deviation : ", format(max(ident$max_margin1_error), digits = 4))
say("max margin 2 deviation : ", format(max(ident$max_margin2_error), digits = 4))
say("")

# ---------------------------------------------------------------------------
# Part C: the article's sample sizes are unchanged
# ---------------------------------------------------------------------------

say("--- Part C: sample sizes of the article ---")

known <- data.frame(
  Test = c("Z-pool", "Z-pool", "Z-pool", "Boschloo", "Boschloo", "Boschloo",
           "Chisq", "Fisher"),
  rho = c(0, 0.3, 0.5, 0, 0.3, 0.5, 0.3, 0.3),
  N_expected = c(144, 142, 140, 144, 142, 140, 142, 150),
  stringsAsFactors = FALSE
)
known$N_now <- NA_integer_
for (k in seq_len(nrow(known))) {
  known$N_now[k] <- ss2BinaryExact(0.54, 0.54, 0.25, 0.25,
                                   known$rho[k], known$rho[k],
                                   1, 0.025, 0.1, known$Test[k])[["N"]]
  say("  ", known$Test[k], " rho=", known$rho[k],
      "  expected=", known$N_expected[k], " now=", known$N_now[k],
      if (known$N_now[k] == known$N_expected[k]) "" else "   *** MISMATCH ***")
}
write.csv(known, file.path(out_dir, "rcpp_samplesize.csv"), row.names = FALSE)
say("")
say("designs unchanged : ", sum(known$N_now == known$N_expected),
    " of ", nrow(known))
say("")

# ---------------------------------------------------------------------------
# Part D: timing
# ---------------------------------------------------------------------------

say("--- Part D: time to build one probability mass matrix ---")

# Repeat the fast calls so that the timer resolution does not dominate
reps_for <- function(N) if (N <= 60) 20L else if (N <= 120) 10L else 5L

tim <- data.frame()
for (N in c(50, 100, 150, 200)) {
  reps <- reps_for(N)
  gam <- to_gamma(0.54, 0.54, 0.3); xi <- 0.54 + gam * (0.54 - 0.54)
  y1 <- rep(0:N, each = N + 1); y2 <- rep(0:N, times = N + 1)

  t_c <- system.time(
    for (i in seq_len(reps)) {
      outer(0:N, 0:N, function(x, y) dbibinom(N, x, y, 0.54, 0.54, 0.3))
    }
  )[["elapsed"]] / reps

  t_r <- system.time(
    for (i in seq_len(reps)) g_r(N, y1, y2, xi, gam)
  )[["elapsed"]] / reps

  tim <- rbind(tim, data.frame(
    N = N, reps = reps,
    seconds_R = signif(t_r, 3), seconds_cpp = signif(t_c, 3),
    speedup = if (t_c > 0) round(t_r / t_c, 1) else NA_real_
  ))
}
write.csv(tim, file.path(out_dir, "rcpp_timing.csv"), row.names = FALSE)
print(tim)
capture.output(print(tim), file = log_con, append = TRUE)
say("")

say("--- End to end: one sample size search ---")
t_ss <- system.time(
  N_ss <- ss2BinaryExact(0.40, 0.35, 0.25, 0.20, 0.5, 0.5,
                         1, 0.025, 0.2, "Boschloo")[["N"]]
)[["elapsed"]]
say("usage example design, Boschloo: N = ", N_ss, " in ", round(t_ss, 1), " s")
say("(expected N = 368)")

close(log_con)
