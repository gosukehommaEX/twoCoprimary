# Verification and benchmark of an algebraic reformulation of the co-primary
# power computation inside power2BinaryExact().
#
# The current implementation materialises two K by K matrices, where K is the
# number of cells in the rejection region RR. Since K grows like n ^ 2, the cost
# grows like n ^ 4 in both time and memory.
#
# Writing A for the row indices and C for the column indices of the cells in RR,
# the quantity being summed is
#
#   sum over l, m of pmass1[A[l], A[m]] pmass2[C[l], C[m]]
#
# which is the same as
#
#   sum over (a, c) in RR and (b, d) in RR of pmass1[a, b] pmass2[c, d]
#     = sum over a, b of pmass1[a, b] (RR pmass2 RR') [a, b]
#     = sum(pmass1 * (RR %*% pmass2 %*% t(RR)))
#
# The right hand side is two matrix products of order n + 1, evaluated by BLAS
# in O(n ^ 3) time and O(n ^ 2) memory.
#
# Run from the package root:
#   source("dev/bench_exact_power.R")
#
# Writes dev/out/exact_power_accuracy.csv, dev/out/exact_power_timing.csv and
# dev/out/exact_power.log

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "exact_power.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_con)
}

say("twoCoprimary version: ", as.character(utils::packageVersion("twoCoprimary")))
say("R version: ", R.version.string)
say("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# Current implementation, transcribed from power2BinaryExact()
power_current <- function(RR, pmass1, pmass2) {
  A <- row(RR)[RR]
  C <- col(RR)[RR]
  sum(t(pmass1[A, ])[A, ] * t(pmass2[C, ])[C, ])
}

# Proposed reformulation
power_matmul <- function(RR, pmass1, pmass2) {
  RRn <- RR * 1
  sum(pmass1 * (RRn %*% pmass2 %*% t(RRn)))
}

# Bivariate binomial probability mass matrix for one group
build_pmass <- function(n, p1, p2, rho) {
  outer(0:n, 0:n, function(x, y) dbibinom(n, x, y, p1, p2, rho))
}

# Pick a correlation a fixed fraction of the way into the admissible interval
pick_rho <- function(p1, p2, frac) {
  b <- corrbound2Binary(p1, p2)
  b[1] + frac * (b[2] - b[1])
}

tests <- c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")

# ---------------------------------------------------------------------------
# Part A: numerical agreement between the two expressions
# ---------------------------------------------------------------------------

settings <- expand.grid(
  n = c(15, 25, 40),
  Test = tests,
  stringsAsFactors = FALSE
)

acc <- data.frame()
for (k in seq_len(nrow(settings))) {
  n <- settings$n[k]
  Test <- settings$Test[k]

  p11 <- 0.55; p12 <- 0.65; p21 <- 0.30; p22 <- 0.40
  rho1 <- pick_rho(p11, p12, 0.6)
  rho2 <- pick_rho(p21, p22, 0.6)

  RR <- rr1Binary(n, n, 0.025, Test)
  pmass1 <- build_pmass(n, p11, p12, rho1)
  pmass2 <- build_pmass(n, p21, p22, rho2)

  v_cur <- power_current(RR, pmass1, pmass2)
  v_new <- power_matmul(RR, pmass1, pmass2)

  # Independent check against the exported function
  v_pkg <- power2BinaryExact(n, n, p11, p12, p21, p22, rho1, rho2,
                             0.025, Test)$powerCoprimary

  acc <- rbind(acc, data.frame(
    n = n, Test = Test, K = sum(RR),
    current = v_cur, matmul = v_new, package = v_pkg,
    abs_diff_new_vs_current = abs(v_new - v_cur),
    abs_diff_current_vs_package = abs(v_cur - v_pkg)
  ))
}

write.csv(acc, file.path(out_dir, "exact_power_accuracy.csv"), row.names = FALSE)

say("--- Part A: numerical agreement ---")
say("max |matmul - current|  = ", format(max(acc$abs_diff_new_vs_current), digits = 3))
say("max |current - package| = ", format(max(acc$abs_diff_current_vs_package), digits = 3))
say("")
print(acc)
capture.output(print(acc), file = log_con, append = TRUE)
say("")

# ---------------------------------------------------------------------------
# Part B: timing as a function of n
# ---------------------------------------------------------------------------
# The current expression allocates two K by K double matrices. Guard against
# exhausting memory by skipping it once K exceeds k_cap.

k_cap <- 6000
ns <- c(20, 40, 60, 80, 100, 150, 200)

tim <- data.frame()
for (n in ns) {
  for (Test in c("Fisher", "Boschloo")) {

    p11 <- 0.55; p12 <- 0.65; p21 <- 0.30; p22 <- 0.40
    rho1 <- pick_rho(p11, p12, 0.6)
    rho2 <- pick_rho(p21, p22, 0.6)

    t_rr <- system.time(RR <- rr1Binary(n, n, 0.025, Test))[["elapsed"]]
    t_pm <- system.time({
      pmass1 <- build_pmass(n, p11, p12, rho1)
      pmass2 <- build_pmass(n, p21, p22, rho2)
    })[["elapsed"]]

    K <- sum(RR)
    t_cur <- NA_real_
    if (K <= k_cap) {
      t_cur <- system.time(power_current(RR, pmass1, pmass2))[["elapsed"]]
    }
    t_new <- system.time(power_matmul(RR, pmass1, pmass2))[["elapsed"]]

    tim <- rbind(tim, data.frame(
      n = n, Test = Test, K = K,
      sec_rr1Binary = t_rr,
      sec_pmass = t_pm,
      sec_power_current = t_cur,
      sec_power_matmul = t_new,
      gb_current_alloc = 2 * K ^ 2 * 8 / 1024 ^ 3
    ))
    rm(RR, pmass1, pmass2)
    gc()
  }
}

write.csv(tim, file.path(out_dir, "exact_power_timing.csv"), row.names = FALSE)

say("--- Part B: timing (seconds) ---")
print(tim)
capture.output(print(tim), file = log_con, append = TRUE)
say("")
say("NA in sec_power_current means K exceeded the guard of ", k_cap,
    " cells and the current expression was not run.")

close(log_con)
