# Does the tie handling in rr1Binary() change the rejection region?
#
# rr1Binary() orders the outcome cells from the most extreme to the least
# extreme and then takes a plain cumulative sum of the null cell probabilities
# along that order. Cells that share the same value of the ordering statistic
# form a tie group, and the tail event "at least as extreme as the observed
# cell" must include every member of that group. A plain cumulative sum instead
# gives the earlier members of a tie group a smaller tail probability than the
# later ones, and the order within a group is whatever order() happens to
# return.
#
# The consequence has a direction. The p-values currently assigned to the
# earlier members of a tie group are too small, so the rejection region can be
# too large and the actual size can exceed alpha.
#
# bbssr handles this explicitly: every member of a tie group receives the
# cumulative sum evaluated at the last member of that group. This script
# reimplements rr1Binary() with that single change and compares the two
# rejection regions. Only "Z-pool" and "Boschloo" are affected, since the other
# three tests read their p-values off a distribution and never order the cells.
#
# Run from the package root:
#   source("dev/check_tie_handling.R")
#
# Writes dev/out/tie_manuscript.csv, dev/out/tie_small.csv and
# dev/out/tie_handling.log

library(twoCoprimary)
library(fpCompare)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "tie_handling.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_con)
}

say("twoCoprimary version: ", as.character(utils::packageVersion("twoCoprimary")))
say("R version: ", R.version.string)
say("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# Position of the last member of the tie group of each element of a sorted
# vector. Relative tolerance, following bbssr, because exact p-values can be
# far below the absolute tolerance of fpCompare.
tie_last <- function(x) {
  n <- length(x)
  if (n <= 1) return(seq_len(n))
  tol <- 1e-10
  ref <- pmax(abs(x[-n]), abs(x[-1]))
  new_grp <- abs(diff(x)) > tol * ref
  grp_end <- c(which(new_grp), n)
  grp_size <- diff(c(0L, grp_end))
  rep(as.integer(grp_end), times = grp_size)
}

# Rejection regions under both tie conventions. Everything except the single
# line marked below is transcribed from rr1Binary().
rr_variants <- function(n1, n2, alpha, Test, grid = 100) {

  if (Test == "Z-pool") {
    Z_ij <- outer(n2 * (0:n1), n1 * (0:n2), "-") / (n1 * n2) /
      sqrt(outer(0:n1, 0:n2, "+") / (n1 * n2) *
             (1 - outer(0:n1, 0:n2, "+") / (n1 + n2)))
    Z_ij[is.na(Z_ij)] <- 0
    sel <- Z_ij %>>% 0
    stat <- Z_ij[sel]
    ord <- order(stat, decreasing = TRUE)
    stat_sorted <- stat[ord]
    base <- Z_ij
  } else if (Test == "Boschloo") {
    p_fisher <- outer(0:n1, 0:n2, function(i, j) {
      phyper(i - 1, n1, n2, i + j, lower.tail = FALSE)
    })
    p_max <- min(p_fisher[(outer(n2 * (0:n1), n1 * (0:n2), "-") / (n1 * n2)) %<=% 0])
    sel <- p_fisher %<<% p_max
    stat <- p_fisher[sel]
    ord <- order(stat, decreasing = FALSE)
    stat_sorted <- stat[ord]
    base <- p_fisher
  } else {
    stop("Test must be Z-pool or Boschloo")
  }

  i <- (row(base)[sel] - 1)[ord]
  j <- (col(base)[sel] - 1)[ord]

  uniq_i <- sort(unique(i))
  uniq_j <- sort(unique(j))
  theta <- seq(0, 1, length.out = grid)
  dbinom_i <- sapply(theta, function(t) dbinom(uniq_i, n1, t))
  dbinom_j <- sapply(theta, function(t) dbinom(uniq_j, n2, t))
  P_H0 <- dbinom_i[match(i, uniq_i), , drop = FALSE] *
    dbinom_j[match(j, uniq_j), , drop = FALSE]

  cum <- apply(P_H0, 2, cumsum)
  if (is.null(dim(cum))) cum <- matrix(cum, nrow = length(i))

  idx_last <- tie_last(stat_sorted)

  p_cur <- apply(cum, 1, max)                                # current
  p_fix <- apply(cum[idx_last, , drop = FALSE], 1, max)      # tie-aware

  mk <- function(p) {
    m <- array(1, dim(base))
    m[cbind(i + 1, j + 1)] <- p
    m %<<% alpha
  }

  list(
    RR_cur = mk(p_cur),
    RR_fix = mk(p_fix),
    n_pos = length(i),
    n_tied_cells = sum(idx_last != seq_along(idx_last)),
    max_tie_size = max(table(idx_last))
  )
}

# Actual size: maximum over the nuisance parameter of the null probability of
# the rejection region
actual_size <- function(RR, n1, n2, grid = 2001) {
  theta <- seq(0, 1, length.out = grid)
  RRn <- RR * 1
  vals <- vapply(theta, function(t) {
    as.numeric(dbinom(0:n1, n1, t) %*% RRn %*% dbinom(0:n2, n2, t))
  }, numeric(1))
  max(vals)
}

scan_designs <- function(designs, label) {
  res <- data.frame()
  for (k in seq_len(nrow(designs))) {
    n1 <- designs$n1[k]; n2 <- designs$n2[k]
    alpha <- designs$alpha[k]; Test <- designs$Test[k]
    v <- rr_variants(n1, n2, alpha, Test)
    n_diff <- sum(v$RR_cur != v$RR_fix)
    res <- rbind(res, data.frame(
      n1 = n1, n2 = n2, alpha = alpha, Test = Test,
      n_positive_cells = v$n_pos,
      n_cells_in_a_tie = v$n_tied_cells,
      max_tie_size = v$max_tie_size,
      n_cells_RR_differ = n_diff,
      size_current = if (n_diff > 0) actual_size(v$RR_cur, n1, n2) else NA_real_,
      size_tiefix  = if (n_diff > 0) actual_size(v$RR_fix, n1, n2) else NA_real_
    ))
    if (k %% 20 == 0) say("  ", label, ": ", k, " / ", nrow(designs), " done")
  }
  res
}

# ---------------------------------------------------------------------------
# Part A: the range of sample sizes the manuscript actually reaches
# ---------------------------------------------------------------------------
# Validation uses p11 = p12 = 0.54, p21 = p22 = 0.25, alpha = 0.025,
# 1 - beta = 0.90, with r = 1 and r = 2. The usage example uses
# p11 = 0.40, p12 = 0.35, p21 = 0.25, p22 = 0.20, alpha = 0.025,
# 1 - beta = 0.80, r = 1. Scanning every integer n over a band that comfortably
# contains those designs settles the question without running the sample size
# search: if no rejection region in the band changes, no number in the article
# can change.

d_r1 <- expand.grid(n2 = 20:150, alpha = 0.025,
                    Test = c("Z-pool", "Boschloo"), stringsAsFactors = FALSE)
d_r1$n1 <- d_r1$n2

d_r2 <- expand.grid(n2 = 20:100, alpha = 0.025,
                    Test = c("Z-pool", "Boschloo"), stringsAsFactors = FALSE)
d_r2$n1 <- 2 * d_r2$n2

designs_ms <- rbind(d_r1[, c("n1", "n2", "alpha", "Test")],
                    d_r2[, c("n1", "n2", "alpha", "Test")])

say("--- Part A: sample sizes relevant to the manuscript (",
    nrow(designs_ms), " designs) ---")
res_ms <- scan_designs(designs_ms, "Part A")
write.csv(res_ms, file.path(out_dir, "tie_manuscript.csv"), row.names = FALSE)

say("designs with at least one tied cell : ", sum(res_ms$n_cells_in_a_tie > 0))
say("designs where the rejection region differs : ",
    sum(res_ms$n_cells_RR_differ > 0))
if (any(res_ms$n_cells_RR_differ > 0)) {
  bad <- res_ms[res_ms$n_cells_RR_differ > 0, ]
  say("")
  say("Designs whose rejection region changes:")
  print(bad)
  capture.output(print(bad), file = log_con, append = TRUE)
  say("max actual size, current  : ", format(max(bad$size_current), digits = 6))
  say("max actual size, tie-fixed: ", format(max(bad$size_tiefix), digits = 6))
} else {
  say("")
  say("No rejection region in the scanned band changes.")
  say("The manuscript numbers therefore cannot move under the tie fix.")
}
say("")

# ---------------------------------------------------------------------------
# Part B: small samples, where ties are most likely to matter
# ---------------------------------------------------------------------------

d_small <- expand.grid(n2 = 5:40, alpha = c(0.025, 0.05),
                       Test = c("Z-pool", "Boschloo"), stringsAsFactors = FALSE)
d_small$n1 <- d_small$n2
d_small_unbal <- d_small
d_small_unbal$n1 <- 2 * d_small_unbal$n2
designs_small <- rbind(d_small[, c("n1", "n2", "alpha", "Test")],
                       d_small_unbal[, c("n1", "n2", "alpha", "Test")])

say("--- Part B: small samples (", nrow(designs_small), " designs) ---")
res_small <- scan_designs(designs_small, "Part B")
write.csv(res_small, file.path(out_dir, "tie_small.csv"), row.names = FALSE)

say("designs with at least one tied cell : ", sum(res_small$n_cells_in_a_tie > 0))
say("designs where the rejection region differs : ",
    sum(res_small$n_cells_RR_differ > 0))
if (any(res_small$n_cells_RR_differ > 0)) {
  bad <- res_small[res_small$n_cells_RR_differ > 0, ]
  say("")
  say("Worst cases by actual size of the current rejection region:")
  bad <- bad[order(-bad$size_current), ]
  print(utils::head(bad, 15))
  capture.output(print(utils::head(bad, 15)), file = log_con, append = TRUE)
  say("")
  say("designs where the current region exceeds alpha : ",
      sum(bad$size_current > bad$alpha))
  say("designs where the tie-fixed region exceeds alpha : ",
      sum(bad$size_tiefix > bad$alpha))
}

close(log_con)
