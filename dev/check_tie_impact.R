# Does the tie fix move the sample sizes reported in the manuscript, and which
# tie convention agrees with independent software?
#
# check_tie_handling.R established that the rejection region changes in 51 of
# the 424 designs in the band the manuscript reaches, including n = 58 and
# n = 59 per group, which is close to where the validation designs land. A
# changed rejection region does not by itself move a sample size: the search
# stops at the first n whose power reaches the target, and a one cell change
# shifts the power by the probability of that cell under the alternative. This
# script settles the question by running the search itself under both
# conventions.
#
# Part B compares the two conventions against the Exact package, which computes
# the tail probability of the ordering statistic directly and therefore includes
# whole tie groups by construction. Whichever convention agrees with it is the
# standard Z-pooled and Boschloo test.
#
# Run from the package root:
#   source("dev/check_tie_impact.R")
#
# Writes dev/out/tie_impact_samplesize.csv, dev/out/tie_impact_external.csv and
# dev/out/tie_impact.log

library(twoCoprimary)
library(fpCompare)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "tie_impact.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_con)
}

say("twoCoprimary version: ", as.character(utils::packageVersion("twoCoprimary")))
say("R version: ", R.version.string)
say("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# Reinstalling a package that is already attached leaves the running session
# with a stale lazy-load database. Fail early with a clear instruction rather
# than part way through the search.
probe <- tryCatch({
  dbibinom(4, 0:2, 0:2, 0.4, 0.5, 0.1)
  TRUE
}, error = function(e) {
  message("Cannot call dbibinom(): ", conditionMessage(e))
  message("Restart the R session (Session > Restart R, or Ctrl+Shift+F10) and ",
          "run this script again.")
  FALSE
})
if (!probe) {
  close(log_con)
  stop("stale package state; restart R and re-run")
}

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

# p-value matrices under both tie conventions
pvals_both <- function(n1, n2, Test, grid = 100) {

  if (Test == "Z-pool") {
    Z_ij <- outer(n2 * (0:n1), n1 * (0:n2), "-") / (n1 * n2) /
      sqrt(outer(0:n1, 0:n2, "+") / (n1 * n2) *
             (1 - outer(0:n1, 0:n2, "+") / (n1 + n2)))
    Z_ij[is.na(Z_ij)] <- 0
    sel <- Z_ij %>>% 0
    stat <- Z_ij[sel]
    ord <- order(stat, decreasing = TRUE)
    base <- Z_ij
  } else if (Test == "Boschloo") {
    p_fisher <- outer(0:n1, 0:n2, function(i, j) {
      phyper(i - 1, n1, n2, i + j, lower.tail = FALSE)
    })
    p_max <- min(p_fisher[(outer(n2 * (0:n1), n1 * (0:n2), "-") / (n1 * n2)) %<=% 0])
    sel <- p_fisher %<<% p_max
    stat <- p_fisher[sel]
    ord <- order(stat, decreasing = FALSE)
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

  idx_last <- tie_last(stat[ord])
  mk <- function(p) {
    m <- array(1, dim(base))
    m[cbind(i + 1, j + 1)] <- p
    m
  }
  list(current = mk(apply(cum, 1, max)),
       tiefix  = mk(apply(cum[idx_last, , drop = FALSE], 1, max)))
}

# Co-primary power for a given rejection region, using the matrix product form
power_given_rr <- function(RR, n1, n2, p11, p12, p21, p22, rho1, rho2) {
  pmass1 <- outer(0:n1, 0:n1, function(x, y) dbibinom(n1, x, y, p11, p12, rho1))
  pmass2 <- outer(0:n2, 0:n2, function(x, y) dbibinom(n2, x, y, p21, p22, rho2))
  RRn <- RR * 1
  sum(pmass1 * (RRn %*% pmass2 %*% t(RRn)))
}

# Sequential search, mirroring .ss_sequential_search()
ss_search <- function(initial_n2, r, target_power, Test, alpha, convention,
                      p11, p12, p21, p22, rho1, rho2) {
  pow <- function(n1, n2) {
    RR <- pvals_both(n1, n2, Test)[[convention]] %<<% alpha
    power_given_rr(RR, n1, n2, p11, p12, p21, p22, rho1, rho2)
  }
  n2 <- initial_n2
  n1 <- ceiling(r * n2)
  power <- pow(n1, n2)
  if (power %>=% target_power) {
    while (power %>=% target_power && n2 > 1) {
      n2 <- n2 - 1
      n1 <- ceiling(r * n2)
      power <- pow(n1, n2)
    }
    n2 <- n2 + 1
  } else {
    while (power %<<% target_power) {
      n2 <- n2 + 1
      n1 <- ceiling(r * n2)
      power <- pow(n1, n2)
    }
  }
  n1 <- ceiling(r * n2)
  list(n1 = n1, n2 = n2, N = n1 + n2, power = pow(n1, n2))
}

# ---------------------------------------------------------------------------
# Part A: the designs the manuscript actually reports
# ---------------------------------------------------------------------------
# Validation, Homma and Yoshida (2025) Table 4:
#   p11 = p12 = 0.54, p21 = p22 = 0.25, alpha = 0.025, beta = 0.1, r = 1 and 2
# Usage example:
#   p11 = 0.40, p12 = 0.35, p21 = 0.25, p22 = 0.20, rho = 0.5,
#   alpha = 0.025, beta = 0.2, r = 1
# Only Z-pool and Boschloo can change; the other three tests never order the
# cells and are therefore untouched.

designs <- rbind(
  expand.grid(source = "validation", p11 = 0.54, p12 = 0.54,
              p21 = 0.25, p22 = 0.25, rho = c(0, 0.3, 0.5, 0.8),
              r = c(1, 2), alpha = 0.025, beta = 0.1,
              Test = c("Z-pool", "Boschloo"), stringsAsFactors = FALSE),
  expand.grid(source = "usage", p11 = 0.40, p12 = 0.35,
              p21 = 0.25, p22 = 0.20, rho = 0.5,
              r = 1, alpha = 0.025, beta = 0.2,
              Test = c("Z-pool", "Boschloo"), stringsAsFactors = FALSE)
)

say("--- Part A: sample sizes under both tie conventions (",
    nrow(designs), " designs) ---")

res <- data.frame()
for (k in seq_len(nrow(designs))) {
  d <- designs[k, ]
  init <- ss2BinaryApprox(d$p11, d$p12, d$p21, d$p22, d$rho, d$rho,
                          d$r, d$alpha, d$beta, "AN")[["n2"]]

  a <- ss_search(init, d$r, 1 - d$beta, d$Test, d$alpha, "current",
                 d$p11, d$p12, d$p21, d$p22, d$rho, d$rho)
  b <- ss_search(init, d$r, 1 - d$beta, d$Test, d$alpha, "tiefix",
                 d$p11, d$p12, d$p21, d$p22, d$rho, d$rho)

  # Cross-check the current convention against the shipped function
  pkg <- ss2BinaryExact(d$p11, d$p12, d$p21, d$p22, d$rho, d$rho,
                        d$r, d$alpha, d$beta, d$Test)[["N"]]

  res <- rbind(res, data.frame(
    source = d$source, Test = d$Test, rho = d$rho, r = d$r,
    N_current = a$N, N_tiefix = b$N, N_package = pkg,
    n2_current = a$n2, n2_tiefix = b$n2,
    power_current = a$power, power_tiefix = b$power,
    diff_N = b$N - a$N,
    current_matches_package = (a$N == pkg)
  ))
  say("  ", k, " / ", nrow(designs), ": ", d$source, " ", d$Test,
      " rho=", d$rho, " r=", d$r,
      "  N current=", a$N, " tiefix=", b$N, " package=", pkg)
}

write.csv(res, file.path(out_dir, "tie_impact_samplesize.csv"), row.names = FALSE)

say("")
say("designs where the sample size changes : ", sum(res$diff_N != 0), " of ", nrow(res))
if (any(res$diff_N != 0)) {
  print(res[res$diff_N != 0, ])
  capture.output(print(res[res$diff_N != 0, ]), file = log_con, append = TRUE)
}
say("current convention reproduces ss2BinaryExact() in all designs : ",
    all(res$current_matches_package))
say("")

# ---------------------------------------------------------------------------
# Part B: which convention agrees with the Exact package
# ---------------------------------------------------------------------------

if (!requireNamespace("Exact", quietly = TRUE)) {
  say("--- Part B skipped: package 'Exact' is not installed ---")
  say("Install it with install.packages('Exact') and re-run to settle which")
  say("tie convention matches independent software.")
} else {
  say("--- Part B: comparison against Exact ", 
      as.character(utils::packageVersion("Exact")), " ---")

  n1 <- 15; n2 <- 15
  pv <- pvals_both(n1, n2, "Z-pool")

  n_error <- 0L
  first_error <- NULL
  ext <- data.frame()
  for (x1 in 0:n1) {
    for (x2 in 0:n2) {
      tab <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2, byrow = TRUE)
      # ref.pvalue = FALSE keeps Exact on the same fixed grid of 100 nuisance
      # values that rr1Binary() uses. With the default TRUE it refines the
      # maximum on a finer grid and the two are no longer comparable.
      p_ex <- tryCatch(
        Exact::exact.test(tab, alternative = "greater", method = "z-pooled",
                          np.interval = FALSE, npNumbers = 100,
                          ref.pvalue = FALSE, to.plot = FALSE)$p.value,
        error = function(e) {
          if (is.null(first_error)) first_error <<- conditionMessage(e)
          n_error <<- n_error + 1L
          NA_real_
        }
      )
      ext <- rbind(ext, data.frame(
        x1 = x1, x2 = x2,
        p_current = pv$current[x1 + 1, x2 + 1],
        p_tiefix  = pv$tiefix[x1 + 1, x2 + 1],
        p_Exact   = p_ex
      ))
    }
  }
  ext$d_current <- abs(ext$p_current - ext$p_Exact)
  ext$d_tiefix  <- abs(ext$p_tiefix  - ext$p_Exact)
  write.csv(ext, file.path(out_dir, "tie_impact_external.csv"), row.names = FALSE)

  if (n_error > 0) {
    say("Exact::exact.test() failed on ", n_error, " of ", nrow(ext), " cells.")
    say("first error: ", first_error)
    say("If this mentions ExactData, install it with")
    say("  install.packages('ExactData', repos='https://pcalhoun1.github.io/drat/', type='source')")
    say("and re-run.")
    say("")
  }

  ok <- !is.na(ext$p_Exact) & ext$p_Exact < 1
  say("cells compared (p_Exact < 1) : ", sum(ok))
  say("max |current - Exact| : ", format(max(ext$d_current[ok]), digits = 4))
  say("max |tiefix  - Exact| : ", format(max(ext$d_tiefix[ok]),  digits = 4))
  say("cells where current differs from Exact by more than 1e-8 : ",
      sum(ext$d_current[ok] > 1e-8))
  say("cells where tiefix  differs from Exact by more than 1e-8 : ",
      sum(ext$d_tiefix[ok] > 1e-8))
  say("")
  if (sum(ok) > 0) {
    say("The convention with the smaller discrepancy is the standard test.")
    sub <- ext[ok, ]
    worst <- utils::head(sub[order(-sub$d_current), ], 10)
    print(worst)
    capture.output(print(worst), file = log_con, append = TRUE)
  } else {
    say("No cell could be compared. Part B is inconclusive.")
  }
}

close(log_con)
