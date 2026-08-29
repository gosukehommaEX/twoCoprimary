# Does the rr1Binary() tie fix change any published result of
# Homma and Yoshida (2025), Statistical Methods in Medical Research?
#
# The tie defect affects only the two exact unconditional tests, so the parts of
# that article at risk are the Z-pool and Boschloo entries of Table 4 and the
# Z-pool and Boschloo curves of Figures 1 to 3. Table 3 reports the chi-squared
# test only and cannot be affected.
#
# For each design the script computes three things:
#   published  the value printed in the article
#   fixed      the value from the installed package, which now uses the
#              tie-aware convention
#   buggy      the value from a reimplementation of the pre-fix convention
#
# The article is safe when published, fixed and buggy all agree. If published
# and fixed disagree anywhere, that entry of the article is wrong and needs to
# be dealt with.
#
# Part 1 covers Table 4 and takes a few minutes. Part 2 covers the three
# scenario cases behind Figures 1 to 3 and takes considerably longer. Set PARTS
# below to run them separately. Results are written after each part, so
# interrupting during Part 2 still leaves Part 1 on disk.
#
# Run from the package root:
#   source("dev/check_hy2025.R")
#
# Writes dev/out/hy2025_table4.csv, dev/out/hy2025_figures.csv and
# dev/out/hy2025.log

PARTS <- c(1, 2)

library(twoCoprimary)
library(fpCompare)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "hy2025.log"), open = "wt")
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
  rr1Binary(5, 5, 0.025, "Z-pool")
  TRUE
}, error = function(e) {
  message("Cannot call rr1Binary(): ", conditionMessage(e))
  message("Reinstall the package, restart the R session and run this again.")
  FALSE
})
if (!probe) {
  close(log_con)
  stop("stale package state; reinstall, restart R and re-run")
}

# --- the pre-fix convention, reimplemented ---------------------------------

rr_buggy <- function(n1, n2, alpha, Test, grid = 100) {

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
    stop("rr_buggy() only applies to Z-pool and Boschloo")
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

  m <- array(1, dim(base))
  m[cbind(i + 1, j + 1)] <- apply(cum, 1, max)
  m %<<% alpha
}

# --- search under a supplied rejection region rule -------------------------

make_engine <- function(p11, p12, p21, p22, rho1, rho2, alpha, Test) {

  pmass_cache <- new.env(parent = emptyenv())
  pmass_of <- function(n, p1, p2, rho, tag) {
    key <- paste0(tag, "_", n)
    if (!exists(key, envir = pmass_cache, inherits = FALSE)) {
      assign(key, outer(0:n, 0:n, function(x, y) dbibinom(n, x, y, p1, p2, rho)),
             envir = pmass_cache)
    }
    get(key, envir = pmass_cache, inherits = FALSE)
  }

  rr_cache <- new.env(parent = emptyenv())
  rr_of <- function(n1, n2, which) {
    key <- paste0(which, "_", n1, "_", n2)
    if (!exists(key, envir = rr_cache, inherits = FALSE)) {
      val <- if (which == "fixed") rr1Binary(n1, n2, alpha, Test)
             else rr_buggy(n1, n2, alpha, Test)
      assign(key, val, envir = rr_cache)
    }
    get(key, envir = rr_cache, inherits = FALSE)
  }

  power_of <- function(n1, n2, which) {
    RR <- rr_of(n1, n2, which) * 1
    pmass1 <- pmass_of(n1, p11, p12, rho1, "g1")
    pmass2 <- pmass_of(n2, p21, p22, rho2, "g2")
    sum(pmass1 * (RR %*% pmass2 %*% t(RR)))
  }

  list(power = power_of)
}

search_n <- function(engine, which, initial_n2, r, target_power) {
  n2 <- initial_n2
  n1 <- ceiling(r * n2)
  power <- engine$power(n1, n2, which)
  if (power %>=% target_power) {
    while (power %>=% target_power && n2 > 1) {
      n2 <- n2 - 1
      n1 <- ceiling(r * n2)
      power <- engine$power(n1, n2, which)
    }
    n2 <- n2 + 1
  } else {
    while (power %<<% target_power) {
      n2 <- n2 + 1
      n1 <- ceiling(r * n2)
      power <- engine$power(n1, n2, which)
    }
  }
  n1 <- ceiling(r * n2)
  list(n1 = n1, n2 = n2, N = n1 + n2, power = engine$power(n1, n2, which))
}

run_designs <- function(designs, csv_name, label) {
  res <- data.frame()
  for (k in seq_len(nrow(designs))) {
    d <- designs[k, ]
    eng <- make_engine(d$p11, d$p12, d$p21, d$p22, d$rho, d$rho, d$alpha, d$Test)
    init <- ss2BinaryApprox(d$p11, d$p12, d$p21, d$p22, d$rho, d$rho,
                            d$r, d$alpha, d$beta, "AN")[["n2"]]
    a <- search_n(eng, "fixed", init, d$r, 1 - d$beta)
    b <- search_n(eng, "buggy", init, d$r, 1 - d$beta)

    row <- cbind(d, data.frame(
      N_fixed = a$N, N_buggy = b$N,
      power_fixed = round(a$power, 4), power_buggy = round(b$power, 4),
      fixed_vs_buggy = ifelse(a$N == b$N, "same", "DIFFERENT")
    ))
    if (!is.null(d$N_published) && !is.na(d$N_published)) {
      row$fixed_vs_published <- ifelse(a$N == d$N_published, "same", "DIFFERENT")
    }
    res <- rbind(res, row)

    say("  ", k, " / ", nrow(designs), ": ", label, " ", d$Test,
        " alpha=", d$alpha, " r=", d$r, " rho=", d$rho,
        "  fixed=", a$N, " buggy=", b$N,
        if (!is.null(d$N_published) && !is.na(d$N_published))
          paste0(" published=", d$N_published) else "",
        if (a$N == b$N &&
            (is.null(d$N_published) || is.na(d$N_published) || a$N == d$N_published))
          "" else "   *** MISMATCH ***")

    write.csv(res, file.path(out_dir, csv_name), row.names = FALSE)
  }
  res
}

# ---------------------------------------------------------------------------
# Part 1: Table 4
# ---------------------------------------------------------------------------
# p11 = p12 = 0.54, p21 = p22 = 0.25, 1 - beta = 0.9

if (1 %in% PARTS) {

  say("--- Part 1: Table 4 of Homma and Yoshida (2025) ---")

  t4 <- expand.grid(rho = c(0, 0.3, 0.5, 0.8), r = c(1, 2),
                    alpha = c(0.025, 0.05, 0.1),
                    Test = c("Z-pool", "Boschloo"), stringsAsFactors = FALSE)
  t4$p11 <- 0.54; t4$p12 <- 0.54; t4$p21 <- 0.25; t4$p22 <- 0.25; t4$beta <- 0.1

  # Values printed in Table 4, in the order produced by expand.grid above
  t4$N_published <- c(
    144, 142, 140, 134, 180, 180, 177, 168,   # Z-pool  alpha = 0.025, r = 1 then 2
    120, 116, 114, 110, 147, 144, 138, 135,   # Z-pool  alpha = 0.05
    104,  98,  96,  92, 114, 114, 111, 108,   # Z-pool  alpha = 0.1
    144, 142, 140, 134, 162, 159, 156, 150,   # Boschloo alpha = 0.025
    120, 116, 114, 110, 135, 132, 132, 126,   # Boschloo alpha = 0.05
     98,  96,  94,  88, 108, 105, 105,  99    # Boschloo alpha = 0.1
  )
  t4$power_published <- c(
    0.905, 0.903, 0.902, 0.901, 0.901, 0.906, 0.903, 0.908,
    0.903, 0.900, 0.902, 0.904, 0.904, 0.903, 0.902, 0.907,
    0.902, 0.901, 0.905, 0.905, 0.904, 0.909, 0.904, 0.908,
    0.905, 0.903, 0.902, 0.901, 0.905, 0.902, 0.901, 0.903,
    0.903, 0.900, 0.902, 0.904, 0.904, 0.901, 0.907, 0.907,
    0.907, 0.904, 0.902, 0.903, 0.908, 0.902, 0.908, 0.906
  )

  res4 <- run_designs(t4, "hy2025_table4.csv", "Table4")

  say("")
  say("Table 4, entries where the fix changes the sample size : ",
      sum(res4$fixed_vs_buggy != "same"), " of ", nrow(res4))
  say("Table 4, entries where the package disagrees with the article : ",
      sum(res4$fixed_vs_published != "same"), " of ", nrow(res4))
  if (any(res4$fixed_vs_published != "same")) {
    bad <- res4[res4$fixed_vs_published != "same",
                c("Test", "alpha", "r", "rho", "N_published", "N_fixed", "N_buggy")]
    print(bad)
    capture.output(print(bad), file = log_con, append = TRUE)
  }
  say("")
}

# ---------------------------------------------------------------------------
# Part 2: the scenarios behind Figures 1 to 3
# ---------------------------------------------------------------------------
# alpha = 0.025, 1 - beta = 0.8, r = 1 and 2, rho = 0 to 0.8 by 0.1.
# The figures show curves rather than printed numbers, so there is no published
# value to compare against. What matters is whether the fix moves the values.

if (2 %in% PARTS) {

  say("--- Part 2: scenario cases behind Figures 1 to 3 ---")

  cases <- rbind(
    data.frame(case = "A", p11 = 0.4, p21 = 0.2, p12 = 0.4, p22 = 0.2),
    data.frame(case = "A", p11 = 0.4, p21 = 0.2, p12 = 0.3, p22 = 0.1),
    data.frame(case = "B", p11 = 0.7, p21 = 0.4, p12 = 0.7, p22 = 0.4),
    data.frame(case = "B", p11 = 0.7, p21 = 0.4, p12 = 0.5, p22 = 0.2),
    data.frame(case = "C", p11 = 0.8, p21 = 0.5, p12 = 0.8, p22 = 0.6),
    data.frame(case = "C", p11 = 0.8, p21 = 0.5, p12 = 0.8, p22 = 0.4)
  )

  fig <- data.frame()
  for (k in seq_len(nrow(cases))) {
    cs <- cases[k, ]
    b1 <- corrbound2Binary(cs$p11, cs$p12)
    b2 <- corrbound2Binary(cs$p21, cs$p22)
    rho_ok <- seq(0, 0.8, by = 0.1)
    rho_ok <- rho_ok[rho_ok <= min(b1[2], b2[2]) & rho_ok >= max(b1[1], b2[1])]
    if (length(rho_ok) == 0) next
    fig <- rbind(fig, expand.grid(
      case = cs$case, p11 = cs$p11, p12 = cs$p12, p21 = cs$p21, p22 = cs$p22,
      rho = rho_ok, r = c(1, 2), alpha = 0.025, beta = 0.2,
      Test = c("Z-pool", "Boschloo"), stringsAsFactors = FALSE
    ))
  }
  fig$N_published <- NA_integer_

  say("designs to run: ", nrow(fig))
  resf <- run_designs(fig, "hy2025_figures.csv", "Figures")

  say("")
  say("Figures 1 to 3, designs where the fix changes the sample size : ",
      sum(resf$fixed_vs_buggy != "same"), " of ", nrow(resf))
  if (any(resf$fixed_vs_buggy != "same")) {
    bad <- resf[resf$fixed_vs_buggy != "same",
                c("case", "p12", "p22", "Test", "r", "rho", "N_fixed", "N_buggy")]
    print(bad)
    capture.output(print(bad), file = log_con, append = TRUE)
  }
}

close(log_con)
