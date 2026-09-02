# Reproduce the published table of
#   Sozu T, Sugimoto T, Hamasaki T (2011). Sample size determination in
#   superiority clinical trials with multiple co-primary correlated endpoints.
#   Journal of Biopharmaceutical Statistics 21(4), 650-668.
#
# Table 1 gives the sample size per group for two co-primary continuous
# endpoints with a balanced design, one-sided alpha = 0.025, power 0.80 and
# known variance, over five standardized effect sizes and four correlations.
# The same table also gives, in its two rightmost columns, the sample size that
# each endpoint would need on its own, which is what ss1Continuous computes.
#
# The article does not state which software produced the table. The values are
# nevertheless reproduced exactly, so the tolerance here is zero and any
# difference is a defect rather than a numerical difference between platforms.
#
# Run with:  source("dev/reproduce_Sozu_et_al_2011.R")

library(twoCoprimary)
source("dev/reproduce_helpers.R")

reproduce_begin(
  stem = "reproduce_Sozu_et_al_2011",
  article = "Sozu, Sugimoto and Hamasaki (2011), J Biopharm Stat 21(4), 650-668",
  tables = c("Table 1 (two continuous endpoints, n per group)",
             "Table 1 columns E1 and E2 (single endpoint)")
)

ALPHA <- 0.025
BETA <- 0.2
RHOS <- c(0.0, 0.3, 0.5, 0.8)

# ---------------------------------------------------------------------------
# Table 1, as printed on page 656. Sample size PER GROUP, n = n1 = n2.
# Columns: the four correlations, then the size for each effect size alone.
# ---------------------------------------------------------------------------
published <- data.frame(
  d1 = c(0.20, 0.20, 0.20, 0.20, 0.20, 0.25, 0.25, 0.25, 0.25,
         0.30, 0.30, 0.30, 0.35, 0.35, 0.40),
  d2 = c(0.20, 0.25, 0.30, 0.35, 0.40, 0.25, 0.30, 0.35, 0.40,
         0.30, 0.35, 0.40, 0.35, 0.40, 0.40),
  rho0.0 = c(516, 432, 402, 394, 393, 330, 284, 263, 254, 230, 201, 186, 169, 150, 129),
  rho0.3 = c(503, 424, 399, 394, 393, 322, 278, 260, 253, 224, 197, 183, 165, 147, 126),
  rho0.5 = c(490, 417, 397, 393, 393, 314, 272, 257, 253, 218, 192, 181, 160, 143, 123),
  rho0.8 = c(458, 401, 393, 393, 393, 294, 260, 253, 252, 204, 183, 176, 150, 136, 115),
  E1     = c(393, 393, 393, 393, 393, 252, 252, 252, 252, 175, 175, 175, 129, 129,  99),
  E2     = c(393, 252, 175, 129,  99, 252, 175, 129,  99, 175, 129,  99, 129,  99,  99)
)

for (i in seq_len(nrow(published))) {
  for (j in seq_along(RHOS)) {
    rho <- RHOS[j]
    res <- ss2Continuous(delta1 = published$d1[i], delta2 = published$d2[i],
                         sd1 = 1, sd2 = 1, rho = rho, r = 1,
                         alpha = ALPHA, beta = BETA, known_var = TRUE)
    cmp("Sozu 2011 Table 1", 
        sprintf("delta* = (%.2f, %.2f), rho = %.1f",
                published$d1[i], published$d2[i], rho),
        published[[paste0("rho", format(rho, nsmall = 1))]][i],
        res$n2, tol = 0)
  }
}
table_report("Sozu 2011 Table 1",
             "n per group, balanced design, known variance, alpha = 0.025, power 0.80")

# ---------------------------------------------------------------------------
# Columns E1 and E2 of the same table: the size each endpoint needs alone
# ---------------------------------------------------------------------------
for (d in sort(unique(c(published$d1, published$d2)))) {
  want <- unique(c(published$E1[published$d1 == d], published$E2[published$d2 == d]))
  stopifnot(length(want) == 1)
  got <- ss1Continuous(delta = d, sd = 1, r = 1, alpha = ALPHA, beta = BETA)$n2
  cmp("Sozu 2011 Table 1, E1 and E2", sprintf("delta* = %.2f", d), want, got, tol = 0)
}
table_report("Sozu 2011 Table 1, E1 and E2",
             "single endpoint size, the two rightmost columns of the table")

# ---------------------------------------------------------------------------
# The article's own reading of the table, which a reproduction should also
# reproduce: for equal effect sizes the size falls by about 11 per cent between
# no correlation and a correlation of 0.8
# ---------------------------------------------------------------------------
eq <- published$d1 == published$d2
drop_pct <- 100 * (published$rho0.0[eq] - published$rho0.8[eq]) / published$rho0.0[eq]
got_pct <- vapply(which(eq), function(i) {
  a <- ss2Continuous(published$d1[i], published$d2[i], 1, 1, 0.0, 1, ALPHA, BETA)$n2
  b <- ss2Continuous(published$d1[i], published$d2[i], 1, 1, 0.8, 1, ALPHA, BETA)$n2
  100 * (a - b) / a
}, numeric(1))
cmp("Sozu 2011 text", "mean reduction from rho = 0 to rho = 0.8, equal effect sizes",
    round(mean(drop_pct), 1), round(mean(got_pct), 1), tol = 0, quantity = "per cent")
table_report("Sozu 2011 text",
             "the article states this reduction is approximately 11 per cent")

# ---------------------------------------------------------------------------
# The same numbers through design_table, so that the two routes agree
# ---------------------------------------------------------------------------
tab <- design_table(
  param_grid = data.frame(delta1 = published$d1, delta2 = published$d2,
                          sd1 = 1, sd2 = 1),
  rho_values = RHOS, r = 1, alpha = ALPHA, beta = BETA,
  endpoint_type = "continuous"
)
for (i in seq_len(nrow(published))) {
  for (rho in RHOS) {
    cmp("Sozu 2011 Table 1 via design_table",
        sprintf("delta* = (%.2f, %.2f), rho = %.1f",
                published$d1[i], published$d2[i], rho),
        published[[paste0("rho", format(rho, nsmall = 1))]][i],
        tab[[paste0("rho_", format(rho, nsmall = 1))]][i] / 2, tol = 0)
  }
}
table_report("Sozu 2011 Table 1 via design_table",
             "design_table returns the total, which is twice the published per-group size")

reproduce_end()
