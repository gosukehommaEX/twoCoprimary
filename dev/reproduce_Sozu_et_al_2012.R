# Reproduce the published tables of
#   Sozu T, Sugimoto T, Hamasaki T (2012). Sample size determination in clinical
#   trials with multiple co-primary endpoints including mixed continuous and
#   binary variables. Biometrical Journal 54(5), 716-729,
# and of its Supporting Information.
#
# Three tables are covered.
#
#   Table 2, the PREMIER illustration. Every input is stated exactly, so this is
#   the table to read first. Tolerance one subject, the article having been
#   computed in SAS.
#
#   Table 1, twenty four settings over three ratios of the two individual
#   powers. The standardized effect size of the continuous endpoint is printed
#   rounded to three decimal places, and was chosen to make the two individual
#   powers stand in the stated ratio, so the size cannot be reproduced to the
#   subject from the printed input. The tolerance is one per cent of the
#   published size, which is what the third decimal place of the effect size is
#   worth. The E2 column, the size the binary endpoint needs on its own, has
#   exact inputs and is checked exactly.
#
#   Table S5 of the Supporting Information, at proportions between 0.80 and
#   0.99, where the asymptotic methods are least reliable. The columns
#   Normal(CC) and Arcsine(CC) are checked here. The Fisher column of that table
#   is checked in dev/verify_mixed_fisher_option.R, which compares the power at
#   the published sample size rather than re-running a Monte Carlo search.
#
# The article parameterizes the allocation as kappa = n2/n1 and all three tables
# are balanced, so r = 1 is used throughout and the two conventions coincide.
#
# Run with:  source("dev/reproduce_Sozu_et_al_2012.R")

library(twoCoprimary)
source("dev/reproduce_helpers.R")

reproduce_begin(
  stem = "reproduce_Sozu_et_al_2012",
  article = "Sozu, Sugimoto and Hamasaki (2012), Biom J 54(5), 716-729",
  tables = c("Table 2 (PREMIER)", "Table 1", "Table S5, Normal(CC) and Arcsine(CC)")
)

ALPHA <- 0.025
BETA <- 0.2
RHOS <- c(0.0, 0.3, 0.5, 0.8)

# ---------------------------------------------------------------------------
# Table 2, page 726. mTSS as the continuous endpoint and ACR50 as the binary
# one, with the normal approximation method for the binary endpoint.
# ---------------------------------------------------------------------------
t2 <- read.csv(text = "
delta,sd,p1,p2,rho0.0,rho0.3,rho0.5,rho0.8,E1,E2
4.4,19.0,0.59,0.46,346,340,334,323,294,231
4.4,20.0,0.59,0.46,369,363,358,347,326,231
4.4,21.0,0.59,0.46,394,389,384,374,359,231
4.4,22.0,0.59,0.46,422,417,413,404,394,231
")

for (i in seq_len(nrow(t2))) {
  for (rho in RHOS) {
    res <- ss2MixedContinuousBinary(delta = t2$delta[i], sd = t2$sd[i],
                                    p1 = t2$p1[i], p2 = t2$p2[i], rho = rho,
                                    r = 1, alpha = ALPHA, beta = BETA,
                                    Test = "AN")
    cmp("Sozu 2012 Table 2",
        sprintf("sigma = %.1f, rho = %.1f", t2$sd[i], rho),
        t2[[paste0("rho", format(rho, nsmall = 1))]][i], res$n2, tol = 1)
  }
}
table_report("Sozu 2012 Table 2",
             "n per group, PREMIER scenario, AN for the binary endpoint, tolerance one subject")

for (i in seq_len(nrow(t2))) {
  cmp("Sozu 2012 Table 2, E1 and E2", sprintf("E1, sigma = %.1f", t2$sd[i]),
      t2$E1[i],
      ss1Continuous(delta = t2$delta[i], sd = t2$sd[i], r = 1,
                    alpha = ALPHA, beta = BETA)$n2, tol = 1)
}
cmp("Sozu 2012 Table 2, E1 and E2", "E2, p = (0.59, 0.46), AN", t2$E2[1],
    ss1BinaryApprox(p1 = 0.59, p2 = 0.46, r = 1, alpha = ALPHA, beta = BETA,
                    Test = "AN")$n2, tol = 0)
table_report("Sozu 2012 Table 2, E1 and E2", "the size each endpoint needs alone")

# ---------------------------------------------------------------------------
# Table 1, page 724. delta1 is the standardized effect size of the continuous
# endpoint, printed to three decimal places.
# ---------------------------------------------------------------------------
t1 <- read.csv(text = "
ratio,delta,p1,p2,rho0.0,rho0.3,rho0.5,rho0.8,E1,E2
1.0,0.100,0.55,0.50,2055,2016,1980,1904,1565,1565
1.0,0.103,0.65,0.60,1931,1895,1862,1793,1471,1471
1.0,0.112,0.75,0.70,1643,1614,1588,1534,1251,1251
1.0,0.132,0.85,0.80,1189,1171,1154,1122,906,906
1.0,0.201,0.60,0.50,509,499,490,472,388,388
1.0,0.210,0.70,0.60,468,459,451,435,356,356
1.0,0.231,0.80,0.70,385,379,373,361,294,294
1.0,0.281,0.90,0.80,262,258,254,248,199,199
1.5,0.115,0.55,0.50,1819,1789,1760,1702,1183,1565
1.5,0.119,0.65,0.60,1710,1682,1656,1603,1112,1471
1.5,0.129,0.75,0.70,1454,1432,1411,1370,946,1251
1.5,0.151,0.85,0.80,1053,1038,1025,1000,685,906
1.5,0.232,0.60,0.50,451,443,436,422,293,388
1.5,0.242,0.70,0.60,414,407,401,389,270,356
1.5,0.266,0.80,0.70,341,336,331,322,222,294
1.5,0.323,0.90,0.80,232,229,226,221,151,199
3.0,0.160,0.55,0.50,1583,1578,1574,1568,611,1565
3.0,0.165,0.65,0.60,1487,1483,1479,1474,574,1471
3.0,0.179,0.75,0.70,1265,1262,1259,1254,489,1251
3.0,0.211,0.85,0.80,916,914,912,909,354,906
3.0,0.322,0.60,0.50,392,391,390,388,152,388
3.0,0.336,0.70,0.60,360,359,358,357,139,356
3.0,0.370,0.80,0.70,297,296,295,294,115,294
3.0,0.450,0.90,0.80,202,201,201,200,78,199
")

for (i in seq_len(nrow(t1))) {
  for (rho in RHOS) {
    want <- t1[[paste0("rho", format(rho, nsmall = 1))]][i]
    res <- ss2MixedContinuousBinary(delta = t1$delta[i], sd = 1,
                                    p1 = t1$p1[i], p2 = t1$p2[i], rho = rho,
                                    r = 1, alpha = ALPHA, beta = BETA,
                                    Test = "AN")
    cmp("Sozu 2012 Table 1",
        sprintf("c1/c2 = %.1f, delta* = %.3f, p = (%.2f, %.2f), rho = %.1f",
                t1$ratio[i], t1$delta[i], t1$p1[i], t1$p2[i], rho),
        want, res$n2, tol = ceiling(0.01 * want))
  }
}
table_report("Sozu 2012 Table 1",
             "n per group; tolerance one per cent, because the printed effect size is rounded to three decimals")

for (i in seq_len(nrow(t1))) {
  cmp("Sozu 2012 Table 1, E1", sprintf("delta* = %.3f", t1$delta[i]), t1$E1[i],
      ss1Continuous(delta = t1$delta[i], sd = 1, r = 1, alpha = ALPHA,
                    beta = BETA)$n2, tol = ceiling(0.01 * t1$E1[i]))
}
table_report("Sozu 2012 Table 1, E1",
             "the continuous endpoint alone, same rounding caveat")

t1_bin <- unique(t1[, c("p1", "p2", "E2")])
for (i in seq_len(nrow(t1_bin))) {
  cmp("Sozu 2012 Table 1, E2",
      sprintf("p = (%.2f, %.2f), AN", t1_bin$p1[i], t1_bin$p2[i]), t1_bin$E2[i],
      ss1BinaryApprox(p1 = t1_bin$p1[i], p2 = t1_bin$p2[i], r = 1,
                      alpha = ALPHA, beta = BETA, Test = "AN")$n2, tol = 0)
}
table_report("Sozu 2012 Table 1, E2",
             "the binary endpoint alone, exact inputs and therefore an exact comparison")

# ---------------------------------------------------------------------------
# Table S5 of the Supporting Information, Section C
# ---------------------------------------------------------------------------
s5 <- read.csv(text = "
delta,p1,p2,method,Test,rho0.0,rho0.3,rho0.5,rho0.8,E1,E2
0.235,0.99,0.95,Normal(CC),ANc,400,397,395,391,285,333
0.235,0.99,0.95,Arcsine(CC),ASc,376,373,371,367,NA,300
0.397,0.99,0.90,Normal(CC),ANc,143,142,141,139,100,121
0.397,0.99,0.90,Arcsine(CC),ASc,129,128,127,125,NA,103
0.521,0.99,0.85,Normal(CC),ANc,84,83,82,81,58,72
0.521,0.99,0.85,Arcsine(CC),ASc,74,74,73,72,NA,60
0.190,0.95,0.90,Normal(CC),ANc,591,585,579,569,435,474
0.190,0.95,0.90,Arcsine(CC),ASc,584,578,572,562,NA,464
0.335,0.95,0.85,Normal(CC),ANc,195,192,191,187,141,160
0.335,0.95,0.85,Arcsine(CC),ASc,190,187,185,182,NA,153
0.457,0.95,0.80,Normal(CC),ANc,106,105,104,102,76,88
0.457,0.95,0.80,Arcsine(CC),ASc,102,101,100,98,NA,83
", check.names = FALSE)

for (i in seq_len(nrow(s5))) {
  for (rho in RHOS) {
    res <- ss2MixedContinuousBinary(delta = s5$delta[i], sd = 1,
                                    p1 = s5$p1[i], p2 = s5$p2[i], rho = rho,
                                    r = 1, alpha = ALPHA, beta = BETA,
                                    Test = s5$Test[i])
    cmp("Sozu 2012 Table S5",
        sprintf("%-11s p = (%.2f, %.2f), delta* = %.3f, rho = %.1f",
                s5$method[i], s5$p1[i], s5$p2[i], s5$delta[i], rho),
        s5[[paste0("rho", format(rho, nsmall = 1))]][i], res$n2, tol = 1)
  }
}
table_report("Sozu 2012 Table S5",
             "n per group at extreme proportions; the Fisher column of the same table is checked in dev/verify_mixed_fisher_option.R")

for (i in seq_len(nrow(s5))) {
  cmp("Sozu 2012 Table S5, E2",
      sprintf("%-11s p = (%.2f, %.2f)", s5$method[i], s5$p1[i], s5$p2[i]),
      s5$E2[i],
      ss1BinaryApprox(p1 = s5$p1[i], p2 = s5$p2[i], r = 1, alpha = ALPHA,
                      beta = BETA, Test = s5$Test[i])$n2, tol = 0)
}
table_report("Sozu 2012 Table S5, E2",
             "the binary endpoint alone under the continuity corrected methods, at proportions up to 0.99; exact reproduction expected since 1.1.1")

reproduce_end()
