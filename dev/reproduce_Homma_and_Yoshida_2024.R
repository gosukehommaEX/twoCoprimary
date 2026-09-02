# Reproduce the published tables of
#   Homma G, Yoshida T (2024). Sample size calculation in clinical trials with
#   two co-primary endpoints including overdispersed count and continuous
#   outcomes. Pharmaceutical Statistics 23(1), 46-59.
#
# Table 1 (balanced design) and Table 2 (two to one allocation) give, for two
# cases and five correlations, the required sample size in each group and the
# correlation between the two test statistics that the sample size calculation
# uses. The simulated powers in the same tables come from a copula simulation
# outside the scope of the package and are not reproduced here.
#
# The article was implemented in R from the outset and the package implements
# the same formulas, so the tolerance is zero for the sample sizes.
#
# Notation. The article indexes the control group by 0 and the treatment group
# by 1, and writes kappa = n1 / n0. The package calls the same groups n1
# (treatment) and n2 (control) and calls the ratio r. So r0 of the article is
# the argument r2 of the package, r1 is r1, mu0 is mu2, mu1 is mu1, and kappa
# is r. In both cases of the article the treatment benefit is a lower event
# rate and a lower mean, which is the convention the package requires.
#
# The correlation between the test statistics is not returned by the package. It
# is recovered here by inverting the bivariate normal probability that produced
# the co-primary power, which tests the mapping from the outcome correlation to
# the test statistic correlation, equations (11) and (12) of the article.
#
# Run with:  source("dev/reproduce_Homma_and_Yoshida_2024.R")

library(twoCoprimary)
source("dev/reproduce_helpers.R")

reproduce_begin(
  stem = "reproduce_Homma_and_Yoshida_2024",
  article = "Homma and Yoshida (2024), Pharm Stat 23(1), 46-59",
  tables = c("Table 1 (kappa = 1)", "Table 2 (kappa = 2)",
             "the correlation between the test statistics in both tables")
)

ALPHA <- 0.025
BETA <- 0.1
TFOLLOW <- 1

# Case A: control rate 1.25, treatment rate 1, mu = (0, -50), sigma = 250
# Case B: control rate 2,    treatment rate 1, mu = (0, -50), sigma = 75
cases <- list(
  A = list(r_treat = 1.0, r_ctrl = 1.25, sd = 250, nus = c(0.8, 1)),
  B = list(r_treat = 1.0, r_ctrl = 2.00, sd = 75,  nus = c(3, 5))
)
MU_TREAT <- -50
MU_CTRL <- 0

# ---------------------------------------------------------------------------
# Tables 1 and 2 as printed on pages 52 and 53. n0 is the control group, n1 the
# treatment group, gamma the correlation between the two test statistics.
# ---------------------------------------------------------------------------
published <- read.csv(text = "
kappa,case,nu,rho,gamma,n0,n1,N
1,A,0.8,0.0,0.000,935,935,1870
1,A,0.8,0.2,0.200,932,932,1864
1,A,0.8,0.4,0.400,927,927,1854
1,A,0.8,0.6,0.600,920,920,1840
1,A,0.8,0.8,0.800,911,911,1822
1,A,1.0,0.0,0.000,846,846,1692
1,A,1.0,0.2,0.200,842,842,1684
1,A,1.0,0.4,0.400,835,835,1670
1,A,1.0,0.6,0.600,825,825,1650
1,A,1.0,0.8,0.800,811,811,1622
1,B,3.0,0.0,0.000,59,59,118
1,B,3.0,0.2,0.199,58,58,116
1,B,3.0,0.4,0.397,57,57,114
1,B,3.0,0.6,0.596,56,56,112
1,B,3.0,0.8,0.795,54,54,108
1,B,5.0,0.0,0.000,55,55,110
1,B,5.0,0.2,0.198,55,55,110
1,B,5.0,0.4,0.396,54,54,108
1,B,5.0,0.6,0.595,53,53,106
1,B,5.0,0.8,0.793,51,51,102
2,A,0.8,0.0,0.000,692,1384,2076
2,A,0.8,0.2,0.200,690,1380,2070
2,A,0.8,0.4,0.400,686,1372,2058
2,A,0.8,0.6,0.600,680,1360,2040
2,A,0.8,0.8,0.800,673,1346,2019
2,A,1.0,0.0,0.000,626,1252,1878
2,A,1.0,0.2,0.200,623,1246,1869
2,A,1.0,0.4,0.400,617,1234,1851
2,A,1.0,0.6,0.600,610,1220,1830
2,A,1.0,0.8,0.800,599,1198,1797
2,B,3.0,0.0,0.000,42,84,126
2,B,3.0,0.2,0.199,42,84,126
2,B,3.0,0.4,0.397,42,84,126
2,B,3.0,0.6,0.596,41,82,123
2,B,3.0,0.8,0.795,39,78,117
2,B,5.0,0.0,0.000,40,80,120
2,B,5.0,0.2,0.198,40,80,120
2,B,5.0,0.4,0.397,39,78,117
2,B,5.0,0.6,0.595,39,78,117
2,B,5.0,0.8,0.793,38,76,114
")

for (kp in c(1, 2)) {
  block <- published[published$kappa == kp, , drop = FALSE]
  tid <- sprintf("Homma and Yoshida 2024 Table %d", kp)
  for (i in seq_len(nrow(block))) {
    g <- block[i, ]
    cs <- cases[[g$case]]
    lab <- sprintf("case %s, nu = %.1f, rho = %.1f", g$case, g$nu, g$rho)

    res <- ss2MixedCountContinuous(
      r1 = cs$r_treat, r2 = cs$r_ctrl, nu = g$nu, t = TFOLLOW,
      mu1 = MU_TREAT, mu2 = MU_CTRL, sd = cs$sd, r = kp,
      rho1 = g$rho, rho2 = g$rho, alpha = ALPHA, beta = BETA
    )
    cmp(tid, paste(lab, "| n0, the control group"), g$n0, res$n2, tol = 0)
    cmp(tid, paste(lab, "| n1, the treatment group"), g$n1, res$n1, tol = 0)
    cmp(tid, paste(lab, "| N, the total"), g$N, res$N, tol = 0)
  }
  table_report(tid, sprintf("alpha = 0.025, power 0.90, kappa = %d, exact reproduction expected", kp))
}

# ---------------------------------------------------------------------------
# The correlation between the two test statistics, recovered from the power
# ---------------------------------------------------------------------------
for (kp in c(1, 2)) {
  block <- published[published$kappa == kp, , drop = FALSE]
  tid <- sprintf("Homma and Yoshida 2024 Table %d, gamma", kp)
  for (i in seq_len(nrow(block))) {
    g <- block[i, ]
    cs <- cases[[g$case]]
    pw <- power2MixedCountContinuous(
      n1 = g$n1, n2 = g$n0, r1 = cs$r_treat, r2 = cs$r_ctrl, nu = g$nu,
      t = TFOLLOW, mu1 = MU_TREAT, mu2 = MU_CTRL, sd = cs$sd,
      rho1 = g$rho, rho2 = g$rho, alpha = ALPHA
    )
    got <- implied_gamma(pw$powerCount, pw$powerCont, pw$powerCoprimary)
    cmp(tid, sprintf("case %s, nu = %.1f, rho = %.1f", g$case, g$nu, g$rho),
        g$gamma, round(got, 3), tol = 0.001, quantity = "correlation")
  }
  table_report(tid,
               "recovered by inverting the bivariate normal probability; the article prints three decimals")
}

reproduce_end()
