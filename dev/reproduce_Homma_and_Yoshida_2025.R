# Reproduce the published tables of
#   Homma G, Yoshida T (2025). Exact power and sample size in clinical trials
#   with two co-primary binary endpoints. Statistical Methods in Medical
#   Research 34(11), 2183-2201.
#
# Table 4 gives, for the mepolizumab EGPA scenario, the required total sample
# size and the exact power at that size for four tests, three significance
# levels, two allocation ratios and four correlations.
#
# Table 3 gives, for four response probability scenarios, the total sample size
# under the chi-squared test together with the local power of endpoint 1, the
# local type I error rate of endpoint 2 and the overall test size. The sample
# size there is computed under an alternative in which the two endpoints have
# the same probabilities, and the three probabilities are then evaluated at the
# nuisance value p that maximizes the type I error, which the table also prints.
#
# The article was implemented in R and the package implements the same
# calculation, so the tolerance is zero for the sample sizes and one unit in the
# last printed digit for the probabilities.
#
# One published entry is known to be wrong and is expected to differ here.
# Table 4, Z-pool, alpha = 0.05, r = 2, rho = 0.3 prints 144 with a power of
# 0.903. The rejection region behind that number accumulated the null tail
# probability along an arbitrary ordering of tied outcomes, which enlarged the
# region and overstated the power. Under the standard convention, which the
# package has used since version 1.1.0 and which agrees cell for cell with the
# Exact package, the answer is 147. The cell is recorded below as printed, so it
# appears as a difference, and the corrected value is then checked separately.
#
# Running time. The exact calculation grows with the cube of the sample size and
# Table 4 requires 24 sample size searches per test. Restricting ALPHAS to 0.025
# reproduces the block that the R Journal article reproduces and takes about a
# third of the time.
#
# Run with:  source("dev/reproduce_Homma_and_Yoshida_2025.R")

library(twoCoprimary)
source("dev/reproduce_helpers.R")

reproduce_begin(
  stem = "reproduce_Homma_and_Yoshida_2025",
  article = "Homma and Yoshida (2025), Stat Methods Med Res 34(11), 2183-2201",
  tables = c("Table 4 (sample size and exact power)",
             "Table 3 (chi-squared test size)")
)

ALPHAS <- c(0.025, 0.05, 0.1)
BETA <- 0.1
TESTS <- c("Chisq", "Fisher", "Z-pool", "Boschloo")
P11 <- 0.54; P12 <- 0.54; P21 <- 0.25; P22 <- 0.25

# ---------------------------------------------------------------------------
# Table 4, page 14. N is the total sample size and the value in parentheses in
# the article is the exact power at that size.
# ---------------------------------------------------------------------------
t4 <- read.csv(text = "
alpha,r,rho,test,N,power
0.025,1,0.0,Chisq,142,0.902
0.025,1,0.0,Fisher,152,0.901
0.025,1,0.0,Z-pool,144,0.905
0.025,1,0.0,Boschloo,144,0.905
0.025,1,0.3,Chisq,142,0.906
0.025,1,0.3,Fisher,150,0.901
0.025,1,0.3,Z-pool,142,0.903
0.025,1,0.3,Boschloo,142,0.903
0.025,1,0.5,Chisq,140,0.905
0.025,1,0.5,Fisher,150,0.906
0.025,1,0.5,Z-pool,140,0.902
0.025,1,0.5,Boschloo,140,0.902
0.025,1,0.8,Chisq,128,0.900
0.025,1,0.8,Fisher,144,0.901
0.025,1,0.8,Z-pool,134,0.901
0.025,1,0.8,Boschloo,134,0.901
0.025,2,0.0,Chisq,162,0.907
0.025,2,0.0,Fisher,174,0.903
0.025,2,0.0,Z-pool,180,0.901
0.025,2,0.0,Boschloo,162,0.905
0.025,2,0.3,Chisq,159,0.902
0.025,2,0.3,Fisher,174,0.908
0.025,2,0.3,Z-pool,180,0.906
0.025,2,0.3,Boschloo,159,0.902
0.025,2,0.5,Chisq,156,0.901
0.025,2,0.5,Fisher,171,0.904
0.025,2,0.5,Z-pool,177,0.903
0.025,2,0.5,Boschloo,156,0.901
0.025,2,0.8,Chisq,147,0.903
0.025,2,0.8,Fisher,159,0.900
0.025,2,0.8,Z-pool,168,0.908
0.025,2,0.8,Boschloo,150,0.903
0.05,1,0.0,Chisq,116,0.902
0.05,1,0.0,Fisher,132,0.901
0.05,1,0.0,Z-pool,120,0.903
0.05,1,0.0,Boschloo,120,0.903
0.05,1,0.3,Chisq,116,0.907
0.05,1,0.3,Fisher,132,0.906
0.05,1,0.3,Z-pool,116,0.900
0.05,1,0.3,Boschloo,116,0.900
0.05,1,0.5,Chisq,114,0.905
0.05,1,0.5,Fisher,128,0.904
0.05,1,0.5,Z-pool,114,0.902
0.05,1,0.5,Boschloo,114,0.902
0.05,1,0.8,Chisq,110,0.906
0.05,1,0.8,Fisher,118,0.900
0.05,1,0.8,Z-pool,110,0.904
0.05,1,0.8,Boschloo,110,0.904
0.05,2,0.0,Chisq,135,0.907
0.05,2,0.0,Fisher,147,0.901
0.05,2,0.0,Z-pool,147,0.904
0.05,2,0.0,Boschloo,135,0.904
0.05,2,0.3,Chisq,132,0.901
0.05,2,0.3,Fisher,147,0.906
0.05,2,0.3,Z-pool,144,0.903
0.05,2,0.3,Boschloo,132,0.901
0.05,2,0.5,Chisq,129,0.901
0.05,2,0.5,Fisher,144,0.902
0.05,2,0.5,Z-pool,138,0.902
0.05,2,0.5,Boschloo,132,0.907
0.05,2,0.8,Chisq,120,0.901
0.05,2,0.8,Fisher,135,0.907
0.05,2,0.8,Z-pool,135,0.907
0.05,2,0.8,Boschloo,126,0.907
0.1,1,0.0,Chisq,92,0.903
0.1,1,0.0,Fisher,108,0.907
0.1,1,0.0,Z-pool,104,0.902
0.1,1,0.0,Boschloo,98,0.907
0.1,1,0.3,Chisq,88,0.904
0.1,1,0.3,Fisher,106,0.905
0.1,1,0.3,Z-pool,98,0.901
0.1,1,0.3,Boschloo,96,0.904
0.1,1,0.5,Chisq,86,0.904
0.1,1,0.5,Fisher,104,0.902
0.1,1,0.5,Z-pool,96,0.905
0.1,1,0.5,Boschloo,94,0.902
0.1,1,0.8,Chisq,82,0.902
0.1,1,0.8,Fisher,100,0.903
0.1,1,0.8,Z-pool,92,0.905
0.1,1,0.8,Boschloo,88,0.903
0.1,2,0.0,Chisq,108,0.908
0.1,2,0.0,Fisher,120,0.901
0.1,2,0.0,Z-pool,114,0.904
0.1,2,0.0,Boschloo,108,0.908
0.1,2,0.3,Chisq,105,0.905
0.1,2,0.3,Fisher,120,0.906
0.1,2,0.3,Z-pool,114,0.909
0.1,2,0.3,Boschloo,105,0.902
0.1,2,0.5,Chisq,99,0.904
0.1,2,0.5,Fisher,117,0.903
0.1,2,0.5,Z-pool,111,0.904
0.1,2,0.5,Boschloo,105,0.908
0.1,2,0.8,Chisq,93,0.903
0.1,2,0.8,Fisher,108,0.907
0.1,2,0.8,Z-pool,108,0.908
0.1,2,0.8,Boschloo,99,0.906
")

t4 <- t4[t4$alpha %in% ALPHAS, , drop = FALSE]
say(sprintf("Table 4: %d of 96 published entries selected by ALPHAS", nrow(t4)))
say("")

for (i in seq_len(nrow(t4))) {
  g <- t4[i, ]
  lab <- sprintf("%-8s alpha = %.3f, r = %d, rho = %.1f",
                 g$test, g$alpha, g$r, g$rho)
  res <- ss2BinaryExact(p11 = P11, p12 = P12, p21 = P21, p22 = P22,
                        rho1 = g$rho, rho2 = g$rho, r = g$r,
                        alpha = g$alpha, beta = BETA, Test = g$test)
  cmp("Homma and Yoshida 2025 Table 4", paste(lab, "| N"), g$N, res$N, tol = 0)

  pw <- power2BinaryExact(n1 = res$n1, n2 = res$n2, p11 = P11, p12 = P12,
                          p21 = P21, p22 = P22, rho1 = g$rho, rho2 = g$rho,
                          alpha = g$alpha, Test = g$test)
  cmp("Homma and Yoshida 2025 Table 4, power", paste(lab, "| exact power"),
      g$power, round(pw$powerCoprimary, 3), tol = 0.001, quantity = "power")
}
table_report("Homma and Yoshida 2025 Table 4",
             "total sample size; one entry, Z-pool at alpha = 0.05, r = 2, rho = 0.3, is expected to differ")
table_report("Homma and Yoshida 2025 Table 4, power",
             "exact power at the size the package returns, compared with the published three decimals")

# ---------------------------------------------------------------------------
# The single entry that the tie handling correction changes
# ---------------------------------------------------------------------------
if (0.05 %in% ALPHAS) {
  res <- ss2BinaryExact(p11 = P11, p12 = P12, p21 = P21, p22 = P22,
                        rho1 = 0.3, rho2 = 0.3, r = 2, alpha = 0.05,
                        beta = BETA, Test = "Z-pool")
  cmp("Homma and Yoshida 2025 Table 4, corrected entry",
      "Z-pool, alpha = 0.05, r = 2, rho = 0.3", 147, res$N, tol = 0)
  pw <- power2BinaryExact(n1 = res$n1, n2 = res$n2, p11 = P11, p12 = P12,
                          p21 = P21, p22 = P22, rho1 = 0.3, rho2 = 0.3,
                          alpha = 0.05, Test = "Z-pool")
  say(sprintf("  exact power at N = %d is %.4f; the published N = 144 gives %.4f",
              res$N, pw$powerCoprimary,
              power2BinaryExact(n1 = 96, n2 = 48, p11 = P11, p12 = P12,
                                p21 = P21, p22 = P22, rho1 = 0.3, rho2 = 0.3,
                                alpha = 0.05, Test = "Z-pool")$powerCoprimary))
  table_report("Homma and Yoshida 2025 Table 4, corrected entry",
               "the value the standard tie convention gives, which the R Journal article discloses")
}

# ---------------------------------------------------------------------------
# Table 3, page 12. The sample size is obtained under an alternative in which
# the two endpoints have the same probabilities. The three probabilities are
# then evaluated with endpoint 2 set to the nuisance value p in both groups,
# which is the configuration that maximizes the type I error.
# ---------------------------------------------------------------------------
t3 <- read.csv(text = "
p11,p21,p,r,rho,N,power1,tie2,size
0.4,0.2,0.09,1,0.0,212,0.896,0.0263,0.0236
0.4,0.2,0.59,1,0.3,206,0.889,0.0269,0.0263
0.4,0.2,0.50,1,0.5,200,0.880,0.0280,0.0279
0.4,0.2,0.99,2,0.0,243,0.898,0.0415,0.0372
0.4,0.2,0.73,2,0.3,234,0.887,0.0270,0.0264
0.4,0.2,0.50,2,0.5,228,0.879,0.0273,0.0273
0.5,0.3,0.62,1,0.0,242,0.895,0.0262,0.0235
0.5,0.3,0.50,1,0.3,234,0.887,0.0289,0.0282
0.5,0.3,0.50,1,0.5,230,0.879,0.0278,0.0277
0.5,0.3,0.40,1,0.8,220,0.860,0.0242,0.0242
0.5,0.3,0.99,2,0.0,273,0.895,0.0417,0.0374
0.5,0.3,0.82,2,0.3,270,0.890,0.0270,0.0264
0.5,0.3,0.63,2,0.5,264,0.880,0.0258,0.0257
0.5,0.3,0.40,2,0.8,246,0.860,0.0243,0.0243
0.7,0.4,0.20,1,0.0,112,0.896,0.0259,0.0232
0.7,0.4,0.35,1,0.3,110,0.889,0.0265,0.0259
0.7,0.4,0.50,1,0.5,102,0.886,0.0297,0.0296
0.7,0.4,0.98,2,0.0,123,0.899,0.0410,0.0368
0.7,0.4,0.52,2,0.3,117,0.888,0.0289,0.0282
0.7,0.4,0.52,2,0.5,117,0.888,0.0289,0.0288
0.8,0.5,0.50,1,0.0,100,0.895,0.0284,0.0255
0.8,0.5,0.50,1,0.3,100,0.895,0.0284,0.0279
0.8,0.5,0.70,1,0.5,96,0.885,0.0261,0.0261
0.8,0.5,0.98,2,0.0,111,0.900,0.0402,0.0362
0.8,0.5,0.81,2,0.3,108,0.893,0.0305,0.0299
0.8,0.5,0.79,2,0.5,102,0.880,0.0297,0.0296
")

for (i in seq_len(nrow(t3))) {
  g <- t3[i, ]
  lab <- sprintf("p = (%.1f, %.1f), r = %d, rho = %.1f", g$p11, g$p21, g$r, g$rho)

  res <- ss2BinaryExact(p11 = g$p11, p12 = g$p11, p21 = g$p21, p22 = g$p21,
                        rho1 = g$rho, rho2 = g$rho, r = g$r, alpha = 0.025,
                        beta = 0.2, Test = "Chisq")
  cmp("Homma and Yoshida 2025 Table 3", paste(lab, "| N"), g$N, res$N, tol = 0)

  pw_alt <- power2BinaryExact(n1 = res$n1, n2 = res$n2, p11 = g$p11,
                              p12 = g$p11, p21 = g$p21, p22 = g$p21,
                              rho1 = g$rho, rho2 = g$rho, alpha = 0.025,
                              Test = "Chisq")
  cmp("Homma and Yoshida 2025 Table 3", paste(lab, "| local power, endpoint 1"),
      g$power1, round(pw_alt$power1, 3), tol = 0.001, quantity = "power")

  # Endpoint 2 set to the nuisance value in both groups
  pw_null <- power2BinaryExact(n1 = res$n1, n2 = res$n2, p11 = g$p11,
                               p12 = g$p, p21 = g$p21, p22 = g$p,
                               rho1 = g$rho, rho2 = g$rho, alpha = 0.025,
                               Test = "Chisq")
  cmp("Homma and Yoshida 2025 Table 3", paste(lab, "| local type I error, endpoint 2"),
      g$tie2, round(pw_null$power2, 4), tol = 0.0001, quantity = "type I error")
  cmp("Homma and Yoshida 2025 Table 3", paste(lab, "| overall test size"),
      g$size, round(pw_null$powerCoprimary, 4), tol = 0.0001, quantity = "test size")
}
table_report("Homma and Yoshida 2025 Table 3",
             "chi-squared test, alpha = 0.025, power at least 0.80 under equal probabilities on the two endpoints")

reproduce_end()
