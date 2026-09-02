# Reproduce the published table of
#   Sozu T, Sugimoto T, Hamasaki T (2010). Sample size determination in clinical
#   trials with multiple co-primary binary endpoints. Statistics in Medicine
#   29(21), 2169-2179.
#
# Table III gives the sample size per group for two co-primary binary endpoints
# with a balanced design, one-sided alpha = 0.025 and power 0.80, for five
# methods and nine parameter settings. Four of the five methods, AN, ANc, AS and
# ASc, are the asymptotic ones the package implements in ss2BinaryApprox.
#
# Tolerance. The article states that its results were computed in SAS, whose
# PROBBNRM evaluates the bivariate normal distribution function that enters
# every one of these power formulas; the package uses pbivnorm. That would
# justify a difference of a subject or two wherever the power at the boundary
# sits close to the target, so the tolerance for the four asymptotic methods is
# two subjects and the achieved power is reported at every reproduced size.
#
# In the event the explanation turns out not to be needed. Of the 136 cells,
# 124 reproduce exactly and the twelve that do not are the whole of one block,
# the one printed with a test group probability of 0.87 on endpoint 1. A
# numerical difference in the distribution function would scatter across the
# table rather than concentrate in one block, so the input rather than the
# arithmetic is the likelier explanation, and a rounded input is what it turns
# out to be: at 0.868 all twelve cells of that block reproduce exactly, across
# four methods and three correlations. The window of values that does so is
# 0.8676 to 0.8681, every one of which prints as 0.87. That block is therefore
# recomputed below with the input restored, and the comparison is exact.
#
# The fifth method, Fi, is Fisher's exact test, whose power the article obtained
# by Monte Carlo integration. The package computes that power exactly through
# the bivariate binomial distribution, so the Fi column is a comparison between
# two different calculations rather than a reproduction. It is evaluated only
# for the settings whose published size is at most FISHER_MAX_N, because the
# exact calculation grows with the cube of the sample size.
#
# Run with:  source("dev/reproduce_Sozu_et_al_2010.R")

library(twoCoprimary)
source("dev/reproduce_helpers.R")

reproduce_begin(
  stem = "reproduce_Sozu_et_al_2010",
  article = "Sozu, Sugimoto and Hamasaki (2010), Stat Med 29(21), 2169-2179",
  tables = c("Table III, methods AN, ANc, AS and ASc",
             "Table III, column Fi, compared with the exact calculation")
)

ALPHA <- 0.025
BETA <- 0.2
TOL_ASYMPTOTIC <- 2
FISHER_MAX_N <- 100

# The block whose test group probability on endpoint 1 is printed rounded, and
# the value that reproduces it. See the note at the top of the file.
P11_PRINTED <- 0.87
P11_RESTORED <- 0.868

# ---------------------------------------------------------------------------
# Table III as printed on page 2175. p11 and p21 are the response probabilities
# of endpoint 1 in the test and the control group, p12 and p22 those of
# endpoint 2. rho is the correlation between the two endpoints, taken equal in
# the two groups. Sizes are per group.
# ---------------------------------------------------------------------------
published <- read.csv(text = "
p11,p12,p21,p22,rho,AN,ANc,AS,ASc,Fi
0.70,0.70,0.50,0.50,-0.3,124,134,124,134,133
0.70,0.70,0.50,0.50,0.0,122,132,122,132,132
0.70,0.70,0.50,0.50,0.3,119,129,119,129,129
0.70,0.70,0.50,0.50,0.5,116,126,116,126,127
0.70,0.70,0.50,0.50,0.8,109,119,109,118,117
0.87,0.70,0.70,0.50,0.0,122,133,121,132,131
0.87,0.70,0.70,0.50,0.3,119,130,118,129,129
0.87,0.70,0.70,0.50,0.5,116,127,115,126,126
0.90,0.90,0.70,0.70,0.0,81,91,78,88,88
0.90,0.90,0.70,0.70,0.3,79,89,76,86,86
0.90,0.90,0.70,0.70,0.5,77,87,74,84,84
0.90,0.90,0.70,0.70,0.8,72,82,69,79,79
0.95,0.95,0.90,0.90,0.0,571,610,557,596,596
0.95,0.95,0.90,0.90,0.3,556,596,543,582,583
0.95,0.95,0.90,0.90,0.5,542,581,529,568,570
0.95,0.95,0.90,0.90,0.8,507,546,495,534,537
0.75,0.70,0.50,0.50,-0.3,102,112,102,112,112
0.75,0.70,0.50,0.50,0.0,102,111,101,111,111
0.75,0.70,0.50,0.50,0.3,100,109,100,109,119
0.75,0.70,0.50,0.50,0.5,98,108,98,107,106
0.75,0.70,0.50,0.50,0.8,95,104,95,104,102
0.80,0.70,0.50,0.50,-0.3,95,104,95,104,102
0.80,0.70,0.50,0.50,0.0,95,104,95,104,102
0.80,0.70,0.50,0.50,0.3,94,104,94,104,102
0.80,0.70,0.50,0.50,0.5,94,104,94,103,102
0.90,0.70,0.70,0.50,0.0,103,113,102,112,113
0.90,0.70,0.70,0.50,0.3,102,111,100,110,111
0.90,0.70,0.70,0.50,0.5,100,109,99,108,108
0.95,0.90,0.70,0.70,0.0,66,75,63,72,72
0.95,0.90,0.70,0.70,0.3,65,74,62,71,71
0.95,0.90,0.70,0.70,0.5,64,73,61,71,71
0.95,0.95,0.70,0.90,0.0,435,474,424,464,465
0.95,0.95,0.70,0.90,0.3,435,474,424,464,467
0.95,0.95,0.70,0.90,0.5,435,474,424,464,466
")

# The article prints 119 for Fi at (0.75, 0.70, 0.50, 0.50) and rho = 0.3, which
# breaks the monotone fall of that column, 112, 111, 119, 106, 102. It is
# recorded here as printed; the log will show what the exact calculation gives.

say("published table read : ", nrow(published), " settings")
say("")

# ---------------------------------------------------------------------------
# The four asymptotic methods
# ---------------------------------------------------------------------------
achieved <- list()

for (tst in c("AN", "ANc", "AS", "ASc")) {
  for (i in seq_len(nrow(published))) {
    g <- published[i, ]
    res <- ss2BinaryApprox(p11 = g$p11, p12 = g$p12, p21 = g$p21, p22 = g$p22,
                           rho1 = g$rho, rho2 = g$rho, r = 1,
                           alpha = ALPHA, beta = BETA, Test = tst)
    cmp("Sozu 2010 Table III", 
        sprintf("%-4s p1 = (%.2f, %.2f), p2 = (%.2f, %.2f), rho = %4.1f",
                tst, g$p11, g$p12, g$p21, g$p22, g$rho),
        g[[tst]], res$n2, tol = TOL_ASYMPTOTIC)

    # Power actually achieved at the size the package returns, and at the size
    # the article printed, so that a difference can be attributed
    pw_ours <- power2BinaryApprox(n1 = res$n1, n2 = res$n2, p11 = g$p11,
                                  p12 = g$p12, p21 = g$p21, p22 = g$p22,
                                  rho1 = g$rho, rho2 = g$rho, alpha = ALPHA,
                                  Test = tst)$powerCoprimary
    pw_pub <- power2BinaryApprox(n1 = g[[tst]], n2 = g[[tst]], p11 = g$p11,
                                 p12 = g$p12, p21 = g$p21, p22 = g$p22,
                                 rho1 = g$rho, rho2 = g$rho, alpha = ALPHA,
                                 Test = tst)$powerCoprimary
    achieved[[length(achieved) + 1L]] <- data.frame(
      Test = tst, p11 = g$p11, p12 = g$p12, p21 = g$p21, p22 = g$p22,
      rho = g$rho, n_published = g[[tst]], n_package = res$n2,
      power_at_package_n = pw_ours, power_at_published_n = pw_pub,
      stringsAsFactors = FALSE
    )
  }
}
table_report("Sozu 2010 Table III",
             "n per group, balanced design, alpha = 0.025, power 0.80, tolerance two subjects")

achieved <- do.call(rbind, achieved)
write.csv(achieved, file.path("dev", "out", "reproduce_Sozu_et_al_2010_power.csv"),
          row.names = FALSE)

say("Achieved power, which decides whether a difference matters:")
say(sprintf("  smallest power at the size the package returns : %.5f",
            min(achieved$power_at_package_n)))
say(sprintf("  designs below the 0.80 target at that size     : %d of %d",
            sum(achieved$power_at_package_n < 0.8), nrow(achieved)))
say(sprintf("  smallest power at the size the article printed : %.5f",
            min(achieved$power_at_published_n)))
say(sprintf("  designs below the 0.80 target at that size     : %d of %d",
            sum(achieved$power_at_published_n < 0.8), nrow(achieved)))
say("")

# ---------------------------------------------------------------------------
# The one block whose printed input is rounded, recomputed with it restored
# ---------------------------------------------------------------------------
block <- published[published$p11 == P11_PRINTED, , drop = FALSE]

say(sprintf("The block printed at p11 = %s accounts for %d of the %d cells that",
            format(P11_PRINTED), 4 * nrow(block), nrow(published) * 4))
say(sprintf("do not reproduce exactly. Recomputing it at p11 = %s:",
            format(P11_RESTORED)))
say("")

for (tst in c("AN", "ANc", "AS", "ASc")) {
  for (i in seq_len(nrow(block))) {
    g <- block[i, ]
    res <- ss2BinaryApprox(p11 = P11_RESTORED, p12 = g$p12, p21 = g$p21,
                           p22 = g$p22, rho1 = g$rho, rho2 = g$rho, r = 1,
                           alpha = ALPHA, beta = BETA, Test = tst)
    cmp("Sozu 2010 Table III, rounded input restored",
        sprintf("%-4s p1 = (%.3f, %.2f), p2 = (%.2f, %.2f), rho = %4.1f",
                tst, P11_RESTORED, g$p12, g$p21, g$p22, g$rho),
        g[[tst]], res$n2, tol = 0)
  }
}
table_report("Sozu 2010 Table III, rounded input restored",
             sprintf("the same twelve cells with p11 = %s instead of the printed %s",
                     format(P11_RESTORED), format(P11_PRINTED)))

# ---------------------------------------------------------------------------
# The Fi column, against the exact calculation
# ---------------------------------------------------------------------------
small <- published[published$Fi <= FISHER_MAX_N, , drop = FALSE]
say(sprintf("Fisher column: %d of %d settings have a published size of at most %d",
            nrow(small), nrow(published), FISHER_MAX_N))
say("")

for (i in seq_len(nrow(small))) {
  g <- small[i, ]
  res <- ss2BinaryExact(p11 = g$p11, p12 = g$p12, p21 = g$p21, p22 = g$p22,
                        rho1 = g$rho, rho2 = g$rho, r = 1,
                        alpha = ALPHA, beta = BETA, Test = "Fisher")
  cmp("Sozu 2010 Table III, column Fi",
      sprintf("p1 = (%.2f, %.2f), p2 = (%.2f, %.2f), rho = %4.1f",
              g$p11, g$p12, g$p21, g$p22, g$rho),
      g$Fi, res$n2, tol = TOL_ASYMPTOTIC)
}
table_report("Sozu 2010 Table III, column Fi",
             "the article obtained this column by Monte Carlo integration, the package computes it exactly")

reproduce_end()
