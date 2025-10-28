############################################ README ############################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
# The "exactPower2Binary" function aims to calculate exact power for trials 
# with two co-primary binary endpoints.
#
## exactPower2Binary has the following arguments.
# n1:          Sample size for group 1
# n2:          Sample size for group 2
# p11:         True probability of responders in group 1 for the 1st outcome
# p12:         True probability of responders in group 1 for the 2nd outcome
# p21:         True probability of responders in group 2 for the 1st outcome
# p22:         True probability of responders in group 2 for the 2nd outcome
# rho1:        Assumed correlation between two outcomes for group 1 (i.e., corr(X_i11, X_i12))
# rho2:        Assumed correlation between two outcomes for group 2 (i.e., corr(X_i21, X_i22))
# alpha:       One-sided level of significance
# Test:        Statistical testing approach
########################################## How to use ##########################################
# exactPower2Binary(
#   n1 = 100,
#   n2 = 50,
#   p11 = 0.5,
#   p12 = 0.4,
#   p21 = 0.3,
#   p22 = 0.2,
#   rho1 = 0.7,
#   rho2 = 0.7,
#   alpha = 0.025,
#   Test = 'Boschloo'
# )
#    Power1          Power2 Power.coprimary 
# 0.6513156       0.7033195       0.5630546
################################################################################################
source('corrBound.R')
source('rrBinary.R')
source('dBiBin.R')
exactPower2Binary <- function(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test) {
  # Check the condition of rho_{j} \in [L_boundary, U_boundary]
  if(rho1 < corrBound(p11, p12)[1] | rho1 > corrBound(p11, p12)[2]) {
    stop(paste0("rho1 must be within [", round(corrBound(p11, p12)[1], 4), ", ", round(corrBound(p11, p12)[2], 4), "]"))
  }
  if(rho2 < corrBound(p21, p22)[1] | rho2 > corrBound(p21, p22)[2]) {
    stop(paste0("rho2 must be within [", round(corrBound(p21, p22)[1], 4), ", ", round(corrBound(p21, p22)[2], 4), "]"))
  }
  
  # Rejection region
  RR <- rrBinary(n1, n2, alpha, Test)
  
  # Power for the individual endpoint
  power1 <- sum(dbinom(0:n1, n1, p11) * pbinom(rowSums(RR) - 1, n2, p21))
  power2 <- sum(dbinom(0:n1, n1, p12) * pbinom(rowSums(RR) - 1, n2, p22))
  
  # Probability mass functions of bivariate binomial distribution for group j(=1,2)
  pmass1 <- outer(0:n1, 0:n1, function(x, y) dBiBin(n1, x, y, p11, p12, rho1, NULL))
  pmass2 <- outer(0:n2, 0:n2, function(x, y) dBiBin(n2, x, y, p21, p22, rho2, NULL))
  
  # Power for the co-primary endpoint
  powerCoprimary <- sum(
    '*'(
      t(pmass1[row(RR)[RR %>>% 0], ])[row(RR)[RR %>>% 0], ],
      t(pmass2[col(RR)[RR %>>% 0], ])[col(RR)[RR %>>% 0], ]
    )
  )
  
  # Return power
  power <- c(power1, power2, powerCoprimary)
  names(power) <- c('Power1', 'Power2', 'powerCoprimary')
  return(power)
}