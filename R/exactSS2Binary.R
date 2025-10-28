################################################# README #################################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
# The "exactSS2Binary" function aims to calculate the required sample size 
# for trials with two co-primary binary endpoints using exact formula proposed by Homma and Yoshida (2025).
#
## exactSS2Binary has the following arguments.
# p11:         True probability of responders in group 1 for the 1st outcome
# p12:         True probability of responders in group 1 for the 2nd outcome
# p21:         True probability of responders in group 2 for the 1st outcome
# p22:         True probability of responders in group 2 for the 2nd outcome
# rho1:        Assumed correlation between two outcomes for group 1 (i.e., corr(X_i11, X_i12))
# rho2:        Assumed correlation between two outcomes for group 2 (i.e., corr(X_i21, X_i22))
# r:           Allocation ratio to group 1 (i.e., group 1:group 2 = r:1, where r > 0)
# alpha:       One-sided level of significance
# beta:        Target type II error rate
# Test:        Statistical testing approach 
############################################### How to use ###############################################
# exactSS2Binary(
#   p11 = 0.5,
#   p12 = 0.4,
#   p21 = 0.3,
#   p22 = 0.2,
#   rho1 = 0.7,
#   rho2 = 0.7,
#   r = 2,
#   alpha = 0.025,
#   beta = 0.1,
#   Test = 'Boschloo'
# )
##########################################################################################################
source('approxSS2Binary.R')
source('exactPower2Binary.R')
exactSS2Binary <- function(p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test) {
  
  # Step 1: set initial sample size
  n2 <- approxSS2Binary(p11, p12, p21, p22, rho1, rho2, r, alpha, beta, "AN")[["n2"]]
  n1 <- ceiling(r * n2)
  
  # Step 2: calculate power given n1 and n2
  power <- exactPower2Binary(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test)[["powerCoprimary"]]
  
  # Step 3 - 1: n2 = n2 - 1, if power >= 1 - beta
  if(power %>=% (1 - beta)) {
    while(power %>=% (1 - beta)) {
      n2 = n2 - 1
      n1 = ceiling(r * n2)
      power = exactPower2Binary(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test)[["powerCoprimary"]]
    }
    n2 = n2 + 1
  } else {
    # Step 3 - 2: n2 = n2 + 1, if power < 1 - beta
    while(power %<<% (1 - beta)) {
      n2 = n2 + 1
      n1 = ceiling(r * n2)
      power = exactPower2Binary(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test)[["powerCoprimary"]]
    }
  }
  
  # Step 4 (determine the final sample size)
  n = n1 + n2
  
  # Return result
  result <- data.frame(
    p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test, 
    n1, n2, n
  )
  return(result)
}