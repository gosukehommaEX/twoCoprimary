################################################# README #################################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
# The "approxSS2Binary" function aims to calculate the required sample size 
# for trials with two co-primary binary endpoints.
#
## approxSS2Binary has the following arguments.
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
# approxSS2Binary(
#   p11 = 0.5,
#   p12 = 0.4,
#   p21 = 0.3,
#   p22 = 0.2,
#   rho1 = 0.7,
#   rho2 = 0.7,
#   r = 2,
#   alpha = 0.025,
#   beta = 0.1,
#   Test = 'AN'
# )
##########################################################################################################
source('approxSS1Binary.R')
approxSS2Binary <- function(p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test) {
  
  # Find required sample size using iterative approach
  # Step 1: Initial estimate using a formula for single endpoint
  # Set initial sample sizes
  n2_0 <- max(approxSS1Binary(c(p11, p12), c(p21, p22), r, alpha, beta)[["n2"]])
  n2_1 <- max(approxSS1Binary(c(p11, p12), c(p21, p22), r, alpha, 1 - (1 - beta) ^ (1/2))[["n2"]])
  
  # Step 2-4: Iterative refinement using linear extrapolation
  while (n2_1 - n2_0 != 0) {
    # Calculate sample size for group 1
    n1_0 <- ceiling(r * n2_0)
    n1_1 <- ceiling(r * n2_1)
    
    # Calculate power at two candidate values
    power_n2_0 <- approxPower2Binary(n1_0, n2_0, p11, p12, p21, p22, rho1, rho2, alpha, Test)[["powerCoprimary"]]
    power_n2_1 <- approxPower2Binary(n1_1, n2_1, p11, p12, p21, p22, rho1, rho2, alpha, Test)[["powerCoprimary"]]
    
    # Linear extrapolation to find D that achieves target power
    n2_updated <- (n2_0 * (power_n2_1 - (1 - beta)) - n2_1 * (power_n2_0 - (1 - beta))) / (power_n2_1 - power_n2_0)
    # Update values for next iteration
    n2_1 <- n2_0
    n2_0 <- ceiling(n2_updated)
  }
  
  # Final sample size
  n2 <- n2_0
  n1 <- ceiling(r * n2)
  n <- n1 + n2
  
  # Return result
  result <- data.frame(
    p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test, 
    n1, n2, n
  )
  return(result)
}