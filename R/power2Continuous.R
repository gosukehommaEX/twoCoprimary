################################################# README #################################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
# The "power2Continuous" function aims to calculate power for trials 
# with two co-primary continuous endpoints proposed by Sozu et al. (2011).
#
## power2Continuous has the following arguments.
# n1:           Sample size for group 1
# n2:           Sample size for group 2
# delta1:       Mean difference for the 1st endpoint.
# delta2:       Mean difference for the 2nd endpoint.
# sd1:          Common standard deviation for the 1st endpoint.
# sd2:          Common standard deviation for the 2nd endpoint.
# rho:          Common correlation between two outcomes
# alpha:        One-sided level of significance
# known_var:    Logical value (if TRUE, power is calculated under known variance situation,
#                              power is calculated unknown variance, otherwise) 
# nMC:          Number of Monte Carlo simulations (default is 1e+4)
############################################### How to use ###############################################
# power2Continuous(
#   n1 = 490,
#   n2 = 490,
#   delta1 = 0.2,
#   delta2 = 0.2,
#   sd1 = 1,
#   sd2 = 1,
#   rho = 0.5,
#   alpha = 0.025,
#   known_var = TRUE,
#   nMC = 1e+4
# )
##########################################################################################################
library(mvtnorm)
power2Continuous <- function(n1, n2, delta1, delta2, sd1, sd2, rho, alpha, known_var = TRUE, nMC = 1e+4) {

  # Test statistic for each endpoint
  Z1 <- delta1 / (sd1 * sqrt(1 / n1 + 1 / n2))
  Z2 <- delta2 / (sd2 * sqrt(1 / n1 + 1 / n2))
  
  if(known_var) {
    
    nMC <- NA
    
    # Power for the individual endpoint
    power1 <- pnorm(-qnorm(1 - alpha) + Z1)
    power2 <- pnorm(-qnorm(1 - alpha) + Z2)
    
    # Power for the two co-primary endpoints
    powerCoprimary <- pmvnorm(
      lower = c(-Inf, -Inf),
      upper = c(-qnorm(1 - alpha) + Z1, -qnorm(1 - alpha) + Z2),
      mean = c(0, 0),
      corr = matrix(c(1, rho, rho, 1), ncol = 2),
      seed = 123
    )[[1]]
    
  } else {
    
    # Degree of freedom
    nu <- n1 + n2 - 1
    
    # Power for the individual endpoint
    power1 <- 1 - pt(qt(1 - alpha, nu), df = nu, ncp = Z1)
    power2 <- 1 - pt(qt(1 - alpha, nu), df = nu, ncp = Z2)
    
    # Power for the two co-primary endpoints
    Sigma <- matrix(c(sd1 ^ 2, rho * sd1 * sd2, rho * sd1 * sd2, sd2 ^ 2), nrow = 2)
    
    Ws <- rWishart(nMC, df = nu, Sigma = Sigma)
    probs <- numeric(nMC)
    for(i in 1:nMC) {
      Wi <- diag(Ws[, , i])
      ci <- qt(1 - alpha, df = nu) * sqrt(Wi / nu) - c(Z1, Z2)
      probs[i] <- pmvnorm(upper = -ci, sigma = Sigma)
    }
    powerCoprimary <- mean(probs)
  }
  
  # Return result
  result <- data.frame(
    n1, n2, delta1, delta2, sd1, sd2, rho, alpha, known_var, nMC,
    power1, power2, powerCoprimary
  )
  return(result)
}