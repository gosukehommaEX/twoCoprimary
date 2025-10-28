################################################# README #################################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
# The "approxPower2Binary" function aims to calculate power for trials 
# with two co-primary binary endpoints using normal approximiations proposed by Sozu et al. (2010).
#
## approxPower2Binary has the following arguments.
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
############################################### How to use ###############################################
# approxPower2Binary(
#   n1 = 200,
#   n2 = 100,
#   p11 = 0.5,
#   p12 = 0.4,
#   p21 = 0.3,
#   p22 = 0.2,
#   rho1 = 0.7,
#   rho2 = 0.7,
#   alpha = 0.025,
#   Test = 'AN'
# )
##########################################################################################################
library(mvtnorm)
approxPower2Binary <- function(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test) {
  
  # Define common parameters
  nu_1k <- c(p11, p12) * (1 - c(p11, p12))
  nu_2k <- c(p21, p22) * (1 - c(p21, p22))
  
  if(Test == "AN" | Test == "ANc") {
    
    # Pooled proportion
    p_k <- (n1 * c(p11, p12) + n2 * c(p21, p22)) / (n1 + n2)
    
    # Treatment effect
    delta_k <- c(p11, p12) - c(p21, p22)
    
    # Define parameter values
    se_k0 <- sqrt((1 / n1 + 1 / n2) * p_k * (1 - p_k))
    se_k <- sqrt(nu_1k / n1 + nu_2k / n2)
    rho_nml <- (rho1 * sqrt(prod(nu_1k)) / n1 + rho2 * sqrt(prod(nu_2k)) / n2) / prod(se_k)
    
    ## Asymptotic normal method without continuity correction (AN)
    if(Test == "AN") {
      
      # Power for the individual endpoint
      power1and2 <- pnorm((delta_k - se_k0 * qnorm(1 - alpha)) / se_k)
      
      # Power for the two co-primary endpoints
      powerCoprimary <- pmvnorm(
        lower = c(-Inf, -Inf),
        upper = (delta_k - se_k0 * qnorm(1 - alpha)) / se_k,
        mean = c(0, 0),
        corr = matrix(c(1, rho_nml, rho_nml, 1), ncol = 2),
        seed = 123
      )[[1]]
      
    } else if(Test == "ANc") {
      
      ## Asymptotic normal method with continuity correction (ANc)
      # Yates’s continuity correction
      ycc <- (1 / 2) * (1 / n1 + 1 / n2)
      
      # Power for the individual endpoint
      power1and2 <- pnorm((delta_k - se_k0 * qnorm(1 - alpha) - ycc) / se_k)
      
      # Power for the two co-primary endpoints
      powerCoprimary <- pmvnorm(
        lower = c(-Inf, -Inf),
        upper = (delta_k - se_k0 * qnorm(1 - alpha) - ycc) / se_k,
        mean = c(0, 0),
        corr = matrix(c(1, rho_nml, rho_nml, 1), ncol = 2),
        seed = 123
      )[[1]]
    }
    
  } else if(Test == "AS" | Test == "ASc") {
    
    # Standard error
    se <- (1 / 2) * sqrt(1 / n1 + 1 / n2)
    
    if(Test == "AS") {
      
      ## Arcsine method without continuity correction (AS)
      # Define parameter values
      delta_k <- asin(sqrt(c(p11, p12))) - asin(sqrt(c(p21, p22)))
      
      # Power for the individual endpoint
      power1and2 <- pnorm(delta_k / se - qnorm(1 - alpha))
      
      # Correlation between two test statistics
      rho_arc <- (n2 * rho1 + n1 * rho2) / (n1 + n2)
      
      # Power for the two co-primary endpoints
      powerCoprimary <- pmvnorm(
        lower = c(-Inf, -Inf),
        upper = delta_k / se - qnorm(1 - alpha),
        mean = c(0, 0),
        corr = matrix(c(1, rho_arc, rho_arc, 1), ncol = 2)
      )[[1]]
      
    } else if(Test == "ASc") {
      
      ## Arcsine method with continuity correction (ASc)
      # Define parameter values
      c1 <- -1 / (2 * n1)
      c2 <-  1 / (2 * n2)
      delta_k <- asin(sqrt(c(p11, p12) + c1)) - asin(sqrt(c(p21, p22) + c2))
      nu_1k_c <- (c(p11, p12) + c1) * (1 - c(p11, p12) - c1)
      nu_2k_c <- (c(p21, p22) + c2) * (1 - c(p21, p22) - c2)
      se_k <- sqrt(nu_1k / (4 * n1 * nu_1k_c) + nu_2k / (4 * n2 * nu_2k_c))
      
      # Power for the individual endpoint
      power1and2 <- pnorm((delta_k - se * qnorm(1 - alpha)) / se_k)
      
      # Correlation between two test statistics
      rho_arc_c <- '*'(
        1 / prod(se_k),
        '+'(
          rho1 * sqrt(prod(nu_1k)) / (4 * n1 * sqrt(prod(nu_1k_c))),
          rho2 * sqrt(prod(nu_2k)) / (4 * n2 * sqrt(prod(nu_2k_c)))
        )
      )
      
      # Power for the two co-primary endpoints
      powerCoprimary <- pmvnorm(
        lower = c(-Inf, -Inf),
        upper = (delta_k - se * qnorm(1 - alpha)) / se_k,
        mean = c(0, 0),
        corr = matrix(c(1, rho_arc_c, rho_arc_c, 1), ncol = 2)
      )[[1]]
    }
  }
  
  # Return result
  result <- data.frame(
    n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test, 
    power1 = power1and2[1], power2 = power1and2[2], powerCoprimary
  )
  return(result)
}