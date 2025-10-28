################################################# README #################################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
# The "ss1Continuous" function aims to calculate the required sample size 
# for trials with a single continuous endpoint.
#
## ss1Continuous has the following arguments.
# delta:       Mean difference for a continuous endpoint (treatment effect).
# sd:          Common standard deviation for a continuous endpoint
# r:           Allocation ratio to group 1 (i.e., group 1:group 2 = r:1, where r > 0)
# alpha:       One-sided level of significance
# beta:        Target type II error rate
############################################### How to use ###############################################
# ss1Continuous(delta = 0.4, sd = 1, r = 2, alpha = 0.025, beta = 0.1)
##########################################################################################################
ss1Continuous <- function(delta, sd, r, alpha, beta) {
  
  # The required sample size
  n2 <- ceiling((1 + 1 / r) * sd ^ 2 * (qnorm(alpha) + qnorm(beta)) ^ 2 / (delta ^ 2))
  n1 <- ceiling(r * n2)
  n <- n1 + n2
  
  # Return result
  result <- data.frame(delta, sd, r, alpha, beta, n1, n2, n)
  return(result)
  
}