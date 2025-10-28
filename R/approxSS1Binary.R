################################################# README #################################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
# The "approxSS1Binary" function aims to calculate the required sample size 
# for trials with a single binary endpoint.
#
## approxSS1Binary has the following arguments.
# p1:         True probability of responders in group 1
# p2:         True probability of responders in group 2
# p22:         True probability of responders in group 2 for the 2nd outcome
# r:           Allocation ratio to group 1 (i.e., group 1:group 2 = r:1, where r > 0)
# alpha:       One-sided level of significance
# beta:        Target type II error rate
############################################### How to use ###############################################
# approxSS1Binary(p1 = 0.6, p2 = 0.4, r = 2, alpha = 0.025, beta = 0.1)
##########################################################################################################
approxSS1Binary <- function(p1, p2, r, alpha, beta) {
  
  # Pooled proportion
  p <- (r * p1 + p2) / (1 + r)
  
  # Treatment effect
  delta <- p1 - p2
  
  # The required sample size
  n2 <- ceiling(
    (1 + 1 / r) / (delta ^ 2) * (qnorm(alpha) * sqrt(p * (1 - p)) + qnorm(beta) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))) ^ 2
  )
  n1 <- ceiling(r * n2)
  n <- n1 + n2
  
  # Return result
  result <- data.frame(p1, p2, r, alpha, beta, n1, n2, n)
  return(result)
}