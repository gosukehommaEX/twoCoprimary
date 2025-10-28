############################################ README ############################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
# The "gammaBound" function aims to calculate a boundary of a dependence parameter gamma_j for
# the bivariate binomial distribution.
#
# gammaBound has the following arguments.
# p1:          True probability of responders for the 1st outcome
# p1:          True probability of responders for the 2nd outcome
########################################## How to use ##########################################
# gammaBound(p1 = 0.3, p2 = 0.5)
################################################################################################
gammaBound <- function(p1, p2) {
  # Boundary of gamma (see Homma and Yoshida (2025))
  
  # Lower bound
  if(p1 > p2) {
    L_bound <- -p2 / (1 - p1 + p2)
  } else if(p1 < p2) {
    L_bound <- -(1 - p2) / (1 + p1 - p2)
  } else {
    L_bound <- -1 / 2
  }
  
  # Upper bound
  if(p1 > p2) {
    U_bound <- p2 / (p1 - p2)
  } else if(p1 < p2) {
    U_bound <- (1 - p2) / (p2 - p1)
  } else {
    U_bound <- Inf
  }
  
  boundary <- c(L_bound, U_bound)
  names(boundary) <- c('L_bound', 'U_bound')
  # Return output
  return(boundary)
}