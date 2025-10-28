####################################### README #######################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
#
# The "gammaBiBin" function aims to obtain a dependence parameter gamma_j for
# the bivariate binomial distribution.
# For the methodological details, see formulas in Homma and Yoshida (2025).
#
# gammaBiBin has the following arguments;
# p1:          True probability of responders for the 1st outcome
# p2:          True probability of responders for the 2nd outcome
# rho:         The correlation between X_{i,j,1} and X_{i,j,2}
##################################### How to use #####################################
# gammaBiBin(0.6, 0.3, 0.5)
######################################################################################
gammaBiBin <- function(p1, p2, rho) {
  
  gamma <- '/'(
    rho * sqrt(p2 * (1 - p2) / (p1 * (1 - p1))),
    (1 - rho * sqrt(p2 * (1 - p2) / (p1 * (1 - p1))))
  )
  
  return(gamma)
}