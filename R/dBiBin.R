####################################### README #######################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
#
# The "dBiBin" function aims to calculate a probability mass 
# function of the bivariate binomial distribution.
# For the methodological details, see formulas in Homma and Yoshida (2025).
#
# dBiBin has the following arguments;
# N:           Sample size
# y1:          Observed value(s) of the random variable Y1 
# y2:          Observed value(s) of the random variable Y2
# p1:          True probability of responders for the 1st outcome
# p1:          True probability of responders for the 2nd outcome
# rho:         The correlation between X_{i,j,1} and X_{i,j,2}
# (Note)       User can specify multiple values of y1 and y2 
##################################### How to use #####################################
# N = 100
# p1 = 0.3
# p2 = 0.5
# rho = 0.5
# sum(outer(0:N, 0:N, function(x, y) dBiBin(N, x, y, p1, p2, rho, gamma = NULL)))
# [1] 1
######################################################################################
dBiBin <- function(N, y1, y2, p1, p2, rho, gamma = NULL) {
  
  if(!is.null(rho) & !is.null(gamma)) {
    stop("either rho or gamma must be specified")
  }
  
  # Check whether or not y1 and y2 have same length
  if(length(y1) != length(y2)) stop('y1 and y2 should be the same lengths')
  
  if(!is.null(rho)) {
    if(rho < corrBound(p1, p2)[1] | rho > corrBound(p1, p2)[2]) {
      stop(paste0("rho must be within [", round(corrBound(p1, p2)[1], 4), ", ", round(corrBound(p1, p2)[2], 4), "]"))
    }
    # Obtain gamma given rho
    gamma <- gammaBiBin(p1, p2, rho)
  }
  
  if(!is.null(gamma)) {
    if(gamma < gammaBound(p1, p2)[1] | gamma > gammaBound(p1, p2)[2]) {
      stop(paste0("gamma must be within [", round(gammaBound(p1, p2)[1], 4), ", ", round(gammaBound(p1, p2)[2], 4), "]"))
    }
  }
  
  # Set M = {m: m = max(0, y2 - (N - y1)),...,min(y1,y2)}
  m <- Map(':', pmax(0, y2 - (N - y1)), pmin(y1, y2))
  m <- t(sapply(m, `length<-`, max(lengths(m))))
  
  # Define xi
  xi <- p2 + gamma * (p2 - p1)
  
  # g(y2|y1, N, p1, p2, gamma) (see formula in Homma and Yoshida (2025))
  g <- (1 + gamma) ^ (-N) * rowSums(
    '*'(
      choose(y1, m) * choose(N - y1, y2 - m) * (xi + gamma) ^ m,
      (1 - xi) ^ (y1 - m) * xi ^ (y2 - m) * (1 - xi + gamma) ^ (N - y1 - (y2 - m))
    ), 
    na.rm = TRUE
  )
  
  # Pr(Y1 = y2, Y2 = y2|N, p1, p2, gamma) (see formula in Homma and Yoshida (2025))
  return(dbinom(y1, N, p1) * g)
}