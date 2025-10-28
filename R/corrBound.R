############################################ README ############################################
# The code has been tested using R version 4.5.1 (2025-06-13 ucrt).
# The "corrBound" function aims to calculate a boundary of a correlation 
# between X_{i,j,1} and X_{i,j,2}
#
# corrBound has the following arguments.
# p1:          True probability of responders for the 1st outcome
# p1:          True probability of responders for the 2nd outcome
########################################## How to use ##########################################
# corrBound(p1 = 0.3, p2 = 0.5)
################################################################################################
corrBound <- function(p1, p2) {
  # Boundary of rho (see Prentice (1988))
  boundary <- c(
    max(
      -sqrt((p1 * p2) / ((1 - p1) * (1 - p2))),
      -sqrt(((1 - p1) * (1 - p2)) / (p1 * p2))
    ),
    min(
      sqrt((p1 * (1 - p2)) / (p2 * (1 - p1))),
      sqrt((p2 * (1 - p1)) / (p1 * (1 - p2)))
    )
  )
  names(boundary) <- c('L_bound', 'U_bound')
  # Return output
  return(boundary)
}