#' Unified Interface for Two Co-Primary Continuous Endpoints
#'
#' This function provides a unified interface for both power calculation and
#' sample size determination for two co-primary continuous endpoints. Depending
#' on which parameters are provided (sample sizes or power), the function
#' automatically determines whether to calculate power or sample size.
#'
#' @param n1 Sample size for group 1 (treatment group). If NULL, will be calculated.
#' @param n2 Sample size for group 2 (control group). If NULL, will be calculated.
#' @param delta1 Mean difference for the first endpoint
#' @param delta2 Mean difference for the second endpoint
#' @param sd1 Common standard deviation for the first endpoint
#' @param sd2 Common standard deviation for the second endpoint
#' @param rho Common correlation between the two outcomes
#' @param power Target power (1 - beta). If NULL, will be calculated.
#' @param r Allocation ratio (n1/n2). Required when calculating sample size.
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param known_var Logical indicating whether variance is known (TRUE) or
#'   unknown (FALSE). Default is TRUE.
#' @param nMC Number of Monte Carlo simulations when known_var = FALSE.
#'   Default is 10000.
#'
#' @return An object of class "twoCoprimary" containing either:
#'   \itemize{
#'     \item Power calculation results (when n1 and n2 are specified)
#'     \item Sample size calculation results (when power and r are specified)
#'   }
#'
#' @details
#' This function serves as a unified interface similar to \code{power.prop.test()}.
#' The function determines the operation mode based on which parameters are NULL:
#' \itemize{
#'   \item If n1 and n2 are provided and power is NULL: calculates power
#'   \item If power and r are provided and n1/n2 are NULL: calculates sample size
#' }
#'
#' Exactly one of \{(n1, n2), (power, r)\} must be NULL.
#'
#' @examples
#' # Calculate power given sample sizes
#' twoCoprimary2Continuous(
#'   n1 = 100, n2 = 100,
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, alpha = 0.025,
#'   known_var = TRUE
#' )
#'
#' # Calculate sample size given target power
#' twoCoprimary2Continuous(
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, power = 0.8,
#'   r = 1, alpha = 0.025,
#'   known_var = TRUE
#' )
#'
#' @export
twoCoprimary2Continuous <- function(n1 = NULL, n2 = NULL,
                                    delta1, delta2, sd1, sd2, rho,
                                    power = NULL, r = NULL, alpha = 0.025,
                                    known_var = TRUE, nMC = 1e4) {

  # Count NULL parameters
  n_null <- sum(is.null(n1), is.null(n2))
  power_null <- is.null(power)
  r_null <- is.null(r)

  # Validate parameter combinations
  if (n_null == 0 && power_null) {
    # POWER CALCULATION MODE
    # Both n1 and n2 are specified, calculate power
    result <- power2Continuous(
      n1 = n1, n2 = n2,
      delta1 = delta1, delta2 = delta2,
      sd1 = sd1, sd2 = sd2,
      rho = rho, alpha = alpha,
      known_var = known_var, nMC = nMC
    )

  } else if (n_null == 2 && !power_null && !r_null) {
    # SAMPLE SIZE CALCULATION MODE
    # Both n1 and n2 are NULL, power and r are specified
    beta <- 1 - power
    result <- ss2Continuous(
      delta1 = delta1, delta2 = delta2,
      sd1 = sd1, sd2 = sd2,
      rho = rho, r = r,
      alpha = alpha, beta = beta,
      known_var = known_var, nMC = nMC
    )

  } else {
    # INVALID PARAMETER COMBINATION
    stop("Exactly one of {(n1, n2), (power, r)} must be NULL.\n",
         "  - For power calculation: provide n1 and n2, set power = NULL\n",
         "  - For sample size calculation: provide power and r, set n1 = n2 = NULL")
  }

  return(result)
}
