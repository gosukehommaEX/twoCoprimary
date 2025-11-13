#' Unified Interface for Mixed Count and Continuous Co-Primary Endpoints
#'
#' This function provides a unified interface for both power calculation and
#' sample size determination for trials with one count endpoint (modeled by
#' negative binomial distribution) and one continuous endpoint.
#'
#' @param n1 Sample size for group 1 (treatment group). If NULL, will be calculated.
#' @param n2 Sample size for group 2 (control group). If NULL, will be calculated.
#' @param r1 Event rate per unit time for the count endpoint in group 1
#' @param r2 Event rate per unit time for the count endpoint in group 2
#' @param nu Dispersion parameter for the negative binomial distribution
#'   (nu > 0). When nu approaches infinity, the distribution converges to Poisson.
#' @param t Follow-up period (time unit)
#' @param mu1 Mean for the continuous endpoint in group 1
#' @param mu2 Mean for the continuous endpoint in group 2
#' @param sd Common standard deviation for the continuous endpoint
#' @param rho1 Correlation between count and continuous endpoints in group 1
#' @param rho2 Correlation between count and continuous endpoints in group 2
#' @param power Target power (1 - beta). If NULL, will be calculated.
#' @param r Allocation ratio (n1/n2). Required when calculating sample size.
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#'
#' @return An object of class "twoCoprimary" containing either:
#'   \itemize{
#'     \item Power calculation results (when n1 and n2 are specified)
#'     \item Sample size calculation results (when power and r are specified)
#'   }
#'
#' @details
#' This function serves as a unified interface similar to \code{power.prop.test()}.
#' The function determines the operation mode based on which parameters are NULL.
#'
#' Exactly one of \{(n1, n2), (power, r)\} must be NULL.
#'
#' The count endpoint is modeled using a negative binomial distribution to account
#' for overdispersion. The dispersion parameter nu controls the variance:
#' Var = lambda + lambda^2/nu.
#'
#' @examples
#' # Calculate power given sample sizes
#' twoCoprimary2MixedCountContinuous(
#'   n1 = 300, n2 = 300,
#'   r1 = 1.0, r2 = 1.25,
#'   nu = 0.8, t = 1,
#'   mu1 = -50, mu2 = 0, sd = 250,
#'   rho1 = 0.5, rho2 = 0.5,
#'   alpha = 0.025
#' )
#'
#' # Calculate sample size given target power
#' twoCoprimary2MixedCountContinuous(
#'   r1 = 1.0, r2 = 1.25,
#'   nu = 0.8, t = 1,
#'   mu1 = -50, mu2 = 0, sd = 250,
#'   rho1 = 0.5, rho2 = 0.5,
#'   power = 0.8, r = 1,
#'   alpha = 0.025
#' )
#'
#' @export
twoCoprimary2MixedCountContinuous <- function(n1 = NULL, n2 = NULL,
                                              r1, r2, nu, t,
                                              mu1, mu2, sd,
                                              rho1, rho2,
                                              power = NULL, r = NULL,
                                              alpha = 0.025) {

  # Count NULL parameters
  n_null <- sum(is.null(n1), is.null(n2))
  power_null <- is.null(power)
  r_null <- is.null(r)

  # Validate parameter combinations
  if (n_null == 0 && power_null) {
    # POWER CALCULATION MODE
    result <- power2MixedCountContinuous(
      n1 = n1, n2 = n2,
      r1 = r1, r2 = r2,
      nu = nu, t = t,
      mu1 = mu1, mu2 = mu2, sd = sd,
      rho1 = rho1, rho2 = rho2,
      alpha = alpha
    )

  } else if (n_null == 2 && !power_null && !r_null) {
    # SAMPLE SIZE CALCULATION MODE
    beta <- 1 - power
    result <- ss2MixedCountContinuous(
      r1 = r1, r2 = r2,
      nu = nu, t = t,
      mu1 = mu1, mu2 = mu2, sd = sd,
      rho1 = rho1, rho2 = rho2,
      r = r, alpha = alpha, beta = beta
    )

  } else {
    # INVALID PARAMETER COMBINATION
    stop("Exactly one of {(n1, n2), (power, r)} must be NULL.\n",
         "  - For power calculation: provide n1 and n2, set power = NULL\n",
         "  - For sample size calculation: provide power and r, set n1 = n2 = NULL")
  }

  return(result)
}
