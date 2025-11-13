#' Unified Interface for Mixed Continuous and Binary Co-Primary Endpoints
#'
#' This function provides a unified interface for both power calculation and
#' sample size determination for trials with one continuous and one binary
#' co-primary endpoint.
#'
#' @param n1 Sample size for group 1 (treatment group). If NULL, will be calculated.
#' @param n2 Sample size for group 2 (control group). If NULL, will be calculated.
#' @param delta Mean difference for the continuous endpoint
#' @param sd Common standard deviation for the continuous endpoint
#' @param p1 True response probability for the binary endpoint in group 1
#' @param p2 True response probability for the binary endpoint in group 2
#' @param rho Biserial correlation between the continuous endpoint and the
#'   latent continuous variable underlying the binary endpoint
#' @param power Target power (1 - beta). If NULL, will be calculated.
#' @param r Allocation ratio (n1/n2). Required when calculating sample size.
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param Test Test method for the binary endpoint: "AN" (asymptotic normal),
#'   "ANc" (with continuity correction), "AS" (arcsine), "ASc" (arcsine with
#'   continuity correction), or "Fisher" (Fisher's exact test)
#' @param nMC Number of Monte Carlo simulations when Test = "Fisher".
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
#' The function determines the operation mode based on which parameters are NULL.
#'
#' Exactly one of \{(n1, n2), (power, r)\} must be NULL.
#'
#' The biserial correlation rho represents the correlation between the observed
#' continuous endpoint and the latent continuous variable underlying the binary
#' endpoint.
#'
#' @examples
#' # Calculate power given sample sizes
#' twoCoprimary2MixedContinuousBinary(
#'   n1 = 100, n2 = 100,
#'   delta = 0.5, sd = 1,
#'   p1 = 0.6, p2 = 0.4,
#'   rho = 0.5,
#'   alpha = 0.025, Test = "AN"
#' )
#'
#' # Calculate sample size given target power
#' twoCoprimary2MixedContinuousBinary(
#'   delta = 0.5, sd = 1,
#'   p1 = 0.6, p2 = 0.4,
#'   rho = 0.5,
#'   power = 0.8, r = 1,
#'   alpha = 0.025, Test = "AN"
#' )
#'
#' @export
twoCoprimary2MixedContinuousBinary <- function(n1 = NULL, n2 = NULL,
                                               delta, sd, p1, p2, rho,
                                               power = NULL, r = NULL,
                                               alpha = 0.025, Test = "AN",
                                               nMC = 1e4) {

  # Count NULL parameters
  n_null <- sum(is.null(n1), is.null(n2))
  power_null <- is.null(power)
  r_null <- is.null(r)

  # Validate parameter combinations
  if (n_null == 0 && power_null) {
    # POWER CALCULATION MODE
    result <- power2MixedContinuousBinary(
      n1 = n1, n2 = n2,
      delta = delta, sd = sd,
      p1 = p1, p2 = p2,
      rho = rho, alpha = alpha,
      Test = Test, nMC = nMC
    )

  } else if (n_null == 2 && !power_null && !r_null) {
    # SAMPLE SIZE CALCULATION MODE
    beta <- 1 - power
    result <- ss2MixedContinuousBinary(
      delta = delta, sd = sd,
      p1 = p1, p2 = p2,
      rho = rho, r = r,
      alpha = alpha, beta = beta,
      Test = Test, nMC = nMC
    )

  } else {
    # INVALID PARAMETER COMBINATION
    stop("Exactly one of {(n1, n2), (power, r)} must be NULL.\n",
         "  - For power calculation: provide n1 and n2, set power = NULL\n",
         "  - For sample size calculation: provide power and r, set n1 = n2 = NULL")
  }

  return(result)
}
