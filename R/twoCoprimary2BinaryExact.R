#' Unified Interface for Two Co-Primary Binary Endpoints (Exact Methods)
#'
#' This function provides a unified interface for both power calculation and
#' sample size determination for two co-primary binary endpoints using exact
#' inference methods.
#'
#' @param n1 Sample size for group 1 (treatment group). If NULL, will be calculated.
#' @param n2 Sample size for group 2 (control group). If NULL, will be calculated.
#' @param p11 True response probability for endpoint 1 in group 1
#' @param p12 True response probability for endpoint 2 in group 1
#' @param p21 True response probability for endpoint 1 in group 2
#' @param p22 True response probability for endpoint 2 in group 2
#' @param rho1 Correlation between endpoints 1 and 2 in group 1
#' @param rho2 Correlation between endpoints 1 and 2 in group 2
#' @param power Target power (1 - beta). If NULL, will be calculated.
#' @param r Allocation ratio (n1/n2). Required when calculating sample size.
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param Test Test method: "Fisher" (Fisher's exact test), "Chisq"
#'   (Chi-squared test), "Z-pooled" (Z-pooled exact unconditional test),
#'   or "Boschloo" (Boschloo's exact unconditional test)
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
#' Exactly one of {(n1, n2), (power, r)} must be NULL.
#'
#' Note: Exact methods are computationally intensive and may take considerable
#' time, especially for large sample sizes.
#'
#' @examples
#' \donttest{
#' # Calculate power given sample sizes
#' twoCoprimary2BinaryExact(
#'   n1 = 50, n2 = 50,
#'   p11 = 0.5, p12 = 0.4,
#'   p21 = 0.3, p22 = 0.2,
#'   rho1 = 0.5, rho2 = 0.5,
#'   alpha = 0.025, Test = "Fisher"
#' )
#'
#' # Calculate sample size given target power
#' twoCoprimary2BinaryExact(
#'   p11 = 0.5, p12 = 0.4,
#'   p21 = 0.3, p22 = 0.2,
#'   rho1 = 0.5, rho2 = 0.5,
#'   power = 0.8, r = 1,
#'   alpha = 0.025, Test = "Chisq"
#' )
#' }
#'
#' @export
twoCoprimary2BinaryExact <- function(n1 = NULL, n2 = NULL,
                                     p11, p12, p21, p22,
                                     rho1, rho2,
                                     power = NULL, r = NULL,
                                     alpha = 0.025, Test = "Fisher") {

  # Count NULL parameters
  n_null <- sum(is.null(n1), is.null(n2))
  power_null <- is.null(power)
  r_null <- is.null(r)

  # Validate parameter combinations
  if (n_null == 0 && power_null) {
    # POWER CALCULATION MODE
    result <- power2BinaryExact(
      n1 = n1, n2 = n2,
      p11 = p11, p12 = p12,
      p21 = p21, p22 = p22,
      rho1 = rho1, rho2 = rho2,
      alpha = alpha, Test = Test
    )

  } else if (n_null == 2 && !power_null && !r_null) {
    # SAMPLE SIZE CALCULATION MODE
    beta <- 1 - power
    result <- ss2BinaryExact(
      p11 = p11, p12 = p12,
      p21 = p21, p22 = p22,
      rho1 = rho1, rho2 = rho2,
      r = r, alpha = alpha, beta = beta,
      Test = Test
    )

  } else {
    # INVALID PARAMETER COMBINATION
    stop("Exactly one of {(n1, n2), (power, r)} must be NULL.\n",
         "  - For power calculation: provide n1 and n2, set power = NULL\n",
         "  - For sample size calculation: provide power and r, set n1 = n2 = NULL")
  }

  return(result)
}
