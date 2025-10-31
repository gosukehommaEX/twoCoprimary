#' Sample Size Calculation for Two Co-Primary Time-to-Event Endpoints
#'
#' Calculates the required sample size for a two-arm superiority trial with two
#' co-primary time-to-event endpoints using the linear extrapolation approach
#' described in Hamasaki et al. (2013).
#'
#' @param lambdaT1 Hazard rate for the first endpoint in test group (λT1 > 0)
#' @param lambdaT2 Hazard rate for the second endpoint in test group (λT2 > 0)
#' @param lambdaC1 Hazard rate for the first endpoint in control group (λC1 > 0)
#' @param lambdaC2 Hazard rate for the second endpoint in control group (λC2 > 0)
#' @param tau0 Accrual period (recruitment time in years)
#' @param tau Follow-up period after recruitment ends (years)
#' @param rho Correlation between the two endpoints (0 ≤ rho ≤ 1)
#' @param r Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#' @param copula Character string specifying copula type. One of:
#'   \itemize{
#'     \item \code{"clayton"}: Clayton copula (default)
#'     \item \code{"gumbel"}: Gumbel/positive stable copula
#'     \item \code{"frank"}: Frank copula
#'   }
#' @param censor Character string specifying censoring approach. One of:
#'   \itemize{
#'     \item \code{"variance"}: Censoring incorporated into variance only (default)
#'     \item \code{"both"}: Censoring incorporated into both variance and correlation
#'   }
#'
#' @return A data frame with columns for parameters and required sample sizes
#'
#' @details
#' This function uses the linear extrapolation approach to find the required
#' sample size iteratively.
#'
#' **Algorithm:**
#' \enumerate{
#'   \item Initialize with maximum of single endpoint sample sizes
#'   \item Calculate two initial values with target and adjusted power
#'   \item Iteratively refine using linear extrapolation until convergence
#' }
#'
#' @references
#' Hamasaki, T., Sugimoto, T., Evans, S. R., & Sozu, T. (2013). Sample size
#' determination for clinical trials with co-primary outcomes: exponential event
#' times. \emph{Pharmaceutical Statistics}, 12(1), 28-34.
#'
#' @examples
#' ss2Survival(
#'   lambdaT1 = log(2) / 5, lambdaT2 = log(2) / 5,
#'   lambdaC1 = log(2) / 3, lambdaC2 = log(2) / 3,
#'   tau0 = 2, tau = 5,
#'   rho = 0.5, r = 1,
#'   alpha = 0.025, beta = 0.2,
#'   copula = "clayton", censor = "variance"
#' )
#'
#' @export
#' @importFrom stats qnorm
ss2Survival <- function(lambdaT1, lambdaT2, lambdaC1, lambdaC2,
                        tau0, tau, rho, r, alpha, beta,
                        copula = "clayton", censor = "variance") {

  # Input validation
  if (lambdaT1 <= 0 || lambdaT2 <= 0 || lambdaC1 <= 0 || lambdaC2 <= 0) {
    stop("All hazard rates must be positive")
  }
  if (tau0 <= 0) stop("tau0 must be positive")
  if (tau <= tau0) stop("tau must be greater than tau0")
  if (rho < 0 || rho > 1) stop("rho must be between 0 and 1")
  if (r <= 0) stop("r must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)")
  if (beta <= 0 || beta >= 1) stop("beta must be in (0, 1)")
  if (!copula %in% c("clayton", "gumbel", "frank")) {
    stop("copula must be 'clayton', 'gumbel', or 'frank'")
  }
  if (!censor %in% c("variance", "both")) {
    stop("censor must be 'variance' or 'both'")
  }

  # Step 1: Initialize with single endpoint sample sizes
  n2_endpoint1 <- ss1Survival(
    lambdaT = lambdaT1, lambdaC = lambdaC1,
    tau0 = tau0, tau = tau, r = r,
    alpha = alpha, beta = beta, censor = censor
  )[["n2"]]

  n2_endpoint2 <- ss1Survival(
    lambdaT = lambdaT2, lambdaC = lambdaC2,
    tau0 = tau0, tau = tau, r = r,
    alpha = alpha, beta = beta, censor = censor
  )[["n2"]]

  n2_0 <- max(n2_endpoint1, n2_endpoint2)

  # Adjusted target power for bracketing
  beta_adj <- 1 - (1 - beta) ^ (1 / 2)

  n2_endpoint1_adj <- ss1Survival(
    lambdaT = lambdaT1, lambdaC = lambdaC1,
    tau0 = tau0, tau = tau, r = r,
    alpha = alpha, beta = beta_adj, censor = censor
  )[["n2"]]

  n2_endpoint2_adj <- ss1Survival(
    lambdaT = lambdaT2, lambdaC = lambdaC2,
    tau0 = tau0, tau = tau, r = r,
    alpha = alpha, beta = beta_adj, censor = censor
  )[["n2"]]

  n2_1 <- max(n2_endpoint1_adj, n2_endpoint2_adj)

  # Step 2-4: Iterative refinement
  max_iter <- 50
  iter <- 0

  while (n2_1 - n2_0 != 0 && iter < max_iter) {
    iter <- iter + 1

    n1_0 <- ceiling(r * n2_0)
    n1_1 <- ceiling(r * n2_1)

    # Calculate powers
    power_n2_0 <- power2Survival(
      n1 = n1_0, n2 = n2_0,
      lambdaT1 = lambdaT1, lambdaT2 = lambdaT2,
      lambdaC1 = lambdaC1, lambdaC2 = lambdaC2,
      tau0 = tau0, tau = tau, rho = rho, alpha = alpha,
      copula = copula, censor = censor
    )[["powerCoprimary"]]

    power_n2_1 <- power2Survival(
      n1 = n1_1, n2 = n2_1,
      lambdaT1 = lambdaT1, lambdaT2 = lambdaT2,
      lambdaC1 = lambdaC1, lambdaC2 = lambdaC2,
      tau0 = tau0, tau = tau, rho = rho, alpha = alpha,
      copula = copula, censor = censor
    )[["powerCoprimary"]]

    # Linear extrapolation
    n2_updated <- (n2_0 * (power_n2_1 - (1 - beta)) -
                     n2_1 * (power_n2_0 - (1 - beta))) /
      (power_n2_1 - power_n2_0)

    n2_1 <- n2_0
    n2_0 <- ceiling(n2_updated)

    if (n2_0 < 10) n2_0 <- 10
  }

  if (iter >= max_iter) {
    warning("Maximum iterations reached. Results may not have converged.")
  }

  # Final sample sizes
  n2 <- n2_0
  n1 <- ceiling(r * n2)
  n <- n1 + n2

  # Return result
  result <- data.frame(
    lambdaT1 = lambdaT1, lambdaT2 = lambdaT2,
    lambdaC1 = lambdaC1, lambdaC2 = lambdaC2,
    tau0 = tau0, tau = tau, rho = rho, r = r,
    alpha = alpha, beta = beta,
    copula = copula, censor = censor,
    n1 = n1, n2 = n2, n = n
  )
  return(result)
}
