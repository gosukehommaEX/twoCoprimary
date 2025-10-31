#' Sample Size Calculation for Two Co-Primary Endpoints (Count and Continuous)
#'
#' Calculates the required sample size for a two-arm superiority trial with one
#' overdispersed count co-primary endpoint and one continuous co-primary endpoint,
#' as described in Homma and Yoshida (2024).
#'
#' @param r0 Mean rate (events per unit time) for the control group (count endpoint)
#' @param r1 Mean rate (events per unit time) for the treatment group (count endpoint)
#' @param nu Common dispersion parameter for the negative binomial distribution (nu > 0).
#'   Can be a vector to calculate sample sizes for multiple dispersion values
#' @param t Common follow-up time period
#' @param mu0 Mean for the control group (continuous endpoint)
#' @param mu1 Mean for the treatment group (continuous endpoint)
#' @param sigma Common standard deviation for the continuous endpoint
#' @param r Allocation ratio (treatment:control = r:1, where r > 0)
#' @param rho0 Correlation between count and continuous outcomes for control group
#' @param rho1 Correlation between count and continuous outcomes for treatment group
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#'
#' @return A data frame with the following columns:
#'   \item{r0, r1}{Mean rates for count endpoint}
#'   \item{nu}{Dispersion parameter}
#'   \item{t}{Follow-up time}
#'   \item{mu0, mu1}{Means for continuous endpoint}
#'   \item{sigma}{Standard deviation for continuous endpoint}
#'   \item{r}{Allocation ratio}
#'   \item{rho0, rho1}{Correlations}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{n0}{Required sample size for control group}
#'   \item{n1}{Required sample size for treatment group}
#'   \item{N}{Total sample size (n0 + n1)}
#'
#' @details
#' This function uses linear extrapolation algorithm for efficient sample size
#' determination:
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas
#'
#' \strong{Step 2:} Set upper bound using adjusted power
#'
#' \strong{Step 3-4:} Iteratively refine using linear extrapolation:
#' \deqn{n_0^{new} = \frac{n_0^{(0)}(power^{(1)} - (1-\beta)) - n_0^{(1)}(power^{(0)} - (1-\beta))}
#'       {power^{(1)} - power^{(0)}}}
#'
#' The algorithm converges when \eqn{n_0^{(1)} = n_0^{(0)}}. If convergence is not
#' achieved within 100 iterations, a warning is issued.
#'
#' The correlation bounds are automatically checked using \code{\link{corrbound2MixedContinuousCount}}.
#' If the specified correlation is outside the valid range, a warning is issued.
#'
#' @references
#' Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical trials
#' with two co-primary endpoints including overdispersed count and continuous
#' outcomes. \emph{Pharmaceutical Statistics}, 23(1), 46-59.
#'
#' @examples
#' # Example 1: COPD trial scenario (Table 1 in Homma and Yoshida 2024)
#' # with nu = 0.8 and rho = 0.5
#' ssr2MixedContinuousCount(
#'   r0 = 1.25, r1 = 1.0, nu = 0.8, t = 1,
#'   mu0 = 0, mu1 = -50, sigma = 250,
#'   r = 1, rho0 = 0.8, rho1 = 0.8,
#'   alpha = 0.025, beta = 0.1
#' )
#'
#' # Example 2: Multiple dispersion parameters
#' ssr2MixedContinuousCount(
#'   r0 = 1.25, r1 = 1.0, nu = c(0.8, 1.0, 2.0), t = 1,
#'   mu0 = 0, mu1 = -50, sigma = 250,
#'   r = 1, rho0 = 0.6, rho1 = 0.6,
#'   alpha = 0.025, beta = 0.1
#' )
#'
#' # Example 3: Zero correlation (independent endpoints)
#' ssr2MixedContinuousCount(
#'   r0 = 1.5, r1 = 1.0, nu = 1.0, t = 1,
#'   mu0 = 0, mu1 = -40, sigma = 200,
#'   r = 1, rho0 = 0, rho1 = 0,
#'   alpha = 0.025, beta = 0.2
#' )
#'
#' @export
ssr2MixedContinuousCount <- function(r0, r1, nu, t, mu0, mu1, sigma, r,
                                     rho0, rho1, alpha, beta) {

  # Calculate z-scores
  za <- qnorm(alpha)
  zb <- qnorm(beta)
  zc <- qnorm(1 - sqrt(1 - beta))  # For upper bound initialization

  # Initialize results data frame
  results <- data.frame()

  # Loop over dispersion parameters (if multiple values provided)
  for (nu_val in nu) {

    # Check correlation bounds
    lambda0 <- r0 * t
    lambda1 <- r1 * t
    bounds0 <- corrbound2MixedContinuousCount(lambda0, nu_val, mu0, sigma, t)
    bounds1 <- corrbound2MixedContinuousCount(lambda1, nu_val, mu1, sigma, t)

    # Check if correlations are within bounds
    if (rho0 < bounds0[1] | rho0 > bounds0[2]) {
      warning(paste0("rho0 = ", rho0, " is outside valid range [",
                     round(bounds0[1], 4), ", ", round(bounds0[2], 4),
                     "] for nu = ", nu_val))
      next
    }
    if (rho1 < bounds1[1] | rho1 > bounds1[2]) {
      warning(paste0("rho1 = ", rho1, " is outside valid range [",
                     round(bounds1[1], 4), ", ", round(bounds1[2], 4),
                     "] for nu = ", nu_val))
      next
    }

    # Step 1: Calculate sample sizes for individual endpoints

    # Sample size for count endpoint alone
    n0_count <- ss1Count(r0, r1, nu_val, t, r, alpha, beta)[["n0"]]

    # Sample size for continuous endpoint alone
    n0_cont <- ss1Continuous(mu1 - mu0, sigma, r, alpha, beta)[["n2"]]

    # Calculate variance components
    Va <- (1 / t) * (1 / r0 + 1 / (r * r1)) + (1 + r) / (nu_val * r)
    V0 <- Va

    # Step 2: Initialize with minimum of individual sample sizes
    n0_min <- min(n0_count, n0_cont)

    # Calculate upper bound for search using adjusted power
    n0_count_upper <- ceiling(((za * sqrt(V0) + zc * sqrt(Va)) ^ 2) / (log(r1 / r0) ^ 2))
    n0_cont_upper <- ss1Continuous(mu1 - mu0, sigma, r, alpha, 1 - sqrt(1 - beta))[["n2"]]
    n0_max <- max(n0_count_upper, n0_cont_upper)

    # Step 3: Initialize two sample sizes for linear extrapolation
    n0_0 <- n0_min
    n0_1 <- n0_max

    # Step 4: Iterative refinement using linear extrapolation
    max_iter <- 100
    iter <- 0

    while (n0_1 - n0_0 != 0 && iter < max_iter) {
      iter <- iter + 1

      # Calculate sample sizes for treatment group
      n1_0 <- ceiling(r * n0_0)
      n1_1 <- ceiling(r * n0_1)

      # Calculate power at two candidate sample sizes
      power_n0_0 <- power2MixedContinuousCount(n0_0, n1_0, r0, r1, nu_val, t, mu0, mu1,
                                               sigma, rho0, rho1, alpha)[["powerCoprimary"]]
      power_n0_1 <- power2MixedContinuousCount(n0_1, n1_1, r0, r1, nu_val, t, mu0, mu1,
                                               sigma, rho0, rho1, alpha)[["powerCoprimary"]]

      # Check if powers are different enough for extrapolation
      if (abs(power_n0_1 - power_n0_0) < 1e-10) {
        warning("Powers are too similar for extrapolation. Using n0_0.")
        break
      }

      # Linear extrapolation to find sample size that achieves target power
      n0_updated <- (n0_0 * (power_n0_1 - (1 - beta)) - n0_1 * (power_n0_0 - (1 - beta))) /
        (power_n0_1 - power_n0_0)

      # Update values for next iteration
      n0_1 <- n0_0
      n0_0 <- ceiling(n0_updated)

      # Safety check for negative or very small sample sizes
      if (n0_0 < 10) {
        warning("Sample size too small. Setting to minimum of 10.")
        n0_0 <- 10
        break
      }
    }

    if (iter >= max_iter) {
      warning(paste0("Maximum iterations reached for nu = ", nu_val,
                     ". Results may not have converged."))
    }

    # Final sample sizes
    n0 <- n0_0
    n1 <- ceiling(r * n0)
    N <- n0 + n1

    # Append to results
    result_row <- data.frame(
      r0, r1, nu = nu_val, t, mu0, mu1, sigma, r,
      rho0, rho1, alpha, beta, n0, n1, N
    )
    results <- rbind(results, result_row)
  }

  return(results)
}
