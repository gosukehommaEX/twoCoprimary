#' Sample Size Calculation for Two Co-Primary Endpoints (Count and Continuous)
#'
#' Calculates the required sample size for a two-arm superiority trial with one
#' overdispersed count and one continuous co-primary endpoints
#' using the linear extrapolation approach.
#'
#' @param r1 Mean rate (events per unit time) for the treatment group (count endpoint)
#' @param r2 Mean rate (events per unit time) for the control group (count endpoint)
#' @param nu Common dispersion parameter for the negative binomial distribution (nu > 0)
#' @param t Common follow-up time period
#' @param mu1 Mean for group 1 (continuous endpoint)
#' @param mu2 Mean for group 2 (continuous endpoint)
#' @param sd Common standard deviation for the continuous endpoint
#' @param r Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r > 0)
#' @param rho1 Correlation between count and continuous outcomes for treatment group
#' @param rho2 Correlation between count and continuous outcomes for control group
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#'
#' @return A data frame with the following columns:
#'   \item{r1}{Mean rate in group 1 for count endpoint}
#'   \item{r2}{Mean rate in group 2 for count endpoint}
#'   \item{nu}{Dispersion parameter}
#'   \item{t}{Follow-up time}
#'   \item{mu1}{Mean in group 1 for continuous endpoint}
#'   \item{mu2}{Mean in group 2 for continuous endpoint}
#'   \item{sd}{Standard deviation for continuous endpoint}
#'   \item{r}{Allocation ratio}
#'   \item{rho1}{Correlation for group 1}
#'   \item{rho2}{Correlation for group 2}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{N}{Total sample size (n1 + n2)}
#'
#' @details
#' This function uses the linear extrapolation approach to find the required
#' sample size iteratively:
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas.
#' Two initial values are calculated:
#' \itemize{
#'   \item n2_0: Based on target power 1 - beta
#'   \item n2_1: Based on adjusted power 1 - sqrt(1-beta) to provide a bracket
#' }
#'
#' \strong{Step 2-4:} Iteratively refine using linear extrapolation until convergence:
#' \deqn{n_2^{new} = \frac{n_2^{(0)}(power^{(1)} - (1-\beta)) - n_2^{(1)}(power^{(0)} - (1-\beta))}
#'       {power^{(1)} - power^{(0)}}}
#'
#' The algorithm converges when \eqn{n_2^{(1)} = n_2^{(0)}}.
#'
#' The correlation bounds are automatically checked using \code{\link{corrbound2MixedContinuousCount}}.
#' If the specified correlation is outside the valid range, an error is issued.
#'
#' @references
#' Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical trials
#' with two co-primary endpoints including overdispersed count and continuous
#' outcomes. \emph{Pharmaceutical Statistics}, 23(1), 46-59.
#'
#' @examples
#' # Sample size calculation for COPD trial scenario
#' ss2MixedCountContinuous(
#'   r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
#'   mu1 = -50, mu2 = 0, sd = 250,
#'   r = 1, rho1 = 0.5, rho2 = 0.5,
#'   alpha = 0.025, beta = 0.1
#' )
#'
#' # Multiple dispersion parameters
#' ss2MixedCountContinuous(
#'   r1 = 1.0, r2 = 1.25, nu = 1, t = 1,
#'   mu1 = -50, mu2 = 0, sd = 250,
#'   r = 1, rho1 = 0.6, rho2 = 0.6,
#'   alpha = 0.025, beta = 0.1
#' )
#'
#' # Zero correlation (independent endpoints)
#' ss2MixedCountContinuous(
#'   r1 = 1.0, r2 = 1.5, nu = 1, t = 1,
#'   mu1 = -40, mu2 = 0, sd = 200,
#'   r = 1, rho1 = 0, rho2 = 0,
#'   alpha = 0.025, beta = 0.2
#' )
#'
#' @export
ss2MixedCountContinuous <- function(r1, r2, nu, t, mu1, mu2, sd, r,
                                    rho1, rho2, alpha, beta) {

  # Input validation
  if (length(r1) != 1 || length(r2) != 1 || length(nu) != 1 ||
      length(t) != 1 || length(mu1) != 1 || length(mu2) != 1 ||
      length(sd) != 1 || length(r) != 1 || length(rho1) != 1 ||
      length(rho2) != 1 || length(alpha) != 1 || length(beta) != 1) {
    stop("All parameters must be scalar values")
  }
  if (r1 <= 0 || r2 <= 0) {
    stop("r1 and r2 must be positive")
  }
  if (nu <= 0) {
    stop("nu must be positive")
  }
  if (t <= 0) {
    stop("t must be positive")
  }
  if (sd <= 0) {
    stop("sd must be positive")
  }
  if (r <= 0) {
    stop("r must be positive")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }
  if (beta <= 0 || beta >= 1) {
    stop("beta must be in (0, 1)")
  }

  # Step 1: Initialize sample sizes using single endpoint formula
  # Calculate sample size for each endpoint separately, then take the maximum
  n2_0 <- max(
    ss1Count(r1, r2, nu, t, r, alpha, beta)[["n2"]],
    ss1Continuous(-(mu1 - mu2), sd, r, alpha, beta)[["n2"]]
  )

  # For the second initial value, use adjusted target power
  # This provides a bracket for the linear extrapolation
  # The adjusted power 1 - sqrt(1-beta) is typically higher than 1 - beta
  n2_1 <- max(
    ss1Count(r1, r2, nu, t, r, alpha, 1 - (1 - beta) ^ (1 / 2))[["n2"]],
    ss1Continuous(-(mu1 - mu2), sd, r, alpha, 1 - (1 - beta) ^ (1 / 2))[["n2"]]
  )

  # Step 2-4: Iterative refinement using linear extrapolation
  while (n2_1 - n2_0 != 0) {

    # Calculate sample sizes for group 1
    n1_0 <- ceiling(r * n2_0)
    n1_1 <- ceiling(r * n2_1)

    # Calculate power at two candidate sample sizes
    power_n2_0 <- power2MixedCountContinuous(
      n1_0, n2_0, r1, r2, nu, t, mu1, mu2, sd, rho1, rho2, alpha
    )[["powerCoprimary"]]
    power_n2_1 <- power2MixedCountContinuous(
      n1_1, n2_1, r1, r2, nu, t, mu1, mu2, sd, rho1, rho2, alpha
    )[["powerCoprimary"]]

    # Linear extrapolation to find sample size that achieves target power
    # This solves: power = (1 - beta) for n2
    n2_updated <- '/'(
      n2_0 * (power_n2_1 - (1 - beta)) - n2_1 * (power_n2_0 - (1 - beta)),
      power_n2_1 - power_n2_0
    )

    # Update values for next iteration
    n2_1 <- n2_0
    n2_0 <- ceiling(n2_updated)
  }

  # Final sample sizes
  n2 <- n2_0
  n1 <- ceiling(r * n2)
  N <- n1 + n2

  # Return results as a data frame
  result <- data.frame(
    r1, r2, nu, t, mu1, mu2, sd, r, rho1, rho2, alpha, beta,
    n1, n2, N
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
