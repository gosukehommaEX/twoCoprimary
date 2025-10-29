#' Sample Size Calculation for Two Co-Primary Continuous Endpoints (Approximate)
#'
#' Calculates the required sample size for a two-arm superiority trial with two
#' co-primary continuous endpoints using linear extrapolation approach, as
#' described in Hamasaki et al. (2013) and Sozu et al. (2011).
#'
#' @param delta1 Mean difference for the first endpoint
#' @param delta2 Mean difference for the second endpoint
#' @param sd1 Common standard deviation for the first endpoint
#' @param sd2 Common standard deviation for the second endpoint
#' @param rho Common correlation between the two outcomes
#' @param r Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#' @param known_var Logical value indicating whether variance is known (TRUE) or
#'   unknown (FALSE). If TRUE, power is calculated analytically; otherwise,
#'   Monte Carlo simulation is used for unknown variance
#' @param nMC Number of Monte Carlo simulations when known_var = FALSE (default is 10000)
#'
#' @return A data frame with the following columns:
#'   \item{delta1, delta2}{Mean differences}
#'   \item{sd1, sd2}{Standard deviations}
#'   \item{rho}{Correlation}
#'   \item{r}{Allocation ratio}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{known_var}{Variance assumption}
#'   \item{nMC}{Number of Monte Carlo simulations (NA if known_var = TRUE)}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{n}{Total sample size (n1 + n2)}
#'
#' @details
#' This function uses a linear extrapolation algorithm (Hamasaki et al. 2013) to
#' efficiently determine the required sample size that achieves power of at least
#' \eqn{1-\beta} for both co-primary endpoints.
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Initialize two sample sizes: \eqn{n_{2,0}} using target power \eqn{1-\beta}
#'         and \eqn{n_{2,1}} using adjusted power \eqn{1-\sqrt{1-\beta}}
#'   \item Calculate power at both sample sizes
#'   \item Use linear extrapolation to estimate sample size achieving target power:
#'         \deqn{n_{2,new} = \frac{n_{2,0}(power_{n_{2,1}} - (1-\beta)) - n_{2,1}(power_{n_{2,0}} - (1-\beta))}{power_{n_{2,1}} - power_{n_{2,0}}}}
#'   \item Update values and repeat until convergence
#' }
#'
#' For known variance, the standardized test statistics are:
#' \deqn{Z_k = \frac{\delta_k}{\sigma_k \sqrt{1/n_1 + 1/n_2}}}
#' For unknown variance, t-statistics with \eqn{\nu = n_1 + n_2 - 2} degrees of
#' freedom are used, and power is calculated using Monte Carlo simulation following
#' Sozu et al. (2011).
#'
#' @references
#' Hamasaki, T., Sugimoto, T., Evans, S. R., & Sozu, T. (2013). Sample size
#' determination for clinical trials with co-primary outcomes: Exponential event
#' times. \emph{Pharmaceutical Statistics}, 12(1), 28-34.
#'
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2011). Sample size determination in
#' superiority clinical trials with multiple co-primary correlated endpoints.
#' \emph{Journal of Biopharmaceutical Statistics}, 21(4), 650-668.
#'
#' @examples
#' # Sample size calculation with known variance
#' ss2ContinuousApprox(
#'   delta1 = 0.2,
#'   delta2 = 0.2,
#'   sd1 = 1,
#'   sd2 = 1,
#'   rho = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   known_var = TRUE
#' )
#'
#' # Sample size calculation with unequal allocation
#' ss2ContinuousApprox(
#'   delta1 = 0.3,
#'   delta2 = 0.25,
#'   sd1 = 1,
#'   sd2 = 1,
#'   rho = 0.3,
#'   r = 2,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   known_var = TRUE
#' )
#'
#' \donttest{
#' # Sample size calculation with unknown variance
#' ss2ContinuousApprox(
#'   delta1 = 0.5,
#'   delta2 = 0.4,
#'   sd1 = 1,
#'   sd2 = 1,
#'   rho = 0.4,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   known_var = FALSE,
#'   nMC = 10000
#' )
#' }
#'
#' @export
#' @importFrom stats qnorm
ss2ContinuousApprox <- function(delta1, delta2, sd1, sd2, rho, r, alpha, beta,
                                known_var = TRUE, nMC = 1e+4) {

  # Input validation
  if (delta1 <= 0 | delta2 <= 0) {
    stop("delta1 and delta2 must be positive")
  }
  if (sd1 <= 0 | sd2 <= 0) {
    stop("sd1 and sd2 must be positive")
  }
  if (abs(rho) >= 1) {
    stop("rho must be in (-1, 1)")
  }
  if (r <= 0) {
    stop("r must be positive")
  }
  if (alpha <= 0 | alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }
  if (beta <= 0 | beta >= 1) {
    stop("beta must be in (0, 1)")
  }

  # Helper function to calculate sample size for single endpoint
  ss1Continuous <- function(delta, sd, r, alpha, beta) {
    z_alpha <- qnorm(1 - alpha)
    z_beta <- qnorm(1 - beta)
    # n2 = [(z_alpha + z_beta) * sd / delta]^2 * (1 + 1/r)
    n2 <- ceiling(((z_alpha + z_beta) * sd / delta) ^ 2 * (1 + 1 / r))
    return(n2)
  }

  # Step 1: Initialize sample sizes using single endpoint formulas
  # Use the maximum sample size from both endpoints
  n2_0 <- max(
    ss1Continuous(delta1, sd1, r, alpha, beta),
    ss1Continuous(delta2, sd2, r, alpha, beta)
  )

  # For the second initial value, use adjusted target power
  # This provides a bracket for the linear extrapolation
  beta_adj <- 1 - (1 - beta) ^ (1 / 2)
  n2_1 <- max(
    ss1Continuous(delta1, sd1, r, alpha, beta_adj),
    ss1Continuous(delta2, sd2, r, alpha, beta_adj)
  )

  # Step 2-4: Iterative refinement using linear extrapolation
  while (n2_1 - n2_0 != 0) {

    # Calculate sample sizes for group 1
    n1_0 <- ceiling(r * n2_0)
    n1_1 <- ceiling(r * n2_1)

    # Calculate power at two candidate sample sizes
    power_n2_0 <- power2Continuous(n1_0, n2_0, delta1, delta2, sd1, sd2,
                                   rho, alpha, known_var, nMC)$powerCoprimary
    power_n2_1 <- power2Continuous(n1_1, n2_1, delta1, delta2, sd1, sd2,
                                   rho, alpha, known_var, nMC)$powerCoprimary

    # Linear extrapolation to find sample size that achieves target power
    # This solves: power = (1 - beta) for n2
    n2_updated <- (n2_0 * (power_n2_1 - (1 - beta)) - n2_1 * (power_n2_0 - (1 - beta))) /
      (power_n2_1 - power_n2_0)

    # Update values for next iteration
    n2_1 <- n2_0
    n2_0 <- ceiling(n2_updated)
  }

  # Final sample sizes
  n2 <- n2_0
  n1 <- ceiling(r * n2)
  n <- n1 + n2

  # Set nMC to NA for known variance
  if (known_var) {
    nMC <- NA
  }

  # Return results as a data frame
  result <- data.frame(
    delta1, delta2, sd1, sd2, rho, r, alpha, beta, known_var, nMC,
    n1, n2, n
  )
  return(result)
}
