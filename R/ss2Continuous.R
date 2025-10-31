#' Sample Size Calculation for Two Co-Primary Continuous Endpoints
#'
#' Determines the sample size for a two-arm superiority trial with two co-primary
#' continuous endpoints, to achieve a specified power at a given significance level.
#'
#' @param delta1 Mean difference for the first endpoint (group 1 - group 2)
#' @param delta2 Mean difference for the second endpoint (group 1 - group 2)
#' @param sd1 Common standard deviation for the first endpoint
#' @param sd2 Common standard deviation for the second endpoint
#' @param rho Common correlation between the two outcomes
#' @param r Allocation ratio n1/n2 where n1 is sample size for group 1
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Type II error rate (typically 0.1 or 0.2). Power = 1 - beta
#' @param known_var Logical value indicating whether variance is known (TRUE) or
#'   unknown (FALSE). If TRUE, linear extrapolation is used; if FALSE, sequential
#'   search is used for sample size determination
#' @param nMC Number of Monte Carlo simulations when known_var = FALSE (default is 10000)
#'
#' @return A data frame with the following columns:
#'   \item{delta1}{Mean difference for endpoint 1}
#'   \item{delta2}{Mean difference for endpoint 2}
#'   \item{sd1}{Standard deviation for endpoint 1}
#'   \item{sd2}{Standard deviation for endpoint 2}
#'   \item{rho}{Correlation between endpoints}
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
#' This function uses two different approaches depending on the variance assumption:
#'
#' **Known Variance:** Uses linear extrapolation algorithm for efficient sample size
#' determination. The algorithm iteratively refines the sample size estimate:
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas
#'
#' \strong{Step 2:} Set upper bound using adjusted power
#'
#' \strong{Step 3-4:} Iteratively refine using linear extrapolation:
#' \deqn{n_2^{new} = \frac{n_2^{(0)}(power^{(1)} - (1-\beta)) - n_2^{(1)}(power^{(0)} - (1-\beta))}
#'       {power^{(1)} - power^{(0)}}}
#'
#' The algorithm converges when \eqn{n_2^{(1)} = n_2^{(0)}}.
#'
#' **Unknown Variance:** Uses sequential search, as Monte Carlo simulation has random
#' variation that can prevent reliable convergence with linear extrapolation.
#' This approach is more stable for Monte Carlo simulation where power estimates have
#' random variation.
#'
#' For known variance, the standardized test statistics are:
#' \deqn{Z_k = \frac{\delta_k}{\sigma_k \sqrt{1/n_1 + 1/n_2}}}
#' For unknown variance, t-statistics with \eqn{\nu = n_1 + n_2 - 2} degrees of
#' freedom are used, and power is calculated using Monte Carlo simulation following
#' Sozu et al. (2011).
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2011). Sample size determination in
#' superiority clinical trials with multiple co-primary correlated endpoints.
#' \emph{Journal of Biopharmaceutical Statistics}, 21(4), 650-668.
#'
#' @examples
#' # Example 1: Known variance with rho = 0.5 (Table 1 in Sozu et al. 2011)
#' ss2Continuous(
#'   delta1 = 0.2,
#'   delta2 = 0.2,
#'   sd1 = 1,
#'   sd2 = 1,
#'   rho = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   known_var = TRUE
#' )
#'
#' # Example 2: Known variance with unequal allocation
#' ss2Continuous(
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
#' # Example 3: Unknown variance (uses sequential search)
#' ss2Continuous(
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
ss2Continuous <- function(delta1, delta2, sd1, sd2, rho, r, alpha, beta,
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

  if (known_var) {
    # ===== KNOWN VARIANCE: Use linear extrapolation =====

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
    max_iter <- 100
    iter <- 0

    while (n2_1 - n2_0 != 0 && iter < max_iter) {
      iter <- iter + 1

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

    if (iter >= max_iter) {
      warning("Maximum iterations reached. Results may not have converged.")
    }

    # Final sample sizes
    n2 <- n2_0
    n1 <- ceiling(r * n2)
    n <- n1 + n2

    # Set nMC to NA for known variance
    nMC <- NA

  } else {
    # ===== UNKNOWN VARIANCE: Use sequential search =====
    # Monte Carlo simulation has random variation, so linear extrapolation
    # may not converge reliably. Instead, use sequential search.

    # Step 1: Initialize with sample size from single endpoint formulas
    n2 <- max(
      ss1Continuous(delta1, sd1, r, alpha, beta),
      ss1Continuous(delta2, sd2, r, alpha, beta)
    )
    n1 <- ceiling(r * n2)

    # Step 2: Calculate power at initial sample size
    power <- power2Continuous(n1, n2, delta1, delta2, sd1, sd2,
                              rho, alpha, known_var, nMC)$powerCoprimary

    # Step 3: Sequential search - increment n2 by 1 until target power is achieved
    while (power < 1 - beta) {
      n2 <- n2 + 1
      n1 <- ceiling(r * n2)
      power <- power2Continuous(n1, n2, delta1, delta2, sd1, sd2,
                                rho, alpha, known_var, nMC)$powerCoprimary
    }

    # Final sample sizes
    n <- n1 + n2
  }

  # Return results as a data frame
  result <- data.frame(
    delta1, delta2, sd1, sd2, rho, r, alpha, beta, known_var, nMC,
    n1, n2, n
  )
  return(result)
}
