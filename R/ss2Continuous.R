#' Sample Size Calculation for Two Co-Primary Continuous Endpoints (Approximate)
#'
#' Calculates the required sample size for a two-arm superiority trial with two
#' co-primary continuous endpoints using sequential search algorithm.
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
#'   unknown (FALSE). If TRUE, power is calculated analytically; if FALSE,
#'   Monte Carlo simulation is used
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
#'   \item{N}{Total sample size (n1 + n2)}
#'
#' @details
#' This function uses a sequential search algorithm (Homma and Yoshida 2025,
#' Algorithm 1) for both known and unknown variance cases:
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas.
#'
#' \strong{Step 2:} Use sequential search:
#' \itemize{
#'   \item Calculate power at initial sample size
#'   \item If power >= target: decrease n2 until power < target, then add 1 back
#'   \item If power < target: increase n2 until power >= target
#' }
#'
#' \strong{Step 3:} Return final sample sizes.
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
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
#'
#' @examples
#' # Sample size calculation with known variance
#' ss2Continuous(
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
#' # Sample size calculation with unknown variance (uses sequential search)
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
  if (length(delta1) != 1 || length(delta2) != 1 || length(sd1) != 1 ||
      length(sd2) != 1 || length(rho) != 1 || length(r) != 1 ||
      length(alpha) != 1 || length(beta) != 1) {
    stop("All parameters must be scalar values")
  }
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
  if (!is.logical(known_var)) {
    stop("known_var must be logical (TRUE or FALSE)")
  }

  # Step 1: Initialize sample size using single endpoint formulas
  n2_initial <- max(
    ss1Continuous(delta1, sd1, r, alpha, beta)[["n2"]],
    ss1Continuous(delta2, sd2, r, alpha, beta)[["n2"]]
  )

  # Step 2: Sequential search to find minimum sample size
  result_ss <- .ss_sequential_search(
    initial_n2 = n2_initial,
    r = r,
    target_power = 1 - beta,
    power_fun = power2Continuous,
    delta1 = delta1, delta2 = delta2,
    sd1 = sd1, sd2 = sd2,
    rho = rho, alpha = alpha,
    known_var = known_var, nMC = nMC
  )

  n1 <- result_ss$n1
  n2 <- result_ss$n2
  N <- result_ss$N

  # Set nMC to NA for known variance
  if (known_var) {
    nMC <- NA
  }

  # Return results as a data frame
  result <- data.frame(
    delta1, delta2, sd1, sd2, rho, r, alpha, beta, known_var, nMC,
    n1, n2, N
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
