#' Sample Size Calculation for Two Co-Primary Endpoints: One Count and One Continuous
#'
#' Determines the sample size for a two-arm superiority trial with two co-primary
#' endpoints where one is a count (negative binomial) and one is continuous (normal),
#' to achieve a specified power at a given significance level.
#'
#' @param r1 Mean count rate in group 1 for the count endpoint
#' @param r2 Mean count rate in group 2 for the count endpoint
#' @param nu Dispersion parameter for the negative binomial distribution (nu > 0).
#'   Smaller values indicate greater overdispersion
#' @param t Follow-up time period
#' @param mu1 Mean for group 1 for the continuous endpoint
#' @param mu2 Mean for group 2 for the continuous endpoint
#' @param sd Common standard deviation for the continuous endpoint
#' @param r Allocation ratio n1/n2 where n1 is sample size for group 1
#' @param rho1 Correlation between count and continuous endpoints in group 1
#' @param rho2 Correlation between count and continuous endpoints in group 2
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Type II error rate (typically 0.1 or 0.2). Power = 1 - beta
#'
#' @return A data frame with the following columns:
#'   \item{r1, r2}{Count rates}
#'   \item{nu}{Dispersion parameter}
#'   \item{t}{Follow-up time}
#'   \item{mu1, mu2}{Means for continuous endpoint}
#'   \item{sd}{Standard deviation for continuous endpoint}
#'   \item{r}{Allocation ratio}
#'   \item{rho1, rho2}{Correlations}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{N}{Total sample size (n1 + n2)}
#'
#' @details
#' This function implements the sample size calculation for mixed count-continuous
#' co-primary endpoints following the methodology in Homma and Yoshida (2024).
#'
#' The sequential search algorithm (Homma and Yoshida 2025, Algorithm 1) is used:
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
#' \strong{Negative Binomial Distribution:}
#' The count endpoint follows a negative binomial distribution NB(lambda, nu) where:
#' \itemize{
#'   \item lambda = r * t is the mean count
#'   \item nu is the dispersion parameter
#'   \item Variance = lambda + lambda^2 / nu
#' }
#'
#' \strong{Correlation:}
#' The correlations rho1 and rho2 must satisfy feasibility constraints that depend
#' on the parameters. Use \code{\link{corrbound2MixedCountContinuous}} to check
#' valid correlation bounds.
#'
#' @references
#' Homma, G., & Yoshida, T. (2024). Sample size calculation for count and
#' continuous multiple co-primary endpoints. \emph{Pharmaceutical Statistics},
#' 23(3), 372-388.
#'
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
#'
#' @examples
#' # Sample size calculation for count and continuous endpoints
#' ss2MixedCountContinuous(
#'   r1 = 1.0,
#'   r2 = 1.25,
#'   nu = 0.8,
#'   t = 1,
#'   mu1 = -50,
#'   mu2 = 0,
#'   sd = 250,
#'   r = 1,
#'   rho1 = 0.4,
#'   rho2 = 0.4,
#'   alpha = 0.025,
#'   beta = 0.2
#' )
#'
#' # With different dispersion parameter (more overdispersion)
#' ss2MixedCountContinuous(
#'   r1 = 1.0,
#'   r2 = 1.25,
#'   nu = 0.5,
#'   t = 1,
#'   mu1 = -50,
#'   mu2 = 0,
#'   sd = 250,
#'   r = 1,
#'   rho1 = 0.4,
#'   rho2 = 0.4,
#'   alpha = 0.025,
#'   beta = 0.2
#' )
#'
#' @export
ss2MixedCountContinuous <- function(r1, r2, nu, t, mu1, mu2, sd, r, rho1, rho2, alpha, beta) {

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
  if (abs(rho1) >= 1) {
    stop("rho1 must be in (-1, 1)")
  }
  if (abs(rho2) >= 1) {
    stop("rho2 must be in (-1, 1)")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }
  if (beta <= 0 || beta >= 1) {
    stop("beta must be in (0, 1)")
  }

  # Check correlation bounds for group 1
  lambda1 <- r1 * t
  bounds1 <- corrbound2MixedCountContinuous(lambda1, nu, mu1, sd)
  if (rho1 < bounds1[1] || rho1 > bounds1[2]) {
    stop(paste0("rho1 must be within [", round(bounds1[1], 4), ", ",
                round(bounds1[2], 4), "]"))
  }

  # Check correlation bounds for group 2
  lambda2 <- r2 * t
  bounds2 <- corrbound2MixedCountContinuous(lambda2, nu, mu2, sd)
  if (rho2 < bounds2[1] || rho2 > bounds2[2]) {
    stop(paste0("rho2 must be within [", round(bounds2[1], 4), ", ",
                round(bounds2[2], 4), "]"))
  }

  # Step 1: Initialize sample size using single endpoint formulas
  n2_initial <- max(
    ss1Count(r1, r2, nu, t, r, alpha, beta)[["n2"]],
    ss1Continuous(abs(mu1 - mu2), sd, r, alpha, beta)[["n2"]]
  )

  # Step 2: Sequential search to find minimum sample size
  result_ss <- .ss_sequential_search(
    initial_n2 = n2_initial,
    r = r,
    target_power = 1 - beta,
    power_fun = power2MixedCountContinuous,
    r1 = r1, r2 = r2, nu = nu, t = t,
    mu1 = mu1, mu2 = mu2, sd = sd,
    rho1 = rho1, rho2 = rho2, alpha = alpha
  )

  n1 <- result_ss$n1
  n2 <- result_ss$n2
  N <- result_ss$N

  # Return results as a data frame
  result <- data.frame(
    r1, r2, nu, t, mu1, mu2, sd, r, rho1, rho2, alpha, beta,
    n1, n2, N
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
