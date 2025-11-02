#' Sample Size Calculation for a Single Continuous Endpoint
#'
#' Calculates the required sample size for a two-arm superiority trial with a
#' single continuous endpoint using the standard formula for normally distributed
#' outcomes.
#'
#' @param delta Mean difference between treatment groups (treatment effect)
#' @param sd Common standard deviation for the continuous endpoint
#' @param r Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#'
#' @return A data frame with the following columns:
#'   \item{delta}{Mean difference (treatment effect)}
#'   \item{sd}{Common standard deviation}
#'   \item{r}{Allocation ratio}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{N}{Total sample size (n1 + n2)}
#'
#' @details
#' The required sample size for group 2 is calculated using the standard formula:
#' \deqn{n_2 = \left\lceil \frac{(1 + 1/r) \sigma^2 (z_\alpha + z_\beta)^2}{\delta^2} \right\rceil}
#' where \eqn{z_\alpha} and \eqn{z_\beta} are the quantiles of the standard normal
#' distribution corresponding to the one-sided significance level \eqn{\alpha} and
#' type II error rate \eqn{\beta}, respectively. The sample size for group 1 is
#' \eqn{n_1 = \lceil r \times n_2 \rceil}.
#'
#' @examples
#' # Balanced design with 1:1 allocation
#' ss1Continuous(delta = 0.4, sd = 1, r = 1, alpha = 0.025, beta = 0.1)
#'
#' # Unbalanced design with 2:1 allocation
#' ss1Continuous(delta = 0.5, sd = 1.2, r = 2, alpha = 0.025, beta = 0.2)
#'
#' # Large treatment effect
#' ss1Continuous(delta = 0.8, sd = 1, r = 1, alpha = 0.025, beta = 0.1)
#'
#' @export
#' @importFrom stats qnorm
ss1Continuous <- function(delta, sd, r, alpha, beta) {

  # Input validation
  if (length(delta) != 1 || length(sd) != 1 || length(r) != 1 ||
      length(alpha) != 1 || length(beta) != 1) {
    stop("All parameters must be scalar values")
  }
  if (delta <= 0) {
    stop("delta must be positive")
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

  # Calculate z-scores for alpha and beta
  za <- qnorm(alpha)
  zb <- qnorm(beta)

  # Calculate required sample size for group 2
  n2 <- ceiling((1 + 1 / r) / delta ^ 2 * (za + zb) ^ 2 * sd ^ 2)

  # Calculate sample size for group 1
  n1 <- ceiling(r * n2)

  # Total sample size
  N <- n1 + n2

  # Return result as a data frame
  result <- data.frame(delta, sd, r, alpha, beta, n1, n2, N)
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
