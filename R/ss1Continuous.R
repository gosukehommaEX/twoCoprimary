#' Sample Size Calculation for a Single Continuous Endpoint
#'
#' Calculates the required sample size for a two-arm superiority trial with a
#' single continuous endpoint.
#'
#' @param delta Mean difference (group 1 - group 2)
#' @param sd Common standard deviation
#' @param r Allocation ratio (group 1:group 2 = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#'
#' @return A data frame with the following columns:
#'   \item{delta}{Mean difference}
#'   \item{sd}{Standard deviation}
#'   \item{r}{Allocation ratio}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{n}{Total sample size (n1 + n2)}
#'
#' @details
#' The required sample size is calculated using the standard formula for
#' continuous endpoints:
#' \deqn{n_2 = \left\lceil \left(\frac{(z_{1-\alpha} + z_{1-\beta}) \sigma}{\delta}\right)^2
#'       \left(1 + \frac{1}{r}\right) \right\rceil}
#' where \eqn{z_{1-\alpha}} and \eqn{z_{1-\beta}} are standard normal quantiles.
#'
#' @examples
#' # Balanced design with 1:1 allocation
#' ss1Continuous(delta = 0.2, sd = 1, r = 1, alpha = 0.025, beta = 0.1)
#'
#' # Unbalanced design with 2:1 allocation
#' ss1Continuous(delta = 0.3, sd = 1, r = 2, alpha = 0.025, beta = 0.2)
#'
#' # Larger effect size
#' ss1Continuous(delta = 0.5, sd = 1, r = 1, alpha = 0.025, beta = 0.1)
#'
#' @export
#' @importFrom stats qnorm
ss1Continuous <- function(delta, sd, r, alpha, beta) {

  # Calculate z-scores for alpha and beta
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(1 - beta)

  # Calculate sample size for group 2
  # n2 = [(z_alpha + z_beta) * sd / delta]^2 * (1 + 1/r)
  n2 <- ceiling(((z_alpha + z_beta) * sd / delta) ^ 2 * (1 + 1 / r))

  # Calculate sample size for group 1
  n1 <- ceiling(r * n2)

  # Total sample size
  n <- n1 + n2

  # Return result as a data frame
  result <- data.frame(delta, sd, r, alpha, beta, n1, n2, n)
  return(result)
}
