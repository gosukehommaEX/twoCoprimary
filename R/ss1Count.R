#' Sample Size Calculation for a Single Count Endpoint (Negative Binomial)
#'
#' Calculates the required sample size for a two-arm superiority trial with a
#' single overdispersed count endpoint following a negative binomial distribution,
#' as described in Homma and Yoshida (2024).
#'
#' @param r1 Mean rate (events per unit time) for the treatment group
#' @param r2 Mean rate (events per unit time) for the control group
#' @param nu Common dispersion parameter for the negative binomial distribution (nu > 0)
#' @param t Common follow-up time period
#' @param r Allocation ratio (treatment:control = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#'
#' @return A data frame with the following columns:
#'   \item{r1}{Mean rate for treatment group}
#'   \item{r2}{Mean rate for control group}
#'   \item{nu}{Dispersion parameter}
#'   \item{t}{Follow-up time}
#'   \item{r}{Allocation ratio}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{n1}{Required sample size for treatment group}
#'   \item{n2}{Required sample size for control group}
#'   \item{N}{Total sample size (n2 + n1)}
#'
#' @details
#' The test statistic for the negative binomial rate ratio is:
#' \deqn{Z_1 = \frac{\hat{\beta}_1}{\sqrt{Var(\hat{\beta}_1)}}}
#' where \eqn{\hat{\beta}_1 = \log(\bar{Y}_1) - \log(\bar{Y}_2)} and the variance is:
#' \deqn{Var(\hat{\beta}_1) = \frac{1}{n_2}\left[\frac{1}{t}\left(\frac{1}{r_2} +
#'       \frac{1}{r \cdot r_1}\right) + \frac{1+r}{\nu \cdot r}\right]}
#'
#' This is equation (8) in Homma and Yoshida (2024).
#'
#' @references
#' Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical trials
#' with two co-primary endpoints including overdispersed count and continuous
#' outcomes. \emph{Pharmaceutical Statistics}, 23(1), 46-59.
#'
#' @examples
#' # Sample size for count endpoint with nu = 0.8
#' ss1Count(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1, r = 1,
#'          alpha = 0.025, beta = 0.1)
#'
#' # Unbalanced design with 2:1 allocation
#' ss1Count(r1 = 1.0, r2 = 1.5, nu = 1.0, t = 1, r = 2,
#'          alpha = 0.025, beta = 0.2)
#'
#' # Higher dispersion
#' ss1Count(r1 = 1.5, r2 = 2.0, nu = 3.0, t = 1, r = 1,
#'          alpha = 0.025, beta = 0.1)
#'
#' @export
#' @importFrom stats qnorm
ss1Count <- function(r1, r2, nu, t, r, alpha, beta) {

  # Input validation
  if (length(r1) != 1 || length(r2) != 1 || length(nu) != 1 ||
      length(t) != 1 || length(r) != 1 || length(alpha) != 1 || length(beta) != 1) {
    stop("All parameters must be scalar values")
  }
  if (r1 <= 0) {
    stop("r1 must be positive")
  }
  if (r2 <= 0) {
    stop("r2 must be positive")
  }
  if (nu <= 0) {
    stop("nu must be positive")
  }
  if (t <= 0) {
    stop("t must be positive")
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

  # Calculate variance components (equation 8 in Homma and Yoshida 2024)
  Va <- (1 / t) * (1 / r2 + 1 / (r * r1)) + (1 + r) / (nu * r)
  V0 <- Va  # Under H0, variance is the same

  # Calculate log rate ratio (treatment effect)
  beta1 <- log(r1 / r2)

  # Calculate required sample size for control group
  # From equation (14) with single endpoint
  n2 <- ceiling(((za * sqrt(V0) + zb * sqrt(Va)) ^ 2) / (beta1 ^ 2))

  # Calculate sample size for treatment group
  n1 <- ceiling(r * n2)

  # Total sample size
  N <- n1 + n2

  # Return result as a data frame
  result <- data.frame(r1, r2, nu, t, r, alpha, beta, n1, n2, N)
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
