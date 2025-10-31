#' Power Calculation for Two Co-Primary Endpoints (Count and Continuous)
#'
#' Calculates the power for a two-arm superiority trial with one overdispersed
#' count co-primary endpoint and one continuous co-primary endpoint, as described
#' in Homma and Yoshida (2024).
#'
#' @param n0 Sample size for control group
#' @param n1 Sample size for treatment group
#' @param r0 Mean rate (events per unit time) for the control group (count endpoint)
#' @param r1 Mean rate (events per unit time) for the treatment group (count endpoint)
#' @param nu Common dispersion parameter for the negative binomial distribution (nu > 0)
#' @param t Common follow-up time period
#' @param mu0 Mean for the control group (continuous endpoint)
#' @param mu1 Mean for the treatment group (continuous endpoint)
#' @param sigma Common standard deviation for the continuous endpoint
#' @param rho0 Correlation between count and continuous outcomes for control group
#' @param rho1 Correlation between count and continuous outcomes for treatment group
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#'
#' @return A data frame with the following columns:
#'   \item{n0, n1}{Sample sizes}
#'   \item{r0, r1}{Mean rates for count endpoint}
#'   \item{nu}{Dispersion parameter}
#'   \item{t}{Follow-up time}
#'   \item{mu0, mu1}{Means for continuous endpoint}
#'   \item{sigma}{Standard deviation for continuous endpoint}
#'   \item{rho0, rho1}{Correlations}
#'   \item{alpha}{One-sided significance level}
#'   \item{power_count}{Power for the count endpoint alone}
#'   \item{power_cont}{Power for the continuous endpoint alone}
#'   \item{powerCoprimary}{Power for both co-primary endpoints}
#'
#' @details
#' The test statistics are (equation 7 in Homma and Yoshida 2024):
#' \deqn{Z_1 = \frac{\hat{\beta}_1}{\sqrt{Var(\hat{\beta}_1)}}, \quad
#'       Z_2 = \frac{\hat{\delta}}{\sigma\sqrt{(1+\kappa)/(\kappa n_0)}}}
#'
#' The joint distribution of (Z1, Z2) follows an asymptotic bivariate normal
#' distribution with correlation gamma (equation 11):
#' \deqn{\gamma = \sum_{j=0,1} \frac{n_0 \rho_j \sqrt{1+\lambda_j/\nu}}
#'       {n_j \sqrt{\lambda_j V_a} \sqrt{(1+\kappa)/\kappa}}}
#'
#' where \eqn{\lambda_j = r_j \times t}.
#'
#' @references
#' Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical trials
#' with two co-primary endpoints including overdispersed count and continuous
#' outcomes. \emph{Pharmaceutical Statistics}, 23(1), 46-59.
#'
#' @examples
#' # Power calculation with moderate correlation
#' power2MixedContinuousCount(
#'   n0 = 300, n1 = 300,
#'   r0 = 1.25, r1 = 1.0, nu = 0.8, t = 1,
#'   mu0 = 0, mu1 = -50, sigma = 250,
#'   rho0 = 0.5, rho1 = 0.5,
#'   alpha = 0.025
#' )
#'
#' # Power calculation with no correlation
#' power2MixedContinuousCount(
#'   n0 = 350, n1 = 350,
#'   r0 = 1.5, r1 = 1.0, nu = 1.0, t = 1,
#'   mu0 = 0, mu1 = -40, sigma = 200,
#'   rho0 = 0, rho1 = 0,
#'   alpha = 0.025
#' )
#'
#' # Unbalanced design
#' power2MixedContinuousCount(
#'   n0 = 200, n1 = 400,
#'   r0 = 1.25, r1 = 1.0, nu = 1.0, t = 1,
#'   mu0 = 0, mu1 = -50, sigma = 250,
#'   rho0 = 0.6, rho1 = 0.6,
#'   alpha = 0.025
#' )
#'
#' @export
#' @importFrom stats pnorm qnorm
#' @importFrom mvtnorm pmvnorm GenzBretz
power2MixedContinuousCount <- function(n0, n1, r0, r1, nu, t, mu0, mu1, sigma,
                                       rho0, rho1, alpha) {

  # Calculate allocation ratio
  kappa <- n1 / n0

  # Calculate lambda (expected number of events)
  lambda0 <- r0 * t
  lambda1 <- r1 * t

  # Calculate variance components for count endpoint (equation 8)
  Va <- (1 / t) * (1 / r0 + 1 / (kappa * r1)) + (1 + kappa) / (nu * kappa)
  V0 <- Va  # Under H0

  # Calculate treatment effects
  beta1 <- log(r1 / r0)  # Log rate ratio
  delta <- mu1 - mu0      # Mean difference

  # Calculate z-scores
  za <- qnorm(alpha)

  # Calculate test statistics under alternative hypothesis
  Z1 <- sqrt(n0 / V0) * beta1
  Z2 <- delta / (sigma * sqrt((1 + kappa) / (kappa * n0)))

  # Calculate correlation between test statistics (equation 11)
  gamma <- (n0 * rho0 * sqrt(1 + lambda0 / nu)) / (n0 * sqrt(lambda0 * Va) * sqrt((1 + kappa) / kappa)) +
    (n0 * rho1 * sqrt(1 + lambda1 / nu)) / (n1 * sqrt(lambda1 * Va) * sqrt((1 + kappa) / kappa))

  # Calculate power for individual endpoints
  power_count <- pnorm(za - Z1)
  power_cont <- pnorm(za - Z2)

  # Calculate power for co-primary endpoints using bivariate normal distribution
  # (equation 13 in Homma and Yoshida 2024)
  powerCoprimary <- pmvnorm(
    lower = c(-Inf, -Inf),
    upper = c(za * sqrt(V0 / Va) - sqrt(n0) * beta1 / sqrt(Va),
              za - sqrt(n0) * delta / (sigma * sqrt((1 + kappa) / kappa))),
    mean = c(0, 0),
    corr = matrix(c(1, gamma, gamma, 1), ncol = 2),
    algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
    seed = 1
  )[[1]]

  # Return results as a data frame
  result <- data.frame(
    n0, n1, r0, r1, nu, t, mu0, mu1, sigma, rho0, rho1, alpha,
    power_count, power_cont, powerCoprimary
  )
  return(result)
}
