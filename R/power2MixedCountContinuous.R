
#' Power Calculation for Two Co-Primary Endpoints (Count and Continuous)
#'
#' Calculates the power for a two-arm superiority trial with one overdispersed
#' count co-primary endpoint and one continuous co-primary endpoint, as described
#' in Homma and Yoshida (2024).
#'
#' @param n1 Sample size for group 1 (test group)
#' @param n2 Sample size for group 2 (control group)
#' @param r1 Mean rate (events per unit time) for the treatment group (count endpoint)
#' @param r2 Mean rate (events per unit time) for the control group (count endpoint)
#' @param nu Common dispersion parameter for the negative binomial distribution (nu > 0)
#' @param t Common follow-up time period
#' @param mu1 Mean for group 1 (continuous endpoint)
#' @param mu2 Mean for group 2 (continuous endpoint)
#' @param sd Common standard deviation for the continuous endpoint
#' @param rho1 Correlation between count and continuous outcomes for treatment group
#' @param rho2 Correlation between count and continuous outcomes for control group
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#'
#' @return A data frame with the following columns:
#'   \item{n1}{Sample size for group 1}
#'   \item{n2}{Sample size for group 2}
#'   \item{r1}{Mean rate in group 1 for count endpoint}
#'   \item{r2}{Mean rate in group 2 for count endpoint}
#'   \item{nu}{Dispersion parameter}
#'   \item{t}{Follow-up time}
#'   \item{mu1}{Mean in group 1 for continuous endpoint}
#'   \item{mu2}{Mean in group 2 for continuous endpoint}
#'   \item{sd}{Standard deviation for continuous endpoint}
#'   \item{rho1}{Correlation for group 1}
#'   \item{rho2}{Correlation for group 2}
#'   \item{alpha}{One-sided significance level}
#'   \item{powerCount}{Power for the count endpoint alone}
#'   \item{powerCont}{Power for the continuous endpoint alone}
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
#' The correlation bounds are automatically checked using \code{\link{corrbound2MixedCountContinuous}}.
#'
#' @references
#' Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical trials
#' with two co-primary endpoints including overdispersed count and continuous
#' outcomes. \emph{Pharmaceutical Statistics}, 23(1), 46-59.
#'
#' @examples
#' # Power calculation with moderate correlation
#' power2MixedCountContinuous(
#'   n1 = 300,
#'   n2 = 300,
#'   r1 = 1.0,
#'   r2 = 1.25,
#'   nu = 0.8,
#'   t = 1,
#'   mu1 = -50,
#'   mu2 = 0,
#'   sd = 250,
#'   rho1 = 0.5,
#'   rho2 = 0.5,
#'   alpha = 0.025
#' )
#'
#' # Power calculation with no correlation
#' power2MixedCountContinuous(
#'   n1 = 350,
#'   n2 = 350,
#'   r1 = 1.0,
#'   r2 = 1.5,
#'   nu = 1,
#'   t = 1,
#'   mu1 = -40,
#'   mu2 = 0,
#'   sd = 200,
#'   rho1 = 0,
#'   rho2 = 0,
#'   alpha = 0.025
#' )
#'
#' # Unbalanced design
#' power2MixedCountContinuous(
#'   n1 = 400,
#'   n2 = 200,
#'   r1 = 1,
#'   r2 = 1.25,
#'   nu = 1,
#'   t = 1,
#'   mu1 = -50,
#'   mu2 = 0,
#'   sd = 250,
#'   rho1 = 0.6,
#'   rho2 = 0.6,
#'   alpha = 0.025
#' )
#'
#' @export
#' @importFrom stats pnorm qnorm
#' @importFrom pbivnorm pbivnorm
power2MixedCountContinuous <- function(n1, n2, r1, r2, nu, t, mu1, mu2, sd,
                                       rho1, rho2, alpha) {

  # Input validation
  if (length(n1) != 1 || length(n2) != 1) {
    stop("n1 and n2 must be scalar values")
  }
  if (n1 <= 0 || n1 != round(n1)) {
    stop("n1 must be a positive integer")
  }
  if (n2 <= 0 || n2 != round(n2)) {
    stop("n2 must be a positive integer")
  }
  if (length(r1) != 1 || length(r2) != 1 || length(nu) != 1 ||
      length(t) != 1 || length(mu1) != 1 || length(mu2) != 1 ||
      length(sd) != 1 || length(rho1) != 1 || length(rho2) != 1 || length(alpha) != 1) {
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
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }

  # Calculate allocation ratio
  kappa <- n1 / n2

  # Calculate lambda (expected number of events)
  lambda1 <- r1 * t
  lambda2 <- r2 * t

  # Check that rho1 is within valid bounds
  bounds1 <- corrbound2MixedCountContinuous(lambda1, nu, mu1, sd)
  if (rho1 < bounds1[1] | rho1 > bounds1[2]) {
    stop(paste0("rho1 must be within [", round(bounds1[1], 4), ", ",
                round(bounds1[2], 4), "]"))
  }

  # Check that rho2 is within valid bounds
  bounds2 <- corrbound2MixedCountContinuous(lambda2, nu, mu2, sd)
  if (rho2 < bounds2[1] | rho2 > bounds2[2]) {
    stop(paste0("rho2 must be within [", round(bounds2[1], 4), ", ",
                round(bounds2[2], 4), "]"))
  }

  # Calculate variance components for count endpoint (equation 8)
  Va <- (1 / t) * (1 / r2 + 1 / (kappa * r1)) + (1 + kappa) / (nu * kappa)
  V0 <- Va  # Under H0

  # Calculate treatment effects
  delta <- mu1 - mu2
  beta1 <- log(r1 / r2)  # Log rate ratio

  # Standard normal quantiles
  z_alpha <- qnorm(alpha)

  # Calculate test statistics under alternative hypothesis
  Z1 <- sqrt(n2 / V0) * beta1
  Z2 <- delta / (sd * sqrt((1 + kappa) / (kappa * n2)))

  # Critical values
  c_val <- c(
    z_alpha - sqrt(V0) * Z1 / sqrt(Va),
    z_alpha - Z2
  )

  # Calculate correlation between test statistics (equation 11)
  gamma <- '+'(
    '/'(
      n2 * rho2 * sqrt(1 + lambda2 / nu),
      n2 * sqrt(lambda2 * Va) * sqrt((1 + kappa) / kappa)
    ),
    '/'(
      n2 * rho1 * sqrt(1 + lambda1 / nu),
      n1 * sqrt(lambda1 * Va) * sqrt((1 + kappa) / kappa)
    )
  )

  # Calculate power for individual endpoints
  power1and2 <- pnorm(c_val)

  # Calculate power for co-primary endpoints using bivariate normal distribution
  # (equation 13 in Homma and Yoshida 2024)
  powerCoprimary <- pbivnorm(x = c_val[1], y = c_val[2], rho = gamma)

  # Return results as a data frame
  result <- data.frame(
    n1, n2, r1, r2, nu, t, mu1, mu2, sd, rho1, rho2, alpha,
    powerCount = power1and2[1], powerCont = power1and2[2], powerCoprimary
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
