#' Power Calculation for Two Co-Primary Continuous Endpoints
#'
#' Calculates the power for a two-arm superiority trial with two co-primary
#' continuous endpoints, as described in Sozu et al. (2011).
#'
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param delta1 Mean difference for the first endpoint
#' @param delta2 Mean difference for the second endpoint
#' @param sd1 Common standard deviation for the first endpoint
#' @param sd2 Common standard deviation for the second endpoint
#' @param rho Common correlation between the two outcomes
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param known_var Logical value indicating whether variance is known (TRUE) or
#'   unknown (FALSE). If TRUE, power is calculated analytically; otherwise,
#'   Monte Carlo simulation is used for unknown variance
#' @param nMC Number of Monte Carlo simulations when known_var = FALSE (default is 10000)
#'
#' @return A data frame with the following columns:
#'   \item{n1}{Sample size for group 1}
#'   \item{n2}{Sample size for group 2}
#'   \item{delta1}{Mean difference for endpoint 1}
#'   \item{delta2}{Mean difference for endpoint 2}
#'   \item{sd1}{Standard deviation for endpoint 1}
#'   \item{sd2}{Standard deviation for endpoint 2}
#'   \item{rho}{Correlation between endpoints}
#'   \item{alpha}{One-sided significance level}
#'   \item{known_var}{Variance assumption}
#'   \item{nMC}{Number of Monte Carlo simulations (NA if known_var = TRUE)}
#'   \item{power1}{Power for the first endpoint alone}
#'   \item{power2}{Power for the second endpoint alone}
#'   \item{powerCoprimary}{Power for both co-primary endpoints}
#'
#' @details
#' For known variance, the power is calculated using the bivariate normal distribution
#' as described in Sozu et al. (2011). The test statistics are:
#' \deqn{Z_k = \frac{\delta_k}{\sigma_k \sqrt{1/n_1 + 1/n_2}}}
#' for k = 1, 2. The co-primary power is:
#' \deqn{1 - \beta = \Phi_2\left(-z_{1-\alpha} + Z_1, -z_{1-\alpha} + Z_2 \mid \rho\right)}
#' where \eqn{\Phi_2} is the cumulative distribution function of the bivariate
#' standard normal distribution.
#'
#' For unknown variance, Monte Carlo simulation is used with Wishart-distributed
#' variance-covariance matrices to account for variance estimation uncertainty,
#' following equation (6) in Sozu et al. (2011):
#' \deqn{\text{Power} = E_W\left[\Phi_2(-c_1^*\sqrt{w_{11}}, -c_2^*\sqrt{w_{22}} | \rho)\right]}
#' where \eqn{c_k^* = t_{\alpha,\nu}\sqrt{\frac{1}{\nu}} - \frac{Z_k}{\sqrt{w_{kk}}}} and
#' \eqn{W} follows a Wishart distribution with \eqn{\nu = n_1 + n_2 - 2} degrees of freedom.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2011). Sample size determination in
#' superiority clinical trials with multiple co-primary correlated endpoints.
#' \emph{Journal of Biopharmaceutical Statistics}, 21(4), 650-668.
#'
#' @examples
#' # Example 1: Known variance with rho = 0.5 (from Sozu et al. 2011)
#' power2Continuous(
#'   n1 = 417,
#'   n2 = 417,
#'   delta1 = 0.2,
#'   delta2 = 0.25,
#'   sd1 = 1,
#'   sd2 = 1,
#'   rho = 0.5,
#'   alpha = 0.025,
#'   known_var = TRUE
#' )
#'
#' # Example 2: Known variance with zero correlation
#' power2Continuous(
#'   n1 = 350,
#'   n2 = 350,
#'   delta1 = 0.2,
#'   delta2 = 0.2,
#'   sd1 = 1,
#'   sd2 = 1,
#'   rho = 0,
#'   alpha = 0.025,
#'   known_var = TRUE
#' )
#'
#' # Example 3: Unequal allocation
#' power2Continuous(
#'   n1 = 200,
#'   n2 = 100,
#'   delta1 = 0.3,
#'   delta2 = 0.25,
#'   sd1 = 1,
#'   sd2 = 1,
#'   rho = 0.3,
#'   alpha = 0.025,
#'   known_var = TRUE
#' )
#'
#' \donttest{
#' # Example 4: Unknown variance (Monte Carlo simulation)
#' set.seed(12345)
#' power2Continuous(
#'   n1 = 150,
#'   n2 = 150,
#'   delta1 = 0.5,
#'   delta2 = 0.4,
#'   sd1 = 1,
#'   sd2 = 1,
#'   rho = 0.4,
#'   alpha = 0.025,
#'   known_var = FALSE,
#'   nMC = 10000
#' )
#' }
#'
#' @export
#' @importFrom stats qnorm pnorm qt pt rWishart
#' @importFrom mvtnorm pmvnorm GenzBretz
power2Continuous <- function(n1, n2, delta1, delta2, sd1, sd2, rho, alpha,
                             known_var = TRUE, nMC = 1e+4) {

  # Calculate the kappa (allocation ratio)
  kappa <- n1 / n2

  if (known_var) {
    # ===== KNOWN VARIANCE =====
    # Calculate Z-statistics under alternative hypothesis
    Z1 <- delta1 / (sd1 * sqrt(1 / n1 + 1 / n2))
    Z2 <- delta2 / (sd2 * sqrt(1 / n1 + 1 / n2))

    # Critical value for one-sided test
    z_alpha <- qnorm(1 - alpha)

    # Power for individual endpoints
    power1 <- pnorm(-z_alpha + Z1)
    power2 <- pnorm(-z_alpha + Z2)

    # Power for co-primary endpoints using bivariate normal distribution
    powerCoprimary <- pmvnorm(
      lower = c(-Inf, -Inf),
      upper = c(-z_alpha + Z1, -z_alpha + Z2),
      mean = c(0, 0),
      corr = matrix(c(1, rho, rho, 1), ncol = 2),
      algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
      seed = 1
    )[[1]]

  } else {
    # ===== UNKNOWN VARIANCE: Use Monte Carlo simulation =====
    # Degrees of freedom
    nu <- n1 + n2 - 2

    # Critical value for t-distribution
    t_alpha <- qt(1 - alpha, nu)

    # Calculate Z-statistics (same as known variance case)
    Z1 <- delta1 / (sd1 * sqrt(1 / n1 + 1 / n2))
    Z2 <- delta2 / (sd2 * sqrt(1 / n1 + 1 / n2))

    # Generate Wishart matrices
    # Identity matrix scaled by nu
    Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
    W <- rWishart(nMC, df = nu, Sigma = Sigma * nu)

    # Calculate power using Monte Carlo simulation (equation 6 in Sozu et al. 2011)
    # For each Monte Carlo sample, calculate the adjusted critical values
    c1_star <- sapply(1:nMC, function(i) {
      w11 <- W[1, 1, i]
      t_alpha * sqrt(1 / nu) - Z1 / sqrt(w11)
    })

    c2_star <- sapply(1:nMC, function(i) {
      w22 <- W[2, 2, i]
      t_alpha * sqrt(1 / nu) - Z2 / sqrt(w22)
    })

    # Power for individual endpoints
    power1 <- mean(pt(c1_star, nu))
    power2 <- mean(pt(c2_star, nu))

    # Power for co-primary endpoints
    # Use pmvnorm for each Monte Carlo sample and average
    powerCoprimary <- mean(sapply(1:nMC, function(i) {
      pmvnorm(
        lower = c(-Inf, -Inf),
        upper = c(c1_star[i], c2_star[i]),
        mean = c(0, 0),
        corr = matrix(c(1, rho, rho, 1), ncol = 2),
        algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
      )[[1]]
    }))
  }

  # Return results as a data frame
  result <- data.frame(
    n1, n2, delta1, delta2, sd1, sd2, rho, alpha, known_var, nMC,
    power1, power2, powerCoprimary
  )
  return(result)
}
