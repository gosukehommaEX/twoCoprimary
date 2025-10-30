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
#' # Example parameters for comparison across methods
#' n1_ex <- 100
#' n2_ex <- 100
#' delta1_ex <- 0.5
#' delta2_ex <- 0.5
#' sd1_ex <- 1
#' sd2_ex <- 1
#' rho_ex <- 0.3
#' alpha_ex <- 0.025
#'
#' # Power calculation with known variance
#' power2Continuous(
#'   n1 = n1_ex,
#'   n2 = n2_ex,
#'   delta1 = delta1_ex,
#'   delta2 = delta2_ex,
#'   sd1 = sd1_ex,
#'   sd2 = sd2_ex,
#'   rho = rho_ex,
#'   alpha = alpha_ex,
#'   known_var = TRUE
#' )
#'
#' \donttest{
#' # Power calculation with unknown variance (Monte Carlo)
#' power2Continuous(
#'   n1 = n1_ex,
#'   n2 = n2_ex,
#'   delta1 = delta1_ex,
#'   delta2 = delta2_ex,
#'   sd1 = sd1_ex,
#'   sd2 = sd2_ex,
#'   rho = rho_ex,
#'   alpha = alpha_ex,
#'   known_var = FALSE,
#'   nMC = 10000
#' )
#' }
#'
#' @export
#' @importFrom mvtnorm pmvnorm GenzBretz
#' @importFrom stats pnorm pt qt rWishart
power2Continuous <- function(n1, n2, delta1, delta2, sd1, sd2, rho, alpha,
                             known_var = TRUE, nMC = 1e+4) {

  # Calculate test statistics for each endpoint
  Z1 <- delta1 / (sd1 * sqrt(1 / n1 + 1 / n2))
  Z2 <- delta2 / (sd2 * sqrt(1 / n1 + 1 / n2))

  if (known_var) {

    # Set nMC to NA for known variance case
    nMC <- NA

    # Calculate power for individual endpoints using normal distribution
    power1 <- pnorm(-qnorm(1 - alpha) + Z1)
    power2 <- pnorm(-qnorm(1 - alpha) + Z2)

    # Calculate power for co-primary endpoints using bivariate normal distribution
    powerCoprimary <- pmvnorm(
      lower = c(-Inf, -Inf),
      upper = c(-qnorm(1 - alpha) + Z1, -qnorm(1 - alpha) + Z2),
      mean = c(0, 0),
      corr = matrix(c(1, rho, rho, 1), ncol = 2),
      algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
      seed = 1
    )[[1]]

  } else {

    # Calculate degrees of freedom for unknown variance case
    nu <- n1 + n2 - 2

    # Calculate power for individual endpoints using t-distribution
    power1 <- 1 - pt(qt(1 - alpha, nu), df = nu, ncp = Z1)
    power2 <- 1 - pt(qt(1 - alpha, nu), df = nu, ncp = Z2)

    # Define variance-covariance matrix
    Sigma <- matrix(c(sd1 ^ 2, rho * sd1 * sd2, rho * sd1 * sd2, sd2 ^ 2), nrow = 2)

    # Monte Carlo approach following Sozu et al. (2011) equation (6)
    # Generate Wishart random matrices
    Ws <- rWishart(nMC, df = nu, Sigma = Sigma)

    # Calculate power by averaging over simulated variance-covariance matrices
    probs <- numeric(nMC)
    for (i in 1:nMC) {
      # Extract diagonal elements (variances) from the i-th Wishart matrix
      Wi <- diag(Ws[, , i])

      # Calculate critical values adjusted by estimated variances
      ci <- qt(1 - alpha, df = nu) * sqrt(Wi / nu) - c(Z1, Z2)

      # Calculate probability for this iteration
      probs[i] <- pmvnorm(
        upper = -ci,
        mean = c(0, 0),
        corr = matrix(c(1, rho, rho, 1), ncol = 2),
        algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
        seed = 1
      )[[1]]
    }

    # Average over all Monte Carlo iterations
    powerCoprimary <- mean(probs)
  }

  # Return results as a data frame
  result <- data.frame(
    n1, n2, delta1, delta2, sd1, sd2, rho, alpha, known_var, nMC,
    power1, power2, powerCoprimary
  )
  return(result)
}
