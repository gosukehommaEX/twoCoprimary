#' Power Calculation for Two Co-Primary Continuous Endpoints
#'
#' Calculates the power for a two-arm superiority trial with two co-primary
#' continuous endpoints, as described in Sozu et al. (2011).
#'
#' @param n1 Sample size for group 1 (test group)
#' @param n2 Sample size for group 2 (control group)
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
#' correlation matrices of the standardized endpoints to account for variance
#' estimation uncertainty, following equation (6) in Sozu et al. (2011):
#' \deqn{\text{Power} = E_W\left[\Phi_2(-c_1^*\sqrt{w_{11}}, -c_2^*\sqrt{w_{22}} | \rho)\right]}
#' where \eqn{c_k^* = t_{\alpha,\nu}\sqrt{\frac{1}{\nu}} - \frac{Z_k}{\sqrt{w_{kk}}}} and
#' \eqn{W} follows a Wishart distribution with \eqn{\nu = n_1 + n_2 - 2} degrees of
#' freedom and the correlation matrix of the standardized endpoints as its scale
#' matrix, so that \eqn{\sqrt{w_{kk} / \nu}} is the ratio of the estimated to the
#' true standard deviation of endpoint k.
#'
#' The Monte Carlo step draws random numbers, so results for
#' \code{known_var = FALSE} vary between calls unless a seed is set with
#' \code{set.seed} beforehand.
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
#' @importFrom pbivnorm pbivnorm
#' @importFrom stats pnorm pt qt rWishart
#' @importFrom mvtnorm rmvnorm
power2Continuous <- function(n1, n2, delta1, delta2, sd1, sd2, rho, alpha,
                             known_var = TRUE, nMC = 1e+4) {

  # Input validation, matching the other power functions of the package. Without
  # it a negative standard deviation, a correlation outside (-1, 1) or a
  # non-integer sample size are carried through the formula and returned as a
  # number rather than refused.
  if (length(n1) != 1 || length(n2) != 1) {
    stop("n1 and n2 must be scalar values")
  }
  if (n1 <= 0 || n1 != round(n1)) {
    stop("n1 must be a positive integer")
  }
  if (n2 <= 0 || n2 != round(n2)) {
    stop("n2 must be a positive integer")
  }
  if (length(delta1) != 1 || length(delta2) != 1 || length(sd1) != 1 ||
      length(sd2) != 1 || length(rho) != 1 || length(alpha) != 1) {
    stop("All parameters must be scalar values")
  }
  if (sd1 <= 0 || sd2 <= 0) {
    stop("sd1 and sd2 must be positive")
  }
  if (abs(rho) >= 1) {
    stop("rho must be in (-1, 1)")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }
  if (!is.logical(known_var) || length(known_var) != 1 || is.na(known_var)) {
    stop("known_var must be logical (TRUE or FALSE)")
  }
  if (!known_var && (length(nMC) != 1 || !is.finite(nMC) || nMC < 1)) {
    stop("nMC must be a single positive number")
  }

  # Standard normal quantiles
  z_alpha <- qnorm(1 - alpha)

  # Test statistics
  Z <- c(delta1, delta2) / (c(sd1, sd2) * sqrt(1 / n1 + 1 / n2))

  if (known_var) {

    # Set nMC to NA for known variance case
    nMC <- NA

    # Critical values
    c_val <- -z_alpha + Z

    # Calculate power for individual endpoints using normal distribution
    power1and2 <- pnorm(c_val)

    # Calculate power for co-primary endpoints using bivariate normal distribution
    powerCoprimary <- .pbivnorm_safe(c_val[1], c_val[2], rho)

  } else {

    # Calculate degrees of freedom for unknown variance case
    nu <- n1 + n2 - 2

    # With fewer than three patients in total the variance cannot be estimated,
    # so the t-test does not exist and the power is zero. Returning zero rather
    # than failing keeps the sequential search well defined at its lower end.
    if (nu < 1) {
      result <- data.frame(
        n1, n2, delta1, delta2, sd1, sd2, rho, alpha, known_var, nMC,
        power1 = 0, power2 = 0, powerCoprimary = 0
      )
      class(result) <- c("twoCoprimary", "data.frame")
      return(result)
    }

    # Calculate power for individual endpoints using t-distribution
    power1and2 <- 1 - pt(qt(1 - alpha, nu), df = nu, ncp = Z)

    # Correlation matrix of the standardized endpoints. The Wishart matrix must
    # be that of the standardized variables, because Z below is already divided
    # by sd1 and sd2. Drawing it from the variance-covariance matrix instead
    # leaves a factor of sd_k in sqrt(W_kk / nu), which cancels only when
    # sd1 = sd2 = 1.
    Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)

    # Monte Carlo approach following Sozu et al. (2011) equation (6)
    # Generate Wishart random matrices. Only the two diagonal entries are used
    # below. rWishart refuses a single degree of freedom for a two by two scale
    # matrix, although the Wishart matrix is defined there as the outer product
    # of one normal draw, so that case is drawn directly.
    if (nu >= 2) {
      Ws <- rWishart(nMC, df = nu, Sigma = Sigma)
      W11 <- Ws[1, 1, ]
      W22 <- Ws[2, 2, ]
    } else {
      Zs <- rmvnorm(nMC, mean = c(0, 0), sigma = Sigma)
      W11 <- Zs[, 1] ^ 2
      W22 <- Zs[, 2] ^ 2
    }

    # Pre-compute constants
    t_alpha <- qt(1 - alpha, df = nu)
    sqrt_nu_inv <- sqrt(1 / nu)

    # VECTORIZED: Calculate critical values for all iterations at once
    c_val1 <- -t_alpha * sqrt(W11) * sqrt_nu_inv + Z[1]
    c_val2 <- -t_alpha * sqrt(W22) * sqrt_nu_inv + Z[2]

    # VECTORIZED: Use pbivnorm for all iterations at once
    probs <- .pbivnorm_safe(c_val1, c_val2, rho)

    # Average over all Monte Carlo iterations. The two marginal powers above are
    # exact while this one is simulated, so at a small number of replications
    # the simulated value can fall outside the interval the exact marginals
    # allow. It is returned at the nearer end of that interval when it does.
    powerCoprimary <- .clamp_coprimary(mean(probs), power1and2[1],
                                       power1and2[2])
  }

  # Return results as a data frame
  result <- data.frame(
    n1, n2, delta1, delta2, sd1, sd2, rho, alpha, known_var, nMC,
    power1 = power1and2[1], power2 = power1and2[2], powerCoprimary
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
