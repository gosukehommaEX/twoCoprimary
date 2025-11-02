#' Calculate Correlation Bounds Between Count and Continuous Outcomes
#'
#' Computes the lower and upper bounds of the correlation coefficient between
#' an overdispersed count outcome (negative binomial) and a continuous outcome
#' (normal), as described in Homma and Yoshida (2024).
#'
#' @param lambda Mean parameter for the negative binomial distribution (lambda > 0)
#' @param nu Dispersion parameter for the negative binomial distribution (nu > 0)
#' @param mu Mean for the continuous outcome
#' @param sd Standard deviation for the continuous outcome (sd > 0)
#'
#' @return A named numeric vector with two elements:
#'   \item{L_bound}{Lower bound of the correlation}
#'   \item{U_bound}{Upper bound of the correlation}
#'
#' @details
#' The correlation bounds are calculated using the Frechet-Hoeffding bounds for
#' copulas, as described in Trivedi and Zimmer (2007). The negative binomial
#' distribution has mean lambda and variance:
#' \deqn{Var(Y_1) = \lambda + \frac{\lambda^2}{\nu}}
#'
#' The variance of the negative binomial distribution is: Var(Y1) = lambda + lambda^2/nu
#'
#' @references
#' Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical trials
#' with two co-primary endpoints including overdispersed count and continuous
#' outcomes. \emph{Pharmaceutical Statistics}, 23(1), 46-59.
#'
#' Trivedi, P. K., & Zimmer, D. M. (2007). Copula modeling: an introduction for
#' practitioners. \emph{Foundations and Trends in Econometrics}, 1(1), 1-111.
#'
#' @examples
#' # Calculate correlation bounds for NB(1.25, 0.8) and N(0, 250)
#' corrbound2MixedCountContinuous(lambda = 1.25, nu = 0.8, mu = 0, sd = 250)
#'
#' # Higher dispersion parameter
#' corrbound2MixedCountContinuous(lambda = 2.0, nu = 2.0, mu = 50, sd = 200)
#'
#' # Different follow-up time
#' corrbound2MixedCountContinuous(lambda = 1.0 * 2, nu = 1.0, mu = 0, sd = 300)
#'
#' @export
#' @importFrom stats pnbinom pnorm integrate
corrbound2MixedCountContinuous <- function(lambda, nu, mu, sd) {

  # Input validation
  if (length(lambda) != 1 || length(nu) != 1 || length(mu) != 1 || length(sd) != 1) {
    stop("All parameters must be scalar values")
  }
  if (lambda <= 0) {
    stop("lambda must be positive")
  }
  if (nu <= 0) {
    stop("nu must be positive")
  }
  if (sd <= 0) {
    stop("sd must be positive")
  }

  # Calculate variance of count outcome: Var(Y1) = lambda + lambda^2/nu
  var_count <- lambda + lambda ^ 2 / nu
  sd_count <- sqrt(var_count)

  # Calculate CDF of negative binomial distribution
  # P(Y <= y) for y = 0, 1, 2, ..., up to a reasonable upper limit
  y_max <- qnbinom(0.9999, mu = lambda, size = nu)
  Fx <- pnbinom(0:y_max, mu = lambda, size = nu)

  # Remove probabilities equal to 1 (beyond effective support)
  Fx <- Fx[Fx < 1]

  # Calculate lower bound using Frechet-Hoeffding lower bound
  # H^-(u,v) = max(u + v - 1, 0)
  LL_cov <- sum(sapply(seq_along(Fx), function(i) {
    integrate(function(y) {
      pmax(Fx[i] + pnorm(y, mu, sd) - 1, 0) - Fx[i] * pnorm(y, mu, sd)
    }, lower = -Inf, upper = Inf,
    rel.tol = .Machine$double.eps^0.5,
    stop.on.error = FALSE)$value
  }))

  LL_rho <- LL_cov / (sd_count * sd)

  # Calculate upper bound using Frechet-Hoeffding upper bound
  # H^+(u,v) = min(u, v)
  UL_cov <- sum(sapply(seq_along(Fx), function(i) {
    integrate(function(y) {
      pmin(Fx[i], pnorm(y, mu, sd)) - Fx[i] * pnorm(y, mu, sd)
    }, lower = -Inf, upper = Inf,
    rel.tol = .Machine$double.eps^0.5,
    stop.on.error = FALSE)$value
  }))

  UL_rho <- UL_cov / (sd_count * sd)

  # Return bounds
  boundary <- c(LL_rho, UL_rho)
  names(boundary) <- c('L_bound', 'U_bound')

  return(boundary)
}
