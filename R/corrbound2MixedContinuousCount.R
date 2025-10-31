#' Calculate Correlation Bounds for Mixed Count and Continuous Outcomes
#'
#' Computes the lower and upper bounds of the correlation coefficient between
#' a negative binomial count outcome and a continuous outcome, using the
#' Frechet-Hoeffding bounds as described in Homma and Yoshida (2024).
#'
#' @param lambda Expected number of events (lambda = r * t, where r is the rate
#'   and t is the follow-up time)
#' @param nu Dispersion parameter for the negative binomial distribution (nu > 0)
#' @param mu Mean of the continuous outcome
#' @param sd Standard deviation of the continuous outcome
#' @param t Follow-up time period
#'
#' @return A named numeric vector with two elements:
#'   \item{L_bound}{Lower bound of the correlation}
#'   \item{U_bound}{Upper bound of the correlation}
#'
#' @details
#' For a negative binomial count outcome Y1 ~ NB(lambda, nu) and a continuous
#' outcome Y2 ~ N(mu, sd^2), the correlation coefficient rho is bounded using
#' the Frechet-Hoeffding bounds on copulas:
#'
#' **Lower bound (Frechet-Hoeffding lower bound):**
#' \deqn{H^-(u, v) = \max(u + v - 1, 0)}
#'
#' **Upper bound (Frechet-Hoeffding upper bound):**
#' \deqn{H^+(u, v) = \min(u, v)}
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
#' corrbound2MixedContinuousCount(lambda = 1.25, nu = 0.8, mu = 0, sd = 250, t = 1)
#'
#' # Higher dispersion parameter
#' corrbound2MixedContinuousCount(lambda = 2.0, nu = 2.0, mu = 50, sd = 200, t = 1)
#'
#' # Different follow-up time
#' corrbound2MixedContinuousCount(lambda = 1.0 * 2, nu = 1.0, mu = 0, sd = 300, t = 2)
#'
#' @export
#' @importFrom stats pnbinom pnorm integrate qnbinom
corrbound2MixedContinuousCount <- function(lambda, nu, mu, sd, t) {

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
