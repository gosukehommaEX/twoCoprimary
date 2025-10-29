#' Calculate Gamma Parameter Bounds for Bivariate Binomial Distribution
#'
#' Computes the lower and upper bounds of the dependence parameter gamma for
#' the bivariate binomial distribution based on marginal probabilities, as
#' described in Homma and Yoshida (2025).
#'
#' @param p1 True probability of responders for the first outcome (0 < p1 < 1)
#' @param p2 True probability of responders for the second outcome (0 < p2 < 1)
#'
#' @return A named numeric vector with two elements:
#'   \item{L_bound}{Lower bound of gamma}
#'   \item{U_bound}{Upper bound of gamma}
#'
#' @details
#' The dependence parameter gamma in the bivariate binomial distribution is
#' bounded as follows:
#' \itemize{
#'   \item If \eqn{p_1 > p_2}: \eqn{L = -p_2/(1-p_1+p_2)}, \eqn{U = p_2/(p_1-p_2)}
#'   \item If \eqn{p_1 < p_2}: \eqn{L = -(1-p_2)/(1+p_1-p_2)}, \eqn{U = (1-p_2)/(p_2-p_1)}
#'   \item If \eqn{p_1 = p_2}: \eqn{L = -1/2}, \eqn{U = \infty}
#' }
#'
#' @references
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
#'
#' @examples
#' # Calculate gamma bounds for two binary outcomes
#' boundGamma(p1 = 0.3, p2 = 0.5)
#'
#' # When p1 > p2
#' boundGamma(p1 = 0.6, p2 = 0.4)
#'
#' # When p1 = p2, upper bound is infinite
#' boundGamma(p1 = 0.5, p2 = 0.5)
#'
#' @export
boundGamma <- function(p1, p2) {
  # Calculate boundary of gamma (see Homma and Yoshida (2025))

  # Lower bound
  if (p1 > p2) {
    L_bound <- -p2 / (1 - p1 + p2)
  } else if (p1 < p2) {
    L_bound <- -(1 - p2) / (1 + p1 - p2)
  } else {
    L_bound <- -1 / 2
  }

  # Upper bound
  if (p1 > p2) {
    U_bound <- p2 / (p1 - p2)
  } else if (p1 < p2) {
    U_bound <- (1 - p2) / (p2 - p1)
  } else {
    U_bound <- Inf
  }

  boundary <- c(L_bound, U_bound)
  names(boundary) <- c('L_bound', 'U_bound')

  return(boundary)
}
