#' Convert Correlation to Gamma Parameter for Bivariate Binomial Distribution
#'
#' Converts a correlation coefficient to the gamma dependence parameter for the
#' bivariate binomial distribution, as described in Homma and Yoshida (2025).
#'
#' @param p1 True probability of responders for the first outcome (0 < p1 < 1)
#' @param p2 True probability of responders for the second outcome (0 < p2 < 1)
#' @param rho Correlation coefficient between the two binary outcomes
#'
#' @return The gamma dependence parameter for the bivariate binomial distribution
#'
#' @details
#' The relationship between the correlation coefficient rho and the gamma
#' parameter is given by equation (4) in Homma and Yoshida (2025):
#' \deqn{\gamma(\rho, p_1, p_2) = \frac{\rho\sqrt{p_2(1-p_2)/(p_1(1-p_1))}}
#'       {1 - \rho\sqrt{p_2(1-p_2)/(p_1(1-p_1))}}}
#'
#' @references
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
#'
#' @examples
#' # Convert correlation to gamma parameter
#' bibinomGamma(p1 = 0.6, p2 = 0.3, rho = 0.5)
#'
#' # With equal probabilities
#' bibinomGamma(p1 = 0.5, p2 = 0.5, rho = 0.3)
#'
#' # Negative correlation
#' bibinomGamma(p1 = 0.4, p2 = 0.6, rho = -0.2)
#'
#' @export
bibinomGamma <- function(p1, p2, rho) {
  # Calculate gamma from rho using equation (4) in Homma and Yoshida (2025)
  gamma <- '/'(
    rho * sqrt(p2 * (1 - p2) / (p1 * (1 - p1))),
    (1 - rho * sqrt(p2 * (1 - p2) / (p1 * (1 - p1))))
  )

  return(gamma)
}
