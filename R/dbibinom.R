#' Probability Mass Function of Bivariate Binomial Distribution
#'
#' Calculates the probability mass function of the bivariate binomial distribution
#' for given parameters, as described in Homma and Yoshida (2025).
#'
#' @param N Sample size (number of trials)
#' @param y1 Observed value(s) of the first random variable (0 to N)
#' @param y2 Observed value(s) of the second random variable (0 to N)
#' @param p1 True probability of responders for the first outcome (0 < p1 < 1)
#' @param p2 True probability of responders for the second outcome (0 < p2 < 1)
#' @param rho Correlation coefficient between the two binary outcomes (must be NULL if gamma is specified)
#'
#' @return Probability mass function value(s) for the bivariate binomial distribution.
#'   If y1 and y2 are vectors, returns a vector of probabilities.
#'
#' @details
#' The bivariate binomial distribution BiBin(N, p1, p2, gamma) has probability mass
#' function given by equation (3) in Homma and Yoshida (2025):
#' \deqn{P(Y_1 = y_1, Y_2 = y_2) = f(y_1|N, p_1) \times g(y_2|y_1, N, p_1, p_2, \gamma)}
#' where
#' \deqn{g(y_2|y_1, N, p_1, p_2, \gamma) = \frac{1}{(1+\gamma)^N}
#'       \sum_{m \in \mathcal{M}} \binom{y_1}{m} \binom{N-y_1}{y_2-m}
#'       (\xi+\gamma)^m (1-\xi)^{y_1-m} \xi^{y_2-m} (1-\xi+\gamma)^{N-y_1-(y_2-m)}}
#' with \eqn{\xi = p_2 + \gamma(p_2 - p_1)} and
#' \eqn{\mathcal{M} = \{m : m = \max(0, y_2-(N-y_1)), \ldots, \min(y_1, y_2)\}}.
#'
#' @references
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
#'
#' @examples
#' # Calculate single probability mass
#' dbibinom(N = 100, y1 = 30, y2 = 50, p1 = 0.3, p2 = 0.5, rho = 0.5)
#'
#' # Verify that probabilities sum to 1
#' N <- 20
#' p1 <- 0.3
#' p2 <- 0.5
#' rho <- 0.5
#' sum(outer(0:N, 0:N, function(x, y) dbibinom(N, x, y, p1, p2, rho)))
#'
#' @export
dbibinom <- function(N, y1, y2, p1, p2, rho) {

  # Check that y1 and y2 have the same length
  if (length(y1) != length(y2)) {
    stop('y1 and y2 must have the same length')
  }

  # Check that rho is within valid bounds
  bounds <- corrbound2Binary(p1, p2)
  if (rho < bounds[1] | rho > bounds[2]) {
    stop(paste0("rho must be within [", round(bounds[1], 4), ", ",
                round(bounds[2], 4), "]"))
  }
  # Convert rho to gamma
  gamma <- '/'(
    rho * sqrt(p2 * (1 - p2) / (p1 * (1 - p1))),
    (1 - rho * sqrt(p2 * (1 - p2) / (p1 * (1 - p1))))
  )

  # Define the set M = {m: m = max(0, y2 - (N - y1)), ..., min(y1, y2)}
  m <- Map(':', pmax(0, y2 - (N - y1)), pmin(y1, y2))
  m <- t(sapply(m, `length<-`, max(lengths(m))))

  # Define xi (see Homma and Yoshida (2025))
  xi <- p2 + gamma * (p2 - p1)

  # Calculate g(y2|y1, N, p1, p2, gamma)
  g <- (1 + gamma) ^ (-N) * rowSums(
    '*'(
      choose(y1, m) * choose(N - y1, y2 - m) * (xi + gamma) ^ m,
      (1 - xi) ^ (y1 - m) * xi ^ (y2 - m) * (1 - xi + gamma) ^ (N - y1 - (y2 - m))
    ),
    na.rm = TRUE
  )

  # Calculate P(Y1 = y1, Y2 = y2|N, p1, p2, gamma) using equation (3)
  return(dbinom(y1, N, p1) * g)
}
