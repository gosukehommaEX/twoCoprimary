#' Calculate Correlation Bounds for Binary Outcomes
#'
#' Computes the lower and upper bounds of the correlation coefficient between
#' two binary outcomes based on their marginal probabilities, as described in
#' Prentice (1988).
#'
#' @param p1 True probability of responders for the first outcome (0 < p1 < 1)
#' @param p2 True probability of responders for the second outcome (0 < p2 < 1)
#'
#' @return A named numeric vector with two elements:
#'   \item{L_bound}{Lower bound of the correlation}
#'   \item{U_bound}{Upper bound of the correlation}
#'
#' @details
#' For two binary outcomes with marginal probabilities p1 and p2, the correlation
#' coefficient rho is bounded by:
#' \deqn{\rho \in [L(p_1, p_2), U(p_1, p_2)]}
#' where
#' \deqn{L(p_1, p_2) = \max\left\{-\sqrt{\frac{p_1 p_2}{(1-p_1)(1-p_2)}},
#'       -\sqrt{\frac{(1-p_1)(1-p_2)}{p_1 p_2}}\right\}}
#' \deqn{U(p_1, p_2) = \min\left\{\sqrt{\frac{p_1(1-p_2)}{p_2(1-p_1)}},
#'       \sqrt{\frac{p_2(1-p_1)}{p_1(1-p_2)}}\right\}}
#'
#' @references
#' Prentice, R. L. (1988). Correlated binary regression with covariates specific
#' to each binary observation. \emph{Biometrics}, 44(4), 1033-1048.
#'
#' @examples
#' # Calculate correlation bounds for two binary outcomes
#' boundRho(p1 = 0.3, p2 = 0.5)
#'
#' # When probabilities are equal, upper bound is 1
#' boundRho(p1 = 0.4, p2 = 0.4)
#'
#' # When p1 + p2 = 1, lower bound is -1
#' boundRho(p1 = 0.3, p2 = 0.7)
#'
#' @export
boundRho <- function(p1, p2) {
  # Calculate boundary of rho (see Prentice (1988))
  boundary <- c(
    max(
      -sqrt((p1 * p2) / ((1 - p1) * (1 - p2))),
      -sqrt(((1 - p1) * (1 - p2)) / (p1 * p2))
    ),
    min(
      sqrt((p1 * (1 - p2)) / (p2 * (1 - p1))),
      sqrt((p2 * (1 - p1)) / (p1 * (1 - p2)))
    )
  )
  names(boundary) <- c('L_bound', 'U_bound')

  return(boundary)
}
