#' Conditional Probability of the Bivariate Binomial Distribution, in R
#'
#' Internal reference implementation of \code{g(y2 | y1, N, p1, p2, gamma)},
#' equation (3) of Homma and Yoshida (2025). The exported \code{dbibinom} calls
#' the compiled version instead, and this function is retained so that the
#' package tests can check the two against each other.
#'
#' @param N Number of subjects in the group
#' @param y1 Responder counts of the first endpoint
#' @param y2 Responder counts of the second endpoint, of the same length as
#'   \code{y1}
#' @param xi Value of \code{p2 + gamma (p2 - p1)}
#' @param gamma Dependence parameter obtained from the correlation
#'
#' @return A numeric vector of the same length as \code{y1}
#'
#' @keywords internal
#' @noRd
.dbibinom_g_r <- function(N, y1, y2, xi, gamma) {

  # Define the set M = {m: m = max(0, y2 - (N - y1)), ..., min(y1, y2)}
  m <- Map(':', pmax(0, y2 - (N - y1)), pmin(y1, y2))
  m <- matrix(unlist(lapply(m, `length<-`, max(lengths(m)))),
              nrow = length(y1), byrow = TRUE)

  # Calculate g(y2|y1, N, p1, p2, gamma)
  (1 + gamma) ^ (-N) * rowSums(
    '*'(
      choose(y1, m) * choose(N - y1, y2 - m) * (xi + gamma) ^ m,
      (1 - xi) ^ (y1 - m) * xi ^ (y2 - m) * (1 - xi + gamma) ^ (N - y1 - (y2 - m))
    ),
    na.rm = TRUE
  )
}
