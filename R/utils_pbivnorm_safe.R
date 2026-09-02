#' Bivariate Normal Distribution Function, Guarded at the Edges of Its Domain
#'
#' Internal wrapper around \code{pbivnorm} that returns a probability for every
#' argument the power functions of this package can produce.
#'
#' @param x Numeric vector of upper limits for the first coordinate.
#' @param y Numeric vector of upper limits for the second coordinate.
#' @param rho Correlation between the two coordinates, a single number.
#'
#' @return Numeric vector of joint probabilities, of the length of the longer of
#'   \code{x} and \code{y}.
#'
#' @details
#' Three situations arise in the power functions and none of them is handled by
#' \code{pbivnorm} itself. A large standardized effect drives an argument
#' hundreds of standard deviations into a tail, where \code{pbivnorm} returns
#' \code{NaN}; beyond eight standard deviations the univariate normal
#' distribution function is zero or one to within 1e-15, so the arguments are
#' held there. A correlation of plus or minus one, which the Prentice bound
#' produces when the two marginal probabilities are equal, is outside the domain
#' of \code{pbivnorm}; the joint probability is the corresponding
#' Frechet-Hoeffding bound there. And the returned value is confined to the
#' interval those bounds allow, so that the joint probability cannot exceed
#' either marginal by a rounding error nor fall below the Bonferroni bound.
#'
#' @keywords internal
#' @noRd
#' @importFrom pbivnorm pbivnorm
#' @importFrom stats pnorm
.pbivnorm_safe <- function(x, y, rho) {

  if (length(rho) != 1 || !is.finite(rho)) {
    stop("the correlation between the two test statistics is not a finite ",
         "single number")
  }

  n <- max(length(x), length(y))
  px <- rep_len(pnorm(x), n)
  py <- rep_len(pnorm(y), n)
  lo <- pmax(0, px + py - 1)
  hi <- pmin(px, py)

  if (rho >= 1) {
    return(hi)
  }
  if (rho <= -1) {
    return(lo)
  }

  p <- pbivnorm(x = pmin(pmax(x, -8), 8), y = pmin(pmax(y, -8), 8), rho = rho)
  p <- rep_len(p, n)

  # A value that is still not a number after the arguments have been held
  # inside the representable range is replaced by the midpoint of the only
  # interval it can lie in
  not_a_number <- !is.finite(p)
  if (any(not_a_number)) {
    p[not_a_number] <- (lo[not_a_number] + hi[not_a_number]) / 2
  }

  pmin(pmax(p, lo), hi)
}
