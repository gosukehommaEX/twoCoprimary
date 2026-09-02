#' Power for Count and Continuous Co-Primary Endpoints Without Argument Checking
#'
#' Computational core of \code{power2MixedCountContinuous}. The exported
#' function validates its arguments, including the two correlation bounds, and
#' then calls this function.
#'
#' The correlation bounds are obtained from
#' \code{corrbound2MixedCountContinuous}, which integrates over the support of
#' the negative binomial distribution and costs about fifty milliseconds per
#' call. Their value depends on \code{lambda}, \code{nu}, \code{mu} and
#' \code{sd} alone and never on \code{n1} or \code{n2}, so the sequential
#' search of \code{ss2MixedCountContinuous} validates them once at the start
#' and then calls this function at every candidate sample size. Calling the
#' exported function instead recomputed the same two numbers at every
#' iteration, which accounted for essentially the whole running time of the
#' search.
#'
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param r1 Mean rate for the count endpoint in group 1
#' @param r2 Mean rate for the count endpoint in group 2
#' @param nu Dispersion parameter of the negative binomial distribution
#' @param t Follow-up time period
#' @param mu1 Mean for the continuous endpoint in group 1
#' @param mu2 Mean for the continuous endpoint in group 2
#' @param sd Common standard deviation for the continuous endpoint
#' @param rho1 Correlation between the two endpoints in group 1
#' @param rho2 Correlation between the two endpoints in group 2
#' @param alpha One-sided significance level
#'
#' @return A data frame of class \code{twoCoprimary} holding the arguments
#'   together with \code{powerCount}, \code{powerCont} and
#'   \code{powerCoprimary}
#'
#' @keywords internal
#' @noRd
.power2MixedCountContinuous_core <- function(n1, n2, r1, r2, nu, t, mu1, mu2, sd,
                                             rho1, rho2, alpha) {

  # Calculate allocation ratio
  kappa <- n1 / n2

  # Calculate lambda (expected number of events)
  lambda1 <- r1 * t
  lambda2 <- r2 * t

  # Calculate variance components for count endpoint (equation 8)
  Va <- (1 / t) * (1 / r2 + 1 / (kappa * r1)) + (1 + kappa) / (nu * kappa)
  V0 <- Va  # Under H0

  # Calculate treatment effects
  delta <- mu1 - mu2
  beta1 <- log(r1 / r2)  # Log rate ratio

  # Standard normal quantiles
  z_alpha <- qnorm(alpha)

  # Calculate test statistics under alternative hypothesis
  Z1 <- sqrt(n2 / V0) * beta1
  Z2 <- delta / (sd * sqrt((1 + kappa) / (kappa * n2)))

  # Critical values
  c_val <- c(
    z_alpha - sqrt(V0) * Z1 / sqrt(Va),
    z_alpha - Z2
  )

  # Calculate correlation between test statistics (equation 11)
  gamma <- '+'(
    '/'(
      n2 * rho2 * sqrt(1 + lambda2 / nu),
      n2 * sqrt(lambda2 * Va) * sqrt((1 + kappa) / kappa)
    ),
    '/'(
      n2 * rho1 * sqrt(1 + lambda1 / nu),
      n1 * sqrt(lambda1 * Va) * sqrt((1 + kappa) / kappa)
    )
  )

  # Calculate power for individual endpoints
  power1and2 <- pnorm(c_val)

  # Calculate power for co-primary endpoints using bivariate normal distribution
  # (equation 13 in Homma and Yoshida 2024)
  powerCoprimary <- .pbivnorm_safe(c_val[1], c_val[2], gamma)

  # Return results as a data frame
  result <- data.frame(
    n1, n2, r1, r2, nu, t, mu1, mu2, sd, rho1, rho2, alpha,
    powerCount = power1and2[1], powerCont = power1and2[2], powerCoprimary
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
