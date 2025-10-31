#' Power Calculation for Two Co-Primary Binary Endpoints (Approximate)
#'
#' Calculates the power for a two-arm superiority trial with two co-primary
#' binary endpoints using approximate methods, as described in Sozu et al. (2010).
#'
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param p11 True probability of responders in group 1 for the first outcome (0 < p11 < 1)
#' @param p12 True probability of responders in group 1 for the second outcome (0 < p12 < 1)
#' @param p21 True probability of responders in group 2 for the first outcome (0 < p21 < 1)
#' @param p22 True probability of responders in group 2 for the second outcome (0 < p22 < 1)
#' @param rho1 Correlation between the two outcomes for group 1
#' @param rho2 Correlation between the two outcomes for group 2
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param Test Statistical testing method. One of:
#'   \itemize{
#'     \item \code{"AN"}: Asymptotic normal method without continuity correction
#'     \item \code{"ANc"}: Asymptotic normal method with continuity correction
#'     \item \code{"AS"}: Arcsine method without continuity correction
#'     \item \code{"ASc"}: Arcsine method with continuity correction
#'   }
#'
#' @return A data frame with the following columns:
#'   \item{n1}{Sample size for group 1}
#'   \item{n2}{Sample size for group 2}
#'   \item{p11, p12, p21, p22}{Response probabilities}
#'   \item{rho1, rho2}{Correlations}
#'   \item{alpha}{One-sided significance level}
#'   \item{Test}{Testing method used}
#'   \item{power1}{Power for the first endpoint alone}
#'   \item{power2}{Power for the second endpoint alone}
#'   \item{powerCoprimary}{Power for both co-primary endpoints}
#'
#' @details
#' This function calculates the power using asymptotic normal approximation or
#' arcsine transformation methods as described in Sozu et al. (2010).
#'
#' The power for co-primary endpoints is calculated using the bivariate normal
#' distribution, accounting for the correlation between endpoints.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination in
#' clinical trials with multiple co-primary binary endpoints. \emph{Statistics in
#' Medicine}, 29(21), 2169-2179.
#'
#' @examples
#' # Example 1: AN method with rho = 0.8 (from Sozu et al. 2010)
#' power2BinaryApprox(
#'   n1 = 507,
#'   n2 = 507,
#'   p11 = 0.95,
#'   p12 = 0.95,
#'   p21 = 0.9,
#'   p22 = 0.9,
#'   rho1 = 0.8,
#'   rho2 = 0.8,
#'   alpha = 0.025,
#'   Test = 'AN'
#' )
#'
#' # Example 2: Balanced design with arcsine method
#' power2BinaryApprox(
#'   n1 = 150,
#'   n2 = 150,
#'   p11 = 0.6,
#'   p12 = 0.5,
#'   p21 = 0.4,
#'   p22 = 0.3,
#'   rho1 = 0.5,
#'   rho2 = 0.5,
#'   alpha = 0.025,
#'   Test = 'AS'
#' )
#'
#' # Example 3: With continuity correction
#' power2BinaryApprox(
#'   n1 = 180,
#'   n2 = 120,
#'   p11 = 0.55,
#'   p12 = 0.45,
#'   p21 = 0.35,
#'   p22 = 0.25,
#'   rho1 = 0.6,
#'   rho2 = 0.6,
#'   alpha = 0.025,
#'   Test = 'ANc'
#' )
#'
#' # Example 4: Zero correlation
#' power2BinaryApprox(
#'   n1 = 250,
#'   n2 = 250,
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0,
#'   rho2 = 0,
#'   alpha = 0.025,
#'   Test = 'AN'
#' )
#'
#' @export
#' @importFrom stats qnorm pnorm
#' @importFrom mvtnorm pmvnorm GenzBretz
power2BinaryApprox <- function(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test) {

  # Check that rho1 is within valid bounds
  bounds1 <- corrbound2Binary(p11, p12)
  if (rho1 < bounds1[1] | rho1 > bounds1[2]) {
    stop(paste0("rho1 must be within [", round(bounds1[1], 4), ", ",
                round(bounds1[2], 4), "]"))
  }

  # Check that rho2 is within valid bounds
  bounds2 <- corrbound2Binary(p21, p22)
  if (rho2 < bounds2[1] | rho2 > bounds2[2]) {
    stop(paste0("rho2 must be within [", round(bounds2[1], 4), ", ",
                round(bounds2[2], 4), "]"))
  }

  # Calculate kappa (allocation ratio)
  kappa <- n1 / n2

  # Standard normal quantile
  z_alpha <- qnorm(1 - alpha)

  if (Test == 'AN') {
    # ===== Asymptotic Normal Method without Continuity Correction =====

    # Calculate variances for each group and endpoint
    nu_1k <- c(p11 * (1 - p11), p12 * (1 - p12))
    nu_2k <- c(p21 * (1 - p21), p22 * (1 - p22))

    # Calculate pooled proportions under H0
    p_k <- (kappa * c(p11, p12) + c(p21, p22)) / (1 + kappa)

    # Calculate differences in proportions
    delta_k <- c(p11 - p21, p12 - p22)

    # Calculate standard errors under H0
    se <- sqrt(p_k * (1 - p_k) * (1 / n1 + 1 / n2))

    # Calculate test statistics under H1
    se_k <- sqrt(nu_1k / n1 + nu_2k / n2)

    # Calculate power for individual endpoints
    power1and2 <- pnorm((delta_k - se * z_alpha) / se_k)

    # Calculate correlation between test statistics
    rho_an <- (rho1 * sqrt(prod(nu_1k)) / n1 + rho2 * sqrt(prod(nu_2k)) / n2) /
      prod(se_k)

    # Calculate power for co-primary endpoints using bivariate normal distribution
    powerCoprimary <- pmvnorm(
      lower = c(-Inf, -Inf),
      upper = (delta_k - se * z_alpha) / se_k,
      mean = c(0, 0),
      corr = matrix(c(1, rho_an, rho_an, 1), ncol = 2),
      algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
      seed = 1
    )[[1]]

  } else if (Test == 'AS') {
    # ===== Arcsine Method without Continuity Correction =====

    # Calculate variances for each group and endpoint
    nu_1k <- c(p11 * (1 - p11), p12 * (1 - p12))
    nu_2k <- c(p21 * (1 - p21), p22 * (1 - p22))

    # Calculate differences using arcsine transformation
    delta_k <- asin(sqrt(c(p11, p12))) - asin(sqrt(c(p21, p22)))

    # Calculate standard errors using arcsine transformation
    se <- sqrt(1 / (4 * n1) + 1 / (4 * n2))
    se_k <- sqrt(nu_1k / (4 * n1 * c(p11, p12)) + nu_2k / (4 * n2 * c(p21, p22)))

    # Calculate power for individual endpoints
    power1and2 <- pnorm((delta_k - se * z_alpha) / se_k)

    # Calculate correlation between test statistics
    rho_arc <- (rho1 * sqrt(prod(nu_1k)) / (4 * n1 * sqrt(prod(c(p11, p12)))) +
                  rho2 * sqrt(prod(nu_2k)) / (4 * n2 * sqrt(prod(c(p21, p22))))) /
      prod(se_k)

    # Calculate power for co-primary endpoints
    powerCoprimary <- pmvnorm(
      lower = c(-Inf, -Inf),
      upper = (delta_k - se * z_alpha) / se_k,
      mean = c(0, 0),
      corr = matrix(c(1, rho_arc, rho_arc, 1), ncol = 2),
      algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
      seed = 1
    )[[1]]

  } else if (Test == 'ANc') {
    # ===== Asymptotic Normal Method with Continuity Correction =====

    # Continuity correction terms
    c1 <- 1 / (2 * n1)
    c2 <- 1 / (2 * n2)

    # Calculate variances for each group and endpoint
    nu_1k <- c(p11, p12) * (1 - c(p11, p12))
    nu_2k <- c(p21, p22) * (1 - c(p21, p22))

    # Calculate pooled proportions under H0 with continuity correction
    p_k <- (kappa * (c(p11, p12) + c1) + c(p21, p22) + c2) / (1 + kappa)

    # Calculate differences in proportions with continuity correction
    delta_k <- c(p11, p12) - c(p21, p22) + c1 - c2

    # Calculate standard errors under H0
    se <- sqrt(p_k * (1 - p_k) * (1 / n1 + 1 / n2))

    # Calculate adjusted variances for continuity correction
    nu_1k_c <- (c(p11, p12) + c1) * (1 - c(p11, p12) - c1)
    nu_2k_c <- (c(p21, p22) + c2) * (1 - c(p21, p22) - c2)
    se_k <- sqrt(nu_1k / (4 * n1 * nu_1k_c) + nu_2k / (4 * n2 * nu_2k_c))

    # Calculate power for individual endpoints with continuity correction
    power1and2 <- pnorm((delta_k - se * z_alpha) / se_k)

    # Calculate correlation between test statistics with continuity correction
    rho_an_c <- (rho1 * sqrt(prod(nu_1k)) / (4 * n1 * sqrt(prod(nu_1k_c))) +
                   rho2 * sqrt(prod(nu_2k)) / (4 * n2 * sqrt(prod(nu_2k_c)))) /
      prod(se_k)

    # Calculate power for co-primary endpoints with continuity correction
    powerCoprimary <- pmvnorm(
      lower = c(-Inf, -Inf),
      upper = (delta_k - se * z_alpha) / se_k,
      mean = c(0, 0),
      corr = matrix(c(1, rho_an_c, rho_an_c, 1), ncol = 2),
      algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
      seed = 1
    )[[1]]

  } else {
    # ===== Arcsine Method with Continuity Correction =====

    # Continuity correction terms
    c1 <- 1 / (2 * n1)
    c2 <- 1 / (2 * n2)

    # Calculate variances for each group and endpoint
    nu_1k <- c(p11, p12) * (1 - c(p11, p12))
    nu_2k <- c(p21, p22) * (1 - c(p21, p22))

    # Calculate differences using arcsine transformation with continuity correction
    delta_k <- (1 / (4 * n1) + 1 / (4 * n2)) * (asin(sqrt(c(p11, p12))) - asin(sqrt(c(p21, p22)))) +
      (1 / (2 * n1) - 1 / (2 * n2)) * (asin(sqrt(c(p11, p12) + c1)) - asin(sqrt(c(p21, p22) + c2)))

    # Calculate standard errors using arcsine transformation
    se <- sqrt(1 / (4 * n1) + 1 / (4 * n2))
    vk_prime <- sqrt(nu_1k / (4 * n1 * c(p11, p12)) + nu_2k / (4 * n2 * c(p21, p22)))
    se_k <- ((1 / (4 * n1) + 1 / (4 * n2)) * (asin(sqrt(c(p11, p12))) - asin(sqrt(c(p21, p22)))) +
               (1 / (2 * n1) - 1 / (2 * n2)) * (asin(sqrt(c(p11, p12) + c1)) - asin(sqrt(c(p21, p22) + c2)))) / vk_prime

    # Calculate adjusted variances for continuity correction
    nu_1k_c <- (c(p11, p12) + c1) * (1 - c(p11, p12) - c1)
    nu_2k_c <- (c(p21, p22) + c2) * (1 - c(p21, p22) - c2)
    se_k <- sqrt(nu_1k / (4 * n1 * nu_1k_c) + nu_2k / (4 * n2 * nu_2k_c))

    # Calculate power for individual endpoints with continuity correction
    power1and2 <- pnorm((delta_k - se * qnorm(1 - alpha)) / se_k)

    # Calculate correlation between test statistics with continuity correction
    rho_arc_c <- '*'(
      1 / prod(se_k),
      '+'(
        rho1 * sqrt(prod(nu_1k)) / (4 * n1 * sqrt(prod(nu_1k_c))),
        rho2 * sqrt(prod(nu_2k)) / (4 * n2 * sqrt(prod(nu_2k_c)))
      )
    )

    # Calculate power for co-primary endpoints with continuity correction
    powerCoprimary <- pmvnorm(
      lower = c(-Inf, -Inf),
      upper = (delta_k - se * qnorm(1 - alpha)) / se_k,
      mean = c(0, 0),
      corr = matrix(c(1, rho_arc_c, rho_arc_c, 1), ncol = 2),
      algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
      seed = 1
    )[[1]]
  }

  # Return results as a data frame
  result <- data.frame(
    n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test,
    power1 = power1and2[1], power2 = power1and2[2], powerCoprimary
  )
  return(result)
}
