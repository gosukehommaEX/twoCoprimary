#' Power Calculation for Two Co-Primary Binary Endpoints (Approximate)
#'
#' Calculates the power for a two-arm superiority trial with two co-primary
#' binary endpoints using various asymptotic normal approximation methods, as
#' described in Sozu et al. (2010).
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
#' This function implements four approximate power calculation methods:
#'
#' \strong{Asymptotic Normal (AN):} Uses the standard normal approximation without
#' continuity correction (equations 3-4 in Sozu et al. 2010).
#'
#' \strong{Asymptotic Normal with Continuity Correction (ANc):} Includes Yates's
#' continuity correction (equation 5 in Sozu et al. 2010).
#'
#' \strong{Arcsine (AS):} Uses arcsine transformation without continuity correction
#' (equation 6 in Sozu et al. 2010).
#'
#' \strong{Arcsine with Continuity Correction (ASc):} Arcsine method with continuity
#' correction (equation 7 in Sozu et al. 2010).
#'
#' The correlation between test statistics for the two endpoints is calculated as:
#' \deqn{\rho_{nml} = \frac{\sum_{j=1}^{2} \rho_j \sqrt{\nu_{j1}\nu_{j2}}/n_j}
#'       {se_1 \times se_2}}
#' where \eqn{\nu_{jk} = p_{jk}(1-p_{jk})}.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination in
#' clinical trials with multiple co-primary binary endpoints. \emph{Statistics in
#' Medicine}, 29(21), 2169-2179.
#'
#' @examples
#' # Power calculation using asymptotic normal method
#' power2BinaryApprox(
#'   n1 = 200,
#'   n2 = 100,
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0.7,
#'   rho2 = 0.7,
#'   alpha = 0.025,
#'   Test = 'AN'
#' )
#'
#' # Power calculation with continuity correction
#' power2BinaryApprox(
#'   n1 = 200,
#'   n2 = 100,
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0.7,
#'   rho2 = 0.7,
#'   alpha = 0.025,
#'   Test = 'ANc'
#' )
#'
#' # Power calculation using arcsine method
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
#' @export
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats qnorm
power2BinaryApprox <- function(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test) {

  # Define common variance parameters
  nu_1k <- c(p11, p12) * (1 - c(p11, p12))
  nu_2k <- c(p21, p22) * (1 - c(p21, p22))

  if (Test == "AN" | Test == "ANc") {

    # Calculate pooled proportions for each endpoint
    p_k <- (n1 * c(p11, p12) + n2 * c(p21, p22)) / (n1 + n2)

    # Calculate treatment effects
    delta_k <- c(p11, p12) - c(p21, p22)

    # Calculate standard errors under null and alternative hypotheses
    se_k0 <- sqrt((1 / n1 + 1 / n2) * p_k * (1 - p_k))
    se_k <- sqrt(nu_1k / n1 + nu_2k / n2)

    # Calculate correlation between test statistics (equation 4 in Sozu et al. 2010)
    rho_nml <- (rho1 * sqrt(prod(nu_1k)) / n1 + rho2 * sqrt(prod(nu_2k)) / n2) / prod(se_k)

    if (Test == "AN") {

      # Asymptotic normal method without continuity correction

      # Calculate power for individual endpoints
      power1and2 <- pnorm((delta_k - se_k0 * qnorm(1 - alpha)) / se_k)

      # Calculate power for co-primary endpoints
      powerCoprimary <- pmvnorm(
        lower = c(-Inf, -Inf),
        upper = (delta_k - se_k0 * qnorm(1 - alpha)) / se_k,
        mean = c(0, 0),
        corr = matrix(c(1, rho_nml, rho_nml, 1), ncol = 2),
        algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
      )[[1]]

    } else if (Test == "ANc") {

      # Asymptotic normal method with Yates's continuity correction

      # Calculate continuity correction
      ycc <- (1 / 2) * (1 / n1 + 1 / n2)

      # Calculate power for individual endpoints with continuity correction
      power1and2 <- pnorm((delta_k - se_k0 * qnorm(1 - alpha) - ycc) / se_k)

      # Calculate power for co-primary endpoints with continuity correction
      powerCoprimary <- pmvnorm(
        lower = c(-Inf, -Inf),
        upper = (delta_k - se_k0 * qnorm(1 - alpha) - ycc) / se_k,
        mean = c(0, 0),
        corr = matrix(c(1, rho_nml, rho_nml, 1), ncol = 2),
        algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
      )[[1]]
    }

  } else if (Test == "AS" | Test == "ASc") {

    # Standard error for arcsine transformation
    se <- (1 / 2) * sqrt(1 / n1 + 1 / n2)

    if (Test == "AS") {

      # Arcsine method without continuity correction

      # Apply arcsine transformation to probabilities
      delta_k <- asin(sqrt(c(p11, p12))) - asin(sqrt(c(p21, p22)))

      # Calculate power for individual endpoints
      power1and2 <- pnorm(delta_k / se - qnorm(1 - alpha))

      # Calculate correlation between test statistics for arcsine method
      rho_arc <- (n2 * rho1 + n1 * rho2) / (n1 + n2)

      # Calculate power for co-primary endpoints
      powerCoprimary <- pmvnorm(
        lower = c(-Inf, -Inf),
        upper = delta_k / se - qnorm(1 - alpha),
        mean = c(0, 0),
        corr = matrix(c(1, rho_arc, rho_arc, 1), ncol = 2),
        algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
      )[[1]]

    } else if (Test == "ASc") {

      # Arcsine method with continuity correction

      # Define continuity correction constants
      c1 <- -1 / (2 * n1)
      c2 <-  1 / (2 * n2)

      # Apply arcsine transformation with continuity correction
      delta_k <- asin(sqrt(c(p11, p12) + c1)) - asin(sqrt(c(p21, p22) + c2))

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
        algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
      )[[1]]
    }
  }

  # Return results as a data frame
  result <- data.frame(
    n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test,
    power1 = power1and2[1], power2 = power1and2[2], powerCoprimary
  )
  return(result)
}
