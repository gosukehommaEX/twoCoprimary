#' Power Calculation for Two Co-Primary Time-to-Event Endpoints
#'
#' Calculates the power for a two-arm superiority trial with two co-primary
#' time-to-event endpoints under exponential distributions, as described in
#' Hamasaki et al. (2013).
#'
#' @param n1 Sample size for group 1 (test group)
#' @param n2 Sample size for group 2 (control group)
#' @param lambdaT1 Hazard rate for the first endpoint in test group (λT1 > 0)
#' @param lambdaT2 Hazard rate for the second endpoint in test group (λT2 > 0)
#' @param lambdaC1 Hazard rate for the first endpoint in control group (λC1 > 0)
#' @param lambdaC2 Hazard rate for the second endpoint in control group (λC2 > 0)
#' @param tau0 Accrual period (recruitment time in years)
#' @param tau Follow-up period after recruitment ends (years)
#' @param rho Correlation between the two endpoints (0 ≤ rho ≤ 1)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param copula Character string specifying copula type. One of:
#'   \itemize{
#'     \item \code{"clayton"}: Clayton copula (default)
#'     \item \code{"gumbel"}: Gumbel/positive stable copula
#'     \item \code{"frank"}: Frank copula
#'   }
#' @param censor Character string specifying censoring approach. One of:
#'   \itemize{
#'     \item \code{"variance"}: Censoring incorporated into variance only (default)
#'     \item \code{"both"}: Censoring incorporated into both variance and correlation
#'   }
#'
#' @return A data frame with columns for sample sizes, parameters, and powers
#'
#' @details
#' This function calculates the power for two co-primary time-to-event endpoints
#' using bivariate normal approximation of the log-rank test statistics.
#'
#' @references
#' Hamasaki, T., Sugimoto, T., Evans, S. R., & Sozu, T. (2013). Sample size
#' determination for clinical trials with co-primary outcomes: exponential event
#' times. \emph{Pharmaceutical Statistics}, 12(1), 28-34.
#'
#' @examples
#' power2Survival(
#'   n1 = 300, n2 = 300,
#'   lambdaT1 = log(2) / 5, lambdaT2 = log(2) / 5,
#'   lambdaC1 = log(2) / 3, lambdaC2 = log(2) / 3,
#'   tau0 = 2, tau = 5,
#'   rho = 0.5, alpha = 0.025,
#'   copula = "clayton"
#' )
#'
#' @export
#' @importFrom stats qnorm pnorm
#' @importFrom mvtnorm pmvnorm GenzBretz
power2Survival <- function(n1, n2, lambdaT1, lambdaT2, lambdaC1, lambdaC2,
                           tau0, tau, rho, alpha,
                           copula = "clayton", censor = "variance") {

  # Input validation
  if (n1 <= 0 || n2 <= 0) stop("Sample sizes must be positive")
  if (lambdaT1 <= 0 || lambdaT2 <= 0 || lambdaC1 <= 0 || lambdaC2 <= 0) {
    stop("All hazard rates must be positive")
  }
  if (tau0 <= 0) stop("tau0 must be positive")
  if (tau <= tau0) stop("tau must be greater than tau0")
  if (rho < 0 || rho > 1) stop("rho must be between 0 and 1")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)")
  if (!copula %in% c("clayton", "gumbel", "frank")) {
    stop("copula must be 'clayton', 'gumbel', or 'frank'")
  }
  if (!censor %in% c("variance", "both")) {
    stop("censor must be 'variance' or 'both'")
  }

  # Total sample size and allocation ratio
  N <- n1 + n2
  r <- n1 / n2

  # Calculate psi function (censoring parameter)
  psi <- function(lambda) {
    (1 - exp(-lambda * (tau0 + tau)) - exp(-lambda * tau)) / (lambda * tau0)
  }

  psi_T1 <- psi(lambdaT1)
  psi_T2 <- psi(lambdaT2)
  psi_C1 <- psi(lambdaC1)
  psi_C2 <- psi(lambdaC2)

  # Log hazard ratios
  log_HR1 <- log(lambdaC1 / lambdaT1)
  log_HR2 <- log(lambdaC2 / lambdaT2)

  # Standard normal quantile for alpha
  z_alpha <- qnorm(1 - alpha)

  # ===== ENDPOINT 1 =====
  pi_T1 <- (r * lambdaT1 + (1 - r) * lambdaC1) * psi_T1
  pi_C1 <- (r * lambdaC1 + (1 - r) * lambdaC1) * psi_C1
  var_H0_1 <- 1 / (N * (pi_T1 + pi_C1 / r))
  var_H1_1 <- 1 / (N * lambdaT1 * psi_T1) + 1 / (N * lambdaC1 * psi_C1)
  c1 <- (log_HR1 - z_alpha * sqrt(var_H0_1)) / sqrt(var_H1_1)

  # ===== ENDPOINT 2 =====
  pi_T2 <- (r * lambdaT2 + (1 - r) * lambdaC2) * psi_T2
  pi_C2 <- (r * lambdaC2 + (1 - r) * lambdaC2) * psi_C2
  var_H0_2 <- 1 / (N * (pi_T2 + pi_C2 / r))
  var_H1_2 <- 1 / (N * lambdaT2 * psi_T2) + 1 / (N * lambdaC2 * psi_C2)
  c2 <- (log_HR2 - z_alpha * sqrt(var_H0_2)) / sqrt(var_H1_2)

  # ===== INDIVIDUAL POWERS =====
  power1 <- pnorm(c1)
  power2 <- pnorm(c2)

  # ===== CORRELATION =====
  if (censor == "variance") {
    rho_test <- rho
  } else {
    rho_test <- funCopula(
      rho = rho,
      lambdaT1 = lambdaT1, lambdaT2 = lambdaT2,
      lambdaC1 = lambdaC1, lambdaC2 = lambdaC2,
      tau0 = tau0, tau = tau,
      copula = copula
    )
  }

  # ===== CO-PRIMARY POWER =====
  powerCoprimary <- pmvnorm(
    lower = c(-Inf, -Inf),
    upper = c(c1, c2),
    mean = c(0, 0),
    corr = matrix(c(1, rho_test, rho_test, 1), ncol = 2),
    algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
    seed = 1
  )[[1]]

  # Return results
  result <- data.frame(
    n1 = n1, n2 = n2,
    lambdaT1 = lambdaT1, lambdaT2 = lambdaT2,
    lambdaC1 = lambdaC1, lambdaC2 = lambdaC2,
    tau0 = tau0, tau = tau, rho = rho, alpha = alpha,
    copula = copula, censor = censor,
    power1 = power1, power2 = power2, powerCoprimary = powerCoprimary
  )
  return(result)
}
