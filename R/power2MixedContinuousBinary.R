#' Power Calculation for Two Co-Primary Mixed Endpoints
#'
#' Calculates the power for a two-arm superiority trial with two co-primary
#' endpoints where one is continuous and one is binary, as described in
#' Sozu et al. (2012).
#'
#' @param n1 Sample size for group 1 (test group)
#' @param n2 Sample size for group 2 (control group)
#' @param delta Mean difference for the continuous endpoint (group 1 - group 2)
#' @param sd Common standard deviation for the continuous endpoint
#' @param p1 Probability of response in group 1 for the binary endpoint (0 < p1 < 1)
#' @param p2 Probability of response in group 2 for the binary endpoint (0 < p2 < 1)
#' @param rho Biserial correlation between the latent continuous variable underlying
#'   the binary endpoint and the observed continuous endpoint
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param Test Statistical testing method for the binary endpoint. One of:
#'   \itemize{
#'     \item \code{"AN"}: Asymptotic normal method without continuity correction
#'     \item \code{"ANc"}: Asymptotic normal method with continuity correction
#'     \item \code{"AS"}: Arcsine method without continuity correction
#'     \item \code{"ASc"}: Arcsine method with continuity correction
#'     \item \code{"Fisher"}: Fisher's exact test (Monte Carlo simulation required)
#'   }
#' @param nMC Number of Monte Carlo replications when Test = "Fisher" (default: 10000)
#'
#' @return A data frame with the following columns:
#'   \item{n1}{Sample size for group 1}
#'   \item{n2}{Sample size for group 2}
#'   \item{delta}{Mean difference for continuous endpoint}
#'   \item{sd}{Standard deviation for continuous endpoint}
#'   \item{p1}{Response probability in group 1 for binary endpoint}
#'   \item{p2}{Response probability in group 2 for binary endpoint}
#'   \item{rho}{Biserial correlation}
#'   \item{alpha}{One-sided significance level}
#'   \item{Test}{Testing method used for binary endpoint}
#'   \item{powerCont}{Power for the continuous endpoint alone}
#'   \item{powerBin}{Power for the binary endpoint alone}
#'   \item{powerCoprimary}{Power for both co-primary endpoints}
#'
#' @details
#' This function implements the power calculation for mixed endpoints (one continuous
#' and one binary) as described in Sozu et al. (2012). The method assumes that the
#' binary variable is derived from a latent continuous variable via dichotomization
#' at a threshold point.
#'
#' **For Fisher's exact test**, Monte Carlo simulation is used because exact calculation
#' is computationally intensive. The continuous endpoint is analyzed using t-test,
#' and the binary endpoint uses Fisher's exact test.
#'
#' **For asymptotic methods (AN, ANc, AS, ASc)**, analytical formulas are used based
#' on bivariate normal approximation. The correlation between test statistics depends
#' on the biserial correlation rho and the specific testing method.
#'
#' **Biserial Correlation:**
#' The biserial correlation rho represents the correlation between the latent continuous
#' variable underlying the binary endpoint and the observed continuous endpoint. This is
#' not the same as the point-biserial correlation observed in the data.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination in
#' clinical trials with multiple co-primary endpoints including mixed continuous
#' and binary variables. \emph{Biometrical Journal}, 54(5), 716-729.
#'
#' @examples
#' # Power calculation using asymptotic normal method
#' power2MixedContinuousBinary(
#'   n1 = 100,
#'   n2 = 100,
#'   delta = 0.5,
#'   sd = 1,
#'   p1 = 0.6,
#'   p2 = 0.4,
#'   rho = 0.5,
#'   alpha = 0.025,
#'   Test = 'AN'
#' )
#'
#' \donttest{
#' # Power calculation with Fisher's exact test (computationally intensive)
#' power2MixedContinuousBinary(
#'   n1 = 50,
#'   n2 = 50,
#'   delta = 0.5,
#'   sd = 1,
#'   p1 = 0.6,
#'   p2 = 0.4,
#'   rho = 0.5,
#'   alpha = 0.025,
#'   Test = 'Fisher',
#'   nMC = 5000
#' )
#' }
#'
#' @export
#' @importFrom mvtnorm pmvnorm rmvnorm GenzBretz
#' @importFrom stats qnorm pnorm qt pt dhyper phyper
power2MixedContinuousBinary <- function(n1, n2, delta, sd, p1, p2, rho, alpha, Test, nMC = 10000) {

  # Input validation
  if (length(n1) != 1 || length(n2) != 1) {
    stop("n1 and n2 must be scalar values")
  }
  if (n1 <= 0 || n1 != round(n1)) {
    stop("n1 must be a positive integer")
  }
  if (n2 <= 0 || n2 != round(n2)) {
    stop("n2 must be a positive integer")
  }
  if (length(delta) != 1 || length(sd) != 1 || length(p1) != 1 ||
      length(p2) != 1 || length(rho) != 1 || length(alpha) != 1) {
    stop("All parameters must be scalar values")
  }
  if (delta <= 0) {
    stop("delta must be positive")
  }
  if (sd <= 0) {
    stop("sd must be positive")
  }
  if (p1 <= 0 || p1 >= 1) {
    stop("p1 must be in (0, 1)")
  }
  if (p2 <= 0 || p2 >= 1) {
    stop("p2 must be in (0, 1)")
  }
  if (abs(rho) >= 1) {
    stop("rho must be in (-1, 1)")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }
  if (!Test %in% c("AN", "ANc", "AS", "ASc", "Fisher")) {
    stop("Test must be one of: AN, ANc, AS, ASc, Fisher")
  }

  # ===== FISHER'S EXACT TEST (MONTE CARLO SIMULATION) =====
  if (Test == "Fisher") {

    # For biserial correlation model, use fixed cutoff point (g = 0)
    # and vary the means of latent variables to achieve different response probabilities
    # πT = P(X >= 0) where X ~ N(μT, 1) implies μT = Φ^{-1}(πT)
    muT_binary <- qnorm(p1)  # mean for group 1 latent binary variable
    muC_binary <- qnorm(p2)  # mean for group 2 latent binary variable

    # Dichotomy point for latent variable
    g <- qnorm(1 - p2)

    # Mean vectors for test and control groups
    muT <- c(delta, qnorm(1 - p1))
    muC <- c(0, qnorm(1 - p2))

    # Covariance matrices (latent scale)
    SigmaT <- matrix(c(sd ^ 2, rho * sd, rho * sd, 1), nrow = 2)
    SigmaC <- matrix(c(sd ^ 2, rho * sd, rho * sd, 1), nrow = 2)

    # Monte Carlo simulation (vectorized for efficiency)
    # Generate all latent variables at once
    XT <- rmvnorm(nMC * n1, mean = muT, sigma = SigmaT)
    XC <- rmvnorm(nMC * n2, mean = muC, sigma = SigmaC)

    # Continuous endpoint (first column)
    Y_cont_T <- matrix(XT[, 1], nrow = nMC, ncol = n1, byrow = TRUE)
    Y_cont_C <- matrix(XC[, 1], nrow = nMC, ncol = n2, byrow = TRUE)

    # Calculate t-test statistics (vectorized)
    bar_Y_cont_T <- rowMeans(Y_cont_T)
    bar_Y_cont_C <- rowMeans(Y_cont_C)
    hat_delta <- bar_Y_cont_T - bar_Y_cont_C
    S2 <- '+'(
      rowSums((Y_cont_T - bar_Y_cont_T) ^ 2),
      rowSums((Y_cont_C - bar_Y_cont_C) ^ 2)
    ) / (n1 + n2 - 2)
    T_stat <- hat_delta / sqrt(S2 * (1 / n1 + 1 / n2))
    crt_ttest <- qt(alpha, df = n1 + n2 - 2, lower.tail = FALSE)

    # Binary endpoint (second column, dichotomized at g = 0)
    Y_bin_T <- rowSums(
      matrix(as.numeric(XT[, 2] >= g), nrow = nMC, ncol = n1, byrow = TRUE)
    )
    Y_bin_C <- rowSums(
      matrix(as.numeric(XC[, 2] >= g), nrow = nMC, ncol = n2, byrow = TRUE)
    )

    # Fisher's exact test p-values using hypergeometric distribution
    # P(X >= Y_bin_T | margins fixed) = P(X > Y_bin_T - 1)
    p_bin <- phyper(Y_bin_T - 1, n1, n2, Y_bin_T + Y_bin_C, lower.tail = FALSE)

    # Calculate empirical power
    powerCont <- sum(T_stat > crt_ttest) / nMC
    powerBin <- sum(p_bin < alpha) / nMC
    powerCoprimary <- sum((T_stat > crt_ttest) & (p_bin < alpha)) / nMC

  }

  # ===== ASYMPTOTIC METHODS (AN, ANc, AS, ASc) =====
  # Set nMC to NA for asymptotic methods
  nMC <- NA

  # Calculate allocation ratio
  kappa <- n2 / n1  # κ = n2/n1

  # Standard normal quantiles
  z_alpha <- qnorm(1 - alpha)

  # ===== CONTINUOUS ENDPOINT =====
  # Standardized effect size (equation 5 in Sozu et al. 2012)
  Z_cont <- delta / sd * sqrt(kappa * n1 / (1 + kappa))

  # Standardized critical value (equation 5)
  c_val_cont <- -z_alpha + Z_cont

  # Power for continuous endpoint
  powerCont <- pnorm(c_val_cont)

  # ===== BINARY ENDPOINT =====
  # Calculate xi (density at cutoff point) for biserial correlation
  # Assuming standardized latent variables: mu = 0, sigma = 1
  # Calculate non-response probabilities
  thetaT <- 1 - p1
  thetaC <- 1 - p2

  # Cutoff point: g = Φ^(-1)(1 - π)
  xiT <- dnorm(qnorm(1 - p1))
  xiC <- dnorm(qnorm(1 - p2))

  # Effect size for binary endpoint
  delta_bin <- p1 - p2

  if (Test == "AN" || Test == "ANc") {

    # Variance components for Normal approximation (Table 1 in Supporting Info)
    vk0 <- sqrt((p1 + kappa * p2) * (thetaT + kappa * thetaC) / (kappa * (1 + kappa)))
    vk <- sqrt((kappa * p1 * thetaT + p2 * thetaC) / kappa)

    # Standardized critical value for binary endpoint (Table 1 in Supporting Info)
    c_val_binary <- (sqrt(n1) * delta_bin - vk0 * z_alpha) / vk

    # Correlation between continuous and binary test statistics (Table 2 in Supporting Info)
    # γ = (κ * Corr(YTj1,YTj2) * √(πT2θT2) + Corr(YCj1,YCj2) * √(πC2θC2)) /
    #     (√(1+κ) * √(κπT2θT2 + πC2θC2))
    # where Corr(YTjk,YTjk') = ρ * ξTk' / √(πTk'θTk')
    gamma <- '/'(
      kappa * rho * xiT + rho * xiC,
      sqrt(1 + kappa) * sqrt(kappa * p1 * thetaT + p2 * thetaC)
    )

    if (Test == "ANc") {
      # Apply continuity correction (Table 1 in Supporting Info)
      c_val_binary <- c_val_binary - (1 + kappa) / (2 * kappa * vk * sqrt(n1))
    }

  } else {  # Test == "AS" or "ASc"

    # Arcsine transformation parameter
    v <- 0.5 * sqrt((1 + kappa) / kappa)

    if (Test == "AS") {

      # Standardized critical value for binary endpoint (Table 1 in Supporting Info)
      c_val_binary <- -z_alpha + sqrt(n1) * (asin(sqrt(p1)) - asin(sqrt(p2))) / v

      # Correlation between continuous and binary test statistics (Table 2 in Supporting Info)
      # For Arcsine method (k ≤ km < k'): γ = (κ * Corr + Corr) / (1+κ)
      # where Corr(YTj1,YTj2) = ρ * ξT / √(πTθT)
      gamma <- '+'(
        kappa * rho * xiT / sqrt(p1 * thetaT),
        rho * xiC / sqrt(p2 * thetaC)
      ) / (1 + kappa)

    } else {  # Test == "ASc"

      # Continuity correction terms
      cT <- -1 / (2 * n1)
      cC <- 1 / (2 * kappa * n1)

      # Corrected probabilities
      p2T <- p1 - 1 / (2 * n1)
      thetacT <- thetaT + 1 / (2 * n1)
      p2C <- p2 + 1 / (2 * kappa * n1)
      thetacC <- thetaC - 1 / (2 * kappa * n1)

      # Modified variance parameter (Table 1 in Supporting Info)
      vk_prime <- 0.5 * sqrt(
        '+'(
          p1 * thetaT / ((p1 + cT) * (thetaT - cT)),
          p2 * thetaC / (kappa * (p2 + cC) * (thetaC - cC))
        )
      )

      # Standardized critical value with continuity correction (Table 1)
      c_val_binary <- '-'(
        sqrt(n1) * (asin(sqrt(p1 + cT)) - asin(sqrt(p2 + cC))),
        v * qnorm(1 - alpha)
      ) / vk_prime

      # Correlation between continuous and binary test statistics (Table 2)
      # Complex formula for Arcsine(CC) with k ≤ km < k'
      # Note: Explicit Corr() structure for clarity
      # Note: Paper has typo - second term should use πCk'θCk' not πTk'θTk'
      gamma <- '/'(
        '+'(
          kappa * rho * xiT * sqrt(p2C * thetacC),
          rho * xiC * sqrt(p2T * thetacT)
        ),
        '*'(
          sqrt(1 + kappa),
          sqrt(kappa * p1 * thetaT * p2C * thetacC + p2T * thetacT * p2 * thetaC)
        )
      )
    }
  }

  # Power for binary endpoint
  powerBin <- pnorm(c_val_binary)

  # ===== CO-PRIMARY POWER =====
  # Calculate power for co-primary endpoints using bivariate normal distribution
  powerCoprimary <- pmvnorm(
    lower = c(-Inf, -Inf),
    upper = c(c_val_cont, c_val_binary),
    mean = c(0, 0),
    corr = matrix(c(1, gamma, gamma, 1), ncol = 2),
    algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
    seed = 1
  )[[1]]

  # Return results as a data frame
  result <- data.frame(
    n1, n2, delta, sd, p1, p2, rho, alpha, Test, nMC,
    powerCont, powerBin, powerCoprimary
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
