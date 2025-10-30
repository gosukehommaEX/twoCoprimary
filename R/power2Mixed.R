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
#' **Notation:**
#' - Group 1 corresponds to test group (T), Group 2 corresponds to control group (C)
#' - κ (kappa) = n2/n1: allocation ratio (control/test)
#' - n = n1: sample size in test group
#' - π1, π2: response probabilities in groups 1 and 2
#' - θ1 = 1 - π1, θ2 = 1 - π2: non-response probabilities
#'
#' For the continuous endpoint, the standardized test statistic is:
#' \deqn{Z_{cont} = \frac{\delta}{\sigma} \sqrt{\frac{\kappa n}{1+\kappa}}}
#'
#' For the binary endpoint (AN method), the test statistic uses:
#' \deqn{v_{k0} = \sqrt{\frac{(\pi_1 + \kappa\pi_2)(\theta_1 + \kappa\theta_2)}{\kappa(1+\kappa)}}}
#' \deqn{v_k = \sqrt{\frac{\kappa\pi_1\theta_1 + \pi_2\theta_2}{\kappa}}}
#'
#' The correlation between the two test statistics (continuous and binary) is:
#' \deqn{\gamma = \frac{\kappa \rho \xi_1 + \rho \xi_2}{\sqrt{1+\kappa}\sqrt{\kappa\pi_1\theta_1 + \pi_2\theta_2}}}
#' where ξ_1 = φ(Φ^{-1}(1-π_1)) and ξ_2 = φ(Φ^{-1}(1-π_2)) are the standard normal
#' density values at the dichotomization thresholds. This γ represents the correlation
#' between the standardized test statistics Z*_cont and Z*_bin (Equation 1 in paper).
#'
#' Note: ρ is the biserial correlation between the latent continuous variable
#' underlying the binary endpoint and the observed continuous endpoint, not the
#' point-biserial correlation between the observed binary and continuous variables.
#' The observed correlation is given by Equation (1): Corr(Y_cont, Y_bin) = ρξ/√(πθ),
#' which is always smaller than ρ in absolute value.
#'
#' **Fisher's Exact Test:**
#' When Test = "Fisher", Monte Carlo simulation is used to estimate power. For each
#' replication, latent bivariate normal variables are generated using rmvnorm(), the
#' binary variable is dichotomized, and both t-test (for continuous) and Fisher's
#' exact test (for binary) are performed. The overall power is the proportion of
#' replications where both tests are significant at the α level. Use set.seed()
#' before calling this function for reproducible results with Fisher's method.
#'
#' The implementation uses a fixed cutoff point (g = 0) and varies the means of
#' latent variables to achieve different response probabilities. This is equivalent
#' to using varying cutoff points with fixed means, and is the standard approach
#' in biserial correlation models. For group i with response probability πᵢ,
#' the latent variable has mean μᵢ = Φ⁻¹(πᵢ) and variance 1, with P(X ≥ 0) = Φ(μᵢ) = πᵢ.
#'
#' **Note on Table 2 in Supporting Information:**
#' There is a typographical error in the original paper's Table 2 for Arcsine(CC)
#' method. The second term in the numerator should use πCk'θCk' instead of πTk'θTk'.
#' This implementation uses the corrected formula.
#'
#' **Understanding Key Quantities:**
#' - ρ (rho): Biserial correlation between the continuous endpoint and the latent
#'   continuous variable underlying the binary endpoint (before dichotomization)
#' - γ (gamma): Correlation between the test statistics of the two endpoints
#'   (after dichotomization and standardization)
#' - ξ: Standard normal density at the dichotomization threshold, calculated as
#'   ξ = φ(Φ⁻¹(1-π)) where π is the response probability
#' - c*₁/c*₂: Ratio of standardized critical values for the two endpoints. When
#'   c*₁/c*₂ ≈ 1, the individual powers of the two endpoints are approximately equal
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination in
#' clinical trials with multiple co-primary endpoints including mixed continuous
#' and binary variables. \emph{Biometrical Journal}, 54(5), 716-729.
#'
#' @examples
#' # Example 1: Reproduce Table 2 from Sozu et al. (2012)
#' # mTSS (continuous) and ACR50 (binary) with correlation ρ = 0.0
#' power2Mixed(
#'   n1 = 346,
#'   n2 = 346,
#'   delta = 4.4,
#'   sd = 19.0,
#'   p1 = 0.59,
#'   p2 = 0.46,
#'   rho = 0.0,
#'   alpha = 0.025,
#'   Test = "AN"
#' )
#'
#' # Example 2: Same setting with correlation ρ = 0.8
#' power2Mixed(
#'   n1 = 323,
#'   n2 = 323,
#'   delta = 4.4,
#'   sd = 19.0,
#'   p1 = 0.59,
#'   p2 = 0.46,
#'   rho = 0.8,
#'   alpha = 0.025,
#'   Test = "AN"
#' )
#'
#' # Example 3: Normal approximation with continuity correction
#' power2Mixed(
#'   n1 = 150,
#'   n2 = 150,
#'   delta = 0.3,
#'   sd = 1.0,
#'   p1 = 0.65,
#'   p2 = 0.50,
#'   rho = 0.5,
#'   alpha = 0.025,
#'   Test = "ANc"
#' )
#'
#' # Example 4: Arcsine transformation without continuity correction
#' power2Mixed(
#'   n1 = 180,
#'   n2 = 120,
#'   delta = 4.5,
#'   sd = 18.0,
#'   p1 = 0.70,
#'   p2 = 0.55,
#'   rho = 0.4,
#'   alpha = 0.025,
#'   Test = "AS"
#' )
#'
#' # Example 5: Arcsine transformation with continuity correction
#' power2Mixed(
#'   n1 = 100,
#'   n2 = 100,
#'   delta = 0.4,
#'   sd = 1.2,
#'   p1 = 0.55,
#'   p2 = 0.40,
#'   rho = 0.3,
#'   alpha = 0.025,
#'   Test = "ASc"
#' )
#'
#' \donttest{
#' # Example 6: Fisher's exact test (Monte Carlo simulation)
#' # Useful for small sample sizes or rare events
#' set.seed(12345)  # for reproducibility
#' power2Mixed(
#'   n1 = 50,
#'   n2 = 50,
#'   delta = 0.5,
#'   sd = 1.0,
#'   p1 = 0.4,
#'   p2 = 0.2,
#'   rho = 0.3,
#'   alpha = 0.025,
#'   Test = "Fisher",
#'   nMC = 10000
#' )
#' }
#'
#' @export
#' @importFrom mvtnorm pmvnorm GenzBretz rmvnorm
#' @importFrom stats qnorm pnorm dnorm phyper pt
power2Mixed <- function(n1, n2, delta, sd, p1, p2, rho, alpha, Test, nMC = 10000) {

  # Input validation
  if (n1 <= 0 || n2 <= 0) {
    stop("Sample sizes must be positive")
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
  if (Test == "Fisher" && nMC < 1000) {
    warning("nMC should be at least 1000 for reliable results with Fisher's exact test")
  }

  # Calculate allocation ratio
  kappa <- n2 / n1  # κ = n2/n1

  # Calculate non-response probabilities
  thetaT <- 1 - p1
  thetaC <- 1 - p2

  # ===== FISHER'S EXACT TEST (MONTE CARLO SIMULATION) =====
  if (Test == "Fisher") {

    # For biserial correlation model, use fixed cutoff point (g = 0)
    # and vary the means of latent variables to achieve different response probabilities
    # πT = P(X >= 0) where X ~ N(μT, 1) implies μT = Φ^{-1}(πT)
    muT_binary <- qnorm(p1)  # mean for group 1 latent binary variable
    muC_binary <- qnorm(p2)  # mean for group 2 latent binary variable

    # Fixed cutoff point for dichotomization
    g <- 0

    # Mean vectors for latent variables under alternative hypothesis
    muT <- c(delta, muT_binary)  # continuous has mean delta, binary has mean μT
    muC <- c(0, muC_binary)      # continuous has mean 0, binary has mean μC

    # Covariance matrices (biserial correlation structure)
    SigmaT <- matrix(c(
      sd^2, rho * sd,
      rho * sd, 1
    ), nrow = 2, byrow = TRUE)

    SigmaC <- matrix(c(
      sd^2, rho * sd,
      rho * sd, 1
    ), nrow = 2, byrow = TRUE)

    # Monte Carlo simulation (vectorized for efficiency)
    # Note: Use set.seed() before calling this function for reproducible results

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
    S2 <- (rowSums((Y_cont_T - bar_Y_cont_T)^2) +
             rowSums((Y_cont_C - bar_Y_cont_C)^2)) / (n1 + n2 - 2)
    T_stat <- hat_delta / sqrt(S2 * (1 / n1 + 1 / n2))
    p_cont <- pt(T_stat, df = n1 + n2 - 2, lower.tail = FALSE)

    # Binary endpoint (second column, dichotomized at g = 0)
    Y_bin_T <- rowSums(matrix(as.numeric(XT[, 2] >= g), nrow = nMC, ncol = n1, byrow = TRUE))
    Y_bin_C <- rowSums(matrix(as.numeric(XC[, 2] >= g), nrow = nMC, ncol = n2, byrow = TRUE))

    # Fisher's exact test p-values using hypergeometric distribution
    # P(X >= Y_bin_T | margins fixed) = P(X > Y_bin_T - 1)
    p_bin <- phyper(Y_bin_T - 1, n1, n2, Y_bin_T + Y_bin_C, lower.tail = FALSE)

    # Calculate empirical power
    powerCont <- sum(p_cont < alpha) / nMC
    powerBin <- sum(p_bin < alpha) / nMC
    powerCoprimary <- sum((p_cont < alpha) & (p_bin < alpha)) / nMC

    # Return results
    result <- data.frame(
      n1 = n1,
      n2 = n2,
      delta = delta,
      sd = sd,
      p1 = p1,
      p2 = p2,
      rho = rho,
      alpha = alpha,
      Test = Test,
      nMC = nMC,
      powerCont = powerCont,
      powerBin = powerBin,
      powerCoprimary = powerCoprimary
    )

    return(result)
  }

  # ===== ASYMPTOTIC METHODS (AN, ANc, AS, ASc) =====

  # ===== CONTINUOUS ENDPOINT =====
  # Standardized effect size (equation 5 in Sozu et al. 2012)
  Z_cont <- delta * sqrt(kappa * n1 / (1 + kappa)) / sd

  # Standardized critical value (equation 5)
  c_cont_star <- qnorm(1 - alpha) - Z_cont

  # For pmvnorm, use negative value
  c_cont <- -c_cont_star

  # Power for continuous endpoint: P(Z* > c_cont_star) = P(Z* < -c_cont_star)
  powerCont <- pnorm(c_cont)

  # ===== BINARY ENDPOINT =====
  # Calculate xi (density at cutoff point) for biserial correlation
  # Assuming standardized latent variables: mu = 0, sigma = 1
  # Cutoff point: g = Φ^(-1)(1 - π)
  xiT <- dnorm(qnorm(1 - p1))
  xiC <- dnorm(qnorm(1 - p2))

  # Effect size for binary endpoint
  delta_bin <- p1 - p2

  if (Test == "AN" || Test == "ANc") {

    # Variance components for Normal approximation (Table 1 in Supporting Info)
    vk0 <- sqrt((p1 + kappa * p2) * (thetaT + kappa * thetaC) / (kappa * (1 + kappa)))
    vk <- sqrt((kappa * p1 * thetaT + p2 * thetaC) / kappa)

    # Standardized critical value for binary endpoint (Table 1)
    c_bin_star <- (vk0 * qnorm(1 - alpha) - sqrt(n1) * delta_bin) / vk

    # Correlation between continuous and binary test statistics (Table 2)
    # γ = (κ * Corr(YTj1,YTj2) * √(πT2θT2) + Corr(YCj1,YCj2) * √(πC2θC2)) /
    #     (√(1+κ) * √(κπT2θT2 + πC2θC2))
    # where Corr(YTjk,YTjk') = ρ * ξTk' / √(πTk'θTk')
    gamma <- (kappa * rho * xiT + rho * xiC) /
      (sqrt(1 + kappa) * sqrt(kappa * p1 * thetaT + p2 * thetaC))

    if (Test == "ANc") {
      # Apply continuity correction (Table 1)
      c_bin_star <- c_bin_star + (1 + kappa) / (2 * kappa * vk * sqrt(n1))
    }

    # For pmvnorm, use negative value
    c_bin <- -c_bin_star

  } else {  # Test == "AS" or "ASc"

    # Arcsine transformation parameter
    v <- 0.5 * sqrt((1 + kappa) / kappa)

    if (Test == "AS") {

      # Standardized critical value for binary endpoint (Table 1)
      c_bin_star <- qnorm(1 - alpha) - sqrt(n1) * (asin(sqrt(p1)) - asin(sqrt(p2))) / v

      # Correlation between continuous and binary test statistics (Table 2)
      # For Arcsine method (k ≤ km < k'): γ = (κ * Corr + Corr) / (1+κ)
      # where Corr(YTj1,YTj2) = ρ * ξT / √(πTθT)
      gamma <- (kappa * rho * xiT / sqrt(p1 * thetaT) +
                  rho * xiC / sqrt(p2 * thetaC)) / (1 + kappa)

      # For pmvnorm, use negative value
      c_bin <- -c_bin_star

    } else {  # Test == "ASc"

      # Continuity correction terms
      cT <- -1 / (2 * n1)
      cC <- 1 / (2 * kappa * n1)

      # Corrected probabilities
      p2T <- p1 - 1 / (2 * n1)
      thetacT <- thetaT + 1 / (2 * n1)
      p2C <- p2 + 1 / (2 * kappa * n1)
      thetacC <- thetaC - 1 / (2 * kappa * n1)

      # Modified variance parameter (Table 1)
      vk_prime <- 0.5 * sqrt(
        p1 * thetaT / ((p1 + cT) * (thetaT - cT)) +
          p2 * thetaC / (kappa * (p2 + cC) * (thetaC - cC))
      )

      # Standardized critical value with continuity correction (Table 1)
      c_bin_star <- (v * qnorm(1 - alpha) -
                       sqrt(n1) * (asin(sqrt(p1 + cT)) - asin(sqrt(p2 + cC)))) / vk_prime

      # Correlation between continuous and binary test statistics (Table 2)
      # Complex formula for Arcsine(CC) with k ≤ km < k'
      # Note: Explicit Corr() structure for clarity
      # Note: Paper has typo - second term should use πCk'θCk' not πTk'θTk'
      gamma <- (kappa * rho * xiT / sqrt(p1 * thetaT) * sqrt(p1 * thetaT * p2C * thetacC) +
                  rho * xiC / sqrt(p2 * thetaC) * sqrt(p2T * thetacT * p2 * thetaC)) /
        (sqrt(1 + kappa) * sqrt(kappa * p1 * thetaT * p2C * thetacC +
                                  p2T * thetacT * p2 * thetaC))

      # For pmvnorm, use negative value
      c_bin <- -c_bin_star
    }
  }

  # Power for binary endpoint: P(Z* > c_bin_star) = P(Z* < -c_bin_star)
  powerBin <- pnorm(c_bin)

  # ===== CO-PRIMARY POWER =====
  # Calculate power for co-primary endpoints using bivariate normal distribution
  powerCoprimary <- pmvnorm(
    lower = c(-Inf, -Inf),
    upper = c(c_cont, c_bin),
    mean = c(0, 0),
    corr = matrix(c(1, gamma, gamma, 1), ncol = 2),
    algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
  )[[1]]

  # Return results as a data frame
  result <- data.frame(
    n1 = n1,
    n2 = n2,
    delta = delta,
    sd = sd,
    p1 = p1,
    p2 = p2,
    rho = rho,
    alpha = alpha,
    Test = Test,
    powerCont = powerCont,
    powerBin = powerBin,
    powerCoprimary = powerCoprimary
  )

  return(result)
}
