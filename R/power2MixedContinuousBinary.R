#' Power Calculation for Two Co-Primary Endpoints: One Continuous and One Binary
#'
#' Calculates the power for a two-arm superiority trial with two co-primary
#' endpoints where one is continuous and one is binary, as described in
#' Sozu et al. (2012). This function is specifically designed for mixed
#' continuous-binary endpoint combinations.
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
#' This function implements the power calculation for mixed continuous-binary endpoints
#' as described in Sozu et al. (2012). The method assumes that the binary variable is
#' derived from a latent continuous variable via dichotomization at a threshold point.
#'
#' **Endpoint Types:**
#' \itemize{
#'   \item \strong{Continuous Endpoint}: Analyzed using t-test for comparing means
#'   \item \strong{Binary Endpoint}: Analyzed using one of several methods (AN, ANc, AS, ASc, or Fisher)
#' }
#'
#' **Biserial Correlation Model:**
#' The binary endpoint is assumed to arise from dichotomizing a latent continuous
#' variable at a threshold. The correlation parameter ρ (rho) represents the biserial
#' correlation between the observed continuous endpoint and this latent continuous
#' variable underlying the binary endpoint.
#'
#' **Key Model Assumptions:**
#' - The latent variable underlying the binary endpoint follows a bivariate normal
#'   distribution with the continuous endpoint
#' - For group i with response probability πᵢ, the latent variable has mean μᵢ = Φ⁻¹(πᵢ)
#'   and variance 1
#' - The dichotomization threshold is set at 0, so P(X ≥ 0) = Φ(μᵢ) = πᵢ
#' - This is equivalent to using varying cutoff points with fixed means, and is the
#'   standard approach in biserial correlation models
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
#' power2MixedContinuousBinary(
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
#' power2MixedContinuousBinary(
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
#' power2MixedContinuousBinary(
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
#' power2MixedContinuousBinary(
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
#' power2MixedContinuousBinary(
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
#' power2MixedContinuousBinary(
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
power2MixedContinuousBinary <- function(n1, n2, delta, sd, p1, p2, rho, alpha, Test, nMC = 10000) {

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
  c_cont_star <- sqrt((kappa / (1 + kappa)) * (delta / sd)^2)

  # Standardized critical value for continuous endpoint
  c_cont <- -qnorm(alpha) + c_cont_star

  # Power for continuous endpoint
  powerCont <- pnorm(c_cont)

  # ===== BINARY ENDPOINT =====
  # Standard normal density at dichotomization thresholds
  xiT <- dnorm(qnorm(1 - p1))
  xiC <- dnorm(qnorm(1 - p2))

  # Select method for binary endpoint
  if (Test == "AN") {
    # Asymptotic Normal method without continuity correction (Table 2, equation 8)

    # Variance components (no continuity correction)
    vk <- p1 * thetaT / n1 + p2 * thetaC / n2

    # Standardized effect size (equation 8)
    c_bin_star <- (p1 - p2) / sqrt(vk)

    # Correlation between continuous and binary test statistics (Table 2)
    gamma <- (rho * xiT * sqrt(p1 * thetaT / n1) + rho * xiC * sqrt(p2 * thetaC / n2)) /
      sqrt((1 + kappa) * vk)

    # For pmvnorm, use negative value
    c_bin <- -c_bin_star

  } else if (Test == "ANc") {
    # Asymptotic Normal method with continuity correction (Table 2, equation 9)

    # Continuity corrections
    cT <- 1 / (2 * n1)
    cC <- 1 / (2 * n2)

    # Variance components with continuity correction
    vk <- (p1 + cT) * (thetaT - cT) / n1 + (p2 + cC) * (thetaC - cC) / n2

    # Standardized effect size (equation 9)
    c_bin_star <- ((p1 + cT) - (p2 + cC)) / sqrt(vk)

    # Additional probabilities for correlation calculation
    p1T <- p1 + cT
    p2T <- p2 + 2 * cT
    p2C <- p2 + cC
    thetacT <- thetaT - cT
    thetacC <- thetaC - cC

    # Correlation between continuous and binary test statistics (Table 2)
    # Complex formula for AN(CC)
    gamma <- (kappa * rho * xiT * sqrt(p1T * thetacT / n1) + rho * xiC * sqrt(p2C * thetacC / n2)) /
      (sqrt(1 + kappa) * sqrt(kappa * p1T * thetacT / n1 + p2C * thetacC / n2))

    # For pmvnorm, use negative value
    c_bin <- -c_bin_star

  } else if (Test == "AS") {
    # Arcsine method without continuity correction (Table 2, equation 10)

    # Variance components (no continuity correction)
    vk_prime <- 1 / (4 * n1) + 1 / (4 * n2)

    # Standardized effect size (equation 10)
    c_bin_star <- (asin(sqrt(p1)) - asin(sqrt(p2))) / sqrt(vk_prime)

    # Correlation between continuous and binary test statistics (Table 2)
    gamma <- (kappa * rho * xiT / sqrt(p1 * thetaT) * sqrt(p1 * thetaT / n1) +
                rho * xiC / sqrt(p2 * thetaC) * sqrt(p2 * thetaC / n2)) /
      (sqrt(1 + kappa) * sqrt(kappa * p1 * thetaT / n1 + p2 * thetaC / n2))

    # For pmvnorm, use negative value
    c_bin <- -c_bin_star

  } else if (Test == "ASc") {
    # Arcsine method with continuity correction (Table 2, equation 11)

    # Continuity corrections
    cT <- 1 / (2 * n1)
    cC <- 1 / (2 * n2)

    # Variance components with continuity correction
    vk_prime <- 1 / (4 * n1) + 1 / (4 * n2)

    # Additional probabilities for correlation calculation
    p1T <- p1 + cT
    p2T <- p2 + 2 * cT
    p2C <- p2 + cC
    thetaT <- 1 - p1
    thetacT <- thetaT - cT
    thetacC <- thetaC - cC

    # Standardized effect size (equation 11)
    c_bin_star <- (sqrt(1 / (4 * n1) + 1 / (4 * n2)) * (asin(sqrt(p1)) - asin(sqrt(p2))) +
                     (1 / (2 * n1) - 1 / (2 * n2)) * (asin(sqrt(p1 + cT)) - asin(sqrt(p2 + cC)))) / vk_prime

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

  # Power for binary endpoint: P(Z* > c_bin_star) = P(Z* < -c_bin_star)
  powerBin <- pnorm(c_bin)

  # ===== CO-PRIMARY POWER =====
  # Calculate power for co-primary endpoints using bivariate normal distribution
  powerCoprimary <- pmvnorm(
    lower = c(-Inf, -Inf),
    upper = c(c_cont, c_bin),
    mean = c(0, 0),
    corr = matrix(c(1, gamma, gamma, 1), ncol = 2),
    algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0),
    seed = 1
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
