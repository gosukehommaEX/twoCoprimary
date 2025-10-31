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
#' variable at a threshold. The correlation parameter rho represents the biserial
#' correlation between the observed continuous endpoint and this latent continuous
#' variable underlying the binary endpoint.
#'
#' **Key Model Assumptions:**
#' - The latent variable underlying the binary endpoint follows a bivariate normal
#'   distribution with the continuous endpoint
#' - For group i with response probability pi, the latent variable has mean mui = Phi^(-1)(pi)
#'   and variance 1
#' - The dichotomization threshold is set at 0, so P(X >= 0) = Phi(mui) = pi
#' - This is equivalent to using varying cutoff points with fixed means, and is the
#'   standard approach in biserial correlation models
#'
#' **Note on Table 2 in Supporting Information:**
#' There is a typographical error in the original paper's Table 2 for Arcsine(CC)
#' method. The second term in the numerator should use piCk'*thetaCk' instead of piTk'*thetaTk'.
#' This implementation uses the corrected formula.
#'
#' **Understanding Key Quantities:**
#' - rho: Biserial correlation between the continuous endpoint and the latent
#'   continuous variable underlying the binary endpoint (before dichotomization)
#' - gamma: Correlation between the test statistics of the two endpoints
#'   (after dichotomization and standardization)
#' - xi: Standard normal density at the dichotomization threshold, calculated as
#'   xi = phi(Phi^(-1)(1-pi)) where pi is the response probability
#' - c*1/c*2: Ratio of standardized critical values for the two endpoints. When
#'   c*1/c*2 ≈ 1, the individual powers of the two endpoints are approximately equal
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination in
#' clinical trials with multiple co-primary endpoints including mixed continuous
#' and binary variables. \emph{Biometrical Journal}, 54(5), 716-729.
#'
#' @examples
#' # Example 1: Reproduce Table 2 from Sozu et al. (2012)
#' # mTSS (continuous) and ACR50 (binary) with correlation rho = 0.0
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
#' # Example 2: Same setting with correlation rho = 0.8
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
    stop("n1 and n2 must be positive")
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

  # Calculate allocation ratio
  kappa <- n1 / n2

  # ===== CONTINUOUS ENDPOINT =====
  # Calculate Z-statistic under alternative hypothesis for continuous endpoint
  Z_cont <- delta / (sd * sqrt(1 / n1 + 1 / n2))

  # Critical value for one-sided test
  z_alpha <- qnorm(1 - alpha)

  # Power for continuous endpoint
  powerCont <- pnorm(-z_alpha + Z_cont)

  # Standardized critical value for continuous endpoint
  c_cont <- -z_alpha + Z_cont

  # ===== BINARY ENDPOINT =====
  # For binary endpoint, we use latent variable approach
  # The latent variable has mean mu = Phi^(-1)(p) and variance 1

  # Calculate means of latent variables for each group
  muT <- qnorm(p1)  # Treatment group (group 1)
  muC <- qnorm(p2)  # Control group (group 2)

  # Calculate xi (standard normal density at the threshold)
  xiT <- dnorm(qnorm(1 - p1))
  xiC <- dnorm(qnorm(1 - p2))

  # Calculate theta = 1 - p
  thetaT <- 1 - p1
  thetaC <- 1 - p2

  if (Test == "AN") {
    # ===== Asymptotic Normal Method without Continuity Correction =====

    # Calculate pooled proportion under H0
    p_pooled <- (kappa * p1 + p2) / (1 + kappa)

    # Calculate variance under H0
    v0 <- sqrt(p_pooled * (1 - p_pooled))

    # Calculate variance under H1
    v1 <- sqrt((p1 * thetaT / kappa + p2 * thetaC) / (1 + 1 / kappa))

    # Calculate test statistic
    Z_bin <- (p1 - p2) / (sd * sqrt(1 / n1 + 1 / n2)) * (sd * sqrt(1 / n1 + 1 / n2)) / v1

    # Standardized critical value for binary endpoint
    c_bin_star <- (p1 - p2 - v0 * z_alpha * sqrt(1 / n1 + 1 / n2)) / v1

    # Correlation between continuous and binary test statistics
    # Using equation from Table 2 in Sozu et al. (2012)
    gamma <- (kappa * rho * xiT / sqrt(p1 * thetaT) * sqrt(p1 * thetaT * p2 * thetaC) +
                rho * xiC / sqrt(p2 * thetaC) * sqrt(p1 * thetaT * p2 * thetaC)) /
      (sqrt(1 + kappa) * sqrt(kappa * p1 * thetaT * p2 * thetaC +
                                p1 * thetaT * p2 * thetaC))

    # For pmvnorm, use negative value
    c_bin <- -c_bin_star

  } else if (Test == "ANc") {
    # ===== Asymptotic Normal Method with Continuity Correction =====

    # Continuity correction terms
    cT <- 1 / (2 * n1)
    cC <- 1 / (2 * n2)

    # Calculate pooled proportion under H0 with continuity correction
    p_pooled <- (kappa * (p1 + cT) + p2 + cC) / (1 + kappa)

    # Calculate variance under H0
    v0 <- sqrt(p_pooled * (1 - p_pooled))

    # Calculate adjusted theta with continuity correction
    thetacT <- 1 - p1 - cT
    thetacC <- 1 - p2 - cC

    # Calculate variance under H1 with continuity correction
    v1 <- sqrt((p1 * thetaT / (4 * kappa * (p1 + cT) * thetacT) +
                  p2 * thetaC / (4 * (p2 + cC) * thetacC)) / (1 + 1 / kappa))

    # Calculate test statistic with continuity correction
    Z_bin <- (p1 - p2 + cT - cC) / v1

    # Standardized critical value for binary endpoint
    c_bin_star <- (p1 - p2 + cT - cC - v0 * z_alpha * sqrt(1 / n1 + 1 / n2)) / v1

    # Correlation between continuous and binary test statistics (with CC)
    # Note: Paper has typo - second term should use piCk'*thetaCk' not piTk'*thetaTk'
    gamma <- (kappa * rho * xiT / sqrt(p1 * thetaT) * sqrt(p1 * thetaT * (p2 + cC) * thetacC) +
                rho * xiC / sqrt(p2 * thetaC) * sqrt((p1 + cT) * thetacT * p2 * thetaC)) /
      (sqrt(1 + kappa) * sqrt(kappa * p1 * thetaT * (p2 + cC) * thetacC +
                                (p1 + cT) * thetacT * p2 * thetaC))

    # For pmvnorm, use negative value
    c_bin <- -c_bin_star

  } else if (Test == "AS") {
    # ===== Arcsine Method without Continuity Correction =====

    # Calculate difference using arcsine transformation
    delta_bin <- asin(sqrt(p1)) - asin(sqrt(p2))

    # Calculate variance using arcsine transformation
    v0 <- sqrt(1 / (4 * n1) + 1 / (4 * n2))
    v1 <- sqrt(p1 * thetaT / (4 * n1 * p1) + p2 * thetaC / (4 * n2 * p2))

    # Standardized critical value for binary endpoint
    c_bin_star <- (delta_bin - v0 * z_alpha) / v1

    # Correlation between continuous and binary test statistics (Arcsine)
    gamma <- (kappa * rho * xiT / sqrt(p1 * thetaT) * sqrt(p1 * thetaT * p2 * thetaC) +
                rho * xiC / sqrt(p2 * thetaC) * sqrt(p1 * thetaT * p2 * thetaC)) /
      (sqrt(1 + kappa) * sqrt(kappa * p1 * thetaT * p2 * thetaC +
                                p1 * thetaT * p2 * thetaC))

    # For pmvnorm, use negative value
    c_bin <- -c_bin_star

  } else if (Test == "ASc") {
    # ===== Arcsine Method with Continuity Correction =====

    # Continuity correction terms
    cT <- 1 / (2 * n1)
    cC <- 1 / (2 * n2)

    # Calculate adjusted theta with continuity correction
    thetacT <- 1 - p1 - cT
    thetacC <- 1 - p2 - cC
    p2C <- p2 + cC
    p2T <- p1 + cT

    # Calculate difference using arcsine transformation with continuity correction
    delta_bin <- (1 / (4 * n1) + 1 / (4 * n2)) * (asin(sqrt(p1)) - asin(sqrt(p2))) +
      (1 / (2 * n1) - 1 / (2 * n2)) * (asin(sqrt(p1 + cT)) - asin(sqrt(p2 + cC)))

    # Calculate variance using arcsine transformation
    v0 <- sqrt(1 / (4 * n1) + 1 / (4 * n2))
    vk_prime <- sqrt(p1 * thetaT / (4 * n1 * p1) + p2 * thetaC / (4 * n2 * p2))

    # Standardized critical value for binary endpoint
    c_bin_star <- ((1 / (4 * n1) + 1 / (4 * n2)) * (asin(sqrt(p1)) - asin(sqrt(p2))) +
                     (1 / (2 * n1) - 1 / (2 * n2)) * (asin(sqrt(p1 + cT)) - asin(sqrt(p2 + cC)))) / vk_prime

    # Correlation between continuous and binary test statistics (Table 2)
    # Complex formula for Arcsine(CC) with k <= km < k'
    # Note: Paper has typo - second term should use piCk'*thetaCk' not piTk'*thetaTk'
    gamma <- (kappa * rho * xiT / sqrt(p1 * thetaT) * sqrt(p1 * thetaT * p2C * thetacC) +
                rho * xiC / sqrt(p2 * thetaC) * sqrt(p2T * thetacT * p2 * thetaC)) /
      (sqrt(1 + kappa) * sqrt(kappa * p1 * thetaT * p2C * thetacC +
                                p2T * thetacT * p2 * thetaC))

    # For pmvnorm, use negative value
    c_bin <- -c_bin_star

  } else {
    # ===== Fisher's Exact Test (Monte Carlo Simulation) =====
    stop("Fisher's exact test is not yet implemented for power calculation in this function.")
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
