#' Sample Size Calculation for Two Co-Primary Endpoints: One Continuous and One Binary
#'
#' Determines the sample size for a two-arm superiority trial with two co-primary
#' endpoints where one is continuous and one is binary, to achieve a specified
#' power at a given significance level.
#'
#' @param delta Mean difference for the continuous endpoint (group 1 - group 2)
#' @param sd Common standard deviation for the continuous endpoint
#' @param p1 Probability of response in group 1 for the binary endpoint (0 < p1 < 1)
#' @param p2 Probability of response in group 2 for the binary endpoint (0 < p2 < 1)
#' @param rho Biserial correlation between the latent continuous variable underlying
#'   the binary endpoint and the observed continuous endpoint
#' @param r Allocation ratio n1/n2 where n1 is sample size for group 1
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Type II error rate (typically 0.1 or 0.2). Power = 1 - beta
#' @param Test Statistical testing method for the binary endpoint. One of:
#'   \itemize{
#'     \item \code{"AN"}: Asymptotic normal method without continuity correction
#'     \item \code{"ANc"}: Asymptotic normal method with continuity correction
#'     \item \code{"AS"}: Arcsine method without continuity correction
#'     \item \code{"ASc"}: Arcsine method with continuity correction
#'     \item \code{"Fisher"}: Fisher's exact test (uses sequential search)
#'   }
#' @param nMC Number of Monte Carlo replications when Test = "Fisher" (default: 10000)
#'
#' @return A data frame with the following columns:
#'   \item{delta}{Mean difference for continuous endpoint}
#'   \item{sd}{Standard deviation for continuous endpoint}
#'   \item{p1}{Response probability in group 1 for binary endpoint}
#'   \item{p2}{Response probability in group 2 for binary endpoint}
#'   \item{rho}{Biserial correlation}
#'   \item{r}{Allocation ratio}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{Test}{Testing method used for binary endpoint}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{N}{Total sample size (n1 + n2)}
#'
#' @details
#' This function implements the sample size calculation for mixed continuous-binary
#' co-primary endpoints following the methodology in Sozu et al. (2012).
#'
#' **For Fisher's exact test**, sequential search is used because the saw-tooth
#' nature of exact power makes linear extrapolation inappropriate.
#'
#' **For asymptotic methods (AN, ANc, AS, ASc)**, linear extrapolation is used:
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas.
#' Two initial values are calculated:
#' \itemize{
#'   \item n2_0: Based on target power 1 - beta
#'   \item n2_1: Based on adjusted power 1 - sqrt(1-beta) to provide a bracket
#' }
#'
#' \strong{Step 2-4:} Iteratively refine using linear extrapolation until convergence:
#' \deqn{n_2^{new} = \frac{n_2^{(0)}(power^{(1)} - (1-\beta)) - n_2^{(1)}(power^{(0)} - (1-\beta))}
#'       {power^{(1)} - power^{(0)}}}
#'
#' The algorithm converges when \eqn{n_2^{(1)} = n_2^{(0)}}.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination in
#' clinical trials with multiple co-primary endpoints including mixed continuous
#' and binary variables. \emph{Biometrical Journal}, 54(5), 716-729.
#'
#' @examples
#' # Sample size calculation using asymptotic normal method
#' ss2MixedContinuousBinary(
#'   delta = 0.5,
#'   sd = 1,
#'   p1 = 0.6,
#'   p2 = 0.4,
#'   rho = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = 'AN'
#' )
#'
#' # With continuity correction
#' ss2MixedContinuousBinary(
#'   delta = 0.5,
#'   sd = 1,
#'   p1 = 0.6,
#'   p2 = 0.4,
#'   rho = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = 'ANc'
#' )
#'
#' \donttest{
#' # Fisher's exact test (computationally intensive)
#' ss2MixedContinuousBinary(
#'   delta = 0.5,
#'   sd = 1,
#'   p1 = 0.6,
#'   p2 = 0.4,
#'   rho = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = 'Fisher',
#'   nMC = 5000
#' )
#' }
#'
#' @export
ss2MixedContinuousBinary <- function(delta, sd, p1, p2, rho, r, alpha, beta, Test, nMC = 10000) {

  # Input validation
  if (length(delta) != 1 || length(sd) != 1 || length(p1) != 1 ||
      length(p2) != 1 || length(rho) != 1 || length(r) != 1 ||
      length(alpha) != 1 || length(beta) != 1) {
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
  if (r <= 0) {
    stop("r must be positive")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }
  if (beta <= 0 || beta >= 1) {
    stop("beta must be in (0, 1)")
  }
  if (!Test %in% c("AN", "ANc", "AS", "ASc", "Fisher")) {
    stop("Test must be one of: AN, ANc, AS, ASc, Fisher")
  }

  if (Test == "Fisher") {
    # ===== FISHER'S EXACT TEST: Use sequential search =====
    # Monte Carlo simulation has random variation and exact power has saw-tooth pattern
    # Instead, use sequential search.

    # Step 1: Initialize with sample size from single endpoint formulas
    n2 <- max(
      ss1Continuous(delta, sd, r, alpha, beta)[["n2"]],
      ss1BinaryApprox(p1, p2, r, alpha, beta, Test = "AN")[["n2"]]
    )
    n1 <- ceiling(r * n2)

    # Step 2: Calculate power at initial sample size
    power <- power2MixedContinuousBinary(
      n1, n2, delta, sd, p1, p2, rho, alpha, Test, nMC
    )[["powerCoprimary"]]

    # Step 3: Sequential search - increment n2 by 1 until target power is achieved
    while (power < 1 - beta) {
      n2 <- n2 + 1
      n1 <- ceiling(r * n2)
      power <- power2MixedContinuousBinary(
        n1, n2, delta, sd, p1, p2, rho, alpha, Test, nMC
      )[["powerCoprimary"]]
    }

    # Final sample sizes
    N <- n1 + n2

  } else {
    # ===== ASYMPTOTIC METHODS: Use linear extrapolation =====

    # Step 1: Initialize sample sizes using single endpoint formula
    # Calculate sample size for each endpoint separately, then take the maximum
    n2_0 <- max(
      ss1Continuous(delta, sd, r, alpha, beta)[["n2"]],
      ss1BinaryApprox(p1, p2, r, alpha, beta, Test = Test)[["n2"]]
    )

    # For the second initial value, use adjusted target power
    # This provides a bracket for the linear extrapolation
    # The adjusted power 1 - sqrt(1-beta) is typically higher than 1 - beta
    n2_1 <- max(
      ss1Continuous(delta, sd, r, alpha, 1 - (1 - beta) ^ (1/2))[["n2"]],
      ss1BinaryApprox(p1, p2, r, alpha, 1 - (1 - beta) ^ (1/2), Test = Test)[["n2"]]
    )

    # Step 2-4: Iterative refinement using linear extrapolation
    while (n2_1 - n2_0 != 0) {

      # Calculate sample sizes for group 1
      n1_0 <- ceiling(r * n2_0)
      n1_1 <- ceiling(r * n2_1)

      # Calculate power at two candidate sample sizes
      power_n2_0 <- power2MixedContinuousBinary(
        n1_0, n2_0, delta, sd, p1, p2, rho, alpha, Test, nMC
      )[["powerCoprimary"]]
      power_n2_1 <- power2MixedContinuousBinary(
        n1_1, n2_1, delta, sd, p1, p2, rho, alpha, Test, nMC
      )[["powerCoprimary"]]

      # Linear extrapolation to find sample size that achieves target power
      # This solves: power = (1 - beta) for n2
      n2_updated <- '/'(
        n2_0 * (power_n2_1 - (1 - beta)) - n2_1 * (power_n2_0 - (1 - beta)),
        power_n2_1 - power_n2_0
      )

      # Update values for next iteration
      n2_1 <- n2_0
      n2_0 <- ceiling(n2_updated)
    }

    # Final sample sizes
    n2 <- n2_0
    n1 <- ceiling(r * n2)
    N <- n1 + n2

    # Set nMC to NA for known variance
    nMC <- NA
  }

  # Return results as a data frame
  result <- data.frame(
    delta, sd, p1, p2, rho, r, alpha, beta, Test, nMC,
    n1, n2, N
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
