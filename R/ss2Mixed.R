#' Sample Size Calculation for Two Co-Primary Mixed Endpoints
#'
#' Calculates the required sample size for a two-arm superiority trial with two
#' co-primary endpoints where one is continuous and one is binary, using the
#' linear extrapolation approach.
#'
#' @param delta Mean difference for the continuous endpoint (group 1 - group 2)
#' @param sd Common standard deviation for the continuous endpoint
#' @param p1 Probability of response in group 1 for the binary endpoint (0 < p1 < 1)
#' @param p2 Probability of response in group 2 for the binary endpoint (0 < p2 < 1)
#' @param rho Biserial correlation between the latent continuous variable underlying
#'   the binary endpoint and the observed continuous endpoint
#' @param r Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#' @param Test Statistical testing method for the binary endpoint. One of:
#'   \itemize{
#'     \item \code{"AN"}: Asymptotic normal method without continuity correction
#'     \item \code{"ANc"}: Asymptotic normal method with continuity correction
#'     \item \code{"AS"}: Arcsine method without continuity correction
#'     \item \code{"ASc"}: Arcsine method with continuity correction
#'     \item \code{"Fisher"}: Fisher's exact test (Monte Carlo simulation)
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
#'   \item{n}{Total sample size (n1 + n2)}
#'
#' @details
#' This function uses the linear extrapolation approach to find the required
#' sample size iteratively, as described in Hamasaki et al. (2013) and extended
#' to mixed endpoints following Sozu et al. (2012).
#'
#' **Algorithm:**
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas.
#' The initial value is the maximum of:
#' \itemize{
#'   \item Sample size for continuous endpoint (from ss1Continuous)
#'   \item Sample size for binary endpoint (from ss1BinaryApprox with Test method)
#' }
#' Two initial values are calculated:
#' \itemize{
#'   \item n2_0: Based on target power 1 - β
#'   \item n2_1: Based on adjusted power 1 - √(1-β) to provide a bracket
#' }
#'
#' \strong{Step 2-4:} Iteratively refine using linear extrapolation until convergence:
#' \deqn{n_2^{new} = \frac{n_2^{(0)}(power^{(1)} - (1-\beta)) - n_2^{(1)}(power^{(0)} - (1-\beta))}
#'       {power^{(1)} - power^{(0)}}}
#'
#' The algorithm converges when \eqn{n_2^{(1)} = n_2^{(0)}}.
#'
#' **Fisher's Exact Test:**
#' When Test = "Fisher", Monte Carlo simulation is used in the power calculation.
#' The linear extrapolation algorithm is still appropriate because the overall
#' co-primary power (considering both continuous and binary endpoints) tends to
#' be smoother than the binary endpoint power alone.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination in
#' clinical trials with multiple co-primary endpoints including mixed continuous
#' and binary variables. \emph{Biometrical Journal}, 54(5), 716-729.
#'
#' Hamasaki, T., Sugimoto, T., Evans, S. R., & Sozu, T. (2013). Sample size
#' determination for clinical trials with co-primary outcomes: exponential event
#' times. \emph{Pharmaceutical Statistics}, 12(1), 28-34.
#'
#' @examples
#' # Example 1: Based on PREMIER study (Table 2 in Sozu et al. 2012)
#' # mTSS (continuous) and ACR50 (binary) with ρ = 0.5
#' ss2Mixed(
#'   delta = 4.4,
#'   sd = 19.0,
#'   p1 = 0.59,
#'   p2 = 0.46,
#'   rho = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   Test = "AN"
#' )
#'
#' # Example 2: With continuity correction
#' ss2Mixed(
#'   delta = 5.0,
#'   sd = 20.0,
#'   p1 = 0.65,
#'   p2 = 0.50,
#'   rho = 0.3,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = "ANc"
#' )
#'
#' # Example 3: Arcsine transformation
#' ss2Mixed(
#'   delta = 4.5,
#'   sd = 18.0,
#'   p1 = 0.70,
#'   p2 = 0.55,
#'   rho = 0.4,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = "AS"
#' )
#'
#' # Example 4: Unequal allocation
#' ss2Mixed(
#'   delta = 0.3,
#'   sd = 1.0,
#'   p1 = 0.6,
#'   p2 = 0.4,
#'   rho = 0.5,
#'   r = 2,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = "AS"
#' )
#'
#' \donttest{
#' # Example 5: Fisher's exact test (computationally intensive)
#' set.seed(12345)
#' ss2Mixed(
#'   delta = 0.5,
#'   sd = 1.0,
#'   p1 = 0.4,
#'   p2 = 0.2,
#'   rho = 0.3,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = "Fisher",
#'   nMC = 10000
#' )
#' }
#'
#' @export
ss2Mixed <- function(delta, sd, p1, p2, rho, r, alpha, beta, Test, nMC = 10000) {

  # Input validation
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

  # Step 1: Initialize sample sizes using single endpoint formulas
  # Calculate sample size for continuous endpoint
  n2_cont <- ss1Continuous(delta, sd, r, alpha, beta)[["n2"]]

  # Calculate sample size for binary endpoint using the same Test method
  n2_bin <- ss1BinaryApprox(p1, p2, r, alpha, beta, Test = Test)[["n2"]]

  # Use the maximum as initial value for target power
  n2_0 <- max(n2_cont, n2_bin)

  # For the second initial value, use adjusted target power
  # This provides a bracket for the linear extrapolation
  beta_adj <- 1 - (1 - beta) ^ (1/2)

  n2_cont_adj <- ss1Continuous(delta, sd, r, alpha, beta_adj)[["n2"]]
  n2_bin_adj <- ss1BinaryApprox(p1, p2, r, alpha, beta_adj, Test = Test)[["n2"]]
  n2_1 <- max(n2_cont_adj, n2_bin_adj)

  # Step 2-4: Iterative refinement using linear extrapolation
  while (n2_1 - n2_0 != 0) {

    # Calculate sample sizes for group 1
    n1_0 <- ceiling(r * n2_0)
    n1_1 <- ceiling(r * n2_1)

    # Calculate power at two candidate sample sizes
    power_n2_0 <- power2Mixed(n1_0, n2_0, delta, sd, p1, p2, rho, alpha, Test, nMC)[["powerCoprimary"]]
    power_n2_1 <- power2Mixed(n1_1, n2_1, delta, sd, p1, p2, rho, alpha, Test, nMC)[["powerCoprimary"]]

    # Linear extrapolation to find sample size that achieves target power
    # This solves: power = (1 - beta) for n2
    n2_updated <- (n2_0 * (power_n2_1 - (1 - beta)) - n2_1 * (power_n2_0 - (1 - beta))) /
      (power_n2_1 - power_n2_0)

    # Update values for next iteration
    n2_1 <- n2_0
    n2_0 <- ceiling(n2_updated)
  }

  # Final sample sizes
  n2 <- n2_0
  n1 <- ceiling(r * n2)
  n <- n1 + n2

  # Return result as a data frame
  result <- data.frame(
    delta = delta,
    sd = sd,
    p1 = p1,
    p2 = p2,
    rho = rho,
    r = r,
    alpha = alpha,
    beta = beta,
    Test = Test,
    n1 = n1,
    n2 = n2,
    n = n
  )

  return(result)
}
