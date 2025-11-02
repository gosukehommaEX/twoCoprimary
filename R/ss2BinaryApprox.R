#' Sample Size Calculation for Two Co-Primary Binary Endpoints (Approximate)
#'
#' Calculates the required sample size for a two-arm superiority trial with two
#' co-primary binary endpoints using the linear extrapolation approach.
#'
#' @param p11 True probability of responders in group 1 for the first outcome (0 < p11 < 1)
#' @param p12 True probability of responders in group 1 for the second outcome (0 < p12 < 1)
#' @param p21 True probability of responders in group 2 for the first outcome (0 < p21 < 1)
#' @param p22 True probability of responders in group 2 for the second outcome (0 < p22 < 1)
#' @param rho1 Correlation between the two outcomes for group 1
#' @param rho2 Correlation between the two outcomes for group 2
#' @param r Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#' @param Test Statistical testing method. One of:
#'   \itemize{
#'     \item \code{"AN"}: Asymptotic normal method without continuity correction
#'     \item \code{"ANc"}: Asymptotic normal method with continuity correction
#'     \item \code{"AS"}: Arcsine method without continuity correction
#'     \item \code{"ASc"}: Arcsine method with continuity correction
#'   }
#'
#' @return A data frame with the following columns:
#'   \item{p11, p12, p21, p22}{Response probabilities}
#'   \item{rho1, rho2}{Correlations}
#'   \item{r}{Allocation ratio}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{Test}{Testing method used}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{N}{Total sample size (n1 + n2)}
#'
#' @details
#' This function uses the linear extrapolation approach to find the required
#' sample size iteratively:
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas.
#' The initial values use the same testing method specified in the Test argument
#' to ensure consistency. Two initial values are calculated:
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
#' \strong{Note on Fisher's Exact Test:}
#' Fisher's exact test is not supported in this function because the saw-tooth
#' nature of exact power makes linear extrapolation inappropriate. For Fisher's
#' exact test with co-primary endpoints, use \code{\link{ss2BinaryExact}} instead.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination in
#' clinical trials with multiple co-primary binary endpoints. \emph{Statistics in
#' Medicine}, 29(21), 2169-2179.
#'
#' @examples
#' # Sample size calculation using asymptotic normal method
#' ss2BinaryApprox(
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0.7,
#'   rho2 = 0.7,
#'   r = 2,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = 'AN'
#' )
#'
#' # Balanced design with arcsine method
#' ss2BinaryApprox(
#'   p11 = 0.6,
#'   p12 = 0.5,
#'   p21 = 0.4,
#'   p22 = 0.3,
#'   rho1 = 0.5,
#'   rho2 = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   Test = 'AS'
#' )
#'
#' # With continuity correction
#' ss2BinaryApprox(
#'   p11 = 0.55,
#'   p12 = 0.45,
#'   p21 = 0.35,
#'   p22 = 0.25,
#'   rho1 = 0.6,
#'   rho2 = 0.6,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = 'ANc'
#' )
#'
#' @export
ss2BinaryApprox <- function(p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test) {

  # Input validation
  if (length(p11) != 1 || length(p12) != 1 || length(p21) != 1 ||
      length(p22) != 1 || length(rho1) != 1 || length(rho2) != 1 ||
      length(r) != 1 || length(alpha) != 1 || length(beta) != 1) {
    stop("All parameters must be scalar values")
  }
  if (p11 <= 0 || p11 >= 1 || p12 <= 0 || p12 >= 1 ||
      p21 <= 0 || p21 >= 1 || p22 <= 0 || p22 >= 1) {
    stop("All probabilities must be in (0, 1)")
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
  if (!Test %in% c("AN", "ANc", "AS", "ASc")) {
    stop("Test must be one of: AN, ANc, AS, ASc. For Fisher's exact test, use ss2BinaryExact().")
  }

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

  # Step 1: Initialize sample sizes using single endpoint formulas
  # IMPORTANT: Use the same Test method for initial values to ensure consistency
  # Calculate sample size for each endpoint separately, then take the maximum
  n2_0 <- max(
    ss1BinaryApprox(p11, p21, r, alpha, beta, Test = Test)[["n2"]],
    ss1BinaryApprox(p12, p22, r, alpha, beta, Test = Test)[["n2"]]
  )

  # For the second initial value, use adjusted target power
  # This provides a bracket for the linear extrapolation
  # The adjusted power 1 - sqrt(1-beta) is typically higher than 1 - beta
  n2_1 <- max(
    ss1BinaryApprox(p11, p21, r, alpha, 1 - (1 - beta) ^ (1/2), Test = Test)[["n2"]],
    ss1BinaryApprox(p12, p22, r, alpha, 1 - (1 - beta) ^ (1/2), Test = Test)[["n2"]]
  )

  # Step 2-4: Iterative refinement using linear extrapolation
  while (n2_1 - n2_0 != 0) {

    # Calculate sample sizes for group 1
    n1_0 <- ceiling(r * n2_0)
    n1_1 <- ceiling(r * n2_1)

    # Calculate power at two candidate sample sizes
    power_n2_0 <- power2BinaryApprox(
      n1_0, n2_0, p11, p12, p21, p22, rho1, rho2, alpha, Test
    )[["powerCoprimary"]]
    power_n2_1 <- power2BinaryApprox(
      n1_1, n2_1, p11, p12, p21, p22, rho1, rho2, alpha, Test
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

  # Return result as a data frame
  result <- data.frame(
    p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test,
    n1, n2, N
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
