#' Sample Size Calculation for Two Co-Primary Binary Endpoints (Asymptotic Approximation)
#'
#' Calculates the required sample size for a two-arm superiority trial with two
#' co-primary binary endpoints using asymptotic normal approximation or arcsine
#' transformation, as described in Sozu et al. (2010).
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
#'   * `"AN"`: Asymptotic normal method without continuity correction
#'   * `"ANc"`: Asymptotic normal method with continuity correction
#'   * `"AS"`: Arcsine method without continuity correction
#'   * `"ASc"`: Arcsine method with continuity correction
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
#' This function uses a sequential search algorithm (Homma and Yoshida 2025,
#' Algorithm 1) to find the minimum sample size:
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas.
#'
#' \strong{Step 2:} Use sequential search:
#' \itemize{
#'   \item Calculate power at initial sample size
#'   \item If power >= target: decrease n2 until power < target, then add 1 back
#'   \item If power < target: increase n2 until power >= target
#' }
#'
#' \strong{Step 3:} Return final sample sizes.
#'
#' The asymptotic normal (AN) and arcsine (AS) methods use normal approximation
#' with or without continuity correction. For small sample sizes or extreme
#' probabilities, consider using exact methods via \code{\link{ss2BinaryExact}}.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination in
#' clinical trials with multiple co-primary binary endpoints. \emph{Statistics
#' in Medicine}, 29(21), 2169-2179.
#'
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
#'
#' @examples
#' # Sample size calculation using asymptotic normal method
#' ss2BinaryApprox(
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0.5,
#'   rho2 = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   Test = 'AN'
#' )
#'
#' # With continuity correction
#' ss2BinaryApprox(
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0.5,
#'   rho2 = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   Test = 'ANc'
#' )
#'
#' # Using arcsine transformation
#' ss2BinaryApprox(
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0.5,
#'   rho2 = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   Test = 'AS'
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
    stop("Test must be one of: AN, ANc, AS, ASc. ",
         "For Fisher's exact test, use ss2BinaryExact().")
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

  # Step 1: Initialize sample size using single endpoint formulas
  # Use the same Test method for initial values to ensure consistency
  n2_initial <- max(
    ss1BinaryApprox(p11, p21, r, alpha, beta, Test = Test)[["n2"]],
    ss1BinaryApprox(p12, p22, r, alpha, beta, Test = Test)[["n2"]]
  )

  # Step 2: Sequential search to find minimum sample size
  result_ss <- .ss_sequential_search(
    initial_n2 = n2_initial,
    r = r,
    target_power = 1 - beta,
    power_fun = power2BinaryApprox,
    p11 = p11, p12 = p12, p21 = p21, p22 = p22,
    rho1 = rho1, rho2 = rho2,
    alpha = alpha, Test = Test
  )

  n1 <- result_ss$n1
  n2 <- result_ss$n2
  N <- result_ss$N

  # Return result as a data frame
  result <- data.frame(
    p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test,
    n1, n2, N
  )
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
