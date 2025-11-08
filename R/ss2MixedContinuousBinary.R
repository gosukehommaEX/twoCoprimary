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
#'   * `"AN"`: Asymptotic normal method without continuity correction
#'   * `"ANc"`: Asymptotic normal method with continuity correction
#'   * `"AS"`: Arcsine method without continuity correction
#'   * `"ASc"`: Arcsine method with continuity correction
#'   * `"Fisher"`: Fisher's exact test (uses sequential search)
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
#'   \item{nMC}{Number of Monte Carlo replications (NA if Test != "Fisher")}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{N}{Total sample size (n1 + n2)}
#'
#' @details
#' This function implements the sample size calculation for mixed continuous-binary
#' co-primary endpoints following the methodology in Sozu et al. (2012).
#'
#' The sequential search algorithm (Homma and Yoshida 2025, Algorithm 1) is used
#' for all testing methods:
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
#' \strong{Biserial Correlation:}
#' The biserial correlation rho represents the correlation between the latent continuous
#' variable underlying the binary endpoint and the observed continuous endpoint. This is
#' not the same as the point-biserial correlation observed in the data.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination in
#' clinical trials with multiple co-primary endpoints including mixed continuous
#' and binary variables. \emph{Biometrical Journal}, 54(5), 716-729.
#'
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
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

  # Step 1: Calculate initial sample size
  # Use AN for binary endpoint initial value for consistency
  n2_initial <- max(
    ss1Continuous(delta, sd, r, alpha, beta)[["n2"]],
    ss1BinaryApprox(p1, p2, r, alpha, beta, Test = "AN")[["n2"]]
  )

  # Step 2: Sequential search to find minimum sample size
  result_ss <- .ss_sequential_search(
    initial_n2 = n2_initial,
    r = r,
    target_power = 1 - beta,
    power_fun = power2MixedContinuousBinary,
    delta = delta, sd = sd, p1 = p1, p2 = p2,
    rho = rho, alpha = alpha, Test = Test, nMC = nMC
  )

  n1 <- result_ss$n1
  n2 <- result_ss$n2
  N <- result_ss$N

  # Set nMC to NA for asymptotic methods
  if (Test != "Fisher") {
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
