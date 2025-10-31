#' Sample Size Calculation for Two Co-Primary Binary Endpoints (Approximate)
#'
#' Determines the sample size for a two-arm superiority trial with two co-primary
#' binary endpoints using approximate methods, to achieve a specified power at a
#' given significance level.
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
#'   \item{n}{Total sample size (n1 + n2)}
#'
#' @details
#' This function uses linear extrapolation algorithm for efficient sample size
#' determination. The algorithm iteratively refines the sample size estimate:
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas
#'
#' \strong{Step 2-4:} Iteratively refine using linear extrapolation:
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
#' # Example 1: AN method with rho = 0.8 (from Sozu et al. 2010)
#' ss2BinaryApprox(
#'   p11 = 0.7,
#'   p12 = 0.7,
#'   p21 = 0.5,
#'   p22 = 0.5,
#'   rho1 = 0.8,
#'   rho2 = 0.8,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   Test = 'AN'
#' )
#'
#' # Example 2: Balanced design with arcsine method
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
#' # Example 3: With continuity correction
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
  if (!Test %in% c("AN", "ANc", "AS", "ASc")) {
    stop("Test must be one of: AN, ANc, AS, ASc. For Fisher's exact test, use ss2BinaryExact().")
  }

  # Step 1: Initialize sample sizes using single endpoint formulas
  # IMPORTANT: Use the same Test method for initial values to ensure consistency
  # Calculate sample size for each endpoint separately, then take the maximum
  n2_endpoint1 <- ss1BinaryApprox(p11, p21, r, alpha, beta, Test = Test)[["n2"]]
  n2_endpoint2 <- ss1BinaryApprox(p12, p22, r, alpha, beta, Test = Test)[["n2"]]
  n2_0 <- max(n2_endpoint1, n2_endpoint2)

  # For the second initial value, use adjusted target power
  # This provides a bracket for the linear extrapolation
  # The adjusted power 1 - sqrt(1 - beta) is typically higher than 1 - beta
  n2_endpoint1_adj <- ss1BinaryApprox(p11, p21, r, alpha, 1 - (1 - beta) ^ (1/2), Test = Test)[["n2"]]
  n2_endpoint2_adj <- ss1BinaryApprox(p12, p22, r, alpha, 1 - (1 - beta) ^ (1/2), Test = Test)[["n2"]]
  n2_1 <- max(n2_endpoint1_adj, n2_endpoint2_adj)

  # Step 2-4: Iterative refinement using linear extrapolation
  max_iter <- 100
  iter <- 0

  while (n2_1 - n2_0 != 0 && iter < max_iter) {
    iter <- iter + 1

    # Calculate sample sizes for group 1
    n1_0 <- ceiling(r * n2_0)
    n1_1 <- ceiling(r * n2_1)

    # Calculate power at two candidate sample sizes
    power_n2_0 <- power2BinaryApprox(n1_0, n2_0, p11, p12, p21, p22, rho1, rho2,
                                     alpha, Test)[["powerCoprimary"]]
    power_n2_1 <- power2BinaryApprox(n1_1, n2_1, p11, p12, p21, p22, rho1, rho2,
                                     alpha, Test)[["powerCoprimary"]]

    # Linear extrapolation to find sample size that achieves target power
    # This solves: power = (1 - beta) for n2
    n2_updated <- (n2_0 * (power_n2_1 - (1 - beta)) - n2_1 * (power_n2_0 - (1 - beta))) /
      (power_n2_1 - power_n2_0)

    # Update values for next iteration
    n2_1 <- n2_0
    n2_0 <- ceiling(n2_updated)
  }

  if (iter >= max_iter) {
    warning("Maximum iterations reached. Results may not have converged.")
  }

  # Final sample sizes
  n2 <- n2_0
  n1 <- ceiling(r * n2)
  n <- n1 + n2

  # Return result as a data frame
  result <- data.frame(
    p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test,
    n1, n2, n
  )
  return(result)
}
