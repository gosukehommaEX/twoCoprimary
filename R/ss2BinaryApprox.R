#' Sample Size Calculation for Two Co-Primary Binary Endpoints (Approximate)
#'
#' Calculates the required sample size for a two-arm superiority trial with two
#' co-primary binary endpoints using the linear extrapolation approach described
#' in Hamasaki et al. (2013).
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
#' This function uses the linear extrapolation approach to find the required
#' sample size iteratively:
#'
#' \strong{Step 1:} Initialize with sample sizes from single endpoint formulas
#'
#' \strong{Step 2-4:} Iteratively refine using linear extrapolation until convergence:
#' \deqn{n_2^{new} = \frac{n_2^{(0)}(power^{(1)} - (1-\beta)) - n_2^{(1)}(power^{(0)} - (1-\beta))}
#'       {power^{(1)} - power^{(0)}}}
#'
#' The algorithm converges when \eqn{n_2^{(1)} = n_2^{(0)}}.
#'
#' @references
#' Hamasaki, T., Sugimoto, T., Evans, S. R., & Sozu, T. (2013). Sample size
#' determination for clinical trials with co-primary outcomes: exponential event
#' times. \emph{Pharmaceutical Statistics}, 12(1), 28-34.
#'
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
#' @export
ss2BinaryApprox <- function(p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test) {

  # Step 1: Initialize sample sizes using single endpoint formulas
  # Use the maximum sample size from both endpoints
  n2_0 <- max(ss1BinaryApprox(c(p11, p12), c(p21, p22), r, alpha, beta)[["n2"]])

  # For the second initial value, use adjusted target power
  # This provides a bracket for the linear extrapolation
  n2_1 <- max(ss1BinaryApprox(c(p11, p12), c(p21, p22), r, alpha, 1 - (1 - beta) ^ (1/2))[["n2"]])

  # Step 2-4: Iterative refinement using linear extrapolation
  while (n2_1 - n2_0 != 0) {

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
