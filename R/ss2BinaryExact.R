#' Exact Sample Size Calculation for Two Co-Primary Binary Endpoints
#'
#' Calculates the required sample size for a two-arm superiority trial with two
#' co-primary binary endpoints using exact methods, as described in Homma and
#' Yoshida (2025).
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
#'     \item \code{"Chisq"}: One-sided Pearson chi-squared test
#'     \item \code{"Fisher"}: Fisher exact test
#'     \item \code{"Fisher-midP"}: Fisher mid-p test
#'     \item \code{"Z-pool"}: Z-pooled exact unconditional test
#'     \item \code{"Boschloo"}: Boschloo exact unconditional test
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
#' This function implements the exact sample size calculation algorithm described
#' in Appendix C of Homma and Yoshida (2025):
#'
#' \strong{Step 1:} Initialize with sample size from approximate method (AN)
#'
#' \strong{Step 2:} Calculate exact power at current sample size
#'
#' \strong{Step 3:} Adjust sample size:
#' \itemize{
#'   \item If power >= 1 - beta: decrease n2 until power drops below target, then add 1
#'   \item If power < 1 - beta: increase n2 until power reaches target
#' }
#'
#' \strong{Step 4:} Return final sample size
#'
#' The required sample size is calculated as equation (10):
#' \deqn{N_2 = \arg\min_{N_2 \in \mathbb{Z}} \{power_A(\theta) \geq 1 - \beta\}}
#'
#' @references
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
#'
#' @examples
#' # Exact sample size calculation using Boschloo test
#' ss2BinaryExact(
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0.7,
#'   rho2 = 0.7,
#'   r = 2,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = 'Boschloo'
#' )
#'
#' # Balanced design with Fisher exact test
#' ss2BinaryExact(
#'   p11 = 0.6,
#'   p12 = 0.5,
#'   p21 = 0.4,
#'   p22 = 0.3,
#'   rho1 = 0.5,
#'   rho2 = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   Test = 'Fisher'
#' )
#'
#' \donttest{
#' # Larger sample sizes (computationally intensive)
#' ss2BinaryExact(
#'   p11 = 0.45,
#'   p12 = 0.40,
#'   p21 = 0.30,
#'   p22 = 0.25,
#'   rho1 = 0.6,
#'   rho2 = 0.6,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.1,
#'   Test = 'Chisq'
#' )
#' }
#'
#' @export
#' @import fpCompare
ss2BinaryExact <- function(p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test) {

  # Step 1: Initialize sample size using approximate method (AN)
  # This provides a good starting point for the exact calculation
  n2 <- ss2BinaryApprox(p11, p12, p21, p22, rho1, rho2, r, alpha, beta, "AN")[["n2"]]
  n1 <- ceiling(r * n2)

  # Step 2: Calculate exact power at current sample size
  power <- power2BinaryExact(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test)[["powerCoprimary"]]

  # Step 3: Adjust sample size to achieve target power

  if (power %>=% (1 - beta)) {

    # Step 3-1: If power is too high, decrease n2 until power drops below target
    while (power %>=% (1 - beta)) {
      n2 <- n2 - 1
      n1 <- ceiling(r * n2)
      power <- power2BinaryExact(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test)[["powerCoprimary"]]
    }

    # Add 1 back to get the minimum sample size that achieves target power
    n2 <- n2 + 1

  } else {

    # Step 3-2: If power is too low, increase n2 until power reaches target
    while (power %<<% (1 - beta)) {
      n2 <- n2 + 1
      n1 <- ceiling(r * n2)
      power <- power2BinaryExact(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test)[["powerCoprimary"]]
    }
  }

  # Step 4: Determine final sample sizes
  n1 <- ceiling(r * n2)
  N <- n1 + n2

  # Return result as a data frame
  result <- data.frame(
    p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test,
    n1, n2, N
  )
  return(result)
}
