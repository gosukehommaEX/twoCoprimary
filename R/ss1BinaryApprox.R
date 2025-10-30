#' Sample Size Calculation for a Single Binary Endpoint (Approximate)
#'
#' Calculates the required sample size for a two-arm superiority trial with a
#' single binary endpoint using the asymptotic normal approximation.
#'
#' @param p1 True probability of responders in group 1 (0 < p1 < 1)
#' @param p2 True probability of responders in group 2 (0 < p2 < 1)
#' @param r Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#'
#' @return A data frame with the following columns:
#'   \item{p1}{Probability of responders in group 1}
#'   \item{p2}{Probability of responders in group 2}
#'   \item{r}{Allocation ratio}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{n}{Total sample size (n1 + n2)}
#'
#' @details
#' The required sample size is calculated using the standard asymptotic normal
#' approximation formula for binary endpoints:
#' \deqn{n_2 = \left\lceil \frac{(1 + 1/r)}{(p_1 - p_2)^2}
#'       \left(z_\alpha \sqrt{\bar{p}(1-\bar{p})} +
#'       z_\beta \sqrt{\frac{p_1(1-p_1)}{r} + p_2(1-p_2)} / (1 + 1/r)\right)^2 \right\rceil}
#' where \eqn{\bar{p} = (r \times p_1 + p_2)/(1 + r)} is the pooled proportion,
#' and \eqn{z_\alpha}, \eqn{z_\beta} are standard normal quantiles.
#'
#' This function can handle vector inputs for p1 and p2, returning the maximum
#' required sample size across all endpoint combinations.
#'
#' @examples
#' # Balanced design with 1:1 allocation
#' ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = 1, alpha = 0.025, beta = 0.1)
#'
#' # Unbalanced design with 2:1 allocation
#' ss1BinaryApprox(p1 = 0.5, p2 = 0.3, r = 2, alpha = 0.025, beta = 0.2)
#'
#' # Vector inputs for multiple endpoints
#' ss1BinaryApprox(p1 = c(0.5, 0.4), p2 = c(0.3, 0.2), r = 2, alpha = 0.025, beta = 0.1)
#'
#' @export
#' @importFrom stats qnorm
ss1BinaryApprox <- function(p1, p2, r, alpha, beta) {

  # Calculate pooled proportion
  p <- (r * p1 + p2) / (1 + r)

  # Calculate treatment effect
  delta <- p1 - p2

  # Calculate the required sample size for group 2 using asymptotic normal approximation
  n2 <- ceiling(
    (1 + 1 / r) / (delta ^ 2) *
      (qnorm(1 - alpha) * sqrt(p * (1 - p)) +
         qnorm(1 - beta) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))) ^ 2
  )

  # Calculate the required sample size for group 1
  n1 <- ceiling(r * n2)

  # Calculate total sample size
  n <- n1 + n2

  # Return result as a data frame
  result <- data.frame(p1, p2, r, alpha, beta, n1, n2, n)
  return(result)
}
