#' Sample Size Calculation for a Single Binary Endpoint
#'
#' Calculates the required sample size for a two-arm superiority trial with a
#' single binary endpoint using various statistical testing methods.
#'
#' @param p1 True probability of responders in group 1 (0 < p1 < 1). Must be
#'   greater than \code{p2}, since the design is a superiority trial
#' @param p2 True probability of responders in group 2 (0 < p2 < 1)
#' @param r Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#' @param Test Statistical testing method. One of:
#'   \itemize{
#'     \item \code{"AN"}: Asymptotic normal method without continuity correction (default)
#'     \item \code{"ANc"}: Asymptotic normal method with continuity correction
#'     \item \code{"AS"}: Arcsine transformation without continuity correction
#'     \item \code{"ASc"}: Arcsine transformation with continuity correction
#'     \item \code{"Fisher"}: Fisher's exact test with iterative sample size determination
#'   }
#'
#' @return A data frame with the following columns:
#'   \item{p1}{Probability of responders in group 1}
#'   \item{p2}{Probability of responders in group 2}
#'   \item{r}{Allocation ratio}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{Test}{Testing method used}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{N}{Total sample size (n1 + n2)}
#'
#' @details
#' This function implements sample size calculations for single binary endpoint
#' trials using five different methods.
#'
#' **Important:** This function is designed for a **single binary endpoint**.
#' For co-primary endpoints, use \code{\link{ss2BinaryApprox}} (for approximate
#' methods) or \code{\link{ss2BinaryExact}} (for exact methods).
#'
#' **Notation:**
#' - r = n1/n2: allocation ratio (group 1 to group 2)
#' - kappa = 1/r = n2/n1: inverse allocation ratio
#' - p1, p2: response probabilities
#' - theta1 = 1 - p1, theta2 = 1 - p2: non-response probabilities
#' - delta = p1 - p2: treatment effect
#'
#' All four asymptotic methods return the smallest \eqn{n_2} whose realized
#' group sizes reach the target power under \code{\link{power2BinaryApprox}}
#' with both endpoints set to the same marginal. The closed forms below give the
#' starting value of that search. They are not used as the answer, because they
#' assume that \eqn{n_1} is exactly \eqn{r n_2}, whereas \eqn{n_1} is rounded
#' up to an integer, and because for the two continuity corrected methods the
#' correction depends on the sample size being solved for.
#'
#' **AN (Asymptotic Normal) Method:**
#' Uses the standard normal approximation with pooled variance under H0:
#' \deqn{n_2 = \left\lceil \frac{(1 + \kappa)}{(\pi_1 - \pi_2)^2}
#'       \left(z_{1-\alpha} \sqrt{\bar{\pi}(1-\bar{\pi})} +
#'       z_{1-\beta} \sqrt{\kappa\pi_1\theta_1 + \pi_2\theta_2}\right)^2 / \kappa \right\rceil}
#' where \eqn{\bar{\pi} = (r\pi_1 + \pi_2)/(1 + r)} is the pooled proportion.
#'
#' **ANc Method:**
#' Adds Yates's continuity correction to the AN method.
#'
#' **AS (Arcsine) Method:**
#' Uses the variance-stabilizing arcsine transformation:
#' \deqn{n_2 = \left\lceil \frac{(z_{1-\alpha} + z_{1-\beta})^2 (1 + \kappa)}{4(\sin^{-1}\sqrt{\pi_1} - \sin^{-1}\sqrt{\pi_2})^2} \right\rceil}
#'
#' **ASc Method:**
#' Applies continuity correction to the arcsine method. The correction moves the
#' two arms toward each other, and the variance is adjusted for the corrected
#' proportions following Sozu et al. (2010).
#'
#' **Fisher Method:**
#' Fisher's exact test does not have a closed-form sample size formula. This method:
#' 1. Starts with the AN method's sample size as initial value
#' 2. Incrementally increases n2 by 1
#' 3. Calculates exact power by summing the binomial probabilities over the
#'    rejection region
#' 4. Stops when power is greater than or equal to 1 - beta
#' 5. Steps back down while the target is still met, so that a starting value
#'    that already reached the target does not hide a smaller one
#'
#' Note: Due to the saw-tooth nature of exact power (power does not increase
#' monotonically with sample size), a sequential search approach is used.
#' The incremental approach ensures the minimum sample size that achieves the
#' target power.
#'
#' @references
#' Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination in
#' clinical trials with multiple co-primary binary endpoints. \emph{Statistics in
#' Medicine}, 29(21), 2169-2179.
#'
#' @examples
#' # Balanced design with 1:1 allocation (AN method)
#' ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = 1, alpha = 0.025, beta = 0.1, Test = "AN")
#'
#' # Unbalanced design with 2:1 allocation (ANc method)
#' ss1BinaryApprox(p1 = 0.5, p2 = 0.3, r = 2, alpha = 0.025, beta = 0.2, Test = "ANc")
#'
#' # Arcsine transformation method
#' ss1BinaryApprox(p1 = 0.55, p2 = 0.35, r = 1, alpha = 0.025, beta = 0.1, Test = "AS")
#'
#' # Arcsine with continuity correction
#' ss1BinaryApprox(p1 = 0.65, p2 = 0.45, r = 1, alpha = 0.025, beta = 0.1, Test = "ASc")
#'
#' # Fisher's exact test
#' ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = 2, alpha = 0.025, beta = 0.1, Test = "Fisher")
#'
#' @export
#' @import fpCompare
#' @importFrom stats qnorm dbinom pbinom
ss1BinaryApprox <- function(p1, p2, r, alpha, beta, Test = "AN") {

  # Input validation
  if (length(p1) != 1 || length(p2) != 1 || length(r) != 1 ||
      length(alpha) != 1 || length(beta) != 1) {
    stop("All parameters must be scalar values")
  }
  if (p1 <= 0 || p1 >= 1) {
    stop("p1 must be in (0, 1)")
  }
  if (p2 <= 0 || p2 >= 1) {
    stop("p2 must be in (0, 1)")
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
  if (p1 <= p2) {
    stop("p1 must be greater than p2 for superiority trial")
  }

  # Calculate parameters
  kappa <- 1 / r
  delta <- p1 - p2
  theta1 <- 1 - p1
  theta2 <- 1 - p2
  z_alpha <- qnorm(alpha)
  z_beta <- qnorm(beta)

  if (Test == "AN" | Test == "ANc") {
    # Pooled proportion
    p_pooled <- (r * p1 + p2) / (1 + r)
    v0 <- sqrt(p_pooled * (1 - p_pooled))
    v1 <- sqrt((p1 * theta1 / r + p2 * theta2) / (1 + 1 / r))

    # Closed form starting value (AN method). One subject is the floor: the
    # closed form returns zero when the target power does not exceed the size
    # of the test, and the search below evaluates the power at its start.
    n2 <- max(1, ceiling((1 + 1 / r) / delta ^ 2 *
                           (z_alpha * v0 + z_beta * v1) ^ 2))

  } else if (Test == "AS" | Test == "ASc") {
    # Arcsine transformation
    delta_as <- asin(sqrt(p1)) - asin(sqrt(p2))

    # Closed form starting value (AS method), with the same floor of one
    n2 <- max(1, ceiling((z_alpha + z_beta) ^ 2 / (4 * delta_as ^ 2) *
                           (1 + kappa)))
  }

  if (Test != "Fisher") {
    # The closed forms above assume that n1 is exactly r n2 and, for the two
    # continuity corrected methods, that the correction does not move with the
    # sample size. Neither holds. The group 1 size is rounded up to an integer,
    # so the realized allocation exceeds r whenever r n2 is not an integer, and
    # the correction is a function of the very sample size being solved for.
    # The size returned is therefore the smallest one whose realized (n1, n2)
    # reaches the target power, which is Algorithm 1 of Homma and Yoshida (2025)
    # applied to a single endpoint. The power is taken from
    # power2BinaryApprox with both endpoints set to the same marginal, so that
    # the size and the power of this package always refer to the same formula:
    # for ASc in particular, Sozu et al. (2010) adjust the variance for the
    # corrected proportions, which the closed form above does not do.
    power1_at <- function(m) {
      power2BinaryApprox(n1 = ceiling(r * m), n2 = m,
                         p11 = p1, p12 = p1, p21 = p2, p22 = p2,
                         rho1 = 0, rho2 = 0, alpha = alpha,
                         Test = Test)[["power1"]]
    }

    target_power <- 1 - beta

    while (power1_at(n2) %<<% target_power) {
      n2 <- n2 + 1
    }
    while (n2 > 1 && power1_at(n2 - 1) %>=% target_power) {
      n2 <- n2 - 1
    }

  } else {  # Test == "Fisher"
    # Fisher's exact test requires iterative calculation
    # Start with AN method as initial value
    p_pooled <- (r * p1 + p2) / (1 + r)
    v0 <- sqrt(p_pooled * (1 - p_pooled))
    v1 <- sqrt((p1 * theta1 / r + p2 * theta2) / (1 + 1 / r))
    n2_initial <- max(1, ceiling((1 + 1 / r) / delta ^ 2 *
                                   (z_alpha * v0 + z_beta * v1) ^ 2))

    # Initialize
    n2 <- n2_initial
    n1 <- ceiling(r * n2)
    power <- 0

    # Increment n2 until target power is achieved
    # Note: Due to saw-tooth problem, we increment by 1
    while (power < 1 - beta) {
      # Calculate rejection region using Fisher's exact test
      RR <- rr1Binary(n1, n2, alpha, Test = 'Fisher')

      # Calculate exact power
      # P(reject H0) = sum over all (x1, x2) where we reject
      # For single endpoint: sum_{x1} P(X1=x1) * P(X2 <= critical value | X1=x1)
      power <- sum(dbinom(0:n1, n1, p1) * pbinom(rowSums(RR) - 1, n2, p2))

      # If power is insufficient, increase sample size
      if (power < 1 - beta) {
        n2 <- n2 + 1
        n1 <- ceiling(r * n2)
      }
    }

    # The loop above only moves upward from the approximate starting value, so
    # if that value already reached the target the smaller sizes were never
    # examined. Step back down while the target is still met, which is Step 2a
    # of the sequential search the two endpoint functions use.
    while (n2 > 1) {
      n2_try <- n2 - 1
      n1_try <- ceiling(r * n2_try)
      RR_try <- rr1Binary(n1_try, n2_try, alpha, Test = 'Fisher')
      power_try <- sum(dbinom(0:n1_try, n1_try, p1) *
                         pbinom(rowSums(RR_try) - 1, n2_try, p2))
      if (power_try < 1 - beta) {
        break
      }
      n2 <- n2_try
    }
  }

  # Calculate n1 and total n
  n1 <- ceiling(r * n2)
  N <- n1 + n2

  # Return result as a data frame
  result <- data.frame(p1, p2, r, alpha, beta, Test, n1, n2, N)
  class(result) <- c("twoCoprimary", "data.frame")

  return(result)
}
