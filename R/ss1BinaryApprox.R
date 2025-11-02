#' Sample Size Calculation for a Single Binary Endpoint
#'
#' Calculates the required sample size for a two-arm superiority trial with a
#' single binary endpoint using various statistical testing methods.
#'
#' @param p1 True probability of responders in group 1 (0 < p1 < 1)
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
#' - κ = 1/r = n2/n1: inverse allocation ratio
#' - π1, π2: response probabilities
#' - θ1 = 1 - π1, θ2 = 1 - π2: non-response probabilities
#' - δ = π1 - π2: treatment effect
#'
#' **AN (Asymptotic Normal) Method:**
#' Uses the standard normal approximation with pooled variance under H0:
#' \deqn{n_2 = \left\lceil \frac{(1 + \kappa)}{(\pi_1 - \pi_2)^2}
#'       \left(z_{1-\alpha} \sqrt{\bar{\pi}(1-\bar{\pi})} +
#'       z_{1-\beta} \sqrt{\kappa\pi_1\theta_1 + \pi_2\theta_2}\right)^2 / \kappa \right\rceil}
#' where \eqn{\bar{\pi} = (r\pi_1 + \pi_2)/(1 + r)} is the pooled proportion.
#'
#' **ANc Method:**
#' Adds continuity correction to the AN method. Uses iterative calculation because
#' the correction term depends on sample size. Converges when the difference between
#' successive iterations is ≤ 1.
#'
#' **AS (Arcsine) Method:**
#' Uses the variance-stabilizing arcsine transformation:
#' \deqn{n_2 = \left\lceil \frac{(z_{1-\alpha} + z_{1-\beta})^2}{4(\sin^{-1}\sqrt{\pi_1} - \sin^{-1}\sqrt{\pi_2})^2} \times \frac{1 + \kappa}{\kappa} \right\rceil}
#'
#' **ASc Method:**
#' Applies continuity correction to the arcsine method. Uses iterative procedure
#' with convergence criterion.
#'
#' **Fisher Method:**
#' Fisher's exact test does not have a closed-form sample size formula. This method:
#' 1. Starts with the AN method's sample size as initial value
#' 2. Incrementally increases n2 by 1
#' 3. Calculates exact power using hypergeometric distribution
#' 4. Stops when power ≥ 1 - β
#'
#' Note: Due to the saw-tooth nature of exact power (power does not increase
#' monotonically with sample size), linear extrapolation is inappropriate for
#' Fisher's method. The incremental approach ensures the minimum sample size
#' that achieves the target power.
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
#' # Compare all methods
#' methods <- c("AN", "ANc", "AS", "ASc", "Fisher")
#' do.call(rbind, lapply(methods, function(m) {
#'   ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = 1, alpha = 0.025, beta = 0.1, Test = m)
#' }))
#'
#' @export
#' @importFrom stats qnorm dbinom pbinom
ss1BinaryApprox <- function(p1, p2, r, alpha, beta, Test = "AN") {

  # Input validation
  if (any(p1 <= 0) || any(p1 >= 1)) {
    stop("p1 must be in (0, 1)")
  }
  if (any(p2 <= 0) || any(p2 >= 1)) {
    stop("p2 must be in (0, 1)")
  }
  if (length(p1) != 1 || length(p2) != 1) {
    stop("p1 and p2 must be scalar values. For co-primary endpoints, use ss2BinaryApprox().")
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

  # Calculate basic quantities
  kappa <- 1 / r  # κ = n2/n1
  theta1 <- 1 - p1
  theta2 <- 1 - p2
  delta <- p1 - p2

  # Standard normal quantiles
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(1 - beta)

  if (Test == "AN") {
    # Asymptotic Normal method without continuity correction
    # Pooled proportion under H0
    p_pooled <- (r * p1 + p2) / (1 + r)

    # Variance under H0
    v0 <- sqrt(p_pooled * (1 - p_pooled))

    # Variance under H1
    v1 <- sqrt((p1 * theta1 / r + p2 * theta2) / (1 + 1 / r))

    # Sample size for group 2
    n2 <- ceiling((1 + 1 / r) / delta ^ 2 * (z_alpha * v0 + z_beta * v1) ^ 2)

  } else if (Test == "ANc") {
    # Asymptotic Normal method with continuity correction
    # This requires iterative calculation because the correction depends on n

    # Initial estimate (without correction)
    p_pooled <- (r * p1 + p2) / (1 + r)
    v0 <- sqrt(p_pooled * (1 - p_pooled))
    v1 <- sqrt((p1 * theta1 / r + p2 * theta2) / (1 + 1 / r))
    n2_initial <- ceiling((1 + 1 / r) / delta ^ 2 * (z_alpha * v0 + z_beta * v1) ^ 2)

    # Iterative refinement using while loop
    n2 <- n2_initial
    n2_prev <- n2 + 10  # Initialize to force at least one iteration
    max_iter <- 100
    iter <- 0

    while (abs(n2 - n2_prev) > 1 && iter < max_iter) {
      iter <- iter + 1
      n2_prev <- n2

      # Continuity correction term
      cc <- (1 + kappa) / (2 * kappa * n2)

      # Adjusted effect size
      delta_adj <- delta - cc

      if (delta_adj <= 0) {
        # Correction is too large, use previous value
        n2 <- n2_prev
        break
      }

      # Recalculate n2
      n2 <- ceiling((1 + 1 / r) / delta_adj ^ 2 * (z_alpha * v0 + z_beta * v1) ^ 2)
    }

  } else if (Test == "AS") {
    # Arcsine transformation without continuity correction
    # The arcsine transformation stabilizes the variance

    # Transformed difference
    delta_as <- asin(sqrt(p1)) - asin(sqrt(p2))

    # Sample size (variance is approximately 1/(4n) for arcsine-transformed proportions)
    n2 <- ceiling((z_alpha + z_beta) ^ 2 / (4 * delta_as ^ 2) * (1 + kappa) / kappa)

  } else if (Test == "ASc") {
    # Arcsine transformation with continuity correction

    # Initial estimate
    delta_as <- asin(sqrt(p1)) - asin(sqrt(p2))
    n2_initial <- ceiling((z_alpha + z_beta) ^ 2 / (4 * delta_as ^ 2) * (1 + kappa) / kappa)

    # Iterative refinement using while loop
    n2 <- n2_initial
    n2_prev <- n2 + 10  # Initialize to force at least one iteration
    max_iter <- 100
    iter <- 0

    while (abs(n2 - n2_prev) > 1 && iter < max_iter) {
      iter <- iter + 1
      n2_prev <- n2

      n1_temp <- ceiling(r * n2)

      # Continuity correction terms
      c1 <- -1 / (2 * n1_temp)
      c2 <- 1 / (2 * n2)

      # Corrected values
      p1_c <- p1 + c1
      p2_c <- p2 + c2

      # Check validity
      if (p1_c <= 0 || p1_c >= 1 || p2_c <= 0 || p2_c >= 1) {
        # Correction leads to invalid values, use previous value
        n2 <- n2_prev
        break
      }

      # Adjusted transformed difference
      delta_as_adj <- asin(sqrt(p1_c)) - asin(sqrt(p2_c))

      # Recalculate n2
      n2 <- ceiling((z_alpha + z_beta) ^ 2 / (4 * delta_as_adj ^ 2) * (1 + kappa) / kappa)
    }

  } else {  # Test == "Fisher"
    # Fisher's exact test requires iterative calculation
    # Start with AN method as initial value
    p_pooled <- (r * p1 + p2) / (1 + r)
    v0 <- sqrt(p_pooled * (1 - p_pooled))
    v1 <- sqrt((p1 * theta1 / r + p2 * theta2) / (1 + 1 / r))
    n2_initial <- ceiling((1 + 1 / r) / delta ^ 2 * (z_alpha * v0 + z_beta * v1) ^ 2)

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
  }

  # Calculate n1 and total n
  n1 <- ceiling(r * n2)
  N <- n1 + n2

  # Return result as a data frame
  result <- data.frame(p1, p2, r, alpha, beta, Test, n1, n2, N)
  return(result)
}
