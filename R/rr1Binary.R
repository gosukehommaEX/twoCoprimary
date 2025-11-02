#' Rejection Region for Two-Arm Trials with a Single Binary Endpoint
#'
#' Calculates the rejection region for two-arm trials with a single binary endpoint
#' using various exact statistical tests, as described in Homma and Yoshida (2025).
#'
#' @param n1 Sample size for group 1 (test group)
#' @param n2 Sample size for group 2 (control group)
#' @param alpha One-sided significance level (typically 0.025)
#' @param Test Type of statistical test. One of:
#'   \itemize{
#'     \item \code{"Chisq"}: One-sided Pearson chi-squared test
#'     \item \code{"Fisher"}: Fisher exact test
#'     \item \code{"Fisher-midP"}: Fisher mid-p test
#'     \item \code{"Z-pool"}: Z-pooled exact unconditional test
#'     \item \code{"Boschloo"}: Boschloo exact unconditional test
#'   }
#'
#' @return A logical matrix of dimensions (n1+1) x (n2+1), where TRUE indicates
#'   rejection of the null hypothesis. Rows correspond to the number of responders
#'   in group 1 (0 to n1), and columns correspond to the number of responders in
#'   group 2 (0 to n2).
#'
#' @details
#' This function computes the rejection region for five different one-sided tests:
#'
#' \strong{Chi-squared test:} Uses the asymptotic normal approximation of the
#' chi-squared statistic.
#'
#' \strong{Fisher exact test:} Uses the hypergeometric distribution to calculate
#' exact p-values conditional on the total number of successes.
#'
#' \strong{Fisher mid-p test:} Modification of Fisher's exact test that adds half
#' the probability of the observed outcome to reduce conservatism.
#'
#' \strong{Z-pooled test:} Exact unconditional test that maximizes p-values over
#' all possible values of the nuisance parameter (common success probability under H0).
#'
#' \strong{Boschloo test:} Exact unconditional test similar to Z-pooled but based
#' on Fisher's exact p-values, maximizing over the nuisance parameter.
#'
#' @references
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
#'
#' @examples
#' # Simple example with small sample sizes
#' n1 <- 5
#' n2 <- 5
#' alpha <- 0.025
#' RR <- rr1Binary(n1, n2, alpha, Test = 'Chisq')
#' print(dim(RR))  # Should be (6, 6)
#'
#' # Fisher exact test
#' RR_fisher <- rr1Binary(n1 = 10, n2 = 10, alpha = 0.025, Test = 'Fisher')
#'
#' \donttest{
#' # More computationally intensive: Boschloo test
#' n1 <- 20
#' n2 <- 10
#' alpha <- 0.025
#' RR <- rr1Binary(n1, n2, alpha, Test = 'Boschloo')
#' print(RR)
#' }
#'
#' @export
#' @import fpCompare
#' @importFrom stats pnorm dbinom phyper dhyper
rr1Binary <- function(n1, n2, alpha, Test) {

  # Input validation
  if (length(n1) != 1 || length(n2) != 1) {
    stop("n1 and n2 must be scalar values")
  }
  if (n1 <= 0 || n1 != round(n1)) {
    stop("n1 must be a positive integer")
  }
  if (n2 <= 0 || n2 != round(n2)) {
    stop("n2 must be a positive integer")
  }
  if (length(alpha) != 1) {
    stop("alpha must be a scalar value")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }
  if (!Test %in% c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")) {
    stop("Test must be one of: Chisq, Fisher, Fisher-midP, Z-pool, Boschloo")
  }

  if ((Test == 'Chisq') | (Test == 'Z-pool')) {

    # Calculate test statistics for chi-squared test over all combinations
    # of i (0,...,n1) and j (0,...,n2)
    Z_ij <- '/'(
      outer(n2 * (0:n1), n1 * (0:n2), '-') / (n1 * n2),
      sqrt(outer(0:n1, 0:n2, '+') / (n1 * n2) * (1 - outer(0:n1, 0:n2, '+') / (n1 + n2)))
    )

    # Replace NaN values (when denominator is 0) with 0
    Z_ij[is.na(Z_ij)] <- 0

    if (Test == 'Chisq') {

      # Calculate p-values for chi-squared test using standard normal approximation
      p_val <- pnorm(Z_ij, lower.tail = FALSE)

    } else if (Test == 'Z-pool') {

      # Extract positive test statistics (negative/zero values cannot be significant)
      Z_ij_posi <- Z_ij[Z_ij %>>% 0]
      order_Z_ij_posi <- order(Z_ij_posi, decreasing = TRUE)

      # Get indices for positive test statistics
      i <- (row(Z_ij)[Z_ij %>>% 0] - 1)[order_Z_ij_posi]
      j <- (col(Z_ij)[Z_ij %>>% 0] - 1)[order_Z_ij_posi]

      # Calculate P_H0(X1 = i, X2 = j | theta) for theta in [0, 1]
      # This is the probability under the null hypothesis as a function of theta
      uniq_i <- sort(unique(i))
      dbinom_i <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq_i, n1, theta))
      uniq_j <- sort(unique(j))
      dbinom_j <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq_j, n2, theta))

      # Joint probability for each (i, j) pair across all theta values
      P_H0 <- dbinom_i[match(i, uniq_i), ] * dbinom_j[match(j, uniq_j), ]

      # Calculate p-values by maximizing over theta (nuisance parameter)
      p_ij <- apply(apply(P_H0, 2, cumsum), 1, max)

      # Assign p-values to the matrix (initialize with 1 for non-significant outcomes)
      p_val <- 1 ^ Z_ij
      p_val[cbind(i + 1, j + 1)] <- p_ij
    }

  } else if (Test == 'Fisher-midP') {

    # Calculate p-values for Fisher mid-p test
    # Mid-p = P(X > x) + 0.5 * P(X = x)
    p_val <- outer(0:n1, 0:n2, function(i, j) {
      phyper(i, n1, n2, i + j, lower.tail = FALSE) + 0.5 * dhyper(i, n1, n2, i + j)
    })

  } else if ((Test == 'Fisher') | (Test == 'Boschloo')) {

    # Calculate p-values for Fisher exact test over all combinations
    p_fisher <- outer(0:n1, 0:n2, function(i, j) {
      phyper(i - 1, n1, n2, i + j, lower.tail = FALSE)
    })

    if (Test == 'Fisher') {

      # Use Fisher p-values directly
      p_val <- p_fisher

    } else if (Test == 'Boschloo') {

      # Boschloo test: maximize Fisher p-values over nuisance parameter

      # Extract p-values where p1 - p2 > 0 (only these can be significant)
      p_max_boschloo <- min(p_fisher[(outer(n2 * (0:n1), n1 * (0:n2), '-') / (n1 * n2)) %<=% 0])
      p_fisher_posi <- p_fisher[p_fisher %<<% p_max_boschloo]
      order_p_fisher_posi <- order(p_fisher_posi, decreasing = FALSE)

      # Get indices for positive outcomes
      i <- (row(p_fisher)[p_fisher %<<% p_max_boschloo] - 1)[order_p_fisher_posi]
      j <- (col(p_fisher)[p_fisher %<<% p_max_boschloo] - 1)[order_p_fisher_posi]

      # Calculate P_H0(X1 = i, X2 = j | theta) for theta in [0, 1]
      uniq_i <- sort(unique(i))
      dbinom_i <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq_i, n1, theta))
      uniq_j <- sort(unique(j))
      dbinom_j <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq_j, n2, theta))

      # Joint probability for each (i, j) pair
      P_H0 <- dbinom_i[match(i, uniq_i) - min(match(i, uniq_i)) + 1, ] *
        dbinom_j[match(j, uniq_j), ]

      # Calculate p-values by maximizing over theta
      p_ij <- apply(apply(P_H0, 2, cumsum), 1, max)

      # Assign p-values to the matrix
      p_val <- 1 ^ p_fisher
      p_val[cbind(i + 1, j + 1)] <- p_ij
    }
  }

  # Define rejection region: reject if p-value < alpha
  RR <- (p_val %<<% alpha)

  return(RR)
}
