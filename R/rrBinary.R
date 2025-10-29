#' Rejection Region for Two-Arm Trials with a Single Binary Endpoint
#'
#' Calculates the rejection region for two-arm trials with a single binary endpoint
#' using various exact statistical tests, as described in Homma and Yoshida (2025).
#'
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
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
#' RR <- rrBinary(n1, n2, alpha, Test = 'Chisq')
#' print(dim(RR))  # Should be (6, 6)
#'
#' # Fisher exact test
#' RR_fisher <- rrBinary(n1 = 10, n2 = 10, alpha = 0.025, Test = 'Fisher')
#'
#' \donttest{
#' # More computationally intensive: Boschloo test
#' n1 <- 20
#' n2 <- 10
#' alpha <- 0.025
#' RR <- rrBinary(n1, n2, alpha, Test = 'Boschloo')
#' print(RR)
#' }
#'
#' @export
#' @import fpCompare
#' @importFrom stats pnorm dbinom phyper dhyper
rrBinary <- function(n1, n2, alpha, Test) {

  if ((Test == 'Chisq') | (Test == 'Z-pool')) {

    # Calculate test statistics for chi-squared test over all combinations
    # of i (0,...,n1) and j (0,...,n2)
    Z.ij <- '/'(
      outer(n2 * (0:n1), n1 * (0:n2), '-') / (n1 * n2),
      sqrt(outer(0:n1, 0:n2, '+') / (n1 * n2) * (1 - outer(0:n1, 0:n2, '+') / (n1 + n2)))
    )

    # Replace NaN values (when denominator is 0) with 0
    Z.ij[is.na(Z.ij)] <- 0

    if (Test == 'Chisq') {

      # Calculate p-values for chi-squared test using standard normal approximation
      p.val <- pnorm(Z.ij, lower.tail = FALSE)

    } else if (Test == 'Z-pool') {

      # Extract positive test statistics (negative/zero values cannot be significant)
      Z.ij.posi <- Z.ij[Z.ij %>>% 0]
      order.Z.ij.posi <- order(Z.ij.posi, decreasing = TRUE)

      # Get indices for positive test statistics
      i <- (row(Z.ij)[Z.ij %>>% 0] - 1)[order.Z.ij.posi]
      j <- (col(Z.ij)[Z.ij %>>% 0] - 1)[order.Z.ij.posi]

      # Calculate P_H0(X1 = i, X2 = j | theta) for theta in [0, 1]
      # This is the probability under the null hypothesis as a function of theta
      uniq.i <- sort(unique(i))
      dbinom.i <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.i, n1, theta))
      uniq.j <- sort(unique(j))
      dbinom.j <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.j, n2, theta))

      # Joint probability for each (i, j) pair across all theta values
      P_H0 <- dbinom.i[match(i, uniq.i), ] * dbinom.j[match(j, uniq.j), ]

      # Calculate p-values by maximizing over theta (nuisance parameter)
      p.ij <- apply(apply(P_H0, 2, cumsum), 1, max)

      # Assign p-values to the matrix (initialize with 1 for non-significant outcomes)
      p.val <- 1 ^ Z.ij
      p.val[cbind(i + 1, j + 1)] <- p.ij
    }

  } else if (Test == 'Fisher-midP') {

    # Calculate p-values for Fisher mid-p test
    # Mid-p = P(X > x) + 0.5 * P(X = x)
    p.val <- outer(0:n1, 0:n2, function(i, j) {
      phyper(i, n1, n2, i + j, lower.tail = FALSE) + 0.5 * dhyper(i, n1, n2, i + j)
    })

  } else if ((Test == 'Fisher') | (Test == 'Boschloo')) {

    # Calculate p-values for Fisher exact test over all combinations
    p.fisher <- outer(0:n1, 0:n2, function(i, j) {
      phyper(i - 1, n1, n2, i + j, lower.tail = FALSE)
    })

    if (Test == 'Fisher') {

      # Use Fisher p-values directly
      p.val <- p.fisher

    } else if (Test == 'Boschloo') {

      # Boschloo test: maximize Fisher p-values over nuisance parameter

      # Extract p-values where p1 - p2 > 0 (only these can be significant)
      p.max.boschloo <- min(p.fisher[(outer(n2 * (0:n1), n1 * (0:n2), '-') / (n1 * n2)) %<=% 0])
      p.fisher.posi <- p.fisher[p.fisher %<<% p.max.boschloo]
      order.p.fisher.posi <- order(p.fisher.posi, decreasing = FALSE)

      # Get indices for positive outcomes
      i <- (row(p.fisher)[p.fisher %<<% p.max.boschloo] - 1)[order.p.fisher.posi]
      j <- (col(p.fisher)[p.fisher %<<% p.max.boschloo] - 1)[order.p.fisher.posi]

      # Calculate P_H0(X1 = i, X2 = j | theta) for theta in [0, 1]
      uniq.i <- sort(unique(i))
      dbinom.i <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.i, n1, theta))
      uniq.j <- sort(unique(j))
      dbinom.j <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.j, n2, theta))

      # Joint probability for each (i, j) pair
      P_H0 <- dbinom.i[match(i, uniq.i) - min(match(i, uniq.i)) + 1, ] *
        dbinom.j[match(j, uniq.j), ]

      # Calculate p-values by maximizing over theta
      p.ij <- apply(apply(P_H0, 2, cumsum), 1, max)

      # Assign p-values to the matrix
      p.val <- 1 ^ p.fisher
      p.val[cbind(i + 1, j + 1)] <- p.ij
    }
  }

  # Define rejection region: reject if p-value < alpha
  RR <- (p.val %<<% alpha)

  return(RR)
}
