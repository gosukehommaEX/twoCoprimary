#' Exact Power Calculation for Two Co-Primary Binary Endpoints
#'
#' Calculates the exact power for a two-arm superiority trial with two co-primary
#' binary endpoints using the bivariate binomial distribution, as described in
#' Homma and Yoshida (2025).
#'
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param p11 True probability of responders in group 1 for the first outcome (0 < p11 < 1)
#' @param p12 True probability of responders in group 1 for the second outcome (0 < p12 < 1)
#' @param p21 True probability of responders in group 2 for the first outcome (0 < p21 < 1)
#' @param p22 True probability of responders in group 2 for the second outcome (0 < p22 < 1)
#' @param rho1 Correlation between the two outcomes for group 1
#' @param rho2 Correlation between the two outcomes for group 2
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param Test Statistical testing method. One of:
#'   \itemize{
#'     \item \code{"Chisq"}: One-sided Pearson chi-squared test
#'     \item \code{"Fisher"}: Fisher exact test
#'     \item \code{"Fisher-midP"}: Fisher mid-p test
#'     \item \code{"Z-pool"}: Z-pooled exact unconditional test
#'     \item \code{"Boschloo"}: Boschloo exact unconditional test
#'   }
#'
#' @return A named numeric vector with three elements:
#'   \item{Power1}{Power for the first endpoint alone}
#'   \item{Power2}{Power for the second endpoint alone}
#'   \item{powerCoprimary}{Exact power for both co-primary endpoints}
#'
#' @details
#' This function calculates exact power using equation (9) in Homma and Yoshida (2025):
#' \deqn{power_A(\theta) = \sum_{(a_{1,1},a_{2,1})\in\mathcal{A}_1}
#'       \sum_{(a_{1,2},a_{2,2})\in\mathcal{A}_2} f(a_{1,1}|N_1,p_{1,1}) \times
#'       f(a_{2,1}|N_2,p_{2,1}) \times g(a_{1,2}|a_{1,1},N_1,p_{1,1},p_{1,2},\gamma_1)
#'       \times g(a_{2,2}|a_{2,1},N_2,p_{2,1},p_{2,2},\gamma_2)}
#'
#' where \eqn{\mathcal{A}_k} is the rejection region for endpoint k, and
#' \eqn{(Y_{j,1}, Y_{j,2}) \sim BiBin(N_j, p_{j,1}, p_{j,2}, \gamma_j)} follows
#' the bivariate binomial distribution.
#'
#' The correlation bounds are automatically checked using \code{\link{boundRho}}.
#'
#' @references
#' Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
#' trials with two co-primary binary endpoints. \emph{Statistical Methods in
#' Medical Research}, 34(1), 1-19.
#'
#' @examples
#' # Exact power calculation using Boschloo test
#' power2BinaryExact(
#'   n1 = 100,
#'   n2 = 50,
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0.7,
#'   rho2 = 0.7,
#'   alpha = 0.025,
#'   Test = 'Boschloo'
#' )
#'
#' # Exact power with Fisher exact test
#' power2BinaryExact(
#'   n1 = 80,
#'   n2 = 80,
#'   p11 = 0.6,
#'   p12 = 0.5,
#'   p21 = 0.4,
#'   p22 = 0.3,
#'   rho1 = 0.5,
#'   rho2 = 0.5,
#'   alpha = 0.025,
#'   Test = 'Fisher'
#' )
#'
#' \donttest{
#' # Larger sample sizes (computationally intensive)
#' power2BinaryExact(
#'   n1 = 200,
#'   n2 = 100,
#'   p11 = 0.5,
#'   p12 = 0.4,
#'   p21 = 0.3,
#'   p22 = 0.2,
#'   rho1 = 0.6,
#'   rho2 = 0.6,
#'   alpha = 0.025,
#'   Test = 'Chisq'
#' )
#' }
#'
#' @export
#' @import fpCompare
#' @importFrom stats dbinom pbinom
power2BinaryExact <- function(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test) {

  # Check that rho1 is within valid bounds
  bounds1 <- boundRho(p11, p12)
  if (rho1 < bounds1[1] | rho1 > bounds1[2]) {
    stop(paste0("rho1 must be within [", round(bounds1[1], 4), ", ",
                round(bounds1[2], 4), "]"))
  }

  # Check that rho2 is within valid bounds
  bounds2 <- boundRho(p21, p22)
  if (rho2 < bounds2[1] | rho2 > bounds2[2]) {
    stop(paste0("rho2 must be within [", round(bounds2[1], 4), ", ",
                round(bounds2[2], 4), "]"))
  }

  # Calculate rejection region for the test
  RR <- rrBinary(n1, n2, alpha, Test)

  # Calculate power for individual endpoints using binomial distribution
  # For endpoint k, we sum over all outcomes where the null is rejected
  power1 <- sum(dbinom(0:n1, n1, p11) * pbinom(rowSums(RR) - 1, n2, p21))
  power2 <- sum(dbinom(0:n1, n1, p12) * pbinom(rowSums(RR) - 1, n2, p22))

  # Calculate probability mass functions of bivariate binomial distribution for each group
  # pmass1[i+1, j+1] = P(Y_{1,1} = i, Y_{1,2} = j) for group 1
  pmass1 <- outer(0:n1, 0:n1, function(x, y) dbibinom(n1, x, y, p11, p12, rho1, NULL))

  # pmass2[i+1, j+1] = P(Y_{2,1} = i, Y_{2,2} = j) for group 2
  pmass2 <- outer(0:n2, 0:n2, function(x, y) dbibinom(n2, x, y, p21, p22, rho2, NULL))

  # Calculate power for co-primary endpoints using equation (9)
  # Sum over all outcomes where BOTH endpoints reject the null
  # This requires extracting the joint probabilities where RR[i1, j1] = TRUE
  # and multiplying by the probabilities for endpoint 2
  powerCoprimary <- sum(
    '*'(
      t(pmass1[row(RR)[RR %>>% 0], ])[row(RR)[RR %>>% 0], ],
      t(pmass2[col(RR)[RR %>>% 0], ])[col(RR)[RR %>>% 0], ]
    )
  )

  # Return power as a named vector
  power <- c(power1, power2, powerCoprimary)
  names(power) <- c('Power1', 'Power2', 'powerCoprimary')
  return(power)
}
