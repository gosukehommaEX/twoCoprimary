#' Internal Sequential Search Algorithm for Sample Size Determination
#'
#' This is an internal function used by all ss2* functions to find the minimum
#' sample size that achieves target power using sequential search. The algorithm
#' follows Homma and Yoshida (2025) Algorithm 1.
#'
#' @param initial_n2 Initial sample size for group 2
#' @param r Allocation ratio n1/n2
#' @param target_power Target power (1 - beta)
#' @param power_fun Function to calculate power. Must accept n1, n2 as the first
#'   two arguments and return a list with 'powerCoprimary' element
#' @param ... Additional arguments passed to power_fun
#'
#' @return List with elements:
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{N}{Total sample size (n1 + n2)}
#'
#' @details
#' Algorithm (Homma and Yoshida 2025, Algorithm 1):
#' \itemize{
#'   \item Step 1: Calculate power at initial sample size
#'   \item Step 2a: If power >= target, decrease n2 until power < target, then add 1
#'   \item Step 2b: If power < target, increase n2 until power >= target
#'   \item Step 3: Return final sample sizes
#' }
#'
#' This ensures the minimum sample size that achieves the target power is found,
#' even when the initial estimate is too large.
#'
#' @keywords internal
#' @noRd
#' @import fpCompare
.ss_sequential_search <- function(initial_n2, r, target_power, power_fun, ...) {

  # Start from initial value
  n2 <- initial_n2
  n1 <- ceiling(r * n2)

  # Step 1: Calculate initial power
  power <- power_fun(n1, n2, ...)[["powerCoprimary"]]

  # Step 2: Adjust sample size to achieve target power
  if (power %>=% target_power) {

    # Step 2a: Power is too high or just right
    # Decrease n2 until power drops below target
    while (power %>=% target_power && n2 > 1) {
      n2 <- n2 - 1
      n1 <- ceiling(r * n2)
      power <- power_fun(n1, n2, ...)[["powerCoprimary"]]
    }

    # Add 1 back to get the minimum sample size that achieves target power
    n2 <- n2 + 1

  } else {

    # Step 2b: Power is too low
    # Increase n2 until power reaches target
    while (power %<<% target_power) {
      n2 <- n2 + 1
      n1 <- ceiling(r * n2)
      power <- power_fun(n1, n2, ...)[["powerCoprimary"]]
    }
  }

  # Step 3: Determine final sample sizes
  n1 <- ceiling(r * n2)
  N <- n1 + n2

  return(list(n1 = n1, n2 = n2, N = N))
}
