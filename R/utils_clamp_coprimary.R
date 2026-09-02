#' Hold a Simulated Co-Primary Power Inside the Bounds of Its Marginals
#'
#' Internal helper for the power functions whose co-primary power is estimated
#' by simulation while the two marginal powers are exact.
#'
#' @param joint Estimated probability that both endpoints are significant.
#' @param power1 Exact power for the first endpoint.
#' @param power2 Exact power for the second endpoint.
#'
#' @return A single number in \code{[max(0, power1 + power2 - 1),
#'   min(power1, power2)]}.
#'
#' @details
#' The probability of an intersection lies between the Bonferroni bound and the
#' smaller of the two marginal probabilities. A simulated estimate paired with
#' exact marginals need not respect either, and at a small number of
#' replications it can visibly fail to, so the estimate is returned at the
#' nearer end of the interval when it falls outside. The adjustment is smaller
#' than the simulation error that caused it.
#'
#' @keywords internal
#' @noRd
.clamp_coprimary <- function(joint, power1, power2) {
  min(max(joint, max(0, power1 + power2 - 1)), min(power1, power2))
}
