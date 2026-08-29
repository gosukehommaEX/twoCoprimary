#' Position of the Last Member of Each Tie Group in a Sorted Vector
#'
#' Internal helper used by \code{rr1Binary} for the exact unconditional tests.
#' Given a vector already sorted from the most extreme to the least extreme
#' value of the ordering statistic, the function returns, for each element, the
#' position of the last element of the tie group it belongs to.
#'
#' The tail event of an exact unconditional test is the set of outcomes at least
#' as extreme as the observed one, so all outcomes sharing the same value of the
#' ordering statistic must receive the same tail probability, namely the one
#' accumulated up to the last member of their group. A plain cumulative sum
#' instead gives the earlier members of a group a smaller tail probability, and
#' which member comes first depends on the order that \code{order} happens to
#' return.
#'
#' A relative tolerance is used rather than the absolute tolerance of
#' \code{fpCompare}, because the sequences passed to this function include exact
#' test p-values that can be far smaller than 1e-12, and an absolute tolerance
#' would merge values that are genuinely distinct.
#'
#' @param x A numeric vector sorted in either increasing or decreasing order
#'
#' @return An integer vector of the same length as \code{x} holding the position
#'   of the last member of the tie group of each element
#'
#' @keywords internal
#' @noRd
.tie_last <- function(x) {

  n <- length(x)

  if (n <= 1) {
    return(seq_len(n))
  }

  # Relative tolerance for treating two adjacent values as tied
  tol <- 1e-10
  ref <- pmax(abs(x[-n]), abs(x[-1]))
  new_group <- abs(diff(x)) > tol * ref

  # Position of the last element of each group, repeated over its members
  group_end <- c(which(new_group), n)
  group_size <- diff(c(0L, group_end))

  rep(as.integer(group_end), times = group_size)
}
