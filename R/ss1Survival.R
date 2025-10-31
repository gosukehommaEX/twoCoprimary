#' Sample Size Calculation for a Single Time-to-Event Endpoint
#'
#' Calculates the required sample size for a two-arm superiority trial with a
#' single time-to-event endpoint assuming exponential distributions, as described
#' in Hamasaki et al. (2013).
#'
#' @param lambdaT Hazard rate for the test group (λT > 0)
#' @param lambdaC Hazard rate for the control group (λC > 0)
#' @param tau0 Accrual period (recruitment time in years)
#' @param tau Follow-up period after recruitment ends (years)
#' @param r Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r > 0)
#' @param alpha One-sided significance level (typically 0.025 or 0.05)
#' @param beta Target type II error rate (typically 0.1 or 0.2)
#' @param censor Character string specifying censoring approach. One of:
#'   \itemize{
#'     \item \code{"variance"}: Censoring incorporated into variance only (default)
#'     \item \code{"both"}: Censoring incorporated into both variance and correlation
#'   }
#'
#' @return A data frame with the following columns:
#'   \item{lambdaT}{Hazard rate for test group}
#'   \item{lambdaC}{Hazard rate for control group}
#'   \item{tau0}{Accrual period}
#'   \item{tau}{Follow-up period}
#'   \item{r}{Allocation ratio}
#'   \item{alpha}{One-sided significance level}
#'   \item{beta}{Type II error rate}
#'   \item{censor}{Censoring approach}
#'   \item{n1}{Required sample size for group 1}
#'   \item{n2}{Required sample size for group 2}
#'   \item{n}{Total sample size (n1 + n2)}
#'
#' @details
#' This function calculates sample size for a single time-to-event endpoint
#' using the log-rank test under the exponential distribution assumption.
#'
#' **Exponential Model:**
#' The survival times are assumed to follow exponential distributions with
#' hazard rates λT and λC for test and control groups, respectively.
#'
#' **Censoring:**
#' Participants are uniformly recruited over [0, τ0] and followed until time τ
#' (τ > τ0). The censoring parameter is:
#' \deqn{ψ(λ) = \frac{1 - \exp(-λ(τ0 + τ)) - \exp(-λτ)}{λτ0}}
#'
#' **Sample Size Formula:**
#' When censoring is incorporated into variance only (equation 3 in Hamasaki et al. 2013):
#' \deqn{n_2 = \left\lceil \frac{(1 + 1/r)}{\log^2(λ_C/λ_T)}
#'       \left(z_\alpha \sqrt{\frac{1}{ψ_T + ψ_C/r}} +
#'       z_\beta \sqrt{\frac{1}{ψ_T} + \frac{1}{ψ_C}}\right)^2 \right\rceil}
#'
#' where ψT = rλT ψ(λT) + (1-r)λC ψ(λC) and ψC = rλC ψ(λC) + (1-r)λC ψ(λC).
#'
#' @references
#' Hamasaki, T., Sugimoto, T., Evans, S. R., & Sozu, T. (2013). Sample size
#' determination for clinical trials with co-primary outcomes: exponential event
#' times. \emph{Pharmaceutical Statistics}, 12(1), 28-34.
#'
#' @examples
#' # Example from Hamasaki et al. (2013)
#' ss1Survival(
#'   lambdaT = log(2) / 5,
#'   lambdaC = log(2) / 3,
#'   tau0 = 2,
#'   tau = 5,
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2
#' )
#'
#' @export
#' @importFrom stats qnorm
ss1Survival <- function(lambdaT, lambdaC, tau0, tau, r, alpha, beta,
                        censor = "variance") {

  # Input validation
  if (lambdaT <= 0) stop("lambdaT must be positive")
  if (lambdaC <= 0) stop("lambdaC must be positive")
  if (tau0 <= 0) stop("tau0 must be positive")
  if (tau <= tau0) stop("tau must be greater than tau0")
  if (r <= 0) stop("r must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)")
  if (beta <= 0 || beta >= 1) stop("beta must be in (0, 1)")
  if (!censor %in% c("variance", "both")) stop("censor must be 'variance' or 'both'")

  # Standard normal quantiles
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(1 - beta)

  # Calculate psi function (censoring parameter)
  psi <- function(lambda) {
    (1 - exp(-lambda * (tau0 + tau)) - exp(-lambda * tau)) / (lambda * tau0)
  }

  psi_T <- psi(lambdaT)
  psi_C <- psi(lambdaC)

  # Calculate effective sample sizes under censoring
  pi_T <- (r * lambdaT + (1 - r) * lambdaC) * psi_T
  pi_C <- (r * lambdaC + (1 - r) * lambdaC) * psi_C

  # Log hazard ratio
  log_HR <- log(lambdaC / lambdaT)

  # Variance under H0 and H1
  var_H0 <- 1 / (pi_T + pi_C / r)
  var_H1 <- 1 / (lambdaT * psi_T) + 1 / (lambdaC * psi_C)

  # Calculate required sample size
  n2 <- ceiling(
    ((z_alpha * sqrt(var_H0) + z_beta * sqrt(var_H1)) / log_HR)^2 * (1 + 1 / r)
  )

  n1 <- ceiling(r * n2)
  n <- n1 + n2

  # Return result
  result <- data.frame(
    lambdaT = lambdaT, lambdaC = lambdaC, tau0 = tau0, tau = tau,
    r = r, alpha = alpha, beta = beta, censor = censor,
    n1 = n1, n2 = n2, n = n
  )
  return(result)
}
