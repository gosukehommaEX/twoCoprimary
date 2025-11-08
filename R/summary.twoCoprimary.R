#' Summary Method for twoCoprimary Objects
#'
#' Creates a summary table showing how sample size or power varies across
#' different correlation values for two co-primary endpoints designs. The
#' output format is similar to tables in Sozu et al. (2011), with parameters
#' displayed vertically and correlation values displayed horizontally.
#'
#' @param object An object of class "twoCoprimary" from power or sample size
#'   calculation functions
#' @param rho_values Numeric vector of correlation values to evaluate. Default
#'   is c(0, 0.3, 0.5, 0.8). For binary endpoints with two correlations (rho1
#'   and rho2), both correlations are set to each value in rho_values.
#' @param digits Number of decimal places for rounding numeric values in the
#'   output. Default is 2. Sample sizes are always rounded to integers.
#' @param ... Additional arguments (currently unused)
#'
#' @return A data.frame in wide format with:
#'   \itemize{
#'     \item First column: Parameter names
#'     \item Subsequent columns: Values for each correlation specified in
#'       rho_values
#'   }
#'
#' @details
#' The function recalculates sample size or power across the specified range
#' of correlation values, holding all other parameters constant. This allows
#' comparison of how the correlation between endpoints affects design
#' characteristics.
#'
#' For power calculation results (when n1 and n2 are specified), the function
#' calculates sample sizes needed to achieve the observed power at different
#' correlation values.
#'
#' For sample size calculation results (when power is specified), the function
#' calculates the power achieved by the calculated sample size at different
#' correlation values.
#'
#' The table format follows Sozu et al. (2011) with effect sizes and standard
#' deviations (or proportions, rates, and dispersions for binary/count
#' endpoints) shown for each group, followed by sample sizes and power values
#' for each correlation.
#'
#' @references
#' Sozu, T., Kanou, T., Hamada, C., & Yoshimura, I. (2011). Power and sample
#' size calculations in clinical trials with multiple primary variables.
#' Japanese Journal of Biometrics, 27, 83-96.
#'
#' @examples
#' # Power calculation result - shows sample size at different correlations
#' result_power <- power2Continuous(
#'   n1 = 100, n2 = 100,
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.5, alpha = 0.025,
#'   known_var = TRUE
#' )
#' summary(result_power)
#' summary(result_power, rho_values = c(0, 0.2, 0.4, 0.6, 0.8))
#'
#' # Sample size calculation result - shows power at different correlations
#' result_ss <- ss2Continuous(
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.5, r = 1,
#'   alpha = 0.025, beta = 0.2,
#'   known_var = TRUE
#' )
#' summary(result_ss)
#'
#' @export
summary.twoCoprimary <- function(object, rho_values = c(0, 0.3, 0.5, 0.8),
                                 digits = 2, ...) {

  # Detect endpoint type
  endpoint_type <- detect_endpoint_type_summary(object)

  # Determine calculation mode
  is_power_calc <- "powerCoprimary" %in% names(object) ||
    "powerCont" %in% names(object) ||
    "powerBin" %in% names(object) ||
    "powerCount" %in% names(object)

  # Initialize result list
  result_list <- list()

  # Add parameter rows based on endpoint type
  result_list <- add_parameter_rows(result_list, object, endpoint_type, digits)

  # Calculate sample size or power for each correlation value
  for (rho in rho_values) {
    col_name <- paste0("rho_", format(rho, nsmall = 1))

    if (is_power_calc) {
      # For power calculation: compute sample size at different rho
      ss_result <- calculate_ss_at_rho(object, endpoint_type, rho)
      result_list[[col_name]] <- extract_ss_values(ss_result, endpoint_type,
                                                   object, digits)
    } else {
      # For sample size calculation: compute power at different rho
      power_result <- calculate_power_at_rho(object, endpoint_type, rho)
      result_list[[col_name]] <- extract_power_values(power_result, endpoint_type,
                                                      object, digits)
    }
  }

  # Convert to data.frame
  result_df <- as.data.frame(result_list, stringsAsFactors = FALSE)

  # Set row names as first column
  param_names <- result_df$Parameter
  result_df <- result_df[, -1, drop = FALSE]
  rownames(result_df) <- param_names

  # Add class for potential future methods
  class(result_df) <- c("summary.twoCoprimary", "data.frame")

  return(result_df)
}


# ============================================================================
# Helper function: Detect endpoint type
# ============================================================================

detect_endpoint_type_summary <- function(x) {
  if (all(c("delta1", "delta2", "sd1", "sd2") %in% names(x))) {
    return("continuous")
  } else if (all(c("p11", "p12", "p21", "p22") %in% names(x))) {
    return("binary")
  } else if (all(c("delta", "sd", "p1", "p2") %in% names(x))) {
    return("mixed_cont_binary")
  } else if (all(c("r1", "r2", "nu", "mu1", "mu2") %in% names(x))) {
    return("mixed_count_cont")
  } else {
    stop("Cannot determine endpoint type from object")
  }
}


# ============================================================================
# Helper function: Add parameter rows
# ============================================================================

add_parameter_rows <- function(result_list, x, endpoint_type, digits) {

  if (endpoint_type == "continuous") {
    # Two continuous endpoints
    param_names <- c(
      paste0("delta1 (sd1 = ", round(x$sd1, digits), ")"),
      paste0("delta2 (sd2 = ", round(x$sd2, digits), ")"),
      "n1", "n2", "N", "Power"
    )

    param_values <- c(
      round(x$delta1, digits),
      round(x$delta2, digits),
      NA, NA, NA, NA  # To be filled for each rho
    )

  } else if (endpoint_type == "binary") {
    # Two binary endpoints
    param_names <- c(
      "p1 (Group 1)",
      "p1 (Group 2)",
      "p2 (Group 1)",
      "p2 (Group 2)",
      "n1", "n2", "N", "Power"
    )

    param_values <- c(
      round(x$p11, digits),
      round(x$p21, digits),
      round(x$p12, digits),
      round(x$p22, digits),
      NA, NA, NA, NA  # To be filled for each rho
    )

  } else if (endpoint_type == "mixed_cont_binary") {
    # Mixed continuous and binary
    param_names <- c(
      paste0("delta (sd = ", round(x$sd, digits), ")"),
      "p (Group 1)",
      "p (Group 2)",
      "n1", "n2", "N", "Power"
    )

    param_values <- c(
      round(x$delta, digits),
      round(x$p1, digits),
      round(x$p2, digits),
      NA, NA, NA, NA  # To be filled for each rho
    )

  } else if (endpoint_type == "mixed_count_cont") {
    # Mixed count and continuous
    param_names <- c(
      paste0("r1 (nu = ", round(x$nu, digits), ")"),
      paste0("r2 (nu = ", round(x$nu, digits), ")"),
      paste0("mu1 (sd = ", round(x$sd, digits), ")"),
      paste0("mu2 (sd = ", round(x$sd, digits), ")"),
      "n1", "n2", "N", "Power"
    )

    param_values <- c(
      round(x$r1, digits),
      round(x$r2, digits),
      round(x$mu1, digits),
      round(x$mu2, digits),
      NA, NA, NA, NA  # To be filled for each rho
    )
  }

  result_list$Parameter <- param_names

  return(result_list)
}


# ============================================================================
# Helper function: Calculate sample size at given rho
# ============================================================================

calculate_ss_at_rho <- function(x, endpoint_type, rho) {

  # Calculate power from existing result
  if ("powerCoprimary" %in% names(x)) {
    target_power <- x$powerCoprimary
  } else if ("powerCont" %in% names(x) && "powerBin" %in% names(x)) {
    target_power <- min(x$powerCont, x$powerBin)
  } else if ("powerCount" %in% names(x) && "powerCont" %in% names(x)) {
    target_power <- min(x$powerCount, x$powerCont)
  } else {
    stop("Cannot determine target power from object")
  }

  beta <- 1 - target_power

  if (endpoint_type == "continuous") {
    ss2Continuous(
      delta1 = x$delta1, delta2 = x$delta2,
      sd1 = x$sd1, sd2 = x$sd2,
      rho = rho, r = x$r,
      alpha = x$alpha, beta = beta,
      known_var = ifelse("known_var" %in% names(x), x$known_var, TRUE),
      nMC = ifelse("nMC" %in% names(x) && !is.na(x$nMC), x$nMC, 1000)
    )

  } else if (endpoint_type == "binary") {
    ss2BinaryApprox(
      p11 = x$p11, p12 = x$p12,
      p21 = x$p21, p22 = x$p22,
      rho1 = rho, rho2 = rho,
      r = x$r, alpha = x$alpha, beta = beta,
      Test = x$Test
    )

  } else if (endpoint_type == "mixed_cont_binary") {
    ss2MixedContinuousBinary(
      delta = x$delta, sd = x$sd,
      p1 = x$p1, p2 = x$p2,
      rho = rho, r = x$r,
      alpha = x$alpha, beta = beta,
      Test = x$Test,
      nMC = ifelse("nMC" %in% names(x) && !is.na(x$nMC), x$nMC, 1000)
    )

  } else if (endpoint_type == "mixed_count_cont") {
    ss2MixedCountContinuous(
      r1 = x$r1, r2 = x$r2,
      nu = x$nu, t = x$t,
      mu1 = x$mu1, mu2 = x$mu2, sd = x$sd,
      rho1 = rho, rho2 = rho,
      r = x$r, alpha = x$alpha, beta = beta
    )
  }
}


# ============================================================================
# Helper function: Calculate power at given rho
# ============================================================================

calculate_power_at_rho <- function(x, endpoint_type, rho) {

  if (endpoint_type == "continuous") {
    power2Continuous(
      n1 = x$n1, n2 = x$n2,
      delta1 = x$delta1, delta2 = x$delta2,
      sd1 = x$sd1, sd2 = x$sd2,
      rho = rho, alpha = x$alpha,
      known_var = ifelse("known_var" %in% names(x), x$known_var, TRUE),
      nMC = ifelse("nMC" %in% names(x) && !is.na(x$nMC), x$nMC, 1000)
    )

  } else if (endpoint_type == "binary") {
    power2BinaryApprox(
      n1 = x$n1, n2 = x$n2,
      p11 = x$p11, p12 = x$p12,
      p21 = x$p21, p22 = x$p22,
      rho1 = rho, rho2 = rho,
      alpha = x$alpha,
      Test = x$Test
    )

  } else if (endpoint_type == "mixed_cont_binary") {
    power2MixedContinuousBinary(
      n1 = x$n1, n2 = x$n2,
      delta = x$delta, sd = x$sd,
      p1 = x$p1, p2 = x$p2,
      rho = rho, alpha = x$alpha,
      Test = x$Test,
      nMC = ifelse("nMC" %in% names(x) && !is.na(x$nMC), x$nMC, 1000)
    )

  } else if (endpoint_type == "mixed_count_cont") {
    power2MixedCountContinuous(
      n1 = x$n1, n2 = x$n2,
      r1 = x$r1, r2 = x$r2,
      nu = x$nu, t = x$t,
      mu1 = x$mu1, mu2 = x$mu2, sd = x$sd,
      rho1 = rho, rho2 = rho,
      alpha = x$alpha
    )
  }
}


# ============================================================================
# Helper function: Extract sample size values
# ============================================================================

extract_ss_values <- function(ss_result, endpoint_type, x, digits) {

  # Get parameter values
  if (endpoint_type == "continuous") {
    param_values <- c(
      round(x$delta1, digits),
      round(x$delta2, digits)
    )
  } else if (endpoint_type == "binary") {
    param_values <- c(
      round(x$p11, digits),
      round(x$p21, digits),
      round(x$p12, digits),
      round(x$p22, digits)
    )
  } else if (endpoint_type == "mixed_cont_binary") {
    param_values <- c(
      round(x$delta, digits),
      round(x$p1, digits),
      round(x$p2, digits)
    )
  } else if (endpoint_type == "mixed_count_cont") {
    param_values <- c(
      round(x$r1, digits),
      round(x$r2, digits),
      round(x$mu1, digits),
      round(x$mu2, digits)
    )
  }

  # Extract target power from the original object x
  # (this is the power achieved by the original n1, n2)
  if ("powerCoprimary" %in% names(x)) {
    target_power <- x$powerCoprimary
  } else if ("powerCont" %in% names(x) && "powerBin" %in% names(x)) {
    target_power <- min(x$powerCont, x$powerBin)
  } else if ("powerCount" %in% names(x) && "powerCont" %in% names(x)) {
    target_power <- min(x$powerCount, x$powerCont)
  } else {
    target_power <- NA
  }

  # Combine parameter values with sample sizes and target power
  values <- c(
    param_values,
    ss_result$n1,
    ss_result$n2,
    ss_result$N,
    round(target_power, digits)
  )

  return(values)
}


# ============================================================================
# Helper function: Extract power values
# ============================================================================

extract_power_values <- function(power_result, endpoint_type, x, digits) {

  # Get parameter values
  if (endpoint_type == "continuous") {
    param_values <- c(
      round(x$delta1, digits),
      round(x$delta2, digits)
    )
  } else if (endpoint_type == "binary") {
    param_values <- c(
      round(x$p11, digits),
      round(x$p21, digits),
      round(x$p12, digits),
      round(x$p22, digits)
    )
  } else if (endpoint_type == "mixed_cont_binary") {
    param_values <- c(
      round(x$delta, digits),
      round(x$p1, digits),
      round(x$p2, digits)
    )
  } else if (endpoint_type == "mixed_count_cont") {
    param_values <- c(
      round(x$r1, digits),
      round(x$r2, digits),
      round(x$mu1, digits),
      round(x$mu2, digits)
    )
  }

  # Extract co-primary power
  if ("powerCoprimary" %in% names(power_result)) {
    coprimary_power <- power_result$powerCoprimary
  } else if ("powerCont" %in% names(power_result) &&
             "powerBin" %in% names(power_result)) {
    coprimary_power <- min(power_result$powerCont, power_result$powerBin)
  } else if ("powerCount" %in% names(power_result) &&
             "powerCont" %in% names(power_result)) {
    coprimary_power <- min(power_result$powerCount, power_result$powerCont)
  } else {
    coprimary_power <- NA
  }

  # Combine parameter values with sample sizes (from original object) and power
  values <- c(
    param_values,
    x$n1,
    x$n2,
    x$N,
    round(coprimary_power, digits)
  )

  return(values)
}
