#' Create Design Comparison Table for Two Co-Primary Endpoints
#'
#' Generates a comprehensive table comparing sample sizes or power across
#' different parameter combinations and correlation values. This function is
#' useful for sensitivity analyses and exploring how design parameters affect
#' statistical properties.
#'
#' @param param_grid A data.frame containing parameter combinations. Required
#'   columns depend on endpoint_type and calculation_mode:
#'   \itemize{
#'     \item For continuous endpoints (sample size): delta1, delta2, sd1, sd2
#'     \item For continuous endpoints (power): n1, n2, delta1, delta2, sd1, sd2
#'     \item For binary endpoints (sample size): p11, p12, p21, p22
#'     \item For binary endpoints (power): n1, n2, p11, p12, p21, p22
#'     \item For mixed continuous-binary (sample size): delta, sd, p1, p2
#'     \item For mixed continuous-binary (power): n1, n2, delta, sd, p1, p2
#'     \item For mixed count-continuous (sample size): r1, r2, nu, t, mu1, mu2, sd
#'     \item For mixed count-continuous (power): n1, n2, r1, r2, nu, t, mu1, mu2, sd
#'   }
#' @param rho_values Numeric vector of correlation values to evaluate. Default
#'   is c(0, 0.3, 0.5, 0.8).
#' @param r Allocation ratio (n1/n2). Required for sample size calculation.
#'   Default is 1.
#' @param alpha One-sided significance level. Default is 0.025.
#' @param beta Type II error rate (1 - power). Required for sample size
#'   calculation. Default is 0.2 (power = 0.8).
#' @param endpoint_type Character string specifying endpoint type:
#'   "continuous", "binary", "mixed_cont_binary", or "mixed_count_cont".
#' @param Test Test method for binary endpoints: "AN" (asymptotic normal),
#'   "ANc" (with continuity correction), "AS" (arcsine), or "ASc". Default is
#'   "AN". Only used for binary and mixed_cont_binary endpoints.
#' @param known_var Logical indicating whether variance is known for continuous
#'   endpoints. Default is TRUE.
#' @param nMC Number of Monte Carlo simulations for certain calculations.
#'   Default is 1000.
#' @param output_var Character string specifying which variable to output in
#'   the result columns: "N" (total sample size, default for sample size
#'   calculation) or "powerCoprimary" (co-primary power, default for power
#'   calculation).
#'
#' @return A data.frame of class "twoCoprimary_table" with:
#'   \itemize{
#'     \item Parameter columns (from param_grid)
#'     \item Result columns for each correlation value (rho_0.0, rho_0.3, etc.)
#'   }
#'
#' @details
#' This function performs systematic calculations across all combinations of
#' parameters specified in param_grid and correlation values in rho_values.
#'
#' The calculation mode (sample size vs power) is automatically determined:
#' \itemize{
#'   \item If param_grid contains n1 and n2: calculates power
#'   \item Otherwise: calculates sample size (requires r, alpha, beta)
#' }
#'
#' For binary endpoints with two correlations (rho1, rho2), both are set to
#' the same value from rho_values for each calculation.
#'
#' The output format follows the style of Sozu et al. (2011), with parameters
#' displayed in the leftmost columns and results for each correlation in
#' subsequent columns.
#'
#' @references
#' Sozu, T., Kanou, T., Hamada, C., & Yoshimura, I. (2011). Power and sample
#' size calculations in clinical trials with multiple primary variables.
#' Japanese Journal of Biometrics, 27, 83-96.
#'
#' @examples
#' # Sample size calculation for continuous endpoints
#' param_grid <- expand.grid(
#'   delta1 = c(0.3, 0.5),
#'   delta2 = c(0.1, 0.2, 0.3),
#'   sd1 = c(1.0, 1.5),
#'   sd2 = c(1.0, 1.5)
#' )
#'
#' result <- design_table(
#'   param_grid = param_grid,
#'   rho_values = c(0, 0.3, 0.5, 0.8),
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   endpoint_type = "continuous"
#' )
#' print(result)
#'
#' # Power calculation for continuous endpoints
#' param_grid_power <- expand.grid(
#'   n1 = c(50, 100),
#'   n2 = c(50, 100),
#'   delta1 = 0.5,
#'   delta2 = 0.5,
#'   sd1 = 1.0,
#'   sd2 = 1.0
#' )
#'
#' result_power <- design_table(
#'   param_grid = param_grid_power,
#'   rho_values = c(0, 0.3, 0.5, 0.8),
#'   alpha = 0.025,
#'   endpoint_type = "continuous"
#' )
#' print(result_power)
#'
#' # Binary endpoints
#' param_grid_binary <- expand.grid(
#'   p11 = c(0.6, 0.7),
#'   p12 = c(0.4, 0.5),
#'   p21 = c(0.4, 0.5),
#'   p22 = c(0.2, 0.3)
#' )
#'
#' result_binary <- design_table(
#'   param_grid = param_grid_binary,
#'   rho_values = c(0.3, 0.5, 0.7),
#'   r = 1,
#'   alpha = 0.025,
#'   beta = 0.2,
#'   endpoint_type = "binary",
#'   Test = "AN"
#' )
#' print(result_binary)
#'
#' @export
design_table <- function(param_grid,
                         rho_values = c(0, 0.3, 0.5, 0.8),
                         r = 1,
                         alpha = 0.025,
                         beta = 0.2,
                         endpoint_type = c("continuous", "binary",
                                           "mixed_cont_binary",
                                           "mixed_count_cont"),
                         Test = "AN",
                         known_var = TRUE,
                         nMC = 1000,
                         output_var = NULL) {

  # Match endpoint_type
  endpoint_type <- match.arg(endpoint_type)

  # Validate param_grid
  if (!is.data.frame(param_grid)) {
    stop("param_grid must be a data.frame")
  }

  # Determine calculation mode
  is_power_calc <- all(c("n1", "n2") %in% names(param_grid))

  # Set default output_var if not specified
  if (is.null(output_var)) {
    output_var <- if (is_power_calc) "powerCoprimary" else "N"
  }

  # Validate required columns based on endpoint type and calculation mode
  validate_param_grid(param_grid, endpoint_type, is_power_calc)

  # Initialize result data.frame with parameter columns
  result_df <- param_grid

  # Calculate for each correlation value
  for (rho in rho_values) {
    col_name <- paste0("rho_", format(rho, nsmall = 1))

    # Calculate for each row of param_grid
    col_values <- numeric(nrow(param_grid))

    for (i in 1:nrow(param_grid)) {
      params <- as.list(param_grid[i, ])

      # Check if rho is within valid bounds for this parameter combination
      if (!is_rho_valid(params, rho, endpoint_type)) {
        col_values[i] <- NA
        next
      }

      tryCatch({
        if (is_power_calc) {
          # Power calculation
          result <- calculate_power_design_table(
            params, rho, endpoint_type, alpha, Test, known_var, nMC
          )
        } else {
          # Sample size calculation
          result <- calculate_ss_design_table(
            params, rho, r, alpha, beta, endpoint_type, Test, known_var, nMC
          )
        }

        # Extract the requested output variable
        col_values[i] <- result[[output_var]]
      }, error = function(e) {
        # If calculation fails, return NA
        col_values[i] <<- NA
        warning("Calculation failed for row ", i, " at rho = ", rho,
                ": ", e$message, call. = FALSE)
      })
    }

    result_df[[col_name]] <- col_values
  }

  # Add class for potential S3 methods
  class(result_df) <- c("twoCoprimary_table", "data.frame")

  return(result_df)
}


# ============================================================================
# Helper function: Validate param_grid
# ============================================================================

validate_param_grid <- function(param_grid, endpoint_type, is_power_calc) {

  if (endpoint_type == "continuous") {
    required_cols <- c("delta1", "delta2", "sd1", "sd2")
    if (is_power_calc) {
      required_cols <- c("n1", "n2", required_cols)
    }

  } else if (endpoint_type == "binary") {
    required_cols <- c("p11", "p12", "p21", "p22")
    if (is_power_calc) {
      required_cols <- c("n1", "n2", required_cols)
    }

  } else if (endpoint_type == "mixed_cont_binary") {
    required_cols <- c("delta", "sd", "p1", "p2")
    if (is_power_calc) {
      required_cols <- c("n1", "n2", required_cols)
    }

  } else if (endpoint_type == "mixed_count_cont") {
    required_cols <- c("r1", "r2", "nu", "t", "mu1", "mu2", "sd")
    if (is_power_calc) {
      required_cols <- c("n1", "n2", required_cols)
    }
  }

  missing_cols <- setdiff(required_cols, names(param_grid))
  if (length(missing_cols) > 0) {
    stop("param_grid is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  invisible(TRUE)
}


# ============================================================================
# Helper function: Check if rho is within valid bounds
# ============================================================================

is_rho_valid <- function(params, rho, endpoint_type) {

  if (endpoint_type == "continuous" || endpoint_type == "mixed_cont_binary") {
    # Continuous endpoints have no strict correlation bounds
    # Mixed continuous-binary uses biserial correlation (typically -1 to 1)
    return(TRUE)

  } else if (endpoint_type == "binary") {
    # Calculate correlation bounds for binary endpoints
    bounds1 <- corrbound2Binary(params$p11, params$p12)
    bounds2 <- corrbound2Binary(params$p21, params$p22)

    # Both correlations must be within their respective bounds
    rho_min <- max(bounds1[1], bounds2[1])
    rho_max <- min(bounds1[2], bounds2[2])

    # Check if rho is within valid range (with small tolerance)
    return(rho >= (rho_min - 1e-6) && rho <= (rho_max + 1e-6))

  } else if (endpoint_type == "mixed_count_cont") {
    # Calculate correlation bounds for mixed count-continuous endpoints
    lambda1 <- params$r1 * params$t
    lambda2 <- params$r2 * params$t

    bounds1 <- corrbound2MixedCountContinuous(lambda1, params$nu,
                                              params$mu1, params$sd)
    bounds2 <- corrbound2MixedCountContinuous(lambda2, params$nu,
                                              params$mu2, params$sd)

    # Both correlations must be within their respective bounds
    rho_min <- max(bounds1[1], bounds2[1])
    rho_max <- min(bounds1[2], bounds2[2])

    # Check if rho is within valid range (with small tolerance)
    return(rho >= (rho_min - 1e-6) && rho <= (rho_max + 1e-6))
  }

  # Default: allow the correlation
  return(TRUE)
}


# ============================================================================
# Helper function: Calculate power for design table
# ============================================================================

calculate_power_design_table <- function(params, rho, endpoint_type,
                                         alpha, Test, known_var, nMC) {

  if (endpoint_type == "continuous") {
    power2Continuous(
      n1 = params$n1, n2 = params$n2,
      delta1 = params$delta1, delta2 = params$delta2,
      sd1 = params$sd1, sd2 = params$sd2,
      rho = rho, alpha = alpha,
      known_var = known_var, nMC = nMC
    )

  } else if (endpoint_type == "binary") {
    # Use exact method for Chisq or Fisher, approximate otherwise
    if (Test %in% c("Chisq", "Fisher")) {
      power2BinaryExact(
        n1 = params$n1, n2 = params$n2,
        p11 = params$p11, p12 = params$p12,
        p21 = params$p21, p22 = params$p22,
        rho1 = rho, rho2 = rho,
        alpha = alpha, Test = Test
      )
    } else {
      power2BinaryApprox(
        n1 = params$n1, n2 = params$n2,
        p11 = params$p11, p12 = params$p12,
        p21 = params$p21, p22 = params$p22,
        rho1 = rho, rho2 = rho,
        alpha = alpha, Test = Test
      )
    }

  } else if (endpoint_type == "mixed_cont_binary") {
    power2MixedContinuousBinary(
      n1 = params$n1, n2 = params$n2,
      delta = params$delta, sd = params$sd,
      p1 = params$p1, p2 = params$p2,
      rho = rho, alpha = alpha,
      Test = Test, nMC = nMC
    )

  } else if (endpoint_type == "mixed_count_cont") {
    power2MixedCountContinuous(
      n1 = params$n1, n2 = params$n2,
      r1 = params$r1, r2 = params$r2,
      nu = params$nu, t = params$t,
      mu1 = params$mu1, mu2 = params$mu2, sd = params$sd,
      rho1 = rho, rho2 = rho,
      alpha = alpha
    )
  }
}


# ============================================================================
# Helper function: Calculate sample size for design table
# ============================================================================

calculate_ss_design_table <- function(params, rho, r, alpha, beta,
                                      endpoint_type, Test, known_var, nMC) {

  if (endpoint_type == "continuous") {
    ss2Continuous(
      delta1 = params$delta1, delta2 = params$delta2,
      sd1 = params$sd1, sd2 = params$sd2,
      rho = rho, r = r,
      alpha = alpha, beta = beta,
      known_var = known_var, nMC = nMC
    )

  } else if (endpoint_type == "binary") {
    # Use exact method for Chisq or Fisher, approximate otherwise
    if (Test %in% c("Chisq", "Fisher")) {
      ss2BinaryExact(
        p11 = params$p11, p12 = params$p12,
        p21 = params$p21, p22 = params$p22,
        rho1 = rho, rho2 = rho,
        r = r, alpha = alpha, beta = beta,
        Test = Test
      )
    } else {
      ss2BinaryApprox(
        p11 = params$p11, p12 = params$p12,
        p21 = params$p21, p22 = params$p22,
        rho1 = rho, rho2 = rho,
        r = r, alpha = alpha, beta = beta,
        Test = Test
      )
    }

  } else if (endpoint_type == "mixed_cont_binary") {
    ss2MixedContinuousBinary(
      delta = params$delta, sd = params$sd,
      p1 = params$p1, p2 = params$p2,
      rho = rho, r = r,
      alpha = alpha, beta = beta,
      Test = Test, nMC = nMC
    )

  } else if (endpoint_type == "mixed_count_cont") {
    ss2MixedCountContinuous(
      r1 = params$r1, r2 = params$r2,
      nu = params$nu, t = params$t,
      mu1 = params$mu1, mu2 = params$mu2, sd = params$sd,
      rho1 = rho, rho2 = rho,
      r = r, alpha = alpha, beta = beta
    )
  }
}


# ============================================================================
# Print method for twoCoprimary_table
# ============================================================================

#' Print Method for twoCoprimary_table Objects
#'
#' Provides a clean display of design comparison tables
#'
#' @param x An object of class "twoCoprimary_table"
#' @param ... Additional arguments passed to print.data.frame
#'
#' @return Invisibly returns the original object
#'
#' @export
print.twoCoprimary_table <- function(x, ...) {
  cat("\nDesign Comparison Table for Two Co-Primary Endpoints\n")
  cat("======================================================\n\n")

  # Print as data.frame without row names
  print.data.frame(x, row.names = FALSE, ...)

  invisible(x)
}
