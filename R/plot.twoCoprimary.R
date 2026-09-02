#' Plot Method for twoCoprimary Objects
#'
#' Visualizes power or sample size relationships for two co-primary endpoints
#' designs. The function automatically determines the appropriate plot based on
#' the input object.
#'
#' @param x An object of class "twoCoprimary" from power or sample size
#'   calculation functions
#' @param type Type of plot to generate:
#'   \describe{
#'     \item{"power_curve"}{Power as a function of sample size (default for
#'       power calculation results)}
#'     \item{"sample_size_rho"}{Sample size as a function of correlation
#'       (default for sample size calculation results)}
#'     \item{"effect_contour"}{Contour plot showing combinations of
#'       standardized effect sizes achieving target power. The two axes, and
#'       the two columns of the returned data, are \eqn{\delta_k / \sigma_k}}
#'   }
#' @param n_points Number of points to compute for the curve. Default is 50.
#'   Each point runs a full calculation, so a smaller value is advisable for
#'   objects produced with an exact binary test or with Monte Carlo
#'   integration.
#' @param n_range Sample size range for power_curve plot. If NULL, automatically
#'   determined from the object.
#' @param rho_range Correlation range for sample_size_rho plot. If NULL, the
#'   range is seq(0, 0.9, length.out = n_points) for continuous and mixed
#'   continuous-binary endpoints, and is derived from the Frechet-Hoeffding
#'   bounds of the supplied marginal parameters for binary and mixed
#'   count-continuous endpoints.
#' @param col Line color. Default is "steelblue".
#' @param lwd Line width. Default is 2.
#' @param main Plot title. If NULL, automatically generated.
#' @param xlab X-axis label. If NULL, automatically generated.
#' @param ylab Y-axis label. If NULL, automatically generated.
#' @param show_reference Logical. If TRUE, shows reference lines (e.g., target
#'   power, current values). Default is TRUE.
#' @param ... Additional graphical parameters passed to plot()
#'
#' @details
#' The function creates publication-quality plots to visualize the relationship
#' between design parameters and statistical properties. The plot type is
#' automatically selected based on the input object, but can be overridden
#' using the \code{type} argument.
#'
#' For power calculation results (when n1 and n2 are specified), the default
#' is to show how power changes with sample size.
#'
#' For sample size calculation results (when power is specified), the default
#' is to show how required sample size changes with correlation.
#'
#' The function works with all endpoint types (continuous, binary, mixed)
#' by automatically detecting the appropriate parameters from the input object.
#'
#' @return Invisibly returns the data used to create the plot as a data frame.
#'
#' @examples
#' # Power calculation result
#' result_power <- power2Continuous(
#'   n1 = 100, n2 = 100,
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.5, alpha = 0.025,
#'   known_var = TRUE
#' )
#' plot(result_power)  # Shows power curve
#'
#' # Sample size calculation result
#' result_ss <- ss2Continuous(
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.5, r = 1,
#'   alpha = 0.025, beta = 0.2,
#'   known_var = TRUE
#' )
#' plot(result_ss)  # Shows sample size vs correlation
#'
#' # Custom plot with specified type
#' plot(result_power, type = "power_curve", n_range = c(50, 200))
#'
#' @export
#' @importFrom graphics abline contour grid points text
plot.twoCoprimary <- function(x, type = NULL, n_points = 50,
                              n_range = NULL, rho_range = NULL,
                              col = "steelblue", lwd = 2,
                              main = NULL, xlab = NULL, ylab = NULL,
                              show_reference = TRUE, ...) {

  # Detect endpoint type and calculation mode
  is_power_calc <- "powerCoprimary" %in% names(x)
  is_ss_calc <- "N" %in% names(x)

  # Determine default plot type if not specified
  if (is.null(type)) {
    type <- if (is_power_calc) "power_curve" else "sample_size_rho"
  }

  # Validate type
  type <- match.arg(type, c("power_curve", "sample_size_rho", "effect_contour"))

  # Detect endpoint type for appropriate function calls
  endpoint_type <- detect_endpoint_type(x)

  # Generate appropriate plot
  if (type == "power_curve") {
    plot_data <- plot_power_curve(x, endpoint_type, n_points, n_range,
                                  col, lwd, main, xlab, ylab,
                                  show_reference, ...)
  } else if (type == "sample_size_rho") {
    plot_data <- plot_sample_size_rho(x, endpoint_type, n_points, rho_range,
                                      col, lwd, main, xlab, ylab,
                                      show_reference, ...)
  } else if (type == "effect_contour") {
    plot_data <- plot_effect_contour(x, endpoint_type, n_points,
                                     main, xlab, ylab, ...)
  }

  invisible(plot_data)
}


# ============================================================================
# Helper function: Detect endpoint type
# ============================================================================

detect_endpoint_type <- function(x) {
  if (all(c("delta1", "delta2", "sd1", "sd2") %in% names(x))) {
    return("continuous")
  } else if (all(c("p11", "p12", "p21", "p22") %in% names(x))) {
    return("binary")
  } else if (all(c("delta", "sd", "p1", "p2") %in% names(x))) {
    return("mixed_cont_binary")
  } else if (all(c("r1", "r2", "nu", "mu1", "mu2") %in% names(x))) {
    return("mixed_count_cont")
  } else {
    stop("plot() supports results from the two co-primary endpoint functions ",
         "(ss2*, power2*, twoCoprimary2*). Results from the single endpoint ",
         "helpers ss1Continuous(), ss1Count() and ss1BinaryApprox() do not ",
         "carry the parameters a co-primary plot needs.")
  }
}

# Exact binary test names, used to route plotting helpers to the exact rather
# than the approximate implementation
EXACT_BINARY_TESTS <- c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")


# ============================================================================
# Helper function: Plot power curve
# ============================================================================

plot_power_curve <- function(x, endpoint_type, n_points, n_range,
                             col, lwd, main, xlab, ylab,
                             show_reference, ...) {

  # Determine sample size range
  if (is.null(n_range)) {
    if ("n2" %in% names(x)) {
      # Both object shapes carry n2, so both use a window around it. The lower
      # end is held at two rather than at a fixed floor of ten, which would
      # invert the window for a design of fewer than twenty per group.
      current_n2 <- x$n2
      lo <- max(2, floor(current_n2 * 0.5))
      hi <- max(lo + 1, ceiling(current_n2 * 1.5))
      n_range <- c(lo, hi)
    } else {
      n_range <- c(50, 200)
    }
  }

  # Generate sequence of n2 values
  n2_seq <- seq(n_range[1], n_range[2], length.out = n_points)
  n2_seq <- round(n2_seq)
  n2_seq <- unique(n2_seq)  # Remove duplicates after rounding

  # Calculate allocation ratio
  if ("r" %in% names(x)) {
    r <- x$r
  } else if ("n1" %in% names(x) && "n2" %in% names(x)) {
    r <- x$n1 / x$n2
  } else {
    r <- 1  # Default to balanced
  }

  n1_seq <- ceiling(r * n2_seq)

  # Calculate power for each sample size
  power_seq <- numeric(length(n2_seq))

  for (i in seq_along(n2_seq)) {
    power_result <- calculate_power_for_plot(x, endpoint_type, n1_seq[i], n2_seq[i])
    power_seq[i] <- power_result$powerCoprimary
  }

  # Create plot
  if (is.null(main)) {
    main <- "Power as a Function of Sample Size"
  }
  if (is.null(xlab)) {
    xlab <- if (abs(r - 1) < 0.01) {
      "Sample Size per Group (n)"
    } else {
      bquote("Control Group Sample Size (" ~ n[2] ~ "), r = " ~ .(r))
    }
  }
  if (is.null(ylab)) {
    ylab <- "Power"
  }

  plot(n2_seq, power_seq, type = "l", col = col, lwd = lwd,
       main = main, xlab = xlab, ylab = ylab,
       ylim = c(0, 1), las = 1, ...)

  # Add grid
  grid(col = "lightgray", lty = "dotted")

  # Add reference lines
  if (show_reference) {
    # Horizontal line at target power or current power
    if ("beta" %in% names(x)) {
      target_power <- 1 - x$beta
      abline(h = target_power, col = "red", lty = 2)
      text(n_range[1], target_power, sprintf("Target: %.1f%%", target_power * 100),
           pos = 3, col = "red", cex = 0.8)
    } else if ("powerCoprimary" %in% names(x)) {
      abline(h = x$powerCoprimary, col = "darkgreen", lty = 2)
    }

    # Vertical line at current n2. A power object carries powerCoprimary and
    # a sample size object does not, so the branches are selected on that
    # column rather than on n2, which both shapes have.
    if ("powerCoprimary" %in% names(x)) {
      abline(v = x$n2, col = "darkgreen", lty = 2)
      points(x$n2, x$powerCoprimary, pch = 19, col = "darkgreen", cex = 1.5)
      # Calculate total N
      total_N <- x$n1 + x$n2
      # Create label with allocation ratio and total N
      if (abs(r - 1) < 0.01) {
        # Balanced design
        label_text <- bquote(n == .(x$n2) * ',' ~ N == .(total_N) * ',' ~
                               power == .(sprintf("%.3f", x$powerCoprimary)))
      } else {
        # Unbalanced design
        label_text <- bquote(
          atop(
            n[2] == .(x$n2) * ',' ~ r == .(sprintf("%.2f", r)) * ',',
            N == .(total_N) * ',' ~ power == .(sprintf("%.3f", x$powerCoprimary))
          )
        )
      }
      text(x$n2, x$powerCoprimary, label_text,
           pos = 4, col = "darkgreen", cex = 0.8)
    } else if ("n2" %in% names(x)) {
      current_n2 <- x$n2
      # Find corresponding power
      idx <- which.min(abs(n2_seq - current_n2))
      abline(v = current_n2, col = "red", lty = 2)
      points(current_n2, power_seq[idx], pch = 19, col = "red", cex = 1.5)
      # Calculate total N
      total_N <- if ("N" %in% names(x)) x$N else x$n1 + x$n2
      # Create label with allocation ratio and total N
      if (abs(r - 1) < 0.01) {
        # Balanced design
        label_text <- bquote(n == .(current_n2) * ',' ~ N == .(total_N) * ',' ~
                               power == .(sprintf("%.3f", power_seq[idx])))
      } else {
        # Unbalanced design
        label_text <- bquote(
          atop(
            n[2] == .(current_n2) * ',' ~ r == .(sprintf("%.2f", r)) * ',',
            N == .(total_N) * ',' ~ power == .(sprintf("%.3f", power_seq[idx]))
          )
        )
      }
      text(current_n2, power_seq[idx], label_text,
           pos = 4, col = "red", cex = 0.8)
    }
  }

  # Return data
  plot_data <- data.frame(n1 = n1_seq, n2 = n2_seq, power = power_seq)
  return(plot_data)
}


# ============================================================================
# Helper function: Plot sample size vs correlation
# ============================================================================

plot_sample_size_rho <- function(x, endpoint_type, n_points, rho_range,
                                 col, lwd, main, xlab, ylab,
                                 show_reference, ...) {

  # A power object carries neither r, beta nor N, so derive each from the
  # columns it does have: the realized allocation, the achieved power as the
  # target, and the total of the two group sizes.
  r <- if ("r" %in% names(x)) x$r else x$n1 / x$n2
  beta <- if ("beta" %in% names(x)) {
    x$beta
  } else {
    1 - x$powerCoprimary
  }
  total_N <- if ("N" %in% names(x)) x$N else x$n1 + x$n2

  # (4) Determine correlation range based on endpoint type
  if (is.null(rho_range)) {
    if (endpoint_type == "binary") {
      # Calculate correlation bounds for binary endpoints
      bounds1 <- corrbound2Binary(x$p11, x$p12)
      bounds2 <- corrbound2Binary(x$p21, x$p22)
      rho_min <- max(bounds1[1], bounds2[1])
      rho_max <- min(bounds1[2], bounds2[2])
      # Add small margin
      rho_min <- max(0, rho_min + 0.05)
      rho_max <- min(0.95, rho_max - 0.05)
      rho_range <- seq(rho_min, rho_max, length.out = n_points)
    } else if (endpoint_type == "mixed_count_cont") {
      # Calculate correlation bounds for mixed count-continuous endpoints
      lambda1 <- x$r1 * x$t
      lambda2 <- x$r2 * x$t
      bounds1 <- corrbound2MixedCountContinuous(lambda1, x$nu, x$mu1, x$sd)
      bounds2 <- corrbound2MixedCountContinuous(lambda2, x$nu, x$mu2, x$sd)
      rho_min <- max(bounds1[1], bounds2[1])
      rho_max <- min(bounds1[2], bounds2[2])
      # Add small margin
      rho_min <- max(-0.95, rho_min + 0.05)
      rho_max <- min(0.95, rho_max - 0.05)
      rho_range <- seq(rho_min, rho_max, length.out = n_points)
    } else {
      # Default range for continuous and mixed continuous-binary
      rho_range <- seq(0, 0.9, length.out = n_points)
    }
  }

  # Calculate sample size for each correlation
  n2_seq <- numeric(length(rho_range))

  for (i in seq_along(rho_range)) {
    ss_result <- calculate_sample_size_for_plot(x, endpoint_type, rho_range[i],
                                               r = r, beta = beta)
    n2_seq[i] <- ss_result$n2
  }

  # Create plot
  if (is.null(main)) {
    main <- "Required Sample Size as a Function of Correlation"
  }
  if (is.null(xlab)) {
    xlab <- "Correlation (rho)"
  }
  if (is.null(ylab)) {
    ylab <- if (abs(r - 1) < 0.01) {
      "Sample Size per Group (n)"
    } else {
      # Include allocation ratio in the label
      bquote("Control Group Sample Size (" ~ n[2] ~ "), r = " ~ .(r))
    }
  }

  plot(rho_range, n2_seq, type = "l", col = col, lwd = lwd,
       main = main, xlab = xlab, ylab = ylab, las = 1, ...)

  # Add grid
  grid(col = "lightgray", lty = "dotted")

  # Add reference lines
  if (show_reference && "rho" %in% names(x)) {
    abline(v = x$rho, col = "darkgreen", lty = 2)
    abline(h = x$n2, col = "darkgreen", lty = 2)
    points(x$rho, x$n2, pch = 19, col = "darkgreen", cex = 1.5)
    # Calculate total N and create label
    if (abs(r - 1) < 0.01) {
      # Balanced design
      label_text <- bquote(rho == .(sprintf("%.2f", x$rho)) * ',' ~ n == .(x$n2) * ',' ~
                             N == .(total_N))
    } else {
      # Unbalanced design
      label_text <- bquote(
        atop(
          rho == .(sprintf("%.2f", x$rho)) * ',' ~ n[2] == .(x$n2) * ',',
          r == .(sprintf("%.2f", r)) * ',' ~ N == .(total_N)
        )
      )
    }
    text(x$rho, max(n2_seq) * 0.95, label_text,
         pos = 4, col = "darkgreen", cex = 0.8)
  } else if (show_reference && "rho1" %in% names(x) && "rho2" %in% names(x)) {
    # For binary endpoints with two correlations
    avg_rho <- (x$rho1 + x$rho2) / 2
    abline(v = avg_rho, col = "darkgreen", lty = 2)
    abline(h = x$n2, col = "darkgreen", lty = 2)
    points(avg_rho, x$n2, pch = 19, col = "darkgreen", cex = 1.5)
    # Calculate total N and create label
    if (abs(r - 1) < 0.01) {
      # Balanced design
      label_text <- bquote(rho == .(sprintf("%.2f", avg_rho)) * ',' ~ n == .(x$n2) * ',' ~
                             N == .(total_N))
    } else {
      # Unbalanced design
      label_text <- bquote(
        atop(
          rho == .(sprintf("%.2f", avg_rho)) * ',' ~ n[2] == .(x$n2) * ',',
          r == .(sprintf("%.2f", r)) * ',' ~ N == .(total_N)
        )
      )
    }
    text(avg_rho, max(n2_seq) * 0.95, label_text,
         pos = 4, col = "darkgreen", cex = 0.8)
  }

  # Return data
  plot_data <- data.frame(rho = rho_range, n2 = n2_seq)
  return(plot_data)
}


# ============================================================================
# Helper function: Plot effect size contour
# ============================================================================

plot_effect_contour <- function(x, endpoint_type, n_points,
                                main, xlab, ylab, ...) {

  if (endpoint_type != "continuous") {
    stop("effect_contour plot is currently only available for continuous endpoints")
  }

  # Grid of standardized effect sizes, which is what the axis labels below
  # state. The power is evaluated at the corresponding mean differences, so the
  # contour is the same picture whatever the two standard deviations are. When
  # both are one the two scales coincide.
  delta1_range <- seq(0.2, 1.0, length.out = n_points)
  delta2_range <- seq(0.2, 1.0, length.out = n_points)

  # Create grid
  grid <- expand.grid(delta1 = delta1_range, delta2 = delta2_range)

  # Calculate power for each combination
  grid$power <- numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    temp_x <- x
    temp_x$delta1 <- grid$delta1[i] * x$sd1
    temp_x$delta2 <- grid$delta2[i] * x$sd2
    power_result <- calculate_power_for_plot(temp_x, endpoint_type,
                                             x$n1, x$n2)
    grid$power[i] <- power_result$powerCoprimary
  }

  # Create contour plot
  power_matrix <- matrix(grid$power, nrow = n_points, ncol = n_points)

  if (is.null(main)) {
    main <- "Power Contours for Different Effect Sizes"
  }
  if (is.null(xlab)) {
    # (3) Fix character encoding for subscripts
    xlab <- expression(paste("Effect Size for Endpoint 1 (", delta[1], "/", sigma[1], ")"))
  }
  if (is.null(ylab)) {
    # (3) Fix character encoding for subscripts
    ylab <- expression(paste("Effect Size for Endpoint 2 (", delta[2], "/", sigma[2], ")"))
  }

  contour(delta1_range, delta2_range, power_matrix,
          levels = seq(0.1, 0.9, by = 0.1),
          labels = paste0(seq(10, 90, by = 10), "%"),
          main = main, xlab = xlab, ylab = ylab,
          col = "steelblue", lwd = 1.5, ...)

  # Add current point if available, on the standardized scale of the axes
  if ("delta1" %in% names(x) && "delta2" %in% names(x)) {
    points(x$delta1 / x$sd1, x$delta2 / x$sd2, pch = 19, col = "red", cex = 1.5)
  }

  # Return data
  return(grid)
}


# ============================================================================
# Helper function: Calculate power for plotting
# ============================================================================

calculate_power_for_plot <- function(x, endpoint_type, n1, n2) {

  if (endpoint_type == "continuous") {
    power2Continuous(
      n1 = n1, n2 = n2,
      delta1 = x$delta1, delta2 = x$delta2,
      sd1 = x$sd1, sd2 = x$sd2,
      rho = x$rho, alpha = x$alpha,
      known_var = ifelse("known_var" %in% names(x), x$known_var, TRUE),
      nMC = ifelse("nMC" %in% names(x) && !is.na(x$nMC), x$nMC, 1000)
    )
  } else if (endpoint_type == "binary") {
    # The Test column decides which family the object came from; routing an
    # exact test name into the approximate function would fail.
    if (as.character(x$Test) %in% EXACT_BINARY_TESTS) {
      power2BinaryExact(
        n1 = n1, n2 = n2,
        p11 = x$p11, p12 = x$p12,
        p21 = x$p21, p22 = x$p22,
        rho1 = x$rho1, rho2 = x$rho2,
        alpha = x$alpha,
        Test = as.character(x$Test),
        n_grid = if ("n_grid" %in% names(x)) x$n_grid else 100
      )
    } else {
      power2BinaryApprox(
        n1 = n1, n2 = n2,
        p11 = x$p11, p12 = x$p12,
        p21 = x$p21, p22 = x$p22,
        rho1 = x$rho1, rho2 = x$rho2,
        alpha = x$alpha,
        Test = as.character(x$Test)
      )
    }
  } else if (endpoint_type == "mixed_cont_binary") {
    power2MixedContinuousBinary(
      n1 = n1, n2 = n2,
      delta = x$delta, sd = x$sd,
      p1 = x$p1, p2 = x$p2,
      rho = x$rho, alpha = x$alpha,
      Test = x$Test,
      nMC = ifelse("nMC" %in% names(x) && !is.na(x$nMC), x$nMC, 1000)
    )
  } else if (endpoint_type == "mixed_count_cont") {
    power2MixedCountContinuous(
      n1 = n1, n2 = n2,
      r1 = x$r1, r2 = x$r2,
      nu = x$nu, t = x$t,
      mu1 = x$mu1, mu2 = x$mu2, sd = x$sd,
      rho1 = x$rho1, rho2 = x$rho2,
      alpha = x$alpha
    )
  }
}


# ============================================================================
# Helper function: Calculate sample size for plotting
# ============================================================================

calculate_sample_size_for_plot <- function(x, endpoint_type, rho_value,
                                           r = NULL, beta = NULL) {

  # r and beta are supplied by the caller because a power object has neither
  if (is.null(r)) r <- if ("r" %in% names(x)) x$r else x$n1 / x$n2
  if (is.null(beta)) {
    beta <- if ("beta" %in% names(x)) x$beta else 1 - x$powerCoprimary
  }

  if (endpoint_type == "continuous") {
    ss2Continuous(
      delta1 = x$delta1, delta2 = x$delta2,
      sd1 = x$sd1, sd2 = x$sd2,
      rho = rho_value, r = r,
      alpha = x$alpha, beta = beta,
      known_var = ifelse("known_var" %in% names(x), x$known_var, TRUE),
      nMC = ifelse("nMC" %in% names(x) && !is.na(x$nMC), x$nMC, 1000)
    )
  } else if (endpoint_type == "binary") {
    if (as.character(x$Test) %in% EXACT_BINARY_TESTS) {
      ss2BinaryExact(
        p11 = x$p11, p12 = x$p12,
        p21 = x$p21, p22 = x$p22,
        rho1 = rho_value, rho2 = rho_value,
        r = r, alpha = x$alpha, beta = beta,
        Test = as.character(x$Test),
        n_grid = if ("n_grid" %in% names(x)) x$n_grid else 100
      )
    } else {
      ss2BinaryApprox(
        p11 = x$p11, p12 = x$p12,
        p21 = x$p21, p22 = x$p22,
        rho1 = rho_value, rho2 = rho_value,
        r = r, alpha = x$alpha, beta = beta,
        Test = as.character(x$Test)
      )
    }
  } else if (endpoint_type == "mixed_cont_binary") {
    ss2MixedContinuousBinary(
      delta = x$delta, sd = x$sd,
      p1 = x$p1, p2 = x$p2,
      rho = rho_value, r = r,
      alpha = x$alpha, beta = beta,
      Test = as.character(x$Test),
      nMC = ifelse("nMC" %in% names(x) && !is.na(x$nMC), x$nMC, 1000)
    )
  } else if (endpoint_type == "mixed_count_cont") {
    ss2MixedCountContinuous(
      r1 = x$r1, r2 = x$r2,
      nu = x$nu, t = x$t,
      mu1 = x$mu1, mu2 = x$mu2, sd = x$sd,
      rho1 = rho_value, rho2 = rho_value,
      r = r, alpha = x$alpha, beta = beta
    )
  }
}
