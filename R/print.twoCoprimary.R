#' Print Method for twoCoprimary Objects
#'
#' Provides a clean, formatted display of sample size and power calculation results
#' from the twoCoprimary package.
#'
#' @param x An object of class "twoCoprimary"
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the original object
#'
#' @details
#' This print method provides a formatted output that displays key parameters
#' and results in an easy-to-read format. The specific format adapts to the
#' type of calculation (sample size vs power) and the type of endpoints involved.
#'
#' @export
print.twoCoprimary <- function(x, ...) {

  # Helper function to format numeric values
  # Integers are displayed without decimal points
  # Decimals are displayed with trailing zeros removed
  format_value <- function(val) {
    if (is.na(val)) {
      return("NA")
    }
    # Check if value is effectively an integer
    if (abs(val - round(val)) < 1e-10) {
      return(sprintf("%.0f", val))
    }
    # For decimal values, format with precision and remove trailing zeros
    formatted <- sprintf("%.6f", val)
    formatted <- sub("0+$", "", formatted)
    # Ensure at least one digit after decimal point
    if (grepl("\\.$", formatted)) {
      formatted <- paste0(formatted, "0")
    }
    formatted
  }

  # Helper function to format a vector of values
  format_values <- function(vals) {
    formatted <- sapply(vals, format_value)
    paste(formatted, collapse = ", ")
  }

  # Determine function type based on column names
  is_power <- "powerCoprimary" %in% names(x)
  is_sample_size <- "N" %in% names(x)

  # Determine endpoint types
  has_continuous <- any(c("delta", "delta1", "delta2", "sd", "sd1", "sd2") %in% names(x))
  has_binary <- any(c("p1", "p2", "p11", "p12", "p21", "p22") %in% names(x))
  has_count <- any(c("r1", "r2", "nu") %in% names(x))

  # Build title based on endpoint types
  if (is_power) {
    title <- "Power calculation for "
  } else {
    title <- "Sample size calculation for "
  }

  # Determine endpoint combination
  if (has_continuous && has_binary) {
    endpoint_type <- "mixed continuous and binary co-primary endpoints"
  } else if (has_continuous && has_count) {
    endpoint_type <- "mixed count and continuous co-primary endpoints"
  } else if (has_continuous && ("delta1" %in% names(x) || "delta2" %in% names(x))) {
    endpoint_type <- "two continuous co-primary endpoints"
  } else if (has_binary && ("p11" %in% names(x) || "p12" %in% names(x))) {
    endpoint_type <- "two binary co-primary endpoints"
  } else if (has_continuous) {
    endpoint_type <- "single continuous endpoint"
  } else if (has_binary) {
    endpoint_type <- "single binary endpoint"
  } else if (has_count) {
    endpoint_type <- "single count endpoint"
  } else {
    endpoint_type <- "co-primary endpoints"
  }

  cat("\n", title, endpoint_type, "\n\n", sep = "")

  # Print parameters based on endpoint type
  if ("n1" %in% names(x) && "n2" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "n1", format_value(x$n1)))
    cat(sprintf("%15s = %s\n", "n2", format_value(x$n2)))
  }

  if ("N" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "N", format_value(x$N)))
  }

  # Print endpoint-specific parameters
  if ("delta" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "delta", format_value(x$delta)))
  }
  if ("delta1" %in% names(x) && "delta2" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "delta", format_values(c(x$delta1, x$delta2))))
  }

  if ("sd" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "sd", format_value(x$sd)))
  }
  if ("sd1" %in% names(x) && "sd2" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "sd", format_values(c(x$sd1, x$sd2))))
  }

  if ("p1" %in% names(x) && "p2" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "p", format_values(c(x$p1, x$p2))))
  }
  if ("p11" %in% names(x) && "p12" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "p (group 1)", format_values(c(x$p11, x$p12))))
  }
  if ("p21" %in% names(x) && "p22" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "p (group 2)", format_values(c(x$p21, x$p22))))
  }

  if ("r1" %in% names(x) && "r2" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "rate", format_values(c(x$r1, x$r2))))
  }

  if ("nu" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "nu", format_value(x$nu)))
  }

  if ("t" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "t", format_value(x$t)))
  }

  if ("mu1" %in% names(x) && "mu2" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "mu", format_values(c(x$mu1, x$mu2))))
  }

  if ("rho" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "rho", format_value(x$rho)))
  }
  if ("rho1" %in% names(x) && "rho2" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "rho", format_values(c(x$rho1, x$rho2))))
  }

  if ("r" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "allocation", format_value(x$r)))
  }

  if ("alpha" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "alpha", format_value(x$alpha)))
  }

  if ("beta" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "beta", format_value(x$beta)))
  }

  if ("Test" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "Test", x$Test))
  }

  if ("known_var" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "known_var", x$known_var))
  }

  if ("nMC" %in% names(x) && !is.na(x$nMC)) {
    cat(sprintf("%15s = %s\n", "nMC", format_value(x$nMC)))
  }

  # Print power results
  if ("power1" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "power1", format_value(x$power1)))
  }
  if ("power2" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "power2", format_value(x$power2)))
  }
  if ("powerCont" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "powerCont", format_value(x$powerCont)))
  }
  if ("powerBin" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "powerBin", format_value(x$powerBin)))
  }
  if ("powerCount" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "powerCount", format_value(x$powerCount)))
  }
  if ("powerCoprimary" %in% names(x)) {
    cat(sprintf("%15s = %s\n", "powerCoprimary", format_value(x$powerCoprimary)))
  }

  cat("\n")

  invisible(x)
}
