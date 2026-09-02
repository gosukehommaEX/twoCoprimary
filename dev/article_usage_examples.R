# Explore the design_table() and plot() output that will be added to the
# Usage examples section of the R Journal article.
#
# Run with:  source("dev/article_usage_examples.R")
#
# Writes everything to dev/out/. Nothing is printed that is not also written
# to a file, so the results can be shared by sending the files.

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(out_dir, "usage_examples.log")

if (file.exists(log_path)) file.remove(log_path)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_path, append = TRUE)
}

say("twoCoprimary version : ", as.character(utils::packageVersion("twoCoprimary")))
say("R version            : ", R.version.string)
say("run at               : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# The article renders console output inside a fixed-width verbatim block.
# Existing output lines in the current PDF reach 79 characters, so anything
# wider than about 80 will overflow the text block.
WIDTH_BUDGET <- 77

capture_print <- function(obj, width) {
  old <- options(width = width)
  on.exit(options(old), add = TRUE)
  utils::capture.output(print(obj))
}

report <- function(tag, lines, elapsed) {
  w <- max(nchar(lines))
  path <- file.path(out_dir, paste0("design_table_", tag, ".txt"))
  writeLines(lines, path)
  say("--- ", tag, " ---")
  say("  rows in output   : ", length(lines))
  say("  max line width   : ", w, "   (budget ", WIDTH_BUDGET + 3, ")")
  say("  fits             : ", if (w <= WIDTH_BUDGET + 3) "YES" else "NO")
  say("  elapsed seconds  : ", round(elapsed, 2))
  say("  written to       : ", path)
  say("")
}

run_one <- function(tag, expr) {
  t0 <- proc.time()[["elapsed"]]
  obj <- tryCatch(eval(expr), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - t0
  if (inherits(obj, "error")) {
    say("--- ", tag, " ---")
    say("  ERROR: ", conditionMessage(obj))
    say("")
    return(invisible(NULL))
  }
  report(tag, capture_print(obj, WIDTH_BUDGET), elapsed)
  invisible(obj)
}

# ---------------------------------------------------------------------------
# Option A: binary endpoints, vary the treatment probability on endpoint 1,
# four correlation values. This is the migraine scenario of the article with
# p11 taken over a plausible range.
# ---------------------------------------------------------------------------
grid_a <- expand.grid(
  p11 = c(0.35, 0.40, 0.45),
  p12 = 0.35, p21 = 0.25, p22 = 0.20
)
tab_a <- run_one("optionA_binary_4rho", quote(
  design_table(
    param_grid = grid_a,
    rho_values = c(0, 0.3, 0.5, 0.8),
    r = 1, alpha = 0.025, beta = 0.2,
    endpoint_type = "binary", Test = "AN"
  )
))

# ---------------------------------------------------------------------------
# Option B: same grid, three correlation values. Fallback if Option A is too
# wide for the verbatim block.
# ---------------------------------------------------------------------------
tab_b <- run_one("optionB_binary_3rho", quote(
  design_table(
    param_grid = grid_a,
    rho_values = c(0, 0.5, 0.8),
    r = 1, alpha = 0.025, beta = 0.2,
    endpoint_type = "binary", Test = "AN"
  )
))

# ---------------------------------------------------------------------------
# Option C: binary endpoints varying both treatment probabilities, four
# correlation values. Six rows, which shows the grid behaviour more clearly.
# ---------------------------------------------------------------------------
grid_c <- expand.grid(
  p11 = c(0.35, 0.40, 0.45),
  p12 = c(0.30, 0.35),
  p21 = 0.25, p22 = 0.20
)
tab_c <- run_one("optionC_binary_grid", quote(
  design_table(
    param_grid = grid_c,
    rho_values = c(0, 0.3, 0.5, 0.8),
    r = 1, alpha = 0.025, beta = 0.2,
    endpoint_type = "binary", Test = "AN"
  )
))

# ---------------------------------------------------------------------------
# Option D: continuous endpoints, in case the binary grid reads better as a
# companion to the Alzheimer's example instead.
# ---------------------------------------------------------------------------
grid_d <- expand.grid(
  delta1 = c(0.4, 0.5, 0.6),
  delta2 = 0.5, sd1 = 1, sd2 = 1
)
tab_d <- run_one("optionD_continuous", quote(
  design_table(
    param_grid = grid_d,
    rho_values = c(0, 0.3, 0.5, 0.8),
    r = 1, alpha = 0.025, beta = 0.2,
    endpoint_type = "continuous", known_var = TRUE
  )
))

# ---------------------------------------------------------------------------
# Figures. The article currently shows only type = "sample_size_rho".
# Reviewer 4 asked for figure generating functions that support protocol
# development, so the other two types are candidates for the new subsection.
# Saved at the same size the article uses for its existing figure.
# ---------------------------------------------------------------------------
result_cont <- ss2Continuous(
  delta1 = 0.5, delta2 = 0.5,
  sd1 = 1, sd2 = 1,
  rho = 0.5, r = 1,
  alpha = 0.025, beta = 0.2,
  known_var = TRUE
)

result_bin <- ss2BinaryApprox(
  p11 = 0.40, p12 = 0.35,
  p21 = 0.25, p22 = 0.20,
  rho1 = 0.5, rho2 = 0.5,
  r = 1, alpha = 0.025, beta = 0.2,
  Test = "AN"
)

save_plot <- function(tag, expr) {
  path <- file.path(out_dir, paste0("plot_", tag, ".png"))
  ok <- tryCatch({
    grDevices::png(path, width = 6, height = 4, units = "in", res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)
    eval(expr)
    TRUE
  }, error = function(e) {
    say("  ERROR in plot ", tag, ": ", conditionMessage(e))
    FALSE
  })
  say("plot ", tag, " : ", if (isTRUE(ok)) paste0("written to ", path) else "FAILED")
}

say("--- figures ---")
save_plot("effect_contour_cont", quote(plot(result_cont, type = "effect_contour")))
save_plot("power_curve_cont", quote(plot(result_cont, type = "power_curve")))
save_plot("sample_size_rho_bin", quote(plot(result_bin, type = "sample_size_rho")))
say("")

# ---------------------------------------------------------------------------
# Numbers the surrounding prose will quote, so the text can be checked
# against the data rather than against memory.
# ---------------------------------------------------------------------------
say("--- reference sample sizes (binary AN, r = 1, alpha = 0.025, beta = 0.2) ---")
ref <- do.call(rbind, lapply(c(0, 0.3, 0.5, 0.8), function(rho) {
  do.call(rbind, lapply(c(0.35, 0.40, 0.45), function(p11) {
    res <- ss2BinaryApprox(
      p11 = p11, p12 = 0.35, p21 = 0.25, p22 = 0.20,
      rho1 = rho, rho2 = rho, r = 1, alpha = 0.025, beta = 0.2,
      Test = "AN"
    )
    data.frame(p11 = p11, rho = rho, n1 = res$n1, n2 = res$n2, N = res$N)
  }))
}))
write.csv(ref, file.path(out_dir, "usage_examples_reference_N.csv"), row.names = FALSE)
for (i in seq_len(nrow(ref))) {
  say(sprintf("  p11 = %.2f  rho = %.1f  n1 = %s  n2 = %s  N = %s",
              ref$p11[i], ref$rho[i], format(ref$n1[i]), format(ref$n2[i]),
              format(ref$N[i])))
}
say("")

say("--- correlation bounds used in the article ---")
b1 <- corrbound2Binary(p1 = 0.40, p2 = 0.35)
b2 <- corrbound2Binary(p1 = 0.25, p2 = 0.20)
say(sprintf("  treatment group: [%.4f, %.4f]", b1[1], b1[2]))
say(sprintf("  control group  : [%.4f, %.4f]", b2[1], b2[2]))
say("")

cat("\nDone. See", log_path, "\n")
