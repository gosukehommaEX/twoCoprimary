# Work out the correlation-estimation examples that will be added to the
# Usage examples section of the R Journal article (Reviewer 4's request for an
# example of preliminary data used to derive a correlation input).
#
# Run with:  source("dev/article_correlation_examples.R")
#
# Everything is written to dev/out/. Console output need not be pasted back.

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(out_dir, "correlation_examples.log")
if (file.exists(log_path)) file.remove(log_path)

say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_path, append = TRUE)
}

WIDTH <- 77
capture_at_width <- function(expr, tag) {
  old <- options(width = WIDTH)
  on.exit(options(old), add = TRUE)
  lines <- utils::capture.output(eval(expr))
  path <- file.path(out_dir, paste0("corr_", tag, ".txt"))
  writeLines(lines, path)
  say(sprintf("  %-26s lines = %2d  max width = %2d  fits = %s",
              tag, length(lines), max(nchar(lines)),
              if (max(nchar(lines)) <= WIDTH + 3) "YES" else "NO"))
  invisible(lines)
}

say("twoCoprimary version : ", as.character(utils::packageVersion("twoCoprimary")))
say("run at               : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# ---------------------------------------------------------------------------
# A. Two continuous endpoints: sample correlation from prior individual data.
# Fifteen patients from a hypothetical natural history cohort, change from
# baseline on a cognitive and a functional scale.
# ---------------------------------------------------------------------------
say("--- A. two continuous endpoints ---")
adas <- c(-18.6, -9.3, -24.0, -4.2, -16.5, -6.6, -21.3, -14.4,
          -2.7, -20.7, -11.1, -15.0, -8.4, -22.8, -12.6)
adl <- c(-12.3, -4.2, -15.6, -9.0, -21.3, -15.6, -17.7, -10.8,
         1.5, -24.9, -10.8, -13.2, -14.7, -12.9, -18.3)

rho_hat_cont <- cor(adas, adl)
say(sprintf("  n                 = %d", length(adas)))
say(sprintf("  sd(adas)          = %.3f", sd(adas)))
say(sprintf("  sd(adl)           = %.3f", sd(adl)))
say(sprintf("  Pearson rho_hat   = %.6f  (rounded: %.2f)",
            rho_hat_cont, round(rho_hat_cont, 2)))

ct <- cor.test(adas, adl)
say(sprintf("  95%% CI            = [%.4f, %.4f]", ct$conf.int[1], ct$conf.int[2]))
say("")

say("  sample size at the point estimate and at the confidence limits:")
cont_rows <- list()
for (lab in c("lower CI", "point", "upper CI", "independence")) {
  rv <- switch(lab, "lower CI" = ct$conf.int[1], "point" = rho_hat_cont,
               "upper CI" = ct$conf.int[2], "independence" = 0)
  res <- ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
                       rho = rv, r = 1, alpha = 0.025, beta = 0.2,
                       known_var = TRUE)
  say(sprintf("    %-13s rho = %6.3f   n per group = %3s   N = %3s",
              lab, rv, format(res$n2), format(res$N)))
  cont_rows[[length(cont_rows) + 1]] <- data.frame(
    scenario = lab, rho = rv, n = res$n2, N = res$N)
}
write.csv(do.call(rbind, cont_rows),
          file.path(out_dir, "corr_continuous_sensitivity.csv"), row.names = FALSE)
say("")

# ---------------------------------------------------------------------------
# B. Two binary endpoints: correlation from a joint 2 x 2 table of responses
# in a prior study, in the sense of Prentice (1988).
# ---------------------------------------------------------------------------
say("--- B. two binary endpoints ---")
prior <- matrix(c(31, 17, 11, 61), nrow = 2, byrow = TRUE,
                dimnames = list(c("EP1 response", "EP1 no response"),
                                c("EP2 response", "EP2 no response")))
capture_at_width(quote(print(prior)), "prior_table")

n_prior <- sum(prior)
p1_hat <- sum(prior[1, ]) / n_prior
p2_hat <- sum(prior[, 1]) / n_prior
phi_hat <- prior[1, 1] / n_prior
rho_hat_bin <- (phi_hat - p1_hat * p2_hat) /
  sqrt(p1_hat * (1 - p1_hat) * p2_hat * (1 - p2_hat))

say(sprintf("  n                 = %d", n_prior))
say(sprintf("  p1_hat            = %.4f", p1_hat))
say(sprintf("  p2_hat            = %.4f", p2_hat))
say(sprintf("  phi_hat           = %.4f", phi_hat))
say(sprintf("  rho_hat           = %.6f  (rounded: %.2f)",
            rho_hat_bin, round(rho_hat_bin, 2)))
say("")

say("  Frechet-Hoeffding bounds for the DESIGN probabilities:")
b_trt <- corrbound2Binary(p1 = 0.40, p2 = 0.35)
b_ctl <- corrbound2Binary(p1 = 0.25, p2 = 0.20)
say(sprintf("    treatment (0.40, 0.35) : [%.4f, %.4f]", b_trt[1], b_trt[2]))
say(sprintf("    control   (0.25, 0.20) : [%.4f, %.4f]", b_ctl[1], b_ctl[2]))
say(sprintf("    rho_hat admissible in both : %s",
            rho_hat_bin >= max(b_trt[1], b_ctl[1]) &&
              rho_hat_bin <= min(b_trt[2], b_ctl[2])))
say(sprintf("    binding bound is the control group upper limit : %.4f", b_ctl[2]))
say("")

say("  sample size across a range of plausible correlations:")
bin_rows <- list()
for (rv in c(0, 0.3, round(rho_hat_bin, 2), 0.7, round(b_ctl[2], 2) - 0.01)) {
  res <- ss2BinaryApprox(p11 = 0.40, p12 = 0.35, p21 = 0.25, p22 = 0.20,
                         rho1 = rv, rho2 = rv, r = 1,
                         alpha = 0.025, beta = 0.2, Test = "AN")
  say(sprintf("    rho = %5.2f   n per group = %3s   N = %3s",
              rv, format(res$n2), format(res$N)))
  bin_rows[[length(bin_rows) + 1]] <- data.frame(rho = rv, n = res$n2, N = res$N)
}
write.csv(do.call(rbind, bin_rows),
          file.path(out_dir, "corr_binary_sensitivity.csv"), row.names = FALSE)
say("")

# ---------------------------------------------------------------------------
# C. The two chunks exactly as they would appear in the article, so that the
# rendered console output can be measured before it is pasted in.
# ---------------------------------------------------------------------------
say("--- C. rendered output of the article chunks ---")
capture_at_width(quote(print(round(c(rho_hat = rho_hat_cont), 3))), "cont_rho")
capture_at_width(quote(print(corrbound2Binary(p1 = 0.25, p2 = 0.20))), "bounds_ctl")
capture_at_width(quote(print(ss2BinaryApprox(
  p11 = 0.40, p12 = 0.35, p21 = 0.25, p22 = 0.20,
  rho1 = round(rho_hat_bin, 2), rho2 = round(rho_hat_bin, 2), r = 1,
  alpha = 0.025, beta = 0.2, Test = "AN"))), "ss_at_rho_hat")
say("")

cat("\nDone. See", log_path, "\n")
