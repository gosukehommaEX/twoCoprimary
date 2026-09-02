# Verify the 1.1.1 repairs, and confirm that nothing the article reports has
# moved. Run AFTER Rcpp::compileAttributes(), devtools::document(), Build and
# Install, and a restart of the R session.
#
# Run with:  source("dev/verify_1_1_1_fixes.R")

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(out_dir, "verify_fixes.log")
if (file.exists(log_path)) file.remove(log_path)

say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_path, append = TRUE)
}

pass <- 0L
fail <- 0L
check <- function(id, ok, detail = "") {
  if (isTRUE(ok)) pass <<- pass + 1L else fail <<- fail + 1L
  say(sprintf("[%s] %-34s %s", if (isTRUE(ok)) "PASS" else "FAIL", id, detail))
}

runs <- function(expr) {
  tryCatch({
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    eval(expr)
    TRUE
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
}

say("twoCoprimary version : ", as.character(utils::packageVersion("twoCoprimary")))
say("run at               : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")
check("version.is.1.1.1",
      as.character(utils::packageVersion("twoCoprimary")) == "1.1.1",
      paste("installed:", as.character(utils::packageVersion("twoCoprimary"))))
say("")

# ---------------------------------------------------------------------------
# The calls that failed in 1.1.0
# ---------------------------------------------------------------------------
ss_cont <- ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
                         rho = 0.5, r = 1, alpha = 0.025, beta = 0.2,
                         known_var = TRUE)
pw_cont <- power2Continuous(n1 = 100, n2 = 100, delta1 = 0.5, delta2 = 0.5,
                            sd1 = 1, sd2 = 1, rho = 0.5, alpha = 0.025,
                            known_var = TRUE)
ss_bin_ap <- ss2BinaryApprox(p11 = 0.40, p12 = 0.35, p21 = 0.25, p22 = 0.20,
                             rho1 = 0.5, rho2 = 0.5, r = 1,
                             alpha = 0.025, beta = 0.2, Test = "AN")
ss_mcb <- ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.60, p2 = 0.40,
                                   rho = 0.5, r = 1, alpha = 0.025,
                                   beta = 0.2, Test = "AN")
ss_mcc <- ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                                  mu1 = -50, mu2 = 0, sd = 250,
                                  rho1 = 0.5, rho2 = 0.5, r = 1,
                                  alpha = 0.025, beta = 0.2)
ss_bin_ex <- ss2BinaryExact(p11 = 0.60, p12 = 0.50, p21 = 0.40, p22 = 0.30,
                            rho1 = 0.3, rho2 = 0.3, r = 1,
                            alpha = 0.025, beta = 0.2, Test = "Fisher")
pw_bin_ex <- power2BinaryExact(n1 = 60, n2 = 60, p11 = 0.60, p12 = 0.50,
                               p21 = 0.40, p22 = 0.30, rho1 = 0.3, rho2 = 0.3,
                               alpha = 0.025, Test = "Fisher")

for (nm in c("continuous", "binary_approx", "mixed_cont_binary", "mixed_count_cont")) {
  obj <- switch(nm, continuous = ss_cont, binary_approx = ss_bin_ap,
                mixed_cont_binary = ss_mcb, mixed_count_cont = ss_mcc)
  o <- runs(bquote(plot(.(obj), type = "power_curve")))
  check(paste0("P1.", nm), isTRUE(o), if (isTRUE(o)) "" else o)
}
o <- runs(quote(plot(pw_cont, type = "sample_size_rho")))
check("P2.power.sample_size_rho", isTRUE(o), if (isTRUE(o)) "" else o)

o <- runs(quote(plot(ss_bin_ex, type = "power_curve", n_points = 5)))
check("P3.ss2BinaryExact", isTRUE(o), if (isTRUE(o)) "" else o)
o <- runs(quote(plot(pw_bin_ex, type = "power_curve", n_points = 5)))
check("P4.power2BinaryExact", isTRUE(o), if (isTRUE(o)) "" else o)

o <- runs(quote(plot(ss1Continuous(delta = 0.5, sd = 1, r = 1,
                                   alpha = 0.025, beta = 0.2))))
check("P5.ss1.message", is.character(o) && grepl("single endpoint", o),
      if (is.character(o)) substr(o, 1, 70) else "no error raised")

o <- runs(quote(power2BinaryApprox(n1 = 100, n2 = 100, p11 = 0.60, p12 = 0.50,
                                   p21 = 0.40, p22 = 0.30, rho1 = 0.3,
                                   rho2 = 0.3, alpha = 0.025, Test = "Fisher")))
check("V1.Test.validated", is.character(o) && grepl("Test must be one of", o),
      if (is.character(o)) substr(o, 1, 60) else "no error raised")

o <- runs(quote(ss2MixedCountContinuous(r1 = 1.25, r2 = 1.0, nu = 0.8, t = 1,
                                        mu1 = 0, mu2 = -50, sd = 250,
                                        rho1 = 0.5, rho2 = 0.5, r = 1,
                                        alpha = 0.025, beta = 0.2)))
check("N4.direction.guard", is.character(o) && grepl("r1 must be less", o),
      if (is.character(o)) substr(o, 1, 60) else "no error raised")

o <- tryCatch({
  v <- ss2Continuous(delta1 = 4, delta2 = 4, sd1 = 1, sd2 = 1, rho = 0.5,
                     r = 1, alpha = 0.025, beta = 0.2, known_var = FALSE)
  paste0("n2 = ", format(v$n2))
}, error = function(e) paste0("ERROR: ", conditionMessage(e)))
check("N5.large.effect", !grepl("^ERROR", o), o)
say("")

# ---------------------------------------------------------------------------
# The arcsine corrections
# ---------------------------------------------------------------------------
arcsine_power <- function(p1, p2, n1, n2, alpha) {
  se <- 0.5 * sqrt(1 / n1 + 1 / n2)
  stats::pnorm((asin(sqrt(p1)) - asin(sqrt(p2))) / se - stats::qnorm(1 - alpha))
}
say("--- ss1BinaryApprox, target power 0.90 ---")
rows <- list()
for (tst in c("AN", "AS", "ASc")) {
  for (rv in c(1, 2, 0.5)) {
    res <- ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = rv, alpha = 0.025,
                           beta = 0.1, Test = tst)
    ap <- arcsine_power(0.6, 0.4, res$n1, res$n2, 0.025)
    say(sprintf("  Test = %-4s r = %.1f  n1 = %4s  n2 = %4s  arcsine power = %.4f",
                tst, rv, format(res$n1), format(res$n2), ap))
    rows[[length(rows) + 1]] <- data.frame(Test = tst, r = rv, n1 = res$n1,
                                           n2 = res$n2, power = ap)
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(out_dir, "verify_ss1BinaryApprox.csv"), row.names = FALSE)
as_ok <- all(abs(tab$power[tab$Test == "AS"] - 0.90) < 0.03)
check("N1.arcsine.on.target", as_ok,
      paste("AS achieved power:", paste(sprintf("%.4f", tab$power[tab$Test == "AS"]),
                                        collapse = ", ")))
asc_ok <- all(tab$n2[tab$Test == "ASc"] >= tab$n2[tab$Test == "AS"])
check("N2.ASc.conservative", asc_ok,
      paste("AS n2:", paste(tab$n2[tab$Test == "AS"], collapse = ", "),
            "| ASc n2:", paste(tab$n2[tab$Test == "ASc"], collapse = ", ")))
say("")
say("  NOTE: at r = 1 the AS result must be unchanged from 1.1.0, which gave")
say("        n1 = n2 = 130. ASc at r = 1 gave 122 in 1.1.0 and should now be")
say("        at least 130.")
say("")

# ---------------------------------------------------------------------------
# Nothing the article reports may move. These are the exact grids the
# manuscript uses in its Validation section.
# ---------------------------------------------------------------------------
say("--- article validation tables, recomputed ---")

val_cont <- design_table(
  param_grid = expand.grid(delta1 = c(0.20, 0.25, 0.30, 0.35, 0.40),
                           delta2 = c(0.20, 0.25, 0.30, 0.35, 0.40),
                           sd1 = 1, sd2 = 1),
  rho_values = c(0, 0.3, 0.5, 0.8), r = 1, alpha = 0.025, beta = 0.2,
  endpoint_type = "continuous", known_var = TRUE)
write.csv(val_cont, file.path(out_dir, "verify_val_continuous.csv"), row.names = FALSE)
say("  continuous grid written (Sozu 2011 Table 1)")

val_bin <- do.call(rbind, lapply(c("AN", "ANc", "AS", "ASc"), function(tst) {
  design_table(
    param_grid = tibble::tibble(p11 = c(0.70, 0.87, 0.90, 0.95),
                                p12 = c(0.70, 0.70, 0.90, 0.95),
                                p21 = c(0.50, 0.70, 0.70, 0.90),
                                p22 = c(0.50, 0.50, 0.70, 0.90)),
    rho_values = c(-0.5, -0.3, 0, 0.3, 0.5, 0.8),
    r = 1, alpha = 0.025, beta = 0.2,
    endpoint_type = "binary", Test = tst) |> dplyr::mutate(Test = tst)
}))
write.csv(val_bin, file.path(out_dir, "verify_val_binary_approx.csv"), row.names = FALSE)
say("  binary asymptotic grid written (Sozu 2010 Table III)")

val_exact <- do.call(rbind, lapply(c("Chisq", "Fisher", "Z-pool", "Boschloo"),
                                   function(tst) {
  do.call(rbind, lapply(1:2, function(rv) {
    design_table(
      param_grid = tibble::tibble(p11 = 0.54, p12 = 0.54, p21 = 0.25, p22 = 0.25),
      rho_values = c(0, 0.3, 0.5, 0.8),
      r = rv, alpha = 0.025, beta = 0.1,
      endpoint_type = "binary", Test = tst) |>
      dplyr::mutate(r = rv, Test = tst)
  }))
}))
write.csv(val_exact, file.path(out_dir, "verify_val_binary_exact.csv"), row.names = FALSE)
say("  binary exact grid written (Homma and Yoshida 2025 Table 4)")

val_cb <- design_table(
  param_grid = expand.grid(delta = 4.4, sd = c(19, 20, 21, 22),
                           p1 = 0.59, p2 = 0.46),
  rho_values = c(0, 0.3, 0.5, 0.8), r = 1, alpha = 0.025, beta = 0.2,
  endpoint_type = "mixed_cont_binary", Test = "AN")
write.csv(val_cb, file.path(out_dir, "verify_val_mixed_cb.csv"), row.names = FALSE)
say("  mixed continuous-binary grid written (Sozu 2012 Table 2)")

val_cc <- do.call(rbind, lapply(c(3, 5), function(nu_val) {
  do.call(rbind, lapply(c(0, 0.2, 0.4, 0.6, 0.8), function(rho_val) {
    res <- ss2MixedCountContinuous(r1 = 1, r2 = 2, nu = nu_val, t = 1,
                                   mu1 = -50, mu2 = 0, sd = 75,
                                   rho1 = rho_val, rho2 = rho_val, r = 1,
                                   alpha = 0.025, beta = 0.1)
    data.frame(nu = nu_val, rho = rho_val, n2 = res$n2, N = res$N)
  }))
}))
write.csv(val_cc, file.path(out_dir, "verify_val_mixed_cc.csv"), row.names = FALSE)
say("  mixed count-continuous grid written (Homma and Yoshida 2024 Table 1)")

usage <- design_table(
  param_grid = expand.grid(p11 = c(0.35, 0.40, 0.45), p12 = 0.35,
                           p21 = 0.25, p22 = 0.20),
  rho_values = c(0, 0.3, 0.5, 0.8), r = 1, alpha = 0.025, beta = 0.2,
  endpoint_type = "binary", Test = "AN")
write.csv(usage, file.path(out_dir, "verify_usage_design_table.csv"), row.names = FALSE)
say("  usage example design_table written")
say("")

say("=== summary ===")
say("  PASS : ", pass)
say("  FAIL : ", fail)
cat("\nDone. See", log_path, "\n")
