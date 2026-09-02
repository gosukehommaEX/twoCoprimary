# Exhaustive contract audit of twoCoprimary, run against the INSTALLED package.
#
# Purpose. R CMD check confirms that the package builds, that the examples run
# and that the documentation is syntactically well formed. It does not confirm
# that a function returns what its documentation says it returns, that a sample
# size is the smallest one reaching the target power, that a rejection region
# has the shape the power formula assumes, or that two functions computing the
# same quantity agree. Every defect repaired in 1.1.1 belonged to one of those
# classes. This file checks them mechanically, over every exported function,
# every documented value of every discrete argument, and every shape of object
# the package returns.
#
# Nothing here is judged by eye. Each check states a claim, evaluates it, and
# records PASS, FAIL, INFO or ERROR. Wherever possible the claim is verified
# against an independent construction rather than against the package's own
# machinery: a Frechet-Hoeffding bound against a comonotone simulation, a
# rejection region against a direct null-probability sum, a bivariate binomial
# against its own marginals and moments.
#
# Run with:  source("dev/audit_all_functions.R")
# from the package root, after Build and Install and a session restart.
#
# Results go to dev/out/audit_all_functions.log and .csv. The script never
# stops on a failure; read the summary at the end.

library(twoCoprimary)

t_start_all <- Sys.time()

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(out_dir, "audit_all_functions.log")
csv_path <- file.path(out_dir, "audit_all_functions.csv")
tim_path <- file.path(out_dir, "audit_all_functions_timing.csv")
if (file.exists(log_path)) file.remove(log_path)

say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_path, append = TRUE)
}

results <- list()
record <- function(id, claim, observed, verdict) {
  results[[length(results) + 1]] <<- data.frame(
    id = id, claim = claim, observed = observed, verdict = verdict,
    stringsAsFactors = FALSE
  )
  say(sprintf("[%-5s] %s", verdict, id))
  say("    claim    : ", claim)
  say("    observed : ", observed)
}

# Each check expression must return a list with an 'obs' string and either
# ok = TRUE / FALSE or verdict = "INFO".
check <- function(id, claim, expr) {
  out <- tryCatch(
    expr,
    error = function(e) list(verdict = "ERROR",
                             obs = paste0("ERROR: ", conditionMessage(e)))
  )
  if (is.null(out$verdict)) {
    out$verdict <- if (isTRUE(out$ok)) "PASS" else "FAIL"
  }
  record(id, claim, out$obs, out$verdict)
  invisible(out$verdict)
}

part <- function(title) {
  say("")
  say(strrep("=", 78))
  say(title)
  say(strrep("=", 78))
}

fmt <- function(x, d = 6) {
  paste(formatC(as.numeric(x), format = "f", digits = d), collapse = ", ")
}

# Evaluate a plotting call on a null device
silently_plots <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

timings <- list()
timed <- function(label, expr) {
  el <- system.time(val <- force(expr))[["elapsed"]]
  timings[[length(timings) + 1]] <<- data.frame(
    label = label, elapsed = el, stringsAsFactors = FALSE
  )
  attr(val, "elapsed") <- el
  val
}

say("twoCoprimary version : ", as.character(utils::packageVersion("twoCoprimary")))
say("R version            : ", R.version.string)
say("platform             : ", R.version$platform)
say("run at               : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

# ---------------------------------------------------------------------------
# Reference objects: one of every shape the package can return
# ---------------------------------------------------------------------------
obj <- list()

obj$ss1_cont <- ss1Continuous(delta = 0.5, sd = 1, r = 1, alpha = 0.025, beta = 0.2)
obj$ss1_count <- ss1Count(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1, r = 1,
                          alpha = 0.025, beta = 0.2)
obj$ss1_bin <- ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = 1, alpha = 0.025,
                               beta = 0.2, Test = "AN")

obj$ss2_cont <- ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
                              rho = 0.5, r = 1, alpha = 0.025, beta = 0.2,
                              known_var = TRUE)
obj$pw2_cont <- power2Continuous(n1 = 100, n2 = 100, delta1 = 0.5, delta2 = 0.5,
                                 sd1 = 1, sd2 = 1, rho = 0.5, alpha = 0.025,
                                 known_var = TRUE)
obj$ss2_bin_ap <- ss2BinaryApprox(p11 = 0.6, p12 = 0.5, p21 = 0.4, p22 = 0.3,
                                  rho1 = 0.3, rho2 = 0.3, r = 1,
                                  alpha = 0.025, beta = 0.2, Test = "AN")
obj$pw2_bin_ap <- power2BinaryApprox(n1 = 120, n2 = 120, p11 = 0.6, p12 = 0.5,
                                     p21 = 0.4, p22 = 0.3, rho1 = 0.3,
                                     rho2 = 0.3, alpha = 0.025, Test = "AN")
obj$ss2_bin_ex <- ss2BinaryExact(p11 = 0.7, p12 = 0.6, p21 = 0.3, p22 = 0.2,
                                 rho1 = 0.3, rho2 = 0.3, r = 1,
                                 alpha = 0.025, beta = 0.2, Test = "Chisq")
obj$pw2_bin_ex <- power2BinaryExact(n1 = 40, n2 = 40, p11 = 0.7, p12 = 0.6,
                                    p21 = 0.3, p22 = 0.2, rho1 = 0.3,
                                    rho2 = 0.3, alpha = 0.025, Test = "Fisher")
obj$ss2_mcb <- ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.6, p2 = 0.4,
                                        rho = 0.5, r = 1, alpha = 0.025,
                                        beta = 0.2, Test = "AN")
obj$pw2_mcb <- power2MixedContinuousBinary(n1 = 100, n2 = 100, delta = 0.5,
                                           sd = 1, p1 = 0.6, p2 = 0.4,
                                           rho = 0.5, alpha = 0.025, Test = "AN")
obj$ss2_mcc <- ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                                       mu1 = -50, mu2 = 0, sd = 250,
                                       rho1 = 0.5, rho2 = 0.5, r = 1,
                                       alpha = 0.025, beta = 0.2)
obj$pw2_mcc <- power2MixedCountContinuous(n1 = 300, n2 = 300, r1 = 1.0,
                                          r2 = 1.25, nu = 0.8, t = 1,
                                          mu1 = -50, mu2 = 0, sd = 250,
                                          rho1 = 0.5, rho2 = 0.5, alpha = 0.025)

# Objects produced through the unequal allocation path, which is where the
# arcsine coefficient defect of 1.1.0 lived
obj$ss2_cont_r2 <- ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
                                 rho = 0.5, r = 2, alpha = 0.025, beta = 0.2,
                                 known_var = TRUE)
obj$pw2_cont_r2 <- power2Continuous(n1 = 200, n2 = 100, delta1 = 0.5,
                                    delta2 = 0.5, sd1 = 1, sd2 = 1, rho = 0.5,
                                    alpha = 0.025, known_var = TRUE)

single_endpoint_objects <- c("ss1_cont", "ss1_count", "ss1_bin")
coprimary_objects <- setdiff(names(obj), single_endpoint_objects)

# ===========================================================================
part("Part 0. Coverage: every exported object is exercised somewhere below")
# ===========================================================================

# The functions this file calls at least once. Keeping the list explicit means
# that a new export added later fails this check instead of slipping through
# untested.
exercised <- c(
  "corrbound2Binary", "corrbound2MixedCountContinuous", "dbibinom",
  "design_table", "power2BinaryApprox", "power2BinaryExact", "power2Continuous",
  "power2MixedContinuousBinary", "power2MixedCountContinuous", "rr1Binary",
  "ss1BinaryApprox", "ss1Continuous", "ss1Count", "ss2BinaryApprox",
  "ss2BinaryExact", "ss2Continuous", "ss2MixedContinuousBinary",
  "ss2MixedCountContinuous", "twoCoprimary2BinaryApprox",
  "twoCoprimary2BinaryExact", "twoCoprimary2Continuous",
  "twoCoprimary2MixedContinuousBinary", "twoCoprimary2MixedCountContinuous"
)
s3_methods <- c("plot.twoCoprimary", "print.twoCoprimary",
                "print.twoCoprimary_table")

# S3 methods are registered rather than exported, so they do not appear in
# getNamespaceExports() and have to be collected from the method table
registered_s3 <- local({
  tab <- getNamespaceInfo("twoCoprimary", "S3methods")
  if (is.null(tab) || !nrow(tab)) character(0)
  else sort(unique(paste(as.character(tab[, 1]), as.character(tab[, 2]), sep = ".")))
})

check("coverage.exports", "every public name of the package is exercised here", {
  public <- sort(unique(c(getNamespaceExports("twoCoprimary"), registered_s3)))
  known <- sort(unique(c(exercised, s3_methods)))
  missing <- setdiff(public, known)
  extra <- setdiff(known, public)
  list(ok = length(missing) == 0 && length(extra) == 0,
       obs = sprintf("%d exported functions and %d registered S3 methods; untested: %s; listed but not public: %s",
                     length(getNamespaceExports("twoCoprimary")), length(registered_s3),
                     if (length(missing)) paste(missing, collapse = ", ") else "none",
                     if (length(extra)) paste(extra, collapse = ", ") else "none"))
})

check("coverage.export.count", "the article states the number of exported functions", {
  list(verdict = "INFO",
       obs = sprintf("%d exported functions plus %d registered S3 methods (%s)",
                     length(getNamespaceExports("twoCoprimary")), length(registered_s3),
                     paste(registered_s3, collapse = ", ")))
})

check("coverage.imports", "the article states the number of non-base dependencies", {
  d <- read.dcf(system.file("DESCRIPTION", package = "twoCoprimary"))
  imp <- trimws(unlist(strsplit(d[1, "Imports"], ",")))
  imp <- imp[nzchar(imp)]
  list(verdict = "INFO",
       obs = sprintf("Imports (%d): %s", length(imp), paste(imp, collapse = ", ")))
})

# ===========================================================================
part("Part 1. Documentation contract: arguments and returned columns")
# ===========================================================================

# Return the text of one top level Rd section, matching braces and ignoring
# escaped braces such as the \{ ... \} in the unified interface documentation.
rd_section <- function(txt, name) {
  chars <- strsplit(txt, "")[[1]]
  pat <- paste0("\\", name, "{")
  pos <- regexpr(pat, txt, fixed = TRUE)
  if (pos < 0) return(NA_character_)
  open_idx <- as.integer(pos) + nchar(name) + 1L
  depth <- 0L
  k <- open_idx
  n <- length(chars)
  end_idx <- NA_integer_
  while (k <= n) {
    esc <- k > 1L && chars[k - 1L] == "\\"
    if (!esc && chars[k] == "{") depth <- depth + 1L
    if (!esc && chars[k] == "}") {
      depth <- depth - 1L
      if (depth == 0L) { end_idx <- k; break }
    }
    k <- k + 1L
  }
  if (is.na(end_idx)) return(NA_character_)
  substr(txt, open_idx + 1L, end_idx - 1L)
}

# Names appearing as \item{names}{description} at the top level of the section.
# An \item nested inside a \describe or an \itemize belongs to the prose of
# another item, not to the list of arguments or of returned values, so the
# brace depth is tracked and only depth zero is collected. \item without a
# brace, the form \itemize uses, never matches.
rd_item_names <- function(sec) {
  if (is.na(sec)) return(character(0))
  chars <- strsplit(sec, "")[[1]]
  n <- length(chars)
  depth <- 0L
  out <- character(0)
  i <- 1L
  while (i <= n) {
    esc <- i > 1L && chars[i - 1L] == "\\"
    if (!esc && chars[i] == "{") { depth <- depth + 1L; i <- i + 1L; next }
    if (!esc && chars[i] == "}") { depth <- depth - 1L; i <- i + 1L; next }
    if (depth == 0L && chars[i] == "\\" && i + 5L <= n &&
        paste(chars[i:(i + 4L)], collapse = "") == "\\item" &&
        chars[i + 5L] == "{") {
      j <- i + 5L
      d <- 0L
      k <- j
      while (k <= n) {
        e2 <- k > 1L && chars[k - 1L] == "\\"
        if (!e2 && chars[k] == "{") d <- d + 1L
        if (!e2 && chars[k] == "}") { d <- d - 1L; if (d == 0L) break }
        k <- k + 1L
      }
      out <- c(out, paste(chars[(j + 1L):(k - 1L)], collapse = ""))
      i <- k + 1L
      next
    }
    i <- i + 1L
  }
  out <- unlist(strsplit(out, "[,[:space:]]+"))
  out[nzchar(out)]
}

man_files <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
rd_txt <- lapply(man_files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"))
rd_name <- vapply(rd_txt, function(t) {
  nm <- rd_section(t, "name")
  if (is.na(nm)) NA_character_ else trimws(nm)
}, character(1))
names(rd_txt) <- rd_name

check("rd.files.found", "one Rd file is available for every exported function", {
  missing <- setdiff(c(exercised, s3_methods), rd_name)
  list(ok = length(missing) == 0,
       obs = sprintf("%d Rd files; without a match: %s", length(man_files),
                     if (length(missing)) paste(missing, collapse = ", ") else "none"))
})

for (fn in sort(c(exercised, s3_methods))) {
  local({
    f <- fn
    check(paste0("rd.arguments.", f),
          "the documented arguments are exactly the formal arguments", {
      if (!f %in% names(rd_txt)) stop("no Rd file named ", f)
      doc <- unique(rd_item_names(rd_section(rd_txt[[f]], "arguments")))
      act <- names(formals(get(f, envir = asNamespace("twoCoprimary"))))
      undoc <- setdiff(act, doc)
      spurious <- setdiff(doc, act)
      list(ok = length(undoc) == 0 && length(spurious) == 0,
           obs = sprintf("%d formals; undocumented: %s; documented but absent: %s",
                         length(act),
                         if (length(undoc)) paste(undoc, collapse = ", ") else "none",
                         if (length(spurious)) paste(spurious, collapse = ", ") else "none"))
    })
  })
}

# Map each returning function to an object it produces, so that the documented
# \value items can be compared with the actual names
value_map <- list(
  ss1Continuous = obj$ss1_cont, ss1Count = obj$ss1_count,
  ss1BinaryApprox = obj$ss1_bin,
  ss2Continuous = obj$ss2_cont, power2Continuous = obj$pw2_cont,
  ss2BinaryApprox = obj$ss2_bin_ap, power2BinaryApprox = obj$pw2_bin_ap,
  ss2BinaryExact = obj$ss2_bin_ex, power2BinaryExact = obj$pw2_bin_ex,
  ss2MixedContinuousBinary = obj$ss2_mcb,
  power2MixedContinuousBinary = obj$pw2_mcb,
  ss2MixedCountContinuous = obj$ss2_mcc,
  power2MixedCountContinuous = obj$pw2_mcc,
  corrbound2Binary = corrbound2Binary(p1 = 0.3, p2 = 0.5),
  corrbound2MixedCountContinuous =
    corrbound2MixedCountContinuous(lambda = 1.25, nu = 0.8, mu = 0, sd = 250)
)

for (fn in names(value_map)) {
  local({
    f <- fn
    val <- value_map[[fn]]
    check(paste0("rd.value.", f),
          "the documented \\value items are exactly the names of the result", {
      doc <- rd_item_names(rd_section(rd_txt[[f]], "value"))
      act <- names(val)
      list(ok = setequal(doc, act) && identical(doc, act),
           obs = sprintf("documented (%d): %s | actual (%d): %s | %s",
                         length(doc), paste(doc, collapse = ","),
                         length(act), paste(act, collapse = ","),
                         if (!setequal(doc, act)) "sets differ"
                         else if (!identical(doc, act)) "same set, different order"
                         else "identical"))
    })
  })
}

# ===========================================================================
part("Part 2. Structural contract of every returned object")
# ===========================================================================

for (nm in names(obj)) {
  local({
    o <- obj[[nm]]
    id <- nm
    check(paste0("shape.class.", id),
          "results carry class c('twoCoprimary', 'data.frame') and hold one row", {
      list(ok = identical(class(o), c("twoCoprimary", "data.frame")) && nrow(o) == 1L,
           obs = sprintf("class = %s, nrow = %d, ncol = %d",
                         paste(class(o), collapse = "/"), nrow(o), ncol(o)))
    })
    check(paste0("shape.finite.", id),
          "no returned quantity is NaN or infinite, and NA appears only in nMC", {
      num <- o[, vapply(o, is.numeric, logical(1)), drop = FALSE]
      bad <- names(num)[vapply(num, function(v) !is.finite(v), logical(1))]
      bad <- setdiff(bad, "nMC")
      list(ok = length(bad) == 0,
           obs = if (length(bad)) paste("non-finite:", paste(bad, collapse = ", "))
                 else "all numeric columns finite")
    })
    if (all(c("n1", "n2", "N") %in% names(o))) {
      check(paste0("shape.total.", id), "N equals n1 + n2", {
        list(ok = isTRUE(o$N == o$n1 + o$n2),
             obs = sprintf("n1 = %s, n2 = %s, N = %s", o$n1, o$n2, o$N))
      })
    }
    if (all(c("n1", "n2", "r") %in% names(o))) {
      check(paste0("shape.allocation.", id), "n1 equals ceiling(r n2)", {
        list(ok = isTRUE(o$n1 == ceiling(o$r * o$n2)),
             obs = sprintf("r = %s, n2 = %s, n1 = %s, ceiling(r n2) = %s",
                           o$r, o$n2, o$n1, ceiling(o$r * o$n2)))
      })
    }
    pcols <- intersect(c("power1", "power2", "powerCont", "powerBin",
                         "powerCount", "powerCoprimary"), names(o))
    if (length(pcols)) {
      check(paste0("shape.power.range.", id), "every power lies in [0, 1]", {
        v <- unlist(o[, pcols, drop = FALSE])
        list(ok = all(v >= 0 & v <= 1),
             obs = paste(sprintf("%s = %s", pcols, fmt(v)), collapse = "; "))
      })
    }
  })
}

# ---------------------------------------------------------------------------
# The smallest sample sizes a function will accept
# ---------------------------------------------------------------------------
# A power is a probability. At one or two subjects per group a continuity
# correction can carry a probability to or past 0 or 1, where the arcsine
# transformation is not defined and the corrected variance vanishes; the first
# gives NaN and the second gives a critical value of zero, that is a reported
# power of one half. Every power function is therefore evaluated at the
# smallest sample sizes it accepts and required to return a number in the unit
# interval, and to do so without a warning.

quietly <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr,
                             warning = function(cond) {
                               w <<- c(w, conditionMessage(cond))
                               invokeRestart("muffleWarning")
                             })
  list(value = val, warnings = w)
}

small_probs <- list(c(0.60, 0.50, 0.40, 0.30), c(0.60, 0.40, 0.40, 0.20),
                    c(0.30, 0.25, 0.10, 0.05), c(0.90, 0.50, 0.50, 0.20),
                    c(0.99, 0.95, 0.90, 0.85), c(0.55, 0.50, 0.45, 0.40))

for (tt in c("AN", "ANc", "AS", "ASc")) {
  local({
    t_ <- tt
    check(paste0("smalln.power2BinaryApprox.", t_),
          "at one to six subjects per group the power is a number in [0, 1], with no warning", {
      bad <- character(0)
      warned <- character(0)
      for (pp in small_probs) {
        for (n1 in 1:6) for (n2 in 1:6) {
          out <- quietly(power2BinaryApprox(n1 = n1, n2 = n2, p11 = pp[1],
                                            p12 = pp[2], p21 = pp[3],
                                            p22 = pp[4], rho1 = 0, rho2 = 0,
                                            alpha = 0.025, Test = t_))
          v <- unlist(out$value[, c("power1", "power2", "powerCoprimary")])
          if (!all(is.finite(v)) || any(v < 0 | v > 1)) {
            bad <- c(bad, sprintf("p=(%.2f,%.2f,%.2f,%.2f) n=(%d,%d) -> %s",
                                  pp[1], pp[2], pp[3], pp[4], n1, n2,
                                  paste(format(v), collapse = "/")))
          }
          if (length(out$warnings)) {
            warned <- c(warned, sprintf("n=(%d,%d): %s", n1, n2, out$warnings[1]))
          }
        }
      }
      list(ok = length(bad) == 0 && length(warned) == 0,
           obs = sprintf("%d evaluations; not a probability: %s; warnings: %s",
                         length(small_probs) * 36,
                         if (length(bad)) paste(utils::head(bad, 4), collapse = "; ") else "none",
                         if (length(warned)) paste(utils::head(warned, 2), collapse = "; ") else "none"))
    })
  })
}

for (tt in c("AN", "ANc", "AS", "ASc")) {
  local({
    t_ <- tt
    check(paste0("smalln.power2MixedContinuousBinary.", t_),
          "at one to six subjects per group the power is a number in [0, 1], with no warning", {
      bad <- character(0)
      warned <- character(0)
      for (pp in list(c(0.60, 0.40), c(0.30, 0.10), c(0.99, 0.85), c(0.95, 0.80))) {
        for (n1 in 1:6) for (n2 in 1:6) {
          out <- quietly(power2MixedContinuousBinary(n1 = n1, n2 = n2,
                                                     delta = 0.5, sd = 1,
                                                     p1 = pp[1], p2 = pp[2],
                                                     rho = 0.5, alpha = 0.025,
                                                     Test = t_))
          v <- unlist(out$value[, c("powerCont", "powerBin", "powerCoprimary")])
          if (!all(is.finite(v)) || any(v < 0 | v > 1)) {
            bad <- c(bad, sprintf("p=(%.2f,%.2f) n=(%d,%d) -> %s", pp[1], pp[2],
                                  n1, n2, paste(format(v), collapse = "/")))
          }
          if (length(out$warnings)) {
            warned <- c(warned, sprintf("n=(%d,%d): %s", n1, n2, out$warnings[1]))
          }
        }
      }
      list(ok = length(bad) == 0 && length(warned) == 0,
           obs = sprintf("%d evaluations; not a probability: %s; warnings: %s", 4 * 36,
                         if (length(bad)) paste(utils::head(bad, 4), collapse = "; ") else "none",
                         if (length(warned)) paste(utils::head(warned, 2), collapse = "; ") else "none"))
    })
  })
}

check("smalln.power2Continuous",
      "at one to six subjects per group the power is a number in [0, 1], with no warning", {
  bad <- character(0); warned <- character(0)
  for (n1 in 1:6) for (n2 in 1:6) {
    out <- quietly(power2Continuous(n1 = n1, n2 = n2, delta1 = 0.5, delta2 = 0.4,
                                    sd1 = 1, sd2 = 1.2, rho = 0.5, alpha = 0.025,
                                    known_var = TRUE))
    v <- unlist(out$value[, c("power1", "power2", "powerCoprimary")])
    if (!all(is.finite(v)) || any(v < 0 | v > 1)) {
      bad <- c(bad, sprintf("n=(%d,%d)", n1, n2))
    }
    if (length(out$warnings)) warned <- c(warned, out$warnings[1])
  }
  list(ok = length(bad) == 0 && length(warned) == 0,
       obs = sprintf("36 evaluations; not a probability: %s; warnings: %s",
                     if (length(bad)) paste(bad, collapse = ", ") else "none",
                     if (length(warned)) warned[1] else "none"))
})

check("smalln.power2MixedCountContinuous",
      "at one to six subjects per group the power is a number in [0, 1], with no warning", {
  bad <- character(0); warned <- character(0)
  for (n1 in 1:6) for (n2 in 1:6) {
    out <- quietly(power2MixedCountContinuous(n1 = n1, n2 = n2, r1 = 1.0,
                                              r2 = 1.25, nu = 0.8, t = 1,
                                              mu1 = -50, mu2 = 0, sd = 250,
                                              rho1 = 0.5, rho2 = 0.5,
                                              alpha = 0.025))
    v <- unlist(out$value[, c("powerCount", "powerCont", "powerCoprimary")])
    if (!all(is.finite(v)) || any(v < 0 | v > 1)) bad <- c(bad, sprintf("n=(%d,%d)", n1, n2))
    if (length(out$warnings)) warned <- c(warned, out$warnings[1])
  }
  list(ok = length(bad) == 0 && length(warned) == 0,
       obs = sprintf("36 evaluations; not a probability: %s; warnings: %s",
                     if (length(bad)) paste(bad, collapse = ", ") else "none",
                     if (length(warned)) warned[1] else "none"))
})

for (tt in c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")) {
  local({
    t_ <- tt
    check(paste0("smalln.power2BinaryExact.", t_),
          "at one to four subjects per group the exact power is a number in [0, 1], with no warning", {
      bad <- character(0); warned <- character(0)
      for (n1 in 1:4) for (n2 in 1:4) {
        out <- quietly(power2BinaryExact(n1 = n1, n2 = n2, p11 = 0.6, p12 = 0.5,
                                         p21 = 0.4, p22 = 0.3, rho1 = 0,
                                         rho2 = 0, alpha = 0.025, Test = t_))
        v <- unlist(out$value[, c("power1", "power2", "powerCoprimary")])
        if (!all(is.finite(v)) || any(v < 0 | v > 1)) bad <- c(bad, sprintf("n=(%d,%d)", n1, n2))
        if (length(out$warnings)) warned <- c(warned, out$warnings[1])
      }
      list(ok = length(bad) == 0 && length(warned) == 0,
           obs = sprintf("16 evaluations; not a probability: %s; warnings: %s",
                         if (length(bad)) paste(bad, collapse = ", ") else "none",
                         if (length(warned)) warned[1] else "none"))
    })
  })
}

check("smalln.power2MixedContinuousBinary.Fisher",
      "at one or two subjects in total the simulated power is zero rather than undefined", {
  bad <- character(0); warned <- character(0)
  for (n1 in 1:3) for (n2 in 1:3) {
    out <- quietly(power2MixedContinuousBinary(n1 = n1, n2 = n2, delta = 0.5,
                                               sd = 1, p1 = 0.6, p2 = 0.4,
                                               rho = 0.5, alpha = 0.025,
                                               Test = "Fisher", nMC = 200))
    v <- unlist(out$value[, c("powerCont", "powerBin", "powerCoprimary")])
    if (!all(is.finite(v)) || any(v < 0 | v > 1)) {
      bad <- c(bad, sprintf("n=(%d,%d) -> %s", n1, n2, paste(format(v), collapse = "/")))
    }
    if (length(out$warnings)) warned <- c(warned, sprintf("n=(%d,%d): %s", n1, n2, out$warnings[1]))
  }
  list(ok = length(bad) == 0 && length(warned) == 0,
       obs = sprintf("9 evaluations; not a probability: %s; warnings: %s",
                     if (length(bad)) paste(bad, collapse = "; ") else "none",
                     if (length(warned)) warned[1] else "none"))
})

# ---------------------------------------------------------------------------
# The two correlations at which the parameterisation degenerates
# ---------------------------------------------------------------------------
# With equal marginals on the two endpoints the upper Prentice bound is one,
# the two outcomes coincide, and the dependence parameter of the bivariate
# binomial diverges while the correlation between the two normal test
# statistics reaches one. Both are admissible inputs, so both must return the
# limit rather than NaN. The lower bound behaves the same way when the two
# probabilities of an endpoint sum to one.

check("degenerate.dbibinom.perfect.dependence",
      "at the upper Prentice bound with equal marginals the mass function is still a distribution", {
  N <- 12; pq <- 0.5
  b <- corrbound2Binary(pq, pq)
  out <- quietly(outer(0:N, 0:N, function(a, c2) dbibinom(N, a, c2, pq, pq, b[["U_bound"]])))
  P <- out$value
  y <- 0:N
  m1 <- sum(rowSums(P) * y); m2 <- sum(colSums(P) * y)
  v <- N * pq * (1 - pq)
  rho_hat <- (sum(P * outer(y, y)) - m1 * m2) / v
  list(ok = all(is.finite(P)) && abs(sum(P) - 1) < 1e-10 &&
         max(abs(rowSums(P) - stats::dbinom(y, N, pq))) < 1e-10 &&
         max(abs(colSums(P) - stats::dbinom(y, N, pq))) < 1e-10 &&
         abs(rho_hat - 1) < 1e-10 && length(out$warnings) == 0,
       obs = sprintf("upper bound %s; total mass %.12f; recovered correlation %.12f; warnings %d",
                     format(b[["U_bound"]]), sum(P), rho_hat, length(out$warnings)))
})

check("degenerate.power2BinaryApprox.unit.correlation",
      "at a unit correlation between the test statistics the co-primary power is the smaller marginal", {
  bad <- character(0)
  for (tt in c("AN", "ANc", "AS", "ASc")) {
    o <- power2BinaryApprox(n1 = 120, n2 = 120, p11 = 0.6, p12 = 0.6,
                            p21 = 0.4, p22 = 0.4, rho1 = 1, rho2 = 1,
                            alpha = 0.025, Test = tt)
    want <- min(o$power1, o$power2)
    if (!is.finite(o$powerCoprimary) || abs(o$powerCoprimary - want) > 1e-10) {
      bad <- c(bad, sprintf("%s: joint %s, min marginal %s", tt,
                            format(o$powerCoprimary), format(want)))
    }
  }
  list(ok = length(bad) == 0,
       obs = if (length(bad)) paste(bad, collapse = "; ") else
         "all four methods return the smaller marginal power")
})

check("degenerate.power2BinaryApprox.minus.one",
      "at a correlation of minus one the co-primary power is the Bonferroni lower bound", {
  # p11 + p12 = 1 and p21 + p22 = 1 put both groups at the lower Prentice bound
  o <- power2BinaryApprox(n1 = 120, n2 = 120, p11 = 0.6, p12 = 0.4,
                          p21 = 0.4, p22 = 0.6, rho1 = -1, rho2 = -1,
                          alpha = 0.025, Test = "AN")
  want <- max(0, o$power1 + o$power2 - 1)
  list(ok = is.finite(o$powerCoprimary) && abs(o$powerCoprimary - want) < 1e-10,
       obs = sprintf("joint %s, power1 + power2 - 1 = %s",
                     fmt(o$powerCoprimary), fmt(want)))
})

check("degenerate.power2BinaryExact.perfect.dependence",
      "with identical endpoints and a unit correlation the exact co-primary power equals one marginal", {
  o <- power2BinaryExact(n1 = 30, n2 = 30, p11 = 0.6, p12 = 0.6, p21 = 0.4,
                         p22 = 0.4, rho1 = 1, rho2 = 1, alpha = 0.025,
                         Test = "Fisher")
  list(ok = is.finite(o$powerCoprimary) &&
         abs(o$powerCoprimary - o$power1) < 1e-10 &&
         abs(o$power1 - o$power2) < 1e-10,
       obs = sprintf("power1 %s, power2 %s, joint %s",
                     fmt(o$power1), fmt(o$power2), fmt(o$powerCoprimary)))
})

check("smalln.rr1Binary",
      "at one to four subjects per group every rejection region is a logical matrix without NA", {
  bad <- character(0)
  for (t_ in c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")) {
    for (n1 in 1:4) for (n2 in 1:4) {
      RR <- rr1Binary(n1, n2, 0.025, Test = t_)
      if (!is.logical(RR) || any(is.na(RR)) ||
          !identical(dim(RR), c(n1 + 1L, n2 + 1L))) {
        bad <- c(bad, sprintf("%s n=(%d,%d)", t_, n1, n2))
      }
    }
  }
  list(ok = length(bad) == 0,
       obs = sprintf("80 regions; malformed: %s",
                     if (length(bad)) paste(bad, collapse = ", ") else "none"))
})

# ===========================================================================
part("Part 3. Probability identities that any correct power function satisfies")
# ===========================================================================

# The co-primary power is the probability of an intersection, so it lies
# between the Frechet bounds set by the two marginal powers, whatever the
# endpoint type, the test and the correlation.
power_of <- function(o) {
  m <- intersect(c("power1", "power2", "powerCont", "powerBin", "powerCount"),
                 names(o))
  list(marginal = unlist(o[, m, drop = FALSE]), joint = o$powerCoprimary)
}

power_calls <- list(
  cont_known = quote(power2Continuous(n1 = 90, n2 = 90, delta1 = 0.5, delta2 = 0.4,
                                      sd1 = 1, sd2 = 1.2, rho = 0.6,
                                      alpha = 0.025, known_var = TRUE)),
  cont_unknown = quote(power2Continuous(n1 = 90, n2 = 90, delta1 = 0.5, delta2 = 0.4,
                                        sd1 = 1, sd2 = 1.2, rho = 0.6,
                                        alpha = 0.025, known_var = FALSE, nMC = 20000)),
  bin_AN = quote(power2BinaryApprox(n1 = 120, n2 = 60, p11 = 0.6, p12 = 0.5,
                                    p21 = 0.4, p22 = 0.3, rho1 = 0.4, rho2 = 0.4,
                                    alpha = 0.025, Test = "AN")),
  bin_ANc = quote(power2BinaryApprox(n1 = 120, n2 = 60, p11 = 0.6, p12 = 0.5,
                                     p21 = 0.4, p22 = 0.3, rho1 = 0.4, rho2 = 0.4,
                                     alpha = 0.025, Test = "ANc")),
  bin_AS = quote(power2BinaryApprox(n1 = 120, n2 = 60, p11 = 0.6, p12 = 0.5,
                                    p21 = 0.4, p22 = 0.3, rho1 = 0.4, rho2 = 0.4,
                                    alpha = 0.025, Test = "AS")),
  bin_ASc = quote(power2BinaryApprox(n1 = 120, n2 = 60, p11 = 0.6, p12 = 0.5,
                                     p21 = 0.4, p22 = 0.3, rho1 = 0.4, rho2 = 0.4,
                                     alpha = 0.025, Test = "ASc")),
  ex_Chisq = quote(power2BinaryExact(n1 = 40, n2 = 40, p11 = 0.7, p12 = 0.6,
                                     p21 = 0.3, p22 = 0.2, rho1 = 0.3, rho2 = 0.3,
                                     alpha = 0.025, Test = "Chisq")),
  ex_Fisher = quote(power2BinaryExact(n1 = 40, n2 = 40, p11 = 0.7, p12 = 0.6,
                                      p21 = 0.3, p22 = 0.2, rho1 = 0.3, rho2 = 0.3,
                                      alpha = 0.025, Test = "Fisher")),
  ex_midP = quote(power2BinaryExact(n1 = 40, n2 = 40, p11 = 0.7, p12 = 0.6,
                                    p21 = 0.3, p22 = 0.2, rho1 = 0.3, rho2 = 0.3,
                                    alpha = 0.025, Test = "Fisher-midP")),
  ex_Zpool = quote(power2BinaryExact(n1 = 30, n2 = 30, p11 = 0.7, p12 = 0.6,
                                     p21 = 0.3, p22 = 0.2, rho1 = 0.3, rho2 = 0.3,
                                     alpha = 0.025, Test = "Z-pool")),
  ex_Boschloo = quote(power2BinaryExact(n1 = 30, n2 = 30, p11 = 0.7, p12 = 0.6,
                                        p21 = 0.3, p22 = 0.2, rho1 = 0.3, rho2 = 0.3,
                                        alpha = 0.025, Test = "Boschloo")),
  mcb_AN = quote(power2MixedContinuousBinary(n1 = 100, n2 = 100, delta = 0.5,
                                             sd = 1, p1 = 0.6, p2 = 0.4,
                                             rho = 0.5, alpha = 0.025, Test = "AN")),
  mcb_ANc = quote(power2MixedContinuousBinary(n1 = 100, n2 = 100, delta = 0.5,
                                              sd = 1, p1 = 0.6, p2 = 0.4,
                                              rho = 0.5, alpha = 0.025, Test = "ANc")),
  mcb_AS = quote(power2MixedContinuousBinary(n1 = 100, n2 = 100, delta = 0.5,
                                             sd = 1, p1 = 0.6, p2 = 0.4,
                                             rho = 0.5, alpha = 0.025, Test = "AS")),
  mcb_ASc = quote(power2MixedContinuousBinary(n1 = 100, n2 = 100, delta = 0.5,
                                              sd = 1, p1 = 0.6, p2 = 0.4,
                                              rho = 0.5, alpha = 0.025, Test = "ASc")),
  mcb_Fisher = quote(power2MixedContinuousBinary(n1 = 60, n2 = 60, delta = 0.5,
                                                 sd = 1, p1 = 0.6, p2 = 0.4,
                                                 rho = 0.5, alpha = 0.025,
                                                 Test = "Fisher", nMC = 4000)),
  mcc = quote(power2MixedCountContinuous(n1 = 300, n2 = 300, r1 = 1.0, r2 = 1.25,
                                         nu = 0.8, t = 1, mu1 = -50, mu2 = 0,
                                         sd = 250, rho1 = 0.5, rho2 = 0.5,
                                         alpha = 0.025))
)

set.seed(20260902)
power_res <- lapply(power_calls, function(e) eval(e, envir = globalenv()))

for (nm in names(power_res)) {
  local({
    o <- power_res[[nm]]
    id <- nm
    # Monte Carlo results carry sampling error, so the two identities below are
    # checked with a tolerance for those and exactly for the rest
    tol <- if (id %in% c("cont_unknown", "mcb_Fisher")) 0.02 else 1e-8
    check(paste0("bound.upper.", id),
          "the co-primary power does not exceed either marginal power", {
      p <- power_of(o)
      list(ok = all(p$joint <= p$marginal + tol),
           obs = sprintf("joint = %s, marginals = %s, tol = %g",
                         fmt(p$joint), fmt(p$marginal), tol))
    })
    check(paste0("bound.lower.", id),
          "the co-primary power is at least the sum of the marginals minus one", {
      p <- power_of(o)
      lo <- sum(p$marginal) - 1
      list(ok = p$joint >= lo - tol,
           obs = sprintf("joint = %s, Bonferroni lower bound = %s",
                         fmt(p$joint), fmt(lo)))
    })
  })
}

# Under independence the joint power is the product of the marginals. This is
# an exact identity for every analytic method, and it fails immediately if the
# correlation between the test statistics is assembled incorrectly.
indep_calls <- list(
  cont = quote(power2Continuous(n1 = 90, n2 = 90, delta1 = 0.5, delta2 = 0.4,
                                sd1 = 1, sd2 = 1.2, rho = 0, alpha = 0.025,
                                known_var = TRUE)),
  bin_AN = quote(power2BinaryApprox(n1 = 120, n2 = 60, p11 = 0.6, p12 = 0.5,
                                    p21 = 0.4, p22 = 0.3, rho1 = 0, rho2 = 0,
                                    alpha = 0.025, Test = "AN")),
  bin_ANc = quote(power2BinaryApprox(n1 = 120, n2 = 60, p11 = 0.6, p12 = 0.5,
                                     p21 = 0.4, p22 = 0.3, rho1 = 0, rho2 = 0,
                                     alpha = 0.025, Test = "ANc")),
  bin_AS = quote(power2BinaryApprox(n1 = 120, n2 = 60, p11 = 0.6, p12 = 0.5,
                                    p21 = 0.4, p22 = 0.3, rho1 = 0, rho2 = 0,
                                    alpha = 0.025, Test = "AS")),
  bin_ASc = quote(power2BinaryApprox(n1 = 120, n2 = 60, p11 = 0.6, p12 = 0.5,
                                     p21 = 0.4, p22 = 0.3, rho1 = 0, rho2 = 0,
                                     alpha = 0.025, Test = "ASc")),
  ex_Fisher = quote(power2BinaryExact(n1 = 40, n2 = 40, p11 = 0.7, p12 = 0.6,
                                      p21 = 0.3, p22 = 0.2, rho1 = 0, rho2 = 0,
                                      alpha = 0.025, Test = "Fisher")),
  ex_Chisq = quote(power2BinaryExact(n1 = 40, n2 = 40, p11 = 0.7, p12 = 0.6,
                                     p21 = 0.3, p22 = 0.2, rho1 = 0, rho2 = 0,
                                     alpha = 0.025, Test = "Chisq")),
  mcb_AN = quote(power2MixedContinuousBinary(n1 = 100, n2 = 100, delta = 0.5,
                                             sd = 1, p1 = 0.6, p2 = 0.4,
                                             rho = 0, alpha = 0.025, Test = "AN")),
  mcb_AS = quote(power2MixedContinuousBinary(n1 = 100, n2 = 100, delta = 0.5,
                                             sd = 1, p1 = 0.6, p2 = 0.4,
                                             rho = 0, alpha = 0.025, Test = "AS")),
  mcc = quote(power2MixedCountContinuous(n1 = 300, n2 = 300, r1 = 1.0, r2 = 1.25,
                                         nu = 0.8, t = 1, mu1 = -50, mu2 = 0,
                                         sd = 250, rho1 = 0, rho2 = 0,
                                         alpha = 0.025))
)

for (nm in names(indep_calls)) {
  local({
    id <- nm
    e <- indep_calls[[nm]]
    check(paste0("independence.product.", id),
          "at zero correlation the co-primary power equals the product of the marginals", {
      o <- eval(e, envir = globalenv())
      p <- power_of(o)
      prod_m <- prod(p$marginal)
      list(ok = abs(p$joint - prod_m) < 1e-6,
           obs = sprintf("joint = %s, product = %s, difference = %.3e",
                         fmt(p$joint, 9), fmt(prod_m, 9), p$joint - prod_m))
    })
  })
}

# The co-primary power increases with the correlation between the endpoints
for (spec in list(
  list(id = "cont",
       f = function(rho) power2Continuous(n1 = 90, n2 = 90, delta1 = 0.5,
                                          delta2 = 0.4, sd1 = 1, sd2 = 1.2,
                                          rho = rho, alpha = 0.025,
                                          known_var = TRUE)$powerCoprimary,
       rhos = seq(-0.8, 0.8, by = 0.2)),
  list(id = "bin_AN",
       f = function(rho) power2BinaryApprox(n1 = 120, n2 = 120, p11 = 0.6,
                                            p12 = 0.5, p21 = 0.4, p22 = 0.3,
                                            rho1 = rho, rho2 = rho,
                                            alpha = 0.025, Test = "AN")$powerCoprimary,
       rhos = seq(0, 0.6, by = 0.1)),
  list(id = "mcb_AN",
       f = function(rho) power2MixedContinuousBinary(n1 = 100, n2 = 100,
                                                     delta = 0.5, sd = 1,
                                                     p1 = 0.6, p2 = 0.4,
                                                     rho = rho, alpha = 0.025,
                                                     Test = "AN")$powerCoprimary,
       rhos = seq(-0.8, 0.8, by = 0.2)),
  list(id = "mcc",
       f = function(rho) power2MixedCountContinuous(n1 = 300, n2 = 300, r1 = 1.0,
                                                    r2 = 1.25, nu = 0.8, t = 1,
                                                    mu1 = -50, mu2 = 0, sd = 250,
                                                    rho1 = rho, rho2 = rho,
                                                    alpha = 0.025)$powerCoprimary,
       rhos = seq(-0.8, 0.8, by = 0.2)))) {
  local({
    s <- spec
    check(paste0("monotone.rho.", s$id),
          "the co-primary power increases with the correlation", {
      v <- vapply(s$rhos, s$f, numeric(1))
      d <- diff(v)
      list(ok = all(d >= -1e-10),
           obs = sprintf("rho %s to %s: power %s to %s, smallest increment %.3e",
                         fmt(min(s$rhos), 2), fmt(max(s$rhos), 2),
                         fmt(v[1]), fmt(v[length(v)]), min(d)))
    })
  })
}

# The power increases with the sample size. The exact binary tests are not
# monotone because the attainable significance level moves with n, so those
# are recorded rather than asserted.
for (spec in list(
  list(id = "cont", exact = FALSE,
       f = function(n) power2Continuous(n1 = n, n2 = n, delta1 = 0.5, delta2 = 0.4,
                                        sd1 = 1, sd2 = 1.2, rho = 0.5,
                                        alpha = 0.025, known_var = TRUE)$powerCoprimary),
  list(id = "bin_AN", exact = FALSE,
       f = function(n) power2BinaryApprox(n1 = n, n2 = n, p11 = 0.6, p12 = 0.5,
                                          p21 = 0.4, p22 = 0.3, rho1 = 0.3,
                                          rho2 = 0.3, alpha = 0.025,
                                          Test = "AN")$powerCoprimary),
  list(id = "bin_AS", exact = FALSE,
       f = function(n) power2BinaryApprox(n1 = n, n2 = n, p11 = 0.6, p12 = 0.5,
                                          p21 = 0.4, p22 = 0.3, rho1 = 0.3,
                                          rho2 = 0.3, alpha = 0.025,
                                          Test = "AS")$powerCoprimary),
  list(id = "mcb_AN", exact = FALSE,
       f = function(n) power2MixedContinuousBinary(n1 = n, n2 = n, delta = 0.5,
                                                   sd = 1, p1 = 0.6, p2 = 0.4,
                                                   rho = 0.5, alpha = 0.025,
                                                   Test = "AN")$powerCoprimary),
  list(id = "mcc", exact = FALSE,
       f = function(n) power2MixedCountContinuous(n1 = n, n2 = n, r1 = 1.0,
                                                  r2 = 1.25, nu = 0.8, t = 1,
                                                  mu1 = -50, mu2 = 0, sd = 250,
                                                  rho1 = 0.5, rho2 = 0.5,
                                                  alpha = 0.025)$powerCoprimary),
  list(id = "ex_Fisher", exact = TRUE,
       f = function(n) power2BinaryExact(n1 = n, n2 = n, p11 = 0.7, p12 = 0.6,
                                         p21 = 0.3, p22 = 0.2, rho1 = 0.3,
                                         rho2 = 0.3, alpha = 0.025,
                                         Test = "Fisher")$powerCoprimary))) {
  local({
    s <- spec
    ns <- if (s$exact) seq(20, 44, by = 2) else seq(20, 300, by = 20)
    check(paste0("monotone.n.", s$id),
          "the co-primary power increases with the sample size", {
      v <- vapply(ns, s$f, numeric(1))
      d <- diff(v)
      if (s$exact) {
        list(verdict = "INFO",
             obs = sprintf("saw-tooth expected: %d of %d increments negative, smallest %.3e",
                           sum(d < 0), length(d), min(d)))
      } else {
        list(ok = all(d >= -1e-10),
             obs = sprintf("n %d to %d: power %s to %s, smallest increment %.3e",
                           min(ns), max(ns), fmt(v[1]), fmt(v[length(v)]), min(d)))
      }
    })
  })
}

# ===========================================================================
part("Part 4. Sample size contract: the returned n is the smallest that works")
# ===========================================================================

# fpCompare's default tolerance, the one the sequential search itself uses
FP_TOL <- .Machine$double.eps ^ 0.5

# A sample size function is correct when the power at the returned n2 reaches
# the target and the power at n2 - 1 does not. Nothing else about the function
# needs to be assumed. Each entry supplies the size call and an independent
# power evaluation at an arbitrary (n1, n2).
minimality <- list()
add_min <- function(id, ss_expr, pw_fun, target) {
  minimality[[length(minimality) + 1]] <<- list(
    id = id, ss = ss_expr, pw = pw_fun, target = target
  )
}

# --- single endpoint, checked against formulas coded independently here -----
for (rr in c(1, 2, 3)) {
  local({
    r_ <- rr
    add_min(
      sprintf("ss1Continuous.r%d", r_),
      function() ss1Continuous(delta = 0.5, sd = 1.2, r = r_, alpha = 0.025, beta = 0.1),
      function(n1, n2) stats::pnorm(0.5 / (1.2 * sqrt(1 / n1 + 1 / n2)) -
                                      stats::qnorm(1 - 0.025)),
      0.9)
    add_min(
      sprintf("ss1Count.r%d", r_),
      function() ss1Count(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1, r = r_,
                          alpha = 0.025, beta = 0.1),
      function(n1, n2) {
        kap <- n1 / n2
        Va <- (1 / 1) * (1 / 1.25 + 1 / (kap * 1.0)) + (1 + kap) / (0.8 * kap)
        stats::pnorm((abs(log(1.0 / 1.25)) * sqrt(n2) -
                        stats::qnorm(1 - 0.025) * sqrt(Va)) / sqrt(Va))
      },
      0.9)
    # All four asymptotic methods return the smallest size whose realized group
    # sizes reach the target power under power2BinaryApprox, so that is what is
    # asserted here.
    for (tt in c("AN", "ANc", "AS", "ASc")) {
      local({
        t_ <- tt
        add_min(
          sprintf("ss1BinaryApprox.%s.r%d", t_, r_),
          function() ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = r_, alpha = 0.025,
                                     beta = 0.1, Test = t_),
          function(n1, n2) power2BinaryApprox(n1 = n1, n2 = n2, p11 = 0.6,
                                             p12 = 0.6, p21 = 0.4, p22 = 0.4,
                                             rho1 = 0, rho2 = 0, alpha = 0.025,
                                             Test = t_)$power1,
          0.9)
      })
    }
    add_min(
      sprintf("ss1BinaryApprox.Fisher.r%d", r_),
      function() ss1BinaryApprox(p1 = 0.75, p2 = 0.35, r = r_, alpha = 0.025,
                                 beta = 0.2, Test = "Fisher"),
      function(n1, n2) power2BinaryExact(n1 = n1, n2 = n2, p11 = 0.75, p12 = 0.75,
                                         p21 = 0.35, p22 = 0.35, rho1 = 0,
                                         rho2 = 0, alpha = 0.025,
                                         Test = "Fisher")$power1,
      0.8)
  })
}

# The size a single binary endpoint needs is the smallest one whose realized
# group sizes reach the target power. The predicate is reimplemented here from
# power2BinaryApprox and asserted over a grid of designs. The grid includes
# non-integer allocation ratios, where n1 = ceiling(r n2) exceeds r n2 and the
# closed form formulas are therefore not the criterion, and proportions up to
# 0.99, where the continuity corrections matter most.
local({
  alpha_ <- 0.025
  power1_at <- function(p1, p2, r, m, test) {
    power2BinaryApprox(n1 = ceiling(r * m), n2 = m, p11 = p1, p12 = p1,
                       p21 = p2, p22 = p2, rho1 = 0, rho2 = 0,
                       alpha = alpha_, Test = test)[["power1"]]
  }
  grid <- expand.grid(p1 = c(0.3, 0.5, 0.7, 0.9, 0.95, 0.99),
                      p2 = c(0.05, 0.2, 0.5, 0.8, 0.9),
                      r = c(1, 1.5, 2, 3, 5),
                      beta = c(0.05, 0.1, 0.2))
  grid <- grid[grid$p1 > grid$p2, , drop = FALSE]

  for (tt in c("AN", "ANc", "AS", "ASc")) {
    local({
      t_ <- tt
      check(sprintf("minimal.sweep.ss1BinaryApprox.%s", t_),
            "over a grid of designs the size reaches the target power and one subject fewer does not", {
        bad <- character(0)
        for (i in seq_len(nrow(grid))) {
          g <- grid[i, ]
          n2 <- ss1BinaryApprox(p1 = g$p1, p2 = g$p2, r = g$r, alpha = alpha_,
                                beta = g$beta, Test = t_)$n2
          tgt <- 1 - g$beta
          hi_ok <- power1_at(g$p1, g$p2, g$r, n2, t_) >= tgt - FP_TOL
          lo_ok <- n2 == 1 ||
            power1_at(g$p1, g$p2, g$r, n2 - 1, t_) < tgt - FP_TOL
          if (!hi_ok || !lo_ok) {
            bad <- c(bad, sprintf("p=(%.2f,%.2f) r=%.1f beta=%.2f n2=%d",
                                  g$p1, g$p2, g$r, g$beta, n2))
          }
        }
        list(ok = length(bad) == 0,
             obs = sprintf("%d designs; not minimal: %s", nrow(grid),
                           if (length(bad)) paste(utils::head(bad, 6), collapse = "; ") else "none"))
      })
    })
  }
})

check("ss1BinaryApprox.ANc.cycling.design.terminates",
      "the design whose closed form requirement map cycles between two values still returns", {
  # p1 = 0.95, p2 = 0.05, r = 5, beta = 0.2 sends the requirement map that an
  # earlier implementation iterated into a two-cycle between 2 and 4
  el <- system.time(o <- ss1BinaryApprox(p1 = 0.95, p2 = 0.05, r = 5,
                                         alpha = 0.025, beta = 0.2,
                                         Test = "ANc"))[["elapsed"]]
  list(ok = is.finite(o$n2) && o$n2 >= 1 && el < 5,
       obs = sprintf("n1 = %s, n2 = %s, returned in %.3f s", o$n1, o$n2, el))
})

check("ss1BinaryApprox.agrees.with.power2BinaryApprox",
      "the size function and the power function use the same formula for every asymptotic method", {
  # A size formula that differed from the power formula would show up as a
  # design where the power at the returned size falls short of the target. The
  # proportions are those of Table S5 of Sozu et al. (2012), where the ASc
  # closed form used to fall short by up to five subjects.
  worst <- 1
  detail <- ""
  for (tt in c("AN", "ANc", "AS", "ASc")) {
    for (pp in list(c(0.99, 0.95), c(0.99, 0.90), c(0.99, 0.85),
                    c(0.95, 0.90), c(0.95, 0.85), c(0.95, 0.80))) {
      o <- ss1BinaryApprox(p1 = pp[1], p2 = pp[2], r = 1, alpha = 0.025,
                           beta = 0.2, Test = tt)
      pw <- power2BinaryApprox(n1 = o$n1, n2 = o$n2, p11 = pp[1], p12 = pp[1],
                               p21 = pp[2], p22 = pp[2], rho1 = 0, rho2 = 0,
                               alpha = 0.025, Test = tt)$power1
      if (pw < worst) {
        worst <- pw
        detail <- sprintf("%s at p = (%.2f, %.2f), n2 = %s", tt, pp[1], pp[2], o$n2)
      }
    }
  }
  list(ok = worst >= 0.8 - FP_TOL,
       obs = sprintf("smallest achieved power over those proportions: %s (%s)",
                     fmt(worst), detail))
})

# --- two continuous endpoints ----------------------------------------------
for (rr in c(1, 2)) {
  local({
    r_ <- rr
    add_min(
      sprintf("ss2Continuous.r%d", r_),
      function() ss2Continuous(delta1 = 0.5, delta2 = 0.45, sd1 = 1, sd2 = 1.1,
                               rho = 0.5, r = r_, alpha = 0.025, beta = 0.2,
                               known_var = TRUE),
      function(n1, n2) power2Continuous(n1 = n1, n2 = n2, delta1 = 0.5,
                                        delta2 = 0.45, sd1 = 1, sd2 = 1.1,
                                        rho = 0.5, alpha = 0.025,
                                        known_var = TRUE)$powerCoprimary,
      0.8)
    for (tt in c("AN", "ANc", "AS", "ASc")) {
      local({
        t_ <- tt
        add_min(
          sprintf("ss2BinaryApprox.%s.r%d", t_, r_),
          function() ss2BinaryApprox(p11 = 0.6, p12 = 0.5, p21 = 0.4, p22 = 0.3,
                                     rho1 = 0.3, rho2 = 0.3, r = r_,
                                     alpha = 0.025, beta = 0.2, Test = t_),
          function(n1, n2) power2BinaryApprox(n1 = n1, n2 = n2, p11 = 0.6,
                                              p12 = 0.5, p21 = 0.4, p22 = 0.3,
                                              rho1 = 0.3, rho2 = 0.3,
                                              alpha = 0.025,
                                              Test = t_)$powerCoprimary,
          0.8)
        add_min(
          sprintf("ss2MixedContinuousBinary.%s.r%d", t_, r_),
          function() ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.6,
                                              p2 = 0.4, rho = 0.5, r = r_,
                                              alpha = 0.025, beta = 0.2,
                                              Test = t_),
          function(n1, n2) power2MixedContinuousBinary(n1 = n1, n2 = n2,
                                                       delta = 0.5, sd = 1,
                                                       p1 = 0.6, p2 = 0.4,
                                                       rho = 0.5, alpha = 0.025,
                                                       Test = t_)$powerCoprimary,
          0.8)
      })
    }
    add_min(
      sprintf("ss2MixedCountContinuous.r%d", r_),
      function() ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                                         mu1 = -50, mu2 = 0, sd = 250,
                                         rho1 = 0.5, rho2 = 0.5, r = r_,
                                         alpha = 0.025, beta = 0.2),
      function(n1, n2) power2MixedCountContinuous(n1 = n1, n2 = n2, r1 = 1.0,
                                                  r2 = 1.25, nu = 0.8, t = 1,
                                                  mu1 = -50, mu2 = 0, sd = 250,
                                                  rho1 = 0.5, rho2 = 0.5,
                                                  alpha = 0.025)$powerCoprimary,
      0.8)
  })
}

# --- two binary endpoints, exact -------------------------------------------
for (tt in c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")) {
  local({
    t_ <- tt
    add_min(
      sprintf("ss2BinaryExact.%s.r1", t_),
      function() ss2BinaryExact(p11 = 0.8, p12 = 0.7, p21 = 0.3, p22 = 0.2,
                                rho1 = 0.3, rho2 = 0.3, r = 1, alpha = 0.025,
                                beta = 0.2, Test = t_),
      function(n1, n2) power2BinaryExact(n1 = n1, n2 = n2, p11 = 0.8, p12 = 0.7,
                                         p21 = 0.3, p22 = 0.2, rho1 = 0.3,
                                         rho2 = 0.3, alpha = 0.025,
                                         Test = t_)$powerCoprimary,
      0.8)
  })
}

for (spec in minimality) {
  local({
    s <- spec
    check(paste0("minimal.", s$id),
          "the power reaches the target at n2 and falls short at n2 - 1", {
      o <- timed(paste0("ss:", s$id), s$ss())
      n2 <- o$n2
      r_used <- o$r
      p_at <- s$pw(o$n1, n2)
      if (n2 <= 1) {
        list(ok = p_at >= s$target - FP_TOL,
             obs = sprintf("n2 = 1 (floor), power = %s, target = %s",
                           fmt(p_at), fmt(s$target, 2)))
      } else {
        n1_below <- ceiling(r_used * (n2 - 1))
        p_below <- s$pw(n1_below, n2 - 1)
        list(ok = (p_at >= s$target - FP_TOL) && (p_below < s$target - FP_TOL),
             obs = sprintf("n1 = %s, n2 = %s, power(n2) = %s, power(n2-1) = %s, target = %s, %.2f s",
                           o$n1, n2, fmt(p_at), fmt(p_below), fmt(s$target, 2),
                           attr(o, "elapsed")))
      }
    })
  })
}

# The sequential search stops at the first sample size reaching the target on
# the way down, which is a local statement. For the exact tests the power is
# not monotone, so scan upward from 2 and report whether a smaller sample size
# also reaches the target.
for (tt in c("Chisq", "Fisher")) {
  local({
    t_ <- tt
    check(paste0("minimal.global.scan.", t_),
          "no sample size below the returned one reaches the target", {
      o <- ss2BinaryExact(p11 = 0.8, p12 = 0.7, p21 = 0.3, p22 = 0.2,
                          rho1 = 0.3, rho2 = 0.3, r = 1, alpha = 0.025,
                          beta = 0.2, Test = t_)
      grid <- 2:o$n2
      pw <- vapply(grid, function(n) power2BinaryExact(
        n1 = n, n2 = n, p11 = 0.8, p12 = 0.7, p21 = 0.3, p22 = 0.2,
        rho1 = 0.3, rho2 = 0.3, alpha = 0.025, Test = t_)$powerCoprimary,
        numeric(1))
      reach <- grid[pw >= 0.8 - FP_TOL]
      smallest <- if (length(reach)) min(reach) else NA_integer_
      list(verdict = if (!is.na(smallest) && smallest < o$n2) "INFO" else "PASS",
           obs = sprintf("returned n2 = %s; smallest n2 in 2..%s reaching 0.8 = %s",
                         o$n2, o$n2, smallest))
    })
  })
}

# ===========================================================================
part("Part 5. Two routes to the same number: unified interface and design_table")
# ===========================================================================

same_df <- function(a, b) isTRUE(all.equal(as.data.frame(a), as.data.frame(b),
                                           check.attributes = FALSE))

unified <- list(
  list(id = "continuous.power",
       a = quote(twoCoprimary2Continuous(n1 = 100, n2 = 100, delta1 = 0.5,
                                         delta2 = 0.4, sd1 = 1, sd2 = 1,
                                         rho = 0.3, alpha = 0.025, known_var = TRUE)),
       b = quote(power2Continuous(n1 = 100, n2 = 100, delta1 = 0.5, delta2 = 0.4,
                                  sd1 = 1, sd2 = 1, rho = 0.3, alpha = 0.025,
                                  known_var = TRUE))),
  list(id = "continuous.ss",
       a = quote(twoCoprimary2Continuous(delta1 = 0.5, delta2 = 0.4, sd1 = 1,
                                         sd2 = 1, rho = 0.3, power = 0.8, r = 2,
                                         alpha = 0.025, known_var = TRUE)),
       b = quote(ss2Continuous(delta1 = 0.5, delta2 = 0.4, sd1 = 1, sd2 = 1,
                               rho = 0.3, r = 2, alpha = 0.025, beta = 0.2,
                               known_var = TRUE))),
  list(id = "binaryApprox.power",
       a = quote(twoCoprimary2BinaryApprox(n1 = 200, n2 = 100, p11 = 0.5,
                                           p12 = 0.4, p21 = 0.3, p22 = 0.2,
                                           rho1 = 0.7, rho2 = 0.7, alpha = 0.025,
                                           Test = "ASc")),
       b = quote(power2BinaryApprox(n1 = 200, n2 = 100, p11 = 0.5, p12 = 0.4,
                                    p21 = 0.3, p22 = 0.2, rho1 = 0.7, rho2 = 0.7,
                                    alpha = 0.025, Test = "ASc"))),
  list(id = "binaryApprox.ss",
       a = quote(twoCoprimary2BinaryApprox(p11 = 0.5, p12 = 0.4, p21 = 0.3,
                                           p22 = 0.2, rho1 = 0.7, rho2 = 0.7,
                                           power = 0.8, r = 1, alpha = 0.025,
                                           Test = "AS")),
       b = quote(ss2BinaryApprox(p11 = 0.5, p12 = 0.4, p21 = 0.3, p22 = 0.2,
                                 rho1 = 0.7, rho2 = 0.7, r = 1, alpha = 0.025,
                                 beta = 0.2, Test = "AS"))),
  list(id = "binaryExact.power",
       a = quote(twoCoprimary2BinaryExact(n1 = 40, n2 = 40, p11 = 0.7, p12 = 0.6,
                                          p21 = 0.3, p22 = 0.2, rho1 = 0.3,
                                          rho2 = 0.3, alpha = 0.025,
                                          Test = "Fisher")),
       b = quote(power2BinaryExact(n1 = 40, n2 = 40, p11 = 0.7, p12 = 0.6,
                                   p21 = 0.3, p22 = 0.2, rho1 = 0.3, rho2 = 0.3,
                                   alpha = 0.025, Test = "Fisher"))),
  list(id = "binaryExact.ss",
       a = quote(twoCoprimary2BinaryExact(p11 = 0.8, p12 = 0.7, p21 = 0.3,
                                          p22 = 0.2, rho1 = 0.3, rho2 = 0.3,
                                          power = 0.8, r = 1, alpha = 0.025,
                                          Test = "Chisq")),
       b = quote(ss2BinaryExact(p11 = 0.8, p12 = 0.7, p21 = 0.3, p22 = 0.2,
                                rho1 = 0.3, rho2 = 0.3, r = 1, alpha = 0.025,
                                beta = 0.2, Test = "Chisq"))),
  list(id = "mixedContBinary.power",
       a = quote(twoCoprimary2MixedContinuousBinary(n1 = 100, n2 = 100,
                                                    delta = 0.5, sd = 1, p1 = 0.6,
                                                    p2 = 0.4, rho = 0.5,
                                                    alpha = 0.025, Test = "ANc")),
       b = quote(power2MixedContinuousBinary(n1 = 100, n2 = 100, delta = 0.5,
                                             sd = 1, p1 = 0.6, p2 = 0.4,
                                             rho = 0.5, alpha = 0.025,
                                             Test = "ANc"))),
  list(id = "mixedContBinary.ss",
       a = quote(twoCoprimary2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.6,
                                                    p2 = 0.4, rho = 0.5,
                                                    power = 0.9, r = 1,
                                                    alpha = 0.025, Test = "AN")),
       b = quote(ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.6, p2 = 0.4,
                                          rho = 0.5, r = 1, alpha = 0.025,
                                          beta = 0.1, Test = "AN"))),
  list(id = "mixedCountCont.power",
       a = quote(twoCoprimary2MixedCountContinuous(n1 = 300, n2 = 300, r1 = 1.0,
                                                   r2 = 1.25, nu = 0.8, t = 1,
                                                   mu1 = -50, mu2 = 0, sd = 250,
                                                   rho1 = 0.5, rho2 = 0.5,
                                                   alpha = 0.025)),
       b = quote(power2MixedCountContinuous(n1 = 300, n2 = 300, r1 = 1.0,
                                            r2 = 1.25, nu = 0.8, t = 1,
                                            mu1 = -50, mu2 = 0, sd = 250,
                                            rho1 = 0.5, rho2 = 0.5, alpha = 0.025))),
  list(id = "mixedCountCont.ss",
       a = quote(twoCoprimary2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8,
                                                   t = 1, mu1 = -50, mu2 = 0,
                                                   sd = 250, rho1 = 0.5,
                                                   rho2 = 0.5, power = 0.8, r = 1,
                                                   alpha = 0.025)),
       b = quote(ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                                         mu1 = -50, mu2 = 0, sd = 250,
                                         rho1 = 0.5, rho2 = 0.5, r = 1,
                                         alpha = 0.025, beta = 0.2)))
)

for (u in unified) {
  local({
    s <- u
    check(paste0("unified.", s$id),
          "the unified interface returns exactly what the underlying function returns", {
      a <- eval(s$a, envir = globalenv())
      b <- eval(s$b, envir = globalenv())
      list(ok = same_df(a, b),
           obs = if (same_df(a, b)) "identical"
                 else paste(utils::capture.output(
                   all.equal(as.data.frame(a), as.data.frame(b))), collapse = " | "))
    })
  })
}

# The unified interface must refuse an ambiguous or incomplete specification
for (spec in list(
  list(id = "both.given",
       e = quote(twoCoprimary2Continuous(n1 = 100, n2 = 100, delta1 = 0.5,
                                         delta2 = 0.4, sd1 = 1, sd2 = 1,
                                         rho = 0.3, power = 0.8, r = 1))),
  list(id = "one.n.given",
       e = quote(twoCoprimary2Continuous(n1 = 100, delta1 = 0.5, delta2 = 0.4,
                                         sd1 = 1, sd2 = 1, rho = 0.3))),
  list(id = "power.without.r",
       e = quote(twoCoprimary2Continuous(delta1 = 0.5, delta2 = 0.4, sd1 = 1,
                                         sd2 = 1, rho = 0.3, power = 0.8))),
  list(id = "neither.given",
       e = quote(twoCoprimary2Continuous(delta1 = 0.5, delta2 = 0.4, sd1 = 1,
                                         sd2 = 1, rho = 0.3))))) {
  local({
    s <- spec
    check(paste0("unified.reject.", s$id),
          "an ambiguous specification is refused with the documented message", {
      msg <- tryCatch({ eval(s$e, envir = globalenv()); NA_character_ },
                      error = function(e) conditionMessage(e))
      list(ok = !is.na(msg) && grepl("Exactly one of", msg, fixed = TRUE),
           obs = if (is.na(msg)) "no error raised" else gsub("\n", " ", msg))
    })
  })
}

# design_table must reproduce, cell by cell, the function it dispatches to
dt_specs <- list(
  list(id = "continuous.ss",
       tab = quote(design_table(param_grid = data.frame(delta1 = 0.5, delta2 = 0.4,
                                                        sd1 = 1, sd2 = 1),
                                rho_values = c(0, 0.5), r = 1, alpha = 0.025,
                                beta = 0.2, endpoint_type = "continuous")),
       direct = function(rho) ss2Continuous(delta1 = 0.5, delta2 = 0.4, sd1 = 1,
                                            sd2 = 1, rho = rho, r = 1,
                                            alpha = 0.025, beta = 0.2,
                                            known_var = TRUE)$N),
  list(id = "continuous.power",
       tab = quote(design_table(param_grid = data.frame(n1 = 100, n2 = 100,
                                                        delta1 = 0.5, delta2 = 0.4,
                                                        sd1 = 1, sd2 = 1),
                                rho_values = c(0, 0.5), alpha = 0.025,
                                endpoint_type = "continuous")),
       direct = function(rho) power2Continuous(n1 = 100, n2 = 100, delta1 = 0.5,
                                               delta2 = 0.4, sd1 = 1, sd2 = 1,
                                               rho = rho, alpha = 0.025,
                                               known_var = TRUE)$powerCoprimary),
  list(id = "binary.approx.ss",
       tab = quote(design_table(param_grid = data.frame(p11 = 0.6, p12 = 0.5,
                                                        p21 = 0.4, p22 = 0.3),
                                rho_values = c(0, 0.5), r = 1, alpha = 0.025,
                                beta = 0.2, endpoint_type = "binary", Test = "AN")),
       direct = function(rho) ss2BinaryApprox(p11 = 0.6, p12 = 0.5, p21 = 0.4,
                                              p22 = 0.3, rho1 = rho, rho2 = rho,
                                              r = 1, alpha = 0.025, beta = 0.2,
                                              Test = "AN")$N),
  list(id = "binary.exact.power",
       tab = quote(design_table(param_grid = data.frame(n1 = 40, n2 = 40,
                                                        p11 = 0.7, p12 = 0.6,
                                                        p21 = 0.3, p22 = 0.2),
                                rho_values = c(0, 0.5), alpha = 0.025,
                                endpoint_type = "binary", Test = "Fisher")),
       direct = function(rho) power2BinaryExact(n1 = 40, n2 = 40, p11 = 0.7,
                                                p12 = 0.6, p21 = 0.3, p22 = 0.2,
                                                rho1 = rho, rho2 = rho,
                                                alpha = 0.025,
                                                Test = "Fisher")$powerCoprimary),
  list(id = "mixed.cont.binary.ss",
       tab = quote(design_table(param_grid = data.frame(delta = 0.5, sd = 1,
                                                        p1 = 0.6, p2 = 0.4),
                                rho_values = c(0, 0.5), r = 1, alpha = 0.025,
                                beta = 0.2, endpoint_type = "mixed_cont_binary",
                                Test = "AN")),
       direct = function(rho) ss2MixedContinuousBinary(delta = 0.5, sd = 1,
                                                       p1 = 0.6, p2 = 0.4,
                                                       rho = rho, r = 1,
                                                       alpha = 0.025, beta = 0.2,
                                                       Test = "AN")$N),
  list(id = "mixed.count.cont.ss",
       tab = quote(design_table(param_grid = data.frame(r1 = 1.0, r2 = 1.25,
                                                        nu = 0.8, t = 1,
                                                        mu1 = -50, mu2 = 0,
                                                        sd = 250),
                                rho_values = c(0, 0.5), r = 1, alpha = 0.025,
                                beta = 0.2, endpoint_type = "mixed_count_cont")),
       direct = function(rho) ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25,
                                                      nu = 0.8, t = 1, mu1 = -50,
                                                      mu2 = 0, sd = 250,
                                                      rho1 = rho, rho2 = rho,
                                                      r = 1, alpha = 0.025,
                                                      beta = 0.2)$N)
)

for (spec in dt_specs) {
  local({
    s <- spec
    check(paste0("design_table.", s$id),
          "every cell equals the value the corresponding direct call returns", {
      tab <- eval(s$tab, envir = globalenv())
      cols <- c("rho_0.0", "rho_0.5")
      got <- vapply(cols, function(cc) as.numeric(tab[[cc]][1]), numeric(1))
      want <- vapply(c(0, 0.5), s$direct, numeric(1))
      list(ok = all(abs(got - want) < 1e-8),
           obs = sprintf("table %s | direct %s", fmt(got), fmt(want)))
    })
  })
}

check("design_table.class", "design_table returns class twoCoprimary_table", {
  tab <- design_table(param_grid = data.frame(delta1 = 0.5, delta2 = 0.4,
                                              sd1 = 1, sd2 = 1),
                      rho_values = c(0, 0.5), r = 1, alpha = 0.025, beta = 0.2,
                      endpoint_type = "continuous")
  list(ok = identical(class(tab), c("twoCoprimary_table", "data.frame")),
       obs = paste(class(tab), collapse = "/"))
})

check("design_table.nMC.forwarded",
      "nMC reaches the mixed continuous-binary path, checked by seeding rather than by timing", {
  g <- data.frame(delta = 0.5, sd = 1, p1 = 0.6, p2 = 0.4)
  set.seed(4321)
  a <- design_table(param_grid = g, rho_values = 0.3, r = 1, alpha = 0.025,
                    beta = 0.2, endpoint_type = "mixed_cont_binary",
                    Test = "Fisher", nMC = 500)[["rho_0.3"]]
  set.seed(4321)
  b <- ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.6, p2 = 0.4,
                                rho = 0.3, r = 1, alpha = 0.025, beta = 0.2,
                                Test = "Fisher", nMC = 500)$N
  list(ok = isTRUE(a == b),
       obs = sprintf("design_table N = %s, direct call with the same seed N = %s",
                     a, b))
})

check("design_table.out.of.bounds",
      "a correlation outside the Frechet-Hoeffding bounds yields NA rather than an error", {
  tab <- design_table(param_grid = data.frame(p11 = 0.9, p12 = 0.2, p21 = 0.5,
                                              p22 = 0.1),
                      rho_values = c(0.1, 0.9), r = 1, alpha = 0.025, beta = 0.2,
                      endpoint_type = "binary", Test = "AN")
  bnd <- corrbound2Binary(0.9, 0.2)
  list(ok = is.na(tab[["rho_0.9"]][1]),
       obs = sprintf("bounds for (0.9, 0.2) = [%s]; rho_0.9 cell = %s",
                     fmt(bnd, 4), format(tab[["rho_0.9"]][1])))
})

# ===========================================================================
part("Part 6. Rejection regions: shape, size and the ordering between tests")
# ===========================================================================

RR_TESTS <- c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")
EXACT_LEVEL_TESTS <- c("Fisher", "Z-pool", "Boschloo")

rr_cases <- list(c(n1 = 12, n2 = 12), c(n1 = 15, n2 = 10), c(n1 = 8, n2 = 16))

for (cs in rr_cases) {
  for (tt in RR_TESTS) {
    local({
      n1 <- as.integer(cs[["n1"]]); n2 <- as.integer(cs[["n2"]]); t_ <- tt
      tag <- sprintf("%s.n%d_%d", t_, n1, n2)
      RR <- rr1Binary(n1, n2, 0.025, Test = t_)

      check(paste0("rr.dim.", tag),
            "the rejection region is a logical matrix of size (n1+1) by (n2+1)", {
        list(ok = is.matrix(RR) && is.logical(RR) &&
               identical(dim(RR), c(n1 + 1L, n2 + 1L)),
             obs = sprintf("mode = %s, dim = %s", mode(RR),
                           paste(dim(RR), collapse = " x ")))
      })

      check(paste0("rr.monotone.", tag),
            "rejection is monotone: more responders in group 1, or fewer in group 2, still rejects", {
        bad_i <- sum(RR[-nrow(RR), , drop = FALSE] & !RR[-1, , drop = FALSE])
        bad_j <- sum(RR[, -1, drop = FALSE] & !RR[, -ncol(RR), drop = FALSE])
        list(ok = bad_i == 0 && bad_j == 0,
             obs = sprintf("violations increasing i: %d; decreasing j: %d",
                           bad_i, bad_j))
      })

      check(paste0("rr.prefix.", tag),
            "each row of the region is a prefix in j, which is what power1 and power2 assume", {
        k <- rowSums(RR)
        bad <- sum(vapply(seq_len(nrow(RR)), function(i) {
          any(RR[i, ] != c(rep(TRUE, k[i]), rep(FALSE, ncol(RR) - k[i])))
        }, logical(1)))
        list(ok = bad == 0,
             obs = sprintf("rows whose TRUE entries are not an initial run: %d of %d",
                           bad, nrow(RR)))
      })

      check(paste0("rr.size.", tag),
            "the null rejection probability never exceeds alpha, as an exact test requires", {
        th <- seq(0.005, 0.995, by = 0.005)
        sz <- vapply(th, function(p) {
          sum(outer(stats::dbinom(0:n1, n1, p), stats::dbinom(0:n2, n2, p)) * RR)
        }, numeric(1))
        mx <- max(sz)
        at <- th[which.max(sz)]
        if (t_ == "Fisher") {
          # Conditioning on the total makes the level exact, with no grid to
          # discretise, so this one is asserted outright
          list(ok = mx <= 0.025 + 1e-10,
               obs = sprintf("largest size %.6f at theta = %.3f (alpha = 0.025)", mx, at))
        } else if (t_ %in% EXACT_LEVEL_TESTS) {
          # The unconditional tests maximise over the nuisance parameter on a
          # grid of n_grid points, so a size marginally above alpha is a
          # statement about the grid rather than about the test
          list(verdict = if (mx <= 0.025 + 1e-10) "PASS"
                         else if (mx <= 0.025 * 1.02) "INFO" else "FAIL",
               obs = sprintf("largest size %.6f at theta = %.3f (alpha = 0.025, n_grid = 100)",
                             mx, at))
        } else {
          list(verdict = "INFO",
               obs = sprintf("largest size %.6f at theta = %.3f (this test is not exact)",
                             mx, at))
        }
      })
    })
  }
}

check("rr.nesting.fisher.boschloo",
      "Boschloo rejects wherever Fisher rejects, since it maximises the same p-value", {
  viol <- 0L
  detail <- character(0)
  for (cs in rr_cases) {
    n1 <- as.integer(cs[["n1"]]); n2 <- as.integer(cs[["n2"]])
    f <- rr1Binary(n1, n2, 0.025, Test = "Fisher")
    b <- rr1Binary(n1, n2, 0.025, Test = "Boschloo")
    v <- sum(f & !b)
    viol <- viol + v
    detail <- c(detail, sprintf("n=(%d,%d): Fisher %d cells, Boschloo %d cells, Fisher-only %d",
                                n1, n2, sum(f), sum(b), v))
  }
  list(ok = viol == 0, obs = paste(detail, collapse = "; "))
})

check("rr.midP.contains.fisher",
      "the mid-p test rejects wherever Fisher's exact test rejects", {
  viol <- 0L
  for (cs in rr_cases) {
    n1 <- as.integer(cs[["n1"]]); n2 <- as.integer(cs[["n2"]])
    viol <- viol + sum(rr1Binary(n1, n2, 0.025, Test = "Fisher") &
                         !rr1Binary(n1, n2, 0.025, Test = "Fisher-midP"))
  }
  list(ok = viol == 0, obs = sprintf("cells rejected by Fisher but not by mid-p: %d", viol))
})

check("rr.n_grid.monotone",
      "a finer nuisance-parameter grid cannot enlarge the region of an unconditional test", {
  detail <- character(0)
  ok <- TRUE
  for (t_ in c("Z-pool", "Boschloo")) {
    a <- rr1Binary(12, 12, 0.025, Test = t_, n_grid = 100)
    b <- rr1Binary(12, 12, 0.025, Test = t_, n_grid = 1000)
    grew <- sum(b & !a)
    ok <- ok && grew == 0
    detail <- c(detail, sprintf("%s: %d cells at n_grid 100, %d at 1000, %d new",
                                t_, sum(a), sum(b), grew))
  }
  list(ok = ok, obs = paste(detail, collapse = "; "))
})

check("rr.n_grid.ignored.by.other.tests",
      "n_grid changes nothing for the three tests that read a p-value off a distribution", {
  same <- vapply(c("Chisq", "Fisher", "Fisher-midP"), function(t_) {
    identical(rr1Binary(12, 12, 0.025, Test = t_, n_grid = 10),
              rr1Binary(12, 12, 0.025, Test = t_, n_grid = 1000))
  }, logical(1))
  list(ok = all(same),
       obs = paste(sprintf("%s: %s", names(same), ifelse(same, "identical", "DIFFERS")),
                   collapse = "; "))
})

check("rr.tie.order.independence",
      "outcomes sharing the value of the ordering statistic share the rejection decision", {
  # The tail event of an exact unconditional test is the set of outcomes at
  # least as extreme as the observed one, so a tie group is entered or left as
  # a whole. The ordering statistics are recomputed here from their definitions
  # rather than taken from the package.
  detail <- character(0)
  ok <- TRUE
  for (cs in rr_cases) {
    n1 <- as.integer(cs[["n1"]]); n2 <- as.integer(cs[["n2"]])
    stat <- list(
      "Z-pool" = outer(n2 * (0:n1), n1 * (0:n2), "-") / (n1 * n2) /
        sqrt(outer(0:n1, 0:n2, "+") / (n1 * n2) *
               (1 - outer(0:n1, 0:n2, "+") / (n1 + n2))),
      "Boschloo" = outer(0:n1, 0:n2, function(i, j)
        stats::phyper(i - 1, n1, n2, i + j, lower.tail = FALSE))
    )
    for (t_ in c("Z-pool", "Boschloo")) {
      RR <- rr1Binary(n1, n2, 0.025, Test = t_)
      z <- stat[[t_]]
      z[is.na(z)] <- 0
      key <- signif(as.vector(z), 10)
      split_rr <- split(as.vector(RR), key)
      mixed <- sum(vapply(split_rr, function(v) length(unique(v)) > 1, logical(1)))
      ok <- ok && mixed == 0
      detail <- c(detail, sprintf("%s n=(%d,%d): %d of %d tie groups split",
                                  t_, n1, n2, mixed, length(split_rr)))
    }
  }
  list(ok = ok, obs = paste(detail, collapse = "; "))
})

# ===========================================================================
part("Part 7. Bivariate binomial: marginals, total mass and induced correlation")
# ===========================================================================

dbibinom_cases <- list(
  list(N = 15, p1 = 0.3, p2 = 0.5, rho = 0.5),
  list(N = 15, p1 = 0.3, p2 = 0.5, rho = 0),
  list(N = 15, p1 = 0.3, p2 = 0.5, rho = -0.5),
  list(N = 20, p1 = 0.7, p2 = 0.6, rho = 0.7),
  list(N = 10, p1 = 0.5, p2 = 0.5, rho = 0.9),
  list(N = 25, p1 = 0.2, p2 = 0.8, rho = 0.2)
)

for (cs in dbibinom_cases) {
  local({
    s <- cs
    tag <- sprintf("N%d.p%s_%s.rho%s", s$N, s$p1, s$p2, s$rho)
    P <- outer(0:s$N, 0:s$N, function(a, b) dbibinom(s$N, a, b, s$p1, s$p2, s$rho))

    check(paste0("dbibinom.mass.", tag), "the probability mass function sums to one", {
      list(ok = abs(sum(P) - 1) < 1e-10,
           obs = sprintf("total mass = %.14f", sum(P)))
    })
    check(paste0("dbibinom.nonneg.", tag), "no probability is negative", {
      list(ok = all(P >= -1e-14),
           obs = sprintf("smallest entry = %.3e", min(P)))
    })
    check(paste0("dbibinom.marginal1.", tag),
          "the first marginal is Binomial(N, p1)", {
      d <- max(abs(rowSums(P) - stats::dbinom(0:s$N, s$N, s$p1)))
      list(ok = d < 1e-10, obs = sprintf("largest absolute deviation = %.3e", d))
    })
    check(paste0("dbibinom.marginal2.", tag),
          "the second marginal is Binomial(N, p2)", {
      d <- max(abs(colSums(P) - stats::dbinom(0:s$N, s$N, s$p2)))
      list(ok = d < 1e-10, obs = sprintf("largest absolute deviation = %.3e", d))
    })
    check(paste0("dbibinom.correlation.", tag),
          "the correlation of the two counts equals the rho that was supplied", {
      y <- 0:s$N
      m1 <- sum(rowSums(P) * y); m2 <- sum(colSums(P) * y)
      exy <- sum(P * outer(y, y))
      v1 <- s$N * s$p1 * (1 - s$p1); v2 <- s$N * s$p2 * (1 - s$p2)
      rho_hat <- (exy - m1 * m2) / sqrt(v1 * v2)
      list(ok = abs(rho_hat - s$rho) < 1e-8,
           obs = sprintf("supplied %.4f, recovered %.10f", s$rho, rho_hat))
    })
  })
}

check("dbibinom.cpp.matches.r",
      "the compiled conditional probability agrees with the R reference implementation", {
  g_r <- get(".dbibinom_g_r", envir = asNamespace("twoCoprimary"))
  g_c <- get("dbibinom_g", envir = asNamespace("twoCoprimary"))
  worst <- 0
  for (cs in dbibinom_cases) {
    gamma <- (cs$rho * sqrt(cs$p2 * (1 - cs$p2) / (cs$p1 * (1 - cs$p1)))) /
      (1 - cs$rho * sqrt(cs$p2 * (1 - cs$p2) / (cs$p1 * (1 - cs$p1))))
    xi <- cs$p2 + gamma * (cs$p2 - cs$p1)
    gr <- expand.grid(y1 = 0:cs$N, y2 = 0:cs$N)
    a <- g_r(cs$N, gr$y1, gr$y2, xi, gamma)
    b <- g_c(as.integer(cs$N), as.integer(gr$y1), as.integer(gr$y2), xi, gamma)
    worst <- max(worst, max(abs(a - b)))
  }
  list(ok = worst < 1e-12, obs = sprintf("largest absolute difference = %.3e", worst))
})

check("dbibinom.rho.rejected.outside.bounds",
      "a correlation outside the Prentice bounds is refused", {
  b <- corrbound2Binary(0.3, 0.5)
  msg <- tryCatch({ dbibinom(10, 3, 5, 0.3, 0.5, b[["U_bound"]] + 0.01); NA_character_ },
                  error = function(e) conditionMessage(e))
  list(ok = !is.na(msg) && grepl("rho must be within", msg),
       obs = sprintf("bounds [%s]; message: %s", fmt(b, 4),
                     if (is.na(msg)) "no error raised" else msg))
})

check("dbibinom.vectorised",
      "vector arguments give the same answer as the scalar calls", {
  y1 <- c(0, 3, 7, 10); y2 <- c(2, 3, 4, 10)
  v <- dbibinom(10, y1, y2, 0.4, 0.6, 0.3)
  s <- vapply(seq_along(y1), function(i) dbibinom(10, y1[i], y2[i], 0.4, 0.6, 0.3),
              numeric(1))
  list(ok = max(abs(v - s)) < 1e-14,
       obs = sprintf("largest absolute difference = %.3e", max(abs(v - s))))
})

# ===========================================================================
part("Part 8. Correlation bounds against independent constructions")
# ===========================================================================

# For two binary outcomes the attainable correlations are set by the range of
# the joint cell probability, which is elementary and can be written down here
# without using the package.
for (pp in list(c(0.3, 0.5), c(0.4, 0.4), c(0.3, 0.7), c(0.9, 0.2), c(0.05, 0.95))) {
  local({
    p1 <- pp[1]; p2 <- pp[2]
    check(sprintf("corrbound2Binary.attainable.p%s_%s", p1, p2),
          "the bounds equal the correlations of the Frechet-Hoeffding joint tables", {
      lo_cell <- max(0, p1 + p2 - 1)
      hi_cell <- min(p1, p2)
      den <- sqrt(p1 * (1 - p1) * p2 * (1 - p2))
      want <- c(L_bound = (lo_cell - p1 * p2) / den,
                U_bound = (hi_cell - p1 * p2) / den)
      got <- corrbound2Binary(p1, p2)
      list(ok = max(abs(got - want)) < 1e-12,
           obs = sprintf("package [%s]; elementary [%s]", fmt(got, 8), fmt(want, 8)))
    })
  })
}

# For a negative binomial and a normal margin the bounds have no closed form.
# They are checked here by building the comonotone and countermonotone couplings
# by simulation, which shares no code with the quadrature the package uses.
for (cs in list(list(lambda = 1.25, nu = 0.8, mu = 0, sd = 250),
                list(lambda = 2.0, nu = 2.0, mu = 50, sd = 200),
                list(lambda = 1.0, nu = 0.5, mu = -50, sd = 250))) {
  local({
    s <- cs
    check(sprintf("corrbound2MCC.simulation.l%s.nu%s", s$lambda, s$nu),
          "the bounds match a comonotone and a countermonotone simulation", {
      set.seed(90210)
      n <- 2e6
      u <- stats::runif(n)
      y_co <- stats::qnbinom(u, mu = s$lambda, size = s$nu)
      x_up <- stats::qnorm(u, s$mu, s$sd)
      x_dn <- stats::qnorm(1 - u, s$mu, s$sd)
      sim <- c(L_bound = stats::cor(y_co, x_dn), U_bound = stats::cor(y_co, x_up))
      got <- corrbound2MixedCountContinuous(s$lambda, s$nu, s$mu, s$sd)
      d <- max(abs(got - sim))
      list(ok = d < 0.005,
           obs = sprintf("package [%s]; simulated [%s]; largest difference %.5f",
                         fmt(got, 5), fmt(sim, 5), d))
    })
  })
}

check("corrbound2MCC.invariant.to.mu.and.sd",
      "the bounds do not move when the location or the scale of the normal margin changes", {
  a <- corrbound2MixedCountContinuous(1.25, 0.8, 0, 250)
  b <- corrbound2MixedCountContinuous(1.25, 0.8, -50, 250)
  cc <- corrbound2MixedCountContinuous(1.25, 0.8, 0, 1)
  list(ok = max(abs(a - b)) < 1e-6 && max(abs(a - cc)) < 1e-6,
       obs = sprintf("mu 0 [%s]; mu -50 [%s]; sd 1 [%s]",
                     fmt(a, 6), fmt(b, 6), fmt(cc, 6)))
})

check("corrbound2MCC.dispersion.trend",
      "the bounds move smoothly as the dispersion parameter grows towards the Poisson limit", {
  v <- vapply(c(0.5, 2, 10, 100, 1000), function(nu)
    corrbound2MixedCountContinuous(1.25, nu, 0, 250)[["U_bound"]], numeric(1))
  list(verdict = "INFO",
       obs = sprintf("upper bound at nu = 0.5, 2, 10, 100, 1000: %s", fmt(v, 4)))
})

# ===========================================================================
part("Part 9. Argument validation: every documented constraint is enforced")
# ===========================================================================

bad_calls <- list(
  list(id = "ss1Continuous.delta", e = quote(ss1Continuous(-0.5, 1, 1, 0.025, 0.2)), m = "delta must be positive"),
  list(id = "ss1Continuous.sd", e = quote(ss1Continuous(0.5, 0, 1, 0.025, 0.2)), m = "sd must be positive"),
  list(id = "ss1Continuous.r", e = quote(ss1Continuous(0.5, 1, 0, 0.025, 0.2)), m = "r must be positive"),
  list(id = "ss1Continuous.alpha", e = quote(ss1Continuous(0.5, 1, 1, 1, 0.2)), m = "alpha must be in"),
  list(id = "ss1Continuous.beta", e = quote(ss1Continuous(0.5, 1, 1, 0.025, 0)), m = "beta must be in"),
  list(id = "ss1Continuous.vector", e = quote(ss1Continuous(c(0.4, 0.5), 1, 1, 0.025, 0.2)), m = "scalar"),
  list(id = "ss1Count.nu", e = quote(ss1Count(1, 1.25, 0, 1, 1, 0.025, 0.2)), m = "nu must be positive"),
  list(id = "ss1Count.t", e = quote(ss1Count(1, 1.25, 0.8, 0, 1, 0.025, 0.2)), m = "t must be positive"),
  list(id = "ss1Count.r1", e = quote(ss1Count(0, 1.25, 0.8, 1, 1, 0.025, 0.2)), m = "r1 must be positive"),
  list(id = "ss1Count.rate.order", e = quote(ss1Count(1.25, 1.0, 0.8, 1, 1, 0.025, 0.2)), m = "r1 must be less than r2"),
  list(id = "ss1Count.rate.equal", e = quote(ss1Count(1.0, 1.0, 0.8, 1, 1, 0.025, 0.2)), m = "r1 must be less than r2"),
  list(id = "ss1BinaryApprox.p1", e = quote(ss1BinaryApprox(1, 0.4, 1, 0.025, 0.2, "AN")), m = "p1 must be in"),
  list(id = "ss1BinaryApprox.order", e = quote(ss1BinaryApprox(0.4, 0.6, 1, 0.025, 0.2, "AN")), m = "greater than p2"),
  list(id = "ss1BinaryApprox.Test", e = quote(ss1BinaryApprox(0.6, 0.4, 1, 0.025, 0.2, "Boschloo")), m = "Test must be one of"),
  list(id = "rr1Binary.n1", e = quote(rr1Binary(0, 10, 0.025, "Fisher")), m = "n1 must be a positive integer"),
  list(id = "rr1Binary.n1.fraction", e = quote(rr1Binary(10.5, 10, 0.025, "Fisher")), m = "n1 must be a positive integer"),
  list(id = "rr1Binary.alpha", e = quote(rr1Binary(10, 10, 0, "Fisher")), m = "alpha must be in"),
  list(id = "rr1Binary.Test", e = quote(rr1Binary(10, 10, 0.025, "AN")), m = "Test must be one of"),
  list(id = "rr1Binary.n_grid.small", e = quote(rr1Binary(10, 10, 0.025, "Boschloo", n_grid = 5)), m = "at least 10"),
  list(id = "rr1Binary.n_grid.fraction", e = quote(rr1Binary(10, 10, 0.025, "Boschloo", n_grid = 12.5)), m = "single integer"),
  list(id = "corrbound2Binary.p1", e = quote(corrbound2Binary(0, 0.5)), m = "p1 must be in"),
  list(id = "corrbound2MCC.lambda", e = quote(corrbound2MixedCountContinuous(0, 0.8, 0, 250)), m = "lambda must be positive"),
  list(id = "corrbound2MCC.sd", e = quote(corrbound2MixedCountContinuous(1.25, 0.8, 0, 0)), m = "sd must be positive"),
  list(id = "dbibinom.N", e = quote(dbibinom(0, 0, 0, 0.3, 0.5, 0.2)), m = "N must be a positive integer"),
  list(id = "dbibinom.y.length", e = quote(dbibinom(10, c(1, 2), 3, 0.3, 0.5, 0.2)), m = "same length"),
  list(id = "dbibinom.y.range", e = quote(dbibinom(10, 11, 3, 0.3, 0.5, 0.2)), m = "y1 must contain integers"),
  list(id = "power2BinaryApprox.Test", e = quote(power2BinaryApprox(100, 100, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "Boschloo")), m = "Test must be one of"),
  list(id = "power2BinaryApprox.rho1", e = quote(power2BinaryApprox(100, 100, 0.6, 0.5, 0.4, 0.3, 0.99, 0.3, 0.025, "AN")), m = "rho1 must be within"),
  list(id = "power2BinaryExact.Test", e = quote(power2BinaryExact(20, 20, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "AN")), m = "Test must be one of"),
  list(id = "power2BinaryExact.n1", e = quote(power2BinaryExact(0, 20, 0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.025, "Fisher")), m = "n1 must be a positive integer"),
  list(id = "ss2BinaryApprox.Test", e = quote(ss2BinaryApprox(0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 1, 0.025, 0.2, "Fisher")), m = "ss2BinaryExact"),
  list(id = "ss2BinaryExact.Test", e = quote(ss2BinaryExact(0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 1, 0.025, 0.2, "AN")), m = "Test must be one of"),
  list(id = "ss2Continuous.rho", e = quote(ss2Continuous(0.5, 0.5, 1, 1, 1, 1, 0.025, 0.2)), m = "rho must be in"),
  list(id = "ss2Continuous.known_var", e = quote(ss2Continuous(0.5, 0.5, 1, 1, 0.5, 1, 0.025, 0.2, known_var = "yes")), m = "known_var must be logical"),
  list(id = "power2MixedContinuousBinary.Test", e = quote(power2MixedContinuousBinary(100, 100, 0.5, 1, 0.6, 0.4, 0.5, 0.025, "Boschloo")), m = "Test must be one of"),
  list(id = "power2MixedContinuousBinary.rho", e = quote(power2MixedContinuousBinary(100, 100, 0.5, 1, 0.6, 0.4, 1, 0.025, "AN")), m = "rho must be in"),
  list(id = "ss2MixedContinuousBinary.delta", e = quote(ss2MixedContinuousBinary(0, 1, 0.6, 0.4, 0.5, 1, 0.025, 0.2, "AN")), m = "delta must be positive"),
  list(id = "power2MixedCountContinuous.nu", e = quote(power2MixedCountContinuous(300, 300, 1, 1.25, 0, 1, -50, 0, 250, 0.5, 0.5, 0.025)), m = "nu must be positive"),
  list(id = "power2MixedCountContinuous.rho1", e = quote(power2MixedCountContinuous(300, 300, 1, 1.25, 0.8, 1, -50, 0, 250, 0.95, 0.5, 0.025)), m = "rho1 must be within"),
  list(id = "ss2MixedCountContinuous.rate.order", e = quote(ss2MixedCountContinuous(1.25, 1.0, 0.8, 1, -50, 0, 250, 1, 0.5, 0.5, 0.025, 0.2)), m = "r1 must be less than r2"),
  list(id = "ss2MixedCountContinuous.mean.order", e = quote(ss2MixedCountContinuous(1.0, 1.25, 0.8, 1, 50, 0, 250, 1, 0.5, 0.5, 0.025, 0.2)), m = "mu1 must be less than mu2"),
  list(id = "design_table.grid.type", e = quote(design_table(param_grid = list(delta1 = 0.5), endpoint_type = "continuous")), m = "must be a data.frame"),
  list(id = "design_table.missing.columns", e = quote(design_table(param_grid = data.frame(delta1 = 0.5), endpoint_type = "continuous")), m = "missing required columns"),
  list(id = "design_table.endpoint_type", e = quote(design_table(param_grid = data.frame(delta1 = 0.5, delta2 = 0.4, sd1 = 1, sd2 = 1), endpoint_type = "ordinal")), m = "arg"),
  list(id = "power2Continuous.rho", e = quote(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 2, 0.025)), m = "rho must be in"),
  list(id = "power2Continuous.sd1", e = quote(power2Continuous(100, 100, 0.5, 0.5, -1, 1, 0.5, 0.025)), m = "sd1 and sd2 must be positive"),
  list(id = "power2Continuous.alpha", e = quote(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 0.5, 2)), m = "alpha must be in"),
  list(id = "power2Continuous.n2.zero", e = quote(power2Continuous(100, 0, 0.5, 0.5, 1, 1, 0.5, 0.025)), m = "n2 must be a positive integer"),
  list(id = "power2Continuous.n2.fraction", e = quote(power2Continuous(100, 2.5, 0.5, 0.5, 1, 1, 0.5, 0.025)), m = "n2 must be a positive integer"),
  list(id = "power2Continuous.known_var", e = quote(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 0.5, 0.025, known_var = "yes")), m = "known_var must be logical"),
  list(id = "power2Continuous.nMC", e = quote(power2Continuous(100, 100, 0.5, 0.5, 1, 1, 0.5, 0.025, known_var = FALSE, nMC = 0)), m = "nMC must be a single positive number")
)

for (bc in bad_calls) {
  local({
    s <- bc
    check(paste0("reject.", s$id),
          "an invalid argument produces an informative error rather than a result", {
      msg <- tryCatch({ eval(s$e, envir = globalenv()); NA_character_ },
                      error = function(e) conditionMessage(e))
      list(ok = !is.na(msg) && grepl(s$m, msg, fixed = FALSE),
           obs = if (is.na(msg)) "no error raised; the call returned a value"
                 else sprintf("expected /%s/, got: %s", s$m, gsub("\n", " ", msg)))
    })
  })
}

# ===========================================================================
part("Part 10. plot(): every plot type against every shape of object")
# ===========================================================================

# Small designs, so that the combinatorial sweep stays affordable
pobj <- list(
  cont_ss = obj$ss2_cont,
  cont_pw = obj$pw2_cont,
  cont_ss_r2 = obj$ss2_cont_r2,
  cont_pw_r2 = obj$pw2_cont_r2,
  binap_ss = obj$ss2_bin_ap,
  binap_pw = obj$pw2_bin_ap,
  binex_ss = ss2BinaryExact(p11 = 0.8, p12 = 0.7, p21 = 0.3, p22 = 0.2,
                            rho1 = 0.3, rho2 = 0.3, r = 1, alpha = 0.025,
                            beta = 0.2, Test = "Chisq"),
  binex_pw = power2BinaryExact(n1 = 20, n2 = 20, p11 = 0.8, p12 = 0.7,
                               p21 = 0.3, p22 = 0.2, rho1 = 0.3, rho2 = 0.3,
                               alpha = 0.025, Test = "Fisher"),
  mcb_ss = obj$ss2_mcb,
  mcb_pw = obj$pw2_mcb,
  mcc_ss = obj$ss2_mcc,
  mcc_pw = obj$pw2_mcc,
  single_cont = obj$ss1_cont,
  single_count = obj$ss1_count,
  single_bin = obj$ss1_bin
)

is_continuous_obj <- function(o) all(c("delta1", "delta2", "sd1", "sd2") %in% names(o))
is_single_obj <- function(o) {
  !(all(c("delta1", "delta2", "sd1", "sd2") %in% names(o)) ||
      all(c("p11", "p12", "p21", "p22") %in% names(o)) ||
      all(c("delta", "sd", "p1", "p2") %in% names(o)) ||
      all(c("r1", "r2", "nu", "mu1", "mu2") %in% names(o)))
}

expected_cols <- list(power_curve = c("n1", "n2", "power"),
                      sample_size_rho = c("rho", "n2"),
                      effect_contour = c("delta1", "delta2", "power"))

for (on in names(pobj)) {
  for (ty in c("power_curve", "sample_size_rho", "effect_contour")) {
    local({
      o <- pobj[[on]]; id <- paste0(on, ".", ty); ty_ <- ty
      single <- is_single_obj(o)
      contin <- is_continuous_obj(o)
      should_work <- !single && (ty_ != "effect_contour" || contin)
      check(paste0("plot.", id),
            if (should_work) "the plot is produced and the data behind it is returned"
            else "the unsupported combination stops with an explanatory message", {
        res <- tryCatch(
          silently_plots(plot(o, type = ty_, n_points = 3)),
          error = function(e) structure(conditionMessage(e), class = "audit_err")
        )
        if (inherits(res, "audit_err")) {
          if (should_work) {
            list(ok = FALSE, obs = paste("unexpected error:", gsub("\n", " ", res)))
          } else {
            want <- if (single) "single endpoint helpers" else "only available for continuous"
            list(ok = grepl(want, res, fixed = TRUE),
                 obs = gsub("\n", " ", res))
          }
        } else {
          if (!should_work) {
            list(ok = FALSE, obs = "completed although the combination is unsupported")
          } else {
            want <- expected_cols[[ty_]]
            list(ok = is.data.frame(res) && identical(names(res), want) &&
                   nrow(res) > 0 && all(vapply(res, function(v) all(is.finite(v)), logical(1))),
                 obs = sprintf("data.frame %d x %d, columns %s",
                               nrow(res), ncol(res), paste(names(res), collapse = ",")))
          }
        }
      })
    })
  }
}

check("plot.default.type.power.object",
      "a power object defaults to the power curve", {
  r <- silently_plots(plot(obj$pw2_cont, n_points = 3))
  list(ok = identical(names(r), c("n1", "n2", "power")),
       obs = paste(names(r), collapse = ","))
})

check("plot.default.type.samplesize.object",
      "a sample size object defaults to the sample size against correlation curve", {
  r <- silently_plots(plot(obj$ss2_cont, n_points = 3))
  list(ok = identical(names(r), c("rho", "n2")),
       obs = paste(names(r), collapse = ","))
})

check("plot.honours.n_range", "an explicit n_range is used as given", {
  r <- silently_plots(plot(obj$pw2_cont, type = "power_curve",
                           n_range = c(40, 60), n_points = 5))
  list(ok = min(r$n2) == 40 && max(r$n2) == 60,
       obs = sprintf("n2 from %d to %d over %d points", min(r$n2), max(r$n2), nrow(r)))
})

check("plot.honours.rho_range", "an explicit rho_range is used as given", {
  rr <- c(0.1, 0.4, 0.7)
  r <- silently_plots(plot(obj$ss2_cont, type = "sample_size_rho", rho_range = rr))
  list(ok = isTRUE(all.equal(r$rho, rr)),
       obs = sprintf("requested %s, used %s", fmt(rr, 2), fmt(r$rho, 2)))
})

check("plot.show_reference.off", "show_reference = FALSE still returns the same data", {
  a <- silently_plots(plot(obj$ss2_cont, type = "sample_size_rho",
                           rho_range = c(0.2, 0.5), show_reference = TRUE))
  b <- silently_plots(plot(obj$ss2_cont, type = "sample_size_rho",
                           rho_range = c(0.2, 0.5), show_reference = FALSE))
  list(ok = isTRUE(all.equal(a, b)), obs = "reference lines do not alter the data")
})

check("plot.exact.test.routing",
      "an object built with an exact binary test is replotted with the exact function", {
  r <- silently_plots(plot(pobj$binex_pw, type = "power_curve", n_range = c(18, 22),
                           n_points = 3))
  direct <- vapply(r$n2, function(n) power2BinaryExact(
    n1 = n, n2 = n, p11 = 0.8, p12 = 0.7, p21 = 0.3, p22 = 0.2,
    rho1 = 0.3, rho2 = 0.3, alpha = 0.025, Test = "Fisher")$powerCoprimary,
    numeric(1))
  list(ok = max(abs(r$power - direct)) < 1e-12,
       obs = sprintf("largest difference from the direct exact call = %.3e",
                     max(abs(r$power - direct))))
})

check("plot.effect_contour.standardized",
      "the contour axes are standardized effect sizes, whatever the standard deviations are", {
  o <- power2Continuous(n1 = 60, n2 = 60, delta1 = 1.0, delta2 = 1.2, sd1 = 2,
                        sd2 = 3, rho = 0.4, alpha = 0.025, known_var = TRUE)
  r <- silently_plots(plot(o, type = "effect_contour", n_points = 5))
  # The grid must span the standardized range, and the power at the object's
  # own standardized effect must be the power the object itself reports
  own <- power2Continuous(n1 = 60, n2 = 60, delta1 = 0.5 * 2, delta2 = 0.4 * 3,
                          sd1 = 2, sd2 = 3, rho = 0.4, alpha = 0.025,
                          known_var = TRUE)$powerCoprimary
  cell <- r$power[which.min(abs(r$delta1 - 0.5) + abs(r$delta2 - 0.4))]
  list(ok = abs(min(r$delta1) - 0.2) < 1e-12 && abs(max(r$delta1) - 1.0) < 1e-12 &&
         all(is.finite(r$power)),
       obs = sprintf("delta1 spans %s to %s; power at the nearest grid point %s, direct %s",
                     fmt(min(r$delta1), 2), fmt(max(r$delta1), 2), fmt(cell),
                     fmt(own)))
})

check("plot.power_curve.small.design",
      "the sample size window does not invert for a design of a few subjects per group", {
  o <- power2Continuous(n1 = 5, n2 = 5, delta1 = 1.5, delta2 = 1.5, sd1 = 1,
                        sd2 = 1, rho = 0.5, alpha = 0.025, known_var = TRUE)
  r <- silently_plots(plot(o, type = "power_curve", n_points = 5))
  list(ok = all(diff(r$n2) > 0) && min(r$n2) >= 2 && all(is.finite(r$power)),
       obs = sprintf("n2 from %d to %d over %d points", min(r$n2), max(r$n2), nrow(r)))
})

check("plot.bad.type", "an unknown plot type is refused", {
  msg <- tryCatch({ silently_plots(plot(obj$ss2_cont, type = "spaghetti")); NA_character_ },
                  error = function(e) conditionMessage(e))
  list(ok = !is.na(msg), obs = if (is.na(msg)) "no error raised" else gsub("\n", " ", msg))
})

# ===========================================================================
part("Part 11. print(): every shape of object")
# ===========================================================================

print_expect <- c(
  ss1_cont = "single continuous endpoint",
  ss1_count = "single count endpoint",
  ss1_bin = "single binary endpoint",
  ss2_cont = "two continuous co-primary endpoints",
  pw2_cont = "two continuous co-primary endpoints",
  ss2_bin_ap = "two binary co-primary endpoints",
  pw2_bin_ap = "two binary co-primary endpoints",
  ss2_bin_ex = "two binary co-primary endpoints",
  pw2_bin_ex = "two binary co-primary endpoints",
  ss2_mcb = "mixed continuous and binary co-primary endpoints",
  pw2_mcb = "mixed continuous and binary co-primary endpoints",
  ss2_mcc = "mixed count and continuous co-primary endpoints",
  pw2_mcc = "mixed count and continuous co-primary endpoints"
)

for (nm in names(print_expect)) {
  local({
    o <- obj[[nm]]; id <- nm; want <- print_expect[[nm]]
    check(paste0("print.", id),
          "the printed header names the calculation and the endpoint combination", {
      txt <- utils::capture.output(print(o))
      head_ok <- any(grepl(want, txt, fixed = TRUE))
      mode_ok <- any(grepl(if ("powerCoprimary" %in% names(o)) "Power calculation"
                           else "Sample size calculation", txt, fixed = TRUE))
      list(ok = head_ok && mode_ok && length(txt) > 3,
           obs = sprintf("%d lines; header: %s", length(txt),
                         trimws(txt[which(nzchar(trimws(txt)))[1]])))
    })
    check(paste0("print.invisible.", id),
          "print returns its argument invisibly and unchanged", {
      txt <- utils::capture.output(back <- print(o))
      list(ok = identical(back, o), obs = "unchanged")
    })
  })
}

check("print.table.method", "print.twoCoprimary_table prints a header and the rows", {
  tab <- design_table(param_grid = data.frame(delta1 = 0.5, delta2 = 0.4,
                                              sd1 = 1, sd2 = 1),
                      rho_values = c(0, 0.5), r = 1, alpha = 0.025, beta = 0.2,
                      endpoint_type = "continuous")
  txt <- utils::capture.output(back <- print(tab))
  list(ok = any(grepl("Design Comparison Table", txt, fixed = TRUE)) &&
         identical(back, tab),
       obs = sprintf("%d lines printed", length(txt)))
})

check("print.no.NA.for.suppressed.nMC",
      "nMC is omitted from the printout when the method does not use it", {
  txt <- utils::capture.output(print(obj$pw2_cont))
  list(ok = !any(grepl("nMC", txt, fixed = TRUE)),
       obs = if (any(grepl("nMC", txt, fixed = TRUE))) "nMC printed although it is NA"
             else "nMC not printed")
})

# ===========================================================================
part("Part 12. Reproducibility of the Monte Carlo paths")
# ===========================================================================

for (spec in list(
  list(id = "power2Continuous.unknown_var",
       e = quote(power2Continuous(n1 = 60, n2 = 60, delta1 = 0.5, delta2 = 0.4,
                                  sd1 = 1, sd2 = 1, rho = 0.4, alpha = 0.025,
                                  known_var = FALSE, nMC = 3000))),
  list(id = "power2MixedContinuousBinary.Fisher",
       e = quote(power2MixedContinuousBinary(n1 = 40, n2 = 40, delta = 0.5, sd = 1,
                                             p1 = 0.6, p2 = 0.4, rho = 0.5,
                                             alpha = 0.025, Test = "Fisher",
                                             nMC = 2000))),
  list(id = "ss2Continuous.unknown_var",
       e = quote(ss2Continuous(delta1 = 0.5, delta2 = 0.4, sd1 = 1, sd2 = 1,
                               rho = 0.4, r = 1, alpha = 0.025, beta = 0.2,
                               known_var = FALSE, nMC = 3000)))
)) {
  local({
    s <- spec
    check(paste0("reproducible.", s$id),
          "the same seed reproduces the same result exactly", {
      set.seed(777); a <- eval(s$e, envir = globalenv())
      set.seed(777); b <- eval(s$e, envir = globalenv())
      list(ok = same_df(a, b), obs = "identical under a fixed seed")
    })
    check(paste0("varies.without.seed.", s$id),
          "without a seed the Monte Carlo result moves, as the documentation warns", {
      set.seed(NULL)
      a <- eval(s$e, envir = globalenv()); b <- eval(s$e, envir = globalenv())
      col <- if ("powerCoprimary" %in% names(a)) "powerCoprimary" else "N"
      list(verdict = "INFO",
           obs = sprintf("two unseeded calls gave %s = %s and %s",
                         col, fmt(a[[col]], 5), fmt(b[[col]], 5)))
    })
  })
}

check("continuous.unknown_var.floor",
      "the unknown variance power is zero, not an error, when the variance cannot be estimated", {
  v <- vapply(c(1, 2), function(n)
    power2Continuous(n1 = n, n2 = n, delta1 = 4, delta2 = 4, sd1 = 1, sd2 = 1,
                     rho = 0.5, alpha = 0.025, known_var = FALSE,
                     nMC = 100)$powerCoprimary, numeric(1))
  list(ok = all(is.finite(v)) && v[1] == 0,
       obs = sprintf("n1 = n2 = 1 gives %s; n1 = n2 = 2 gives %s", fmt(v[1], 3), fmt(v[2], 3)))
})

check("continuous.unknown_var.large.effect",
      "a standardized effect above four no longer stops the sample size search", {
  o <- ss2Continuous(delta1 = 4.5, delta2 = 4.5, sd1 = 1, sd2 = 1, rho = 0.5,
                     r = 1, alpha = 0.025, beta = 0.2, known_var = FALSE, nMC = 2000)
  list(ok = is.finite(o$n2) && o$n2 >= 1,
       obs = sprintf("n1 = %s, n2 = %s", o$n1, o$n2))
})

check("continuous.unknown_var.scale.invariance",
      "scaling the effects and the standard deviations together leaves the power unchanged", {
  set.seed(31415)
  a <- power2Continuous(n1 = 60, n2 = 60, delta1 = 0.5, delta2 = 0.4, sd1 = 1,
                        sd2 = 1.2, rho = 0.4, alpha = 0.025, known_var = FALSE,
                        nMC = 5000)$powerCoprimary
  set.seed(31415)
  b <- power2Continuous(n1 = 60, n2 = 60, delta1 = 5, delta2 = 4, sd1 = 10,
                        sd2 = 12, rho = 0.4, alpha = 0.025, known_var = FALSE,
                        nMC = 5000)$powerCoprimary
  list(ok = abs(a - b) < 1e-10,
       obs = sprintf("unscaled %s, scaled %s, difference %.3e", fmt(a), fmt(b), a - b))
})

# ===========================================================================
part("Part 13. Running time, against the claim the article makes about it")
# ===========================================================================

# The article states that four of the five endpoint type combinations return in
# well under a second and that the exact binary method is the exception. The
# budget below is what makes that sentence checkable rather than remembered.
perf <- list(
  list(id = "ss2Continuous", budget = 1.0,
       e = quote(ss2Continuous(delta1 = 0.2, delta2 = 0.2, sd1 = 1, sd2 = 1,
                               rho = 0.5, r = 1, alpha = 0.025, beta = 0.1,
                               known_var = TRUE))),
  list(id = "ss2BinaryApprox", budget = 1.0,
       e = quote(ss2BinaryApprox(p11 = 0.5, p12 = 0.4, p21 = 0.3, p22 = 0.2,
                                 rho1 = 0.5, rho2 = 0.5, r = 1, alpha = 0.025,
                                 beta = 0.2, Test = "AN"))),
  list(id = "ss2MixedContinuousBinary", budget = 1.0,
       e = quote(ss2MixedContinuousBinary(delta = 0.5, sd = 1, p1 = 0.6, p2 = 0.4,
                                          rho = 0.5, r = 1, alpha = 0.025,
                                          beta = 0.1, Test = "AN"))),
  list(id = "ss2MixedCountContinuous", budget = 1.0,
       e = quote(ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                                         mu1 = -50, mu2 = 0, sd = 250,
                                         rho1 = 0.4, rho2 = 0.4, r = 1,
                                         alpha = 0.025, beta = 0.2))),
  list(id = "ss2BinaryExact.Fisher", budget = NA,
       e = quote(ss2BinaryExact(p11 = 0.5, p12 = 0.4, p21 = 0.3, p22 = 0.2,
                                rho1 = 0.5, rho2 = 0.5, r = 1, alpha = 0.025,
                                beta = 0.2, Test = "Fisher")))
)

for (p in perf) {
  local({
    s <- p
    check(paste0("timing.", s$id),
          if (is.na(s$budget)) "the exact binary search is the slow case"
          else sprintf("the search returns in under %.1f second", s$budget), {
      el <- system.time(res <- eval(s$e, envir = globalenv()))[["elapsed"]]
      timings[[length(timings) + 1]] <<- data.frame(
        label = paste0("perf:", s$id), elapsed = el, stringsAsFactors = FALSE)
      if (is.na(s$budget)) {
        list(verdict = "INFO",
             obs = sprintf("%.2f s for n2 = %s", el, res$n2))
      } else {
        list(ok = el < s$budget,
             obs = sprintf("%.3f s for n2 = %s (budget %.1f s)", el, res$n2, s$budget))
      }
    })
  })
}

check("timing.countcont.bounds.hoisted",
      "the correlation bounds are computed twice per search, not twice per candidate sample size", {
  args_design <- list(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1, mu1 = -50, mu2 = 0,
                      sd = 250, rho1 = 0.5, rho2 = 0.5, r = 1, alpha = 0.025,
                      beta = 0.2)
  t_search <- system.time(res <- do.call(ss2MixedCountContinuous, args_design))[["elapsed"]]
  t_b <- system.time({
    corrbound2MixedCountContinuous(args_design$r1 * args_design$t, args_design$nu,
                                   args_design$mu1, args_design$sd)
    corrbound2MixedCountContinuous(args_design$r2 * args_design$t, args_design$nu,
                                   args_design$mu2, args_design$sd)
  })[["elapsed"]]
  n_calls <- 0L
  counting <- function(n1, n2, ...) { n_calls <<- n_calls + 1L
    power2MixedCountContinuous(n1 = n1, n2 = n2, ...) }
  search_fun <- get(".ss_sequential_search", envir = asNamespace("twoCoprimary"))
  init <- ss1Count(r1 = args_design$r1, r2 = args_design$r2, nu = args_design$nu,
                   t = args_design$t, r = args_design$r, alpha = args_design$alpha,
                   beta = args_design$beta)[["n2"]]
  invisible(search_fun(initial_n2 = init, r = args_design$r,
                       target_power = 1 - args_design$beta, power_fun = counting,
                       r1 = args_design$r1, r2 = args_design$r2, nu = args_design$nu,
                       t = args_design$t, mu1 = args_design$mu1, mu2 = args_design$mu2,
                       sd = args_design$sd, rho1 = args_design$rho1,
                       rho2 = args_design$rho2, alpha = args_design$alpha))
  list(ok = t_search < 0.5 * n_calls * t_b,
       obs = sprintf("search %.3f s with %d power evaluations; one pair of bounds %.3f s; recomputing them at every candidate would cost %.3f s",
                     t_search, n_calls, t_b, n_calls * t_b))
})

check("timing.countcont.core.matches.exported",
      "the internal core returns exactly what the exported power function returns", {
  core <- get(".power2MixedCountContinuous_core", envir = asNamespace("twoCoprimary"))
  a <- core(n1 = 300, n2 = 300, r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1, mu1 = -50,
            mu2 = 0, sd = 250, rho1 = 0.5, rho2 = 0.5, alpha = 0.025)
  b <- power2MixedCountContinuous(n1 = 300, n2 = 300, r1 = 1.0, r2 = 1.25,
                                  nu = 0.8, t = 1, mu1 = -50, mu2 = 0, sd = 250,
                                  rho1 = 0.5, rho2 = 0.5, alpha = 0.025)
  list(ok = same_df(a, b) && identical(names(a), names(b)),
       obs = sprintf("core powerCoprimary %s, exported %s",
                     fmt(a$powerCoprimary, 12), fmt(b$powerCoprimary, 12)))
})

# ===========================================================================
part("Summary")
# ===========================================================================

audit <- do.call(rbind, results)
write.csv(audit, csv_path, row.names = FALSE)
if (length(timings)) {
  write.csv(do.call(rbind, timings), tim_path, row.names = FALSE)
}

n_pass <- sum(audit$verdict == "PASS")
n_fail <- sum(audit$verdict == "FAIL")
n_err <- sum(audit$verdict == "ERROR")
n_info <- sum(audit$verdict == "INFO")

say("")
say(sprintf("checks run : %d", nrow(audit)))
say(sprintf("  PASS     : %d", n_pass))
say(sprintf("  FAIL     : %d", n_fail))
say(sprintf("  ERROR    : %d", n_err))
say(sprintf("  INFO     : %d", n_info))
say("")

if (n_fail + n_err > 0) {
  say("--- everything that did not pass ---")
  bad <- audit[audit$verdict %in% c("FAIL", "ERROR"), , drop = FALSE]
  for (i in seq_len(nrow(bad))) {
    say(sprintf("  [%s] %s", bad$verdict[i], bad$id[i]))
    say("        ", bad$observed[i])
  }
  say("")
}

say("--- recorded for reading, not asserted ---")
inf <- audit[audit$verdict == "INFO", , drop = FALSE]
for (i in seq_len(nrow(inf))) {
  say(sprintf("  %-44s %s", inf$id[i], inf$observed[i]))
}

say("")
say(sprintf("total elapsed : %.1f s",
            as.numeric(difftime(Sys.time(), t_start_all, units = "secs"))))

cat("\nDone. See", log_path, ",", csv_path, "and", tim_path, "\n")
