# Shared machinery for the five dev/reproduce_*.R scripts, each of which
# reproduces the published tables of one of the source articles.
#
# The scripts are not a substitute for dev/audit_all_functions.R. That file asks
# whether the package is internally consistent; these ask whether it returns the
# numbers that the peer-reviewed articles it implements actually printed. A
# package can be perfectly self-consistent and still implement the wrong
# formula, and only an external table can catch that.
#
# Every published value is typed in from the article and carried in the script,
# so the comparison is made against the paper rather than against a previous run
# of the package. Each comparison is written out as the published value, the
# computed value and their difference, and a tolerance is stated for every
# table. Where a difference is expected, for instance because the original was
# computed in SAS, the tolerance says so rather than the prose.
#
# This file defines no top level analysis; it is sourced by the scripts.

.rep <- new.env(parent = emptyenv())

say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  if (!is.null(.rep$log)) cat(txt, "\n", sep = "", file = .rep$log, append = TRUE)
}

reproduce_begin <- function(stem, article, tables) {
  .rep$stem <- stem
  out_dir <- file.path("dev", "out")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  .rep$log <- file.path(out_dir, paste0(stem, ".log"))
  .rep$csv <- file.path(out_dir, paste0(stem, ".csv"))
  if (file.exists(.rep$log)) file.remove(.rep$log)
  .rep$rows <- list()
  .rep$t0 <- Sys.time()
  say(strrep("=", 78))
  say("Reproduction of the published tables of")
  say("  ", article)
  say(strrep("=", 78))
  say("tables covered      : ", paste(tables, collapse = "; "))
  say("twoCoprimary version: ", as.character(utils::packageVersion("twoCoprimary")))
  say("R version           : ", R.version.string)
  say("run at              : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  say("")
}

# Record one comparison. published and computed are single numbers; tol is the
# largest difference that still counts as agreement for this table.
cmp <- function(table_id, cell, published, computed, tol = 0,
                quantity = "sample size") {
  d <- suppressWarnings(as.numeric(computed) - as.numeric(published))
  ok <- is.finite(d) && abs(d) <= tol + 1e-9
  .rep$rows[[length(.rep$rows) + 1L]] <- data.frame(
    table = table_id, cell = cell, quantity = quantity,
    published = as.numeric(published), computed = as.numeric(computed),
    difference = d, tolerance = tol,
    verdict = if (!is.finite(d)) "MISSING" else if (ok) "MATCH" else "DIFFERS",
    stringsAsFactors = FALSE
  )
  invisible(ok)
}

# Print the comparisons of one table as published against computed, and give the
# counts and the largest deviation. Nothing here is judged by eye: the verdict
# column comes from the tolerance stated when the values were recorded.
table_report <- function(table_id, note = NULL) {
  all_rows <- do.call(rbind, .rep$rows)
  d <- all_rows[all_rows$table == table_id, , drop = FALSE]
  say(strrep("-", 78))
  say("Table: ", table_id)
  if (!is.null(note)) say("  ", note)
  say(strrep("-", 78))
  say(sprintf("  %-46s %10s %10s %6s", "cell", "published", "computed", "diff"))
  for (i in seq_len(nrow(d))) {
    mark <- if (d$verdict[i] == "MATCH") " " else "*"
    say(sprintf("%s %-46s %10s %10s %6s", mark, d$cell[i],
                format(d$published[i]), format(d$computed[i]),
                if (is.finite(d$difference[i])) format(d$difference[i]) else "NA"))
  }
  fin <- d[is.finite(d$difference), , drop = FALSE]
  say("")
  tol_txt <- if (length(unique(d$tolerance)) == 1L) format(d$tolerance[1]) else
    sprintf("%s to %s", format(min(d$tolerance)), format(max(d$tolerance)))
  say(sprintf("  %d cells: %d match, %d differ, %d missing; tolerance %s; largest deviation %s",
              nrow(d), sum(d$verdict == "MATCH"), sum(d$verdict == "DIFFERS"),
              sum(d$verdict == "MISSING"), tol_txt,
              if (nrow(fin)) format(max(abs(fin$difference))) else "NA"))
  say("")
}

reproduce_end <- function() {
  all_rows <- do.call(rbind, .rep$rows)
  write.csv(all_rows, .rep$csv, row.names = FALSE)
  say(strrep("=", 78))
  say("Summary")
  say(strrep("=", 78))
  for (tb in unique(all_rows$table)) {
    d <- all_rows[all_rows$table == tb, , drop = FALSE]
    fin <- d[is.finite(d$difference), , drop = FALSE]
    say(sprintf("  %-44s %3d cells, %3d match, %3d differ, largest deviation %s",
                tb, nrow(d), sum(d$verdict == "MATCH"),
                sum(d$verdict == "DIFFERS"),
                if (nrow(fin)) format(max(abs(fin$difference))) else "NA"))
  }
  say("")
  say(sprintf("  overall : %d cells, %d match, %d differ, %d missing",
              nrow(all_rows), sum(all_rows$verdict == "MATCH"),
              sum(all_rows$verdict == "DIFFERS"),
              sum(all_rows$verdict == "MISSING")))
  bad <- all_rows[all_rows$verdict != "MATCH", , drop = FALSE]
  if (nrow(bad)) {
    say("")
    say("  cells outside the stated tolerance:")
    for (i in seq_len(nrow(bad))) {
      say(sprintf("    %-30s %-40s published %s, computed %s",
                  bad$table[i], bad$cell[i], format(bad$published[i]),
                  format(bad$computed[i])))
    }
  }
  say("")
  say(sprintf("  elapsed : %.1f s",
              as.numeric(difftime(Sys.time(), .rep$t0, units = "secs"))))
  cat("\nDone. See", .rep$log, "and", .rep$csv, "\n")
}

# Recover the correlation between the two test statistics that a power result
# implies, by inverting the bivariate normal probability. The package does not
# return this quantity, and two of the articles tabulate it.
implied_gamma <- function(marginal1, marginal2, joint) {
  c1 <- stats::qnorm(marginal1)
  c2 <- stats::qnorm(marginal2)
  f <- function(g) pbivnorm::pbivnorm(x = c1, y = c2, rho = g) - joint
  if (f(-0.999) > 0 || f(0.999) < 0) return(NA_real_)
  stats::uniroot(f, c(-0.999, 0.999), tol = 1e-10)$root
}
