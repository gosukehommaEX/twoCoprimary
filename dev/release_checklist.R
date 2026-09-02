# Final consistency check before the CRAN submission.
#
# Part 0  The library that library() attached was built from the sources as they
#         stand now. Every other part of this script, and every other script in
#         dev/, reads the installed build, so a failed or forgotten reinstall
#         makes all of them report on code that is no longer on disk. The
#         version number cannot catch this, because a rebuild during a release
#         carries the same version.
# Part A  The three worked examples printed in README.md still produce the
#         numbers shown there. README is plain markdown, so nothing regenerates
#         it and a stale number would go unnoticed.
# Part B  Every citation that appears in the package documentation matches the
#         reference list, checked mechanically rather than by eye.
# Part C  The version, the NEWS heading and the DESCRIPTION agree.
# Part D  Every script in dev/ is described in dev/README.md and every script
#         described there exists, so that the index cannot drift.
# Part E  Every log in dev/out/ was written by the installed version. A result
#         file from an earlier version is the one way a number quoted in the
#         article or the response letter can be silently stale, so this part
#         names any that are, and they are then regenerated rather than read.
#
# Run from the package root:
#   source("dev/release_checklist.R")
#
# Writes dev/out/release_check.log

library(twoCoprimary)

out_dir <- file.path("dev", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "release_check.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_con)
}

say("twoCoprimary version: ", as.character(utils::packageVersion("twoCoprimary")))
say("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("")

# ---------------------------------------------------------------------------
# Part 0: the attached build is newer than the sources
# ---------------------------------------------------------------------------

say("--- Part 0: the attached build was made from the current sources ---")

lib_path <- find.package("twoCoprimary")
meta_rds <- file.path(lib_path, "Meta", "package.rds")
if (!file.exists(meta_rds)) {
  say("  The attached twoCoprimary is a development load, not an installed ",
      "build.")
  say("  Restart R, install the package, and run this script again.")
  close(log_con)
  stop("development load attached; install the package and restart R")
}

built_at <- file.mtime(meta_rds)
src_files <- c(
  list.files("R", pattern = "\\.R$", full.names = TRUE),
  list.files("src", pattern = "\\.(cpp|h)$", full.names = TRUE),
  list.files("man", pattern = "\\.Rd$", full.names = TRUE),
  c("DESCRIPTION", "NAMESPACE")
)
src_files <- src_files[file.exists(src_files)]
src_time <- file.mtime(src_files)
newer <- src_files[src_time > built_at]

say("  build attached from  : ", lib_path)
say("  built at             : ", format(built_at, "%Y-%m-%d %H:%M:%S"))
say("  newest source file   : ",
    format(max(src_time), "%Y-%m-%d %H:%M:%S"))
say("  sources newer than it: ", length(newer))
if (length(newer)) {
  for (f in newer) say("    ", f)
  say("  *** the attached build is stale; nothing below describes the current ",
      "code ***")
  close(log_con)
  stop("the attached build is older than ", length(newer),
       " source file(s); document, install, restart R and re-run every script")
}
say("")

# ---------------------------------------------------------------------------
# Part A: the README examples
# ---------------------------------------------------------------------------

say("--- Part A: worked examples printed in README.md ---")

readme <- data.frame(
  example = c("two continuous", "two binary, exact Fisher", "count and continuous"),
  n1_readme = c(79, 111, 705),
  n2_readme = c(79, 111, 705),
  N_readme = c(158, 222, 1410),
  stringsAsFactors = FALSE
)

r1 <- ss2Continuous(delta1 = 0.5, delta2 = 0.5, sd1 = 1, sd2 = 1,
                    rho = 0.5, r = 1, alpha = 0.025, beta = 0.2,
                    known_var = TRUE)
r2 <- ss2BinaryExact(p11 = 0.50, p12 = 0.30, p21 = 0.30, p22 = 0.10,
                     rho1 = 0.3, rho2 = 0.3, r = 1,
                     alpha = 0.025, beta = 0.2, Test = "Fisher")
r3 <- ss2MixedCountContinuous(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1,
                              mu1 = -50, mu2 = 0, sd = 250, r = 1,
                              rho1 = 0.5, rho2 = 0.5, alpha = 0.025, beta = 0.2)

readme$n1_now <- c(r1$n1, r2$n1, r3$n1)
readme$n2_now <- c(r1$n2, r2$n2, r3$n2)
readme$N_now <- c(r1$N, r2$N, r3$N)
readme$agrees <- with(readme, n1_now == n1_readme & n2_now == n2_readme &
                              N_now == N_readme)

print(readme[, c("example", "N_readme", "N_now", "agrees")])
capture.output(print(readme[, c("example", "N_readme", "N_now", "agrees")]),
               file = log_con, append = TRUE)
say("")
say("README examples still correct : ", sum(readme$agrees), " of ", nrow(readme))
say("")

# ---------------------------------------------------------------------------
# Part B: citations
# ---------------------------------------------------------------------------

say("--- Part B: citations in the installed help files ---")

expected <- c(
  "Statistics in Medicine}, 29(21), 2169-2179",
  "Journal of Biopharmaceutical Statistics}, 21(4), 650-668",
  "Biometrical Journal}, 54(5), 716-729",
  "Pharmaceutical Statistics}, 23(1), 46-59",
  "Medical Research}, 34(11), 2183-2201"
)
wrong <- c(
  "34(1), 1-19", "2219-2227", "23(3), 368-392", "Kanou",
  "Japanese Journal of Biometrics", "Evans, S. R."
)

# Rd objects keep the markup as separate tokens, so the text has to be
# flattened and its whitespace normalised before any citation string can be
# matched against it.
flatten_rd <- function(x) {
  out <- character(0)
  walk <- function(node) {
    if (is.list(node)) {
      for (el in node) walk(el)
    } else {
      out <<- c(out, as.character(node))
    }
  }
  walk(x)
  gsub("\\s+", " ", paste(out, collapse = ""))
}

# Rd_db reads the installed help database, which is stale until the package has
# been reinstalled after a documentation change. Fail with an instruction rather
# than the raw "installed help is corrupt" error.
db <- tryCatch(tools::Rd_db("twoCoprimary"), error = function(e) {
  message("Cannot read the installed help database: ", conditionMessage(e))
  message("Run devtools::document(), reinstall the package, restart the R ",
          "session, and run this script again.")
  NULL
})
if (is.null(db)) {
  close(log_con)
  stop("stale help database; document, reinstall, restart R and re-run")
}
txt <- paste(vapply(db, flatten_rd, character(1)), collapse = " ")

for (e in expected) {
  hit <- grepl(gsub("\\}", "", e), gsub("[{}]", "", txt), fixed = TRUE)
  say("  present : ", formatC(gsub("\\}", "", e), width = -54), hit,
      if (hit) "" else "   *** MISSING ***")
}
say("")
for (w in wrong) {
  hit <- grepl(w, gsub("[{}]", "", txt), fixed = TRUE)
  say("  absent  : ", formatC(w, width = -54), !hit,
      if (hit) "   *** STILL PRESENT ***" else "")
}
say("")
say("help topics scanned : ", length(db))
say("")

# ---------------------------------------------------------------------------
# Part C: version consistency
# ---------------------------------------------------------------------------

say("--- Part C: version consistency ---")

desc_version <- as.character(utils::packageVersion("twoCoprimary"))
news_path <- system.file("NEWS.md", package = "twoCoprimary")
if (news_path == "") news_path <- "NEWS.md"
news_first <- readLines(news_path, n = 1, warn = FALSE)

say("  DESCRIPTION version : ", desc_version)
say("  NEWS first heading  : ", news_first)
say("  they agree          : ",
    grepl(desc_version, news_first, fixed = TRUE))
say("  NeedsCompilation    : ",
    length(list.files(file.path(find.package("twoCoprimary"), "libs"),
                      recursive = TRUE)) > 0)
say("")

# ---------------------------------------------------------------------------
# Part D: dev/README.md against the contents of dev/
# ---------------------------------------------------------------------------

say("--- Part D: dev/README.md is a complete index of dev/ ---")

readme_path <- file.path("dev", "README.md")
scripts <- sort(setdiff(list.files("dev", pattern = "\\.R$"), character(0)))

if (!file.exists(readme_path)) {
  say("  *** dev/README.md is missing ***")
} else {
  readme <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")
  documented <- unique(regmatches(readme, gregexpr("[A-Za-z0-9_]+\\.R\\b", readme))[[1]])
  undocumented <- setdiff(scripts, documented)
  phantom <- setdiff(documented, scripts)
  say("  scripts in dev/           : ", length(scripts))
  say("  named in README.md        : ", length(intersect(scripts, documented)))
  say("  present but not described : ",
      if (length(undocumented)) paste(undocumented, collapse = ", ") else "none")
  say("  described but absent      : ",
      if (length(phantom)) paste(phantom, collapse = ", ") else "none")
  if (length(undocumented) || length(phantom)) {
    say("  *** dev/README.md is out of date ***")
  }
}
say("")

# ---------------------------------------------------------------------------
# Part E: every result in dev/out/ came from the installed version
# ---------------------------------------------------------------------------

say("--- Part E: results in dev/out/ are current ---")

logs <- list.files(out_dir, pattern = "\\.log$", full.names = TRUE)
logs <- setdiff(logs, file.path(out_dir, "release_check.log"))

stale <- character(0)
unstamped <- character(0)
for (f in logs) {
  head_txt <- readLines(f, n = 8, warn = FALSE)
  hit <- grep("twoCoprimary version", head_txt, value = TRUE)
  if (!length(hit)) {
    unstamped <- c(unstamped, basename(f))
  } else if (!grepl(desc_version, hit[1], fixed = TRUE)) {
    stale <- c(stale, sprintf("%s (%s)", basename(f),
                              trimws(sub(".*version[^0-9]*", "", hit[1]))))
  }
}

say("  logs found                : ", length(logs))
say("  written by ", desc_version, "        : ",
    length(logs) - length(stale) - length(unstamped))
say("  written by another version: ",
    if (length(stale)) paste(stale, collapse = ", ") else "none")
say("  without a version stamp   : ",
    if (length(unstamped)) paste(unstamped, collapse = ", ") else "none")
if (length(stale) || length(unstamped)) {
  say("  *** do not quote a number from those files; delete dev/out/ and rerun ***")
}

close(log_con)
