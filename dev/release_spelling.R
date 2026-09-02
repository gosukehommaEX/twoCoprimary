# Triage the output of spelling::spell_check_package() so that no word has to be
# judged by eye. Every flagged word is placed in exactly one of three buckets.
#
#   A  already in inst/WORDLIST          no action
#   B  occurs only inside math or code   LaTeX markup, not prose
#   C  occurs in prose                   the only bucket that needs judgment
#
# Bucket B exists because spelling reads the Rmd source, where \frac, \geq and
# friends are indistinguishable from words. Bucket C is what a reviewer reads.

pkg <- "."

# Remove code chunks, display math and inline math from a vector of Rmd lines,
# leaving only prose.
strip_math <- function(lines) {
  keep <- character(0)
  in_chunk <- FALSE
  in_display <- FALSE
  for (l in lines) {
    if (grepl("^\\s*```", l)) {
      in_chunk <- !in_chunk
      next
    }
    if (in_chunk) next
    n_dd <- lengths(regmatches(l, gregexpr("\\$\\$", l)))
    if (in_display) {
      if (n_dd > 0) in_display <- FALSE
      next
    }
    if (n_dd %% 2 == 1) {
      in_display <- TRUE
      l <- sub("\\$\\$.*$", "", l)
    } else if (n_dd > 0) {
      l <- gsub("\\$\\$[^$]*\\$\\$", " ", l)
    }
    l <- gsub("\\$[^$]*\\$", " ", l)     # inline math
    l <- gsub("`[^`]*`", " ", l)         # inline code
    keep <- c(keep, l)
  }
  keep
}

prose_corpus <- function() {
  rmd <- list.files(file.path(pkg, "vignettes"), pattern = "\\.Rmd$",
                    full.names = TRUE)
  txt <- unlist(lapply(rmd, function(f)
    strip_math(readLines(f, warn = FALSE, encoding = "UTF-8"))))
  other <- c(file.path(pkg, c("README.md", "NEWS.md", "DESCRIPTION")),
             list.files(file.path(pkg, "man"), pattern = "\\.Rd$",
                        full.names = TRUE))
  other <- other[file.exists(other)]
  for (f in other) {
    txt <- c(txt, readLines(f, warn = FALSE, encoding = "UTF-8"))
  }
  txt
}

wl_path <- file.path(pkg, "inst", "WORDLIST")
wordlist <- if (file.exists(wl_path)) {
  readLines(wl_path, warn = FALSE, encoding = "UTF-8")
} else {
  character(0)
}

flagged <- spelling::spell_check_package(pkg)
words <- flagged$word
corpus <- prose_corpus()

in_prose <- vapply(words, function(w) {
  pat <- paste0("(^|[^A-Za-z\\\\])", gsub("([.\\\\|()\\[{^$*+?])", "\\\\\\1", w),
                "([^A-Za-z]|$)")
  any(grepl(pat, corpus, perl = TRUE))
}, logical(1))

bucket <- ifelse(words %in% wordlist, "A",
                 ifelse(in_prose, "C", "B"))

cat("flagged words :", length(words), "\n")
cat("  A already in WORDLIST :", sum(bucket == "A"), "\n")
cat("  B math or code only   :", sum(bucket == "B"), "\n")
cat("  C prose, needs review :", sum(bucket == "C"), "\n\n")

if (any(bucket == "B")) {
  cat("--- B: LaTeX markup, safe to accept ---\n")
  cat(strwrap(paste(sort(words[bucket == "B"]), collapse = ", "), width = 76),
      sep = "\n")
  cat("\n")
}

if (any(bucket == "C")) {
  cat("--- C: appears in prose, read these ---\n")
  for (w in sort(words[bucket == "C"])) {
    hit <- grep(paste0("(^|[^A-Za-z\\\\])",
                       gsub("([.\\\\|()\\[{^$*+?])", "\\\\\\1", w),
                       "([^A-Za-z]|$)"), corpus, perl = TRUE, value = TRUE)[1]
    cat(sprintf("  %-18s %s\n", w, substr(trimws(hit), 1, 55)))
  }
} else {
  cat("--- C is empty: no prose word needs review ---\n")
}
