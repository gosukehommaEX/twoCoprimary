# Audit the claims that the prose of this project makes.
#
# Every other script in dev/ checks the code. This one checks the English, and
# it exists because the code was exhaustively checked while the prose was not,
# so errors kept surfacing in the vignettes, NEWS, README, cran-comments and the
# manuscript long after the code was settled.
#
# Five parts. A to C decide by themselves and fail. D and E cannot be decided by
# a machine and print a list for a person to read; they exist so that nothing is
# missed, not so that everything is judged.
#
#   Part A  Every path a prose file names exists.
#   Part B  Every function name a prose file names is defined by the package.
#   Part C  Every count a prose file states matches the thing it counts.
#   Part D  Sentences shared by two or more files, listed so that a correction
#           made in one place is made in all of them. This is the check that
#           would have caught the n_grid range and the SAS explanation, each of
#           which was corrected in the manuscript and left wrong in a vignette.
#   Part E  Sentences carrying a quantifier or a number, listed per file, to be
#           read against the table or log that produces the value.
#
# Run from the package root with:  source("dev/audit_prose_claims.R")
# The manuscript and the response letter are picked up automatically when the
# working copy sits beside 04_Manuscript and 05_Response.

# ---------------------------------------------------------------- setup ----

out_dir <- "dev/out"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
log_path <- file.path(out_dir, "audit_prose_claims.log")
con <- file(log_path, open = "wt", encoding = "UTF-8")
say <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  writeLines(line, con)
}

n_fail <- 0L
fail <- function(...) {
  n_fail <<- n_fail + 1L
  say("  FAIL  ", ...)
}

repo_root <- normalizePath("../..", mustWork = FALSE)

prose_files <- c(
  list.files("vignettes", "\\.Rmd$", full.names = TRUE),
  list.files("man", "\\.Rd$", full.names = TRUE),
  "NEWS.md", "README.md", "cran-comments.md", "DESCRIPTION",
  file.path(repo_root, "04_Manuscript", "twoCoprimary.Rmd"),
  file.path(repo_root, "05_Response", "response-to-reviewers.tex")
)
prose_files <- prose_files[file.exists(prose_files)]

read_prose <- function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

# Prose only. Code chunks state what the code does and are checked by running it.
strip_code <- function(txt) {
  txt <- gsub("```\\{[^\n]*\n.*?\n```", " ", txt)
  txt <- gsub("```.*?```", " ", txt)
  txt <- gsub("\\\\begin\\{verbatim\\}.*?\\\\end\\{verbatim\\}", " ", txt)
  txt <- gsub("\\\\examples\\{.*", " ", txt)
  txt
}

say("==============================================================================")
say("Claims made by the prose of this project")
say("==============================================================================")
say("run at   : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("files    : ", length(prose_files))
say("")

# --------------------------------------------------- Part A: paths named ----

say("--- Part A: every path named in prose exists ---")
path_pat <- "(?:dev|tests|inst|R|src|man|vignettes|tools|scripts)/[A-Za-z0-9_./-]+"
n_paths <- 0L
for (f in prose_files) {
  hits <- unique(regmatches(read_prose(f), gregexpr(path_pat, read_prose(f), perl = TRUE))[[1]])
  for (h in hits) {
    p <- sub("[.,;:)]+$", "", h)
    n_paths <- n_paths + 1L
    if (!file.exists(p) && !dir.exists(p)) fail(basename(f), " names ", p, ", which does not exist")
  }
}
say("  ", n_paths, " path references checked")
say("")

# ----------------------------------------------- Part B: functions named ----

say("--- Part B: every function named in prose is defined here or is base R ---")
exported <- gsub('"', "", regmatches(
  readLines("NAMESPACE", warn = FALSE),
  regexpr("(?<=export\\().*(?=\\))", readLines("NAMESPACE", warn = FALSE), perl = TRUE)
))
exported <- exported[nzchar(exported)]
internal <- unlist(lapply(list.files("R", "\\.R$", full.names = TRUE), function(f) {
  l <- readLines(f, warn = FALSE)
  m <- regmatches(l, regexpr("^\\.?[A-Za-z_][A-Za-z0-9._]*(?=\\s*<-\\s*function)", l, perl = TRUE))
  m
}))
defined <- unique(c(exported, internal))

is_known <- function(nm) {
  nm %in% defined ||
    any(vapply(search(), function(env) exists(nm, where = env, inherits = FALSE), logical(1)))
}

n_names <- 0L
for (f in prose_files) {
  txt <- strip_code(read_prose(f))
  nms <- unique(c(
    regmatches(txt, gregexpr("(?<=`)[A-Za-z_][A-Za-z0-9._]*(?=\\(\\)`)", txt, perl = TRUE))[[1]],
    regmatches(txt, gregexpr("(?<=\\\\code\\{)[A-Za-z_][A-Za-z0-9._]*(?=\\(\\)\\})", txt, perl = TRUE))[[1]],
    regmatches(txt, gregexpr("(?<=\\\\texttt\\{)[A-Za-z_][A-Za-z0-9._]*(?=\\(\\)\\})", txt, perl = TRUE))[[1]]
  ))
  for (nm in nms) {
    n_names <- n_names + 1L
    if (!is_known(nm)) fail(basename(f), " names ", nm, "(), which is not defined")
  }
}
say("  ", n_names, " function references checked")
say("")

# ------------------------------------------------------ Part C: counts ------

say("--- Part C: counts stated in prose against the thing they count ---")

news <- read_prose("NEWS.md")
sec_111 <- sub("(?s)# twoCoprimary 1\\.1\\.0.*", "", news, perl = TRUE)
bug_bullets <- lengths(regmatches(
  sub("(?s)## Performance.*", "", sec_111, perl = TRUE),
  gregexpr("\n\\* ", sub("(?s)## Performance.*", "", sec_111, perl = TRUE))
))

check_count <- function(label, stated, actual) {
  if (is.na(stated)) {
    say("  skip  ", label, " (not stated)")
  } else if (stated != actual) {
    fail(label, ": prose says ", stated, ", actual is ", actual)
  } else {
    say("  ok    ", label, " = ", actual)
  }
}

stated_in <- function(path, pattern) {
  if (!file.exists(path)) return(NA_integer_)
  m <- regmatches(read_prose(path), regexpr(pattern, read_prose(path), perl = TRUE))
  if (!length(m)) return(NA_integer_)
  as.integer(regmatches(m, regexpr("[0-9]+", m)))
}

check_count("defects fixed in 1.1.1, cran-comments vs NEWS bullets",
            stated_in("cran-comments.md", "fixes \\d+ defects"), bug_bullets)
check_count("exported functions, manuscript vs NAMESPACE",
            stated_in(file.path(repo_root, "04_Manuscript", "twoCoprimary.Rmd"),
                      "\\d+ exported functions"),
            length(exported))
check_count("vignettes, manuscript vs vignettes/",
            stated_in(file.path(repo_root, "04_Manuscript", "twoCoprimary.Rmd"),
                      "\\d+ vignettes"),
            length(list.files("vignettes", "\\.Rmd$")))

# Counts that come from a reproduction log rather than from a directory.
cell_counts <- list(
  "Homma and Yoshida 2025 Table 4" = "reproduce_Homma_and_Yoshida_2025.csv",
  "Sozu 2012 Table S5"             = "reproduce_Sozu_et_al_2012.csv",
  "Sozu 2010 Table III"            = "reproduce_Sozu_et_al_2010.csv"
)
for (tab in names(cell_counts)) {
  f <- file.path(out_dir, cell_counts[[tab]])
  if (!file.exists(f)) {
    say("  skip  ", tab, " (run the reproduction script first)")
    next
  }
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  say("  info  ", tab, ": ", sum(d$table == tab), " cells, ",
      sum(d$table == tab & d$diff != 0), " differ")
}
say("")

# ---------------------------------------- Part D: sentences shared by files ----

say("--- Part D: sentences that appear in more than one file ---")
say("  A correction applied to one copy and not the other is how the same error")
say("  survived in the manuscript and a vignette. Read each group and confirm the")
say("  copies still agree.")
say("")

sentences_of <- function(f) {
  txt <- strip_code(read_prose(f))
  txt <- gsub("\\s+", " ", txt)
  s <- unlist(strsplit(txt, "(?<=[.]) (?=[A-Z])", perl = TRUE))
  s <- trimws(s)
  s[nchar(s) >= 60]
}
norm <- function(s) tolower(gsub("[^a-z0-9 ]", "", tolower(s)))

all_s <- do.call(rbind, lapply(prose_files, function(f) {
  s <- sentences_of(f)
  if (!length(s)) return(NULL)
  data.frame(file = basename(f), sentence = s, key = norm(s), stringsAsFactors = FALSE)
}))
dup_keys <- unique(all_s$key[duplicated(all_s$key)])
if (!length(dup_keys)) {
  say("  none")
} else {
  for (k in dup_keys) {
    grp <- all_s[all_s$key == k, ]
    if (length(unique(grp$file)) < 2) next
    say("  * ", paste(unique(grp$file), collapse = " + "))
    say("    ", substr(grp$sentence[1], 1, 200))
  }
}
say("")

# ------------------------------ Part E: sentences to be read against a source ----

say("--- Part E: sentences carrying a number or a quantifier ---")
say("  Read each against the table or the log that produces the value. This list")
say("  is generated, so it cannot omit a sentence; the judgement stays human.")
say("")

quant <- paste0("\\b(all|every|always|never|none|only|exactly|identical|unchanged|",
                "smallest|largest|each|reproduces?|agrees?|no more than|at most|at least)\\b")
numeric_claim <- "(?<![A-Za-z0-9._])[0-9]+(\\.[0-9]+)?(?![A-Za-z0-9._])"

for (f in prose_files) {
  s <- sentences_of(f)
  keep <- grepl(quant, s, perl = TRUE, ignore.case = TRUE) |
    grepl(numeric_claim, s, perl = TRUE)
  s <- s[keep]
  if (!length(s)) next
  say("  ", basename(f), "  (", length(s), " sentences)")
  for (x in s) say("      - ", substr(x, 1, 200))
  say("")
}

# ---------------------------------------------------------------- summary ----

say("==============================================================================")
say("Summary")
say("==============================================================================")
say("  automatic failures (Parts A to C) : ", n_fail)
say("  Parts D and E are lists for a person to read, not verdicts.")
say("")
say("  log written to ", log_path)
close(con)

if (n_fail > 0L) stop("audit_prose_claims: ", n_fail, " failure(s); see ", log_path)
invisible(NULL)
