# Audit the claims that the prose of this project makes.
#
# Every other script in dev/ checks the code. This one checks the English, and
# it exists because the code was swept exhaustively while the prose was not, so
# errors kept surfacing in the vignettes, NEWS, README, cran-comments and the
# manuscript long after the code was settled.
#
#   Part A  Every path a prose file names exists.
#   Part B  Every function name a prose file names is defined here or is base R.
#   Part C  Every count a prose file states matches the thing it counts.
#   Part D  Sentences shared by two or more files, so that a correction made in
#           one place is made in all of them. This is the check that would have
#           caught the n_grid range and the SAS explanation, each corrected in
#           the manuscript and left wrong in a vignette.
#   Part E  Sentences carrying a number or a quantifier, listed per file, to be
#           read against the table or log that produces the value.
#
# Parts A to C decide by themselves and stop the script. Parts D and E are lists
# for a person to read; they exist so that no sentence is missed, not so that
# every sentence is judged.
#
# Run from the package root with:  source("dev/audit_prose_claims.R")
# The manuscript and the response letter are picked up when the working copy
# sits beside 04_Manuscript and 05_Response.

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
manuscript <- file.path(repo_root, "04_Manuscript", "twoCoprimary.Rmd")
response <- file.path(repo_root, "05_Response", "response-to-reviewers.tex")

prose_files <- c(
  list.files("vignettes", "\\.Rmd$", full.names = TRUE),
  list.files("man", "\\.Rd$", full.names = TRUE),
  "NEWS.md", "README.md", "cran-comments.md", "DESCRIPTION",
  manuscript, response
)
prose_files <- prose_files[file.exists(prose_files)]

read_prose <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

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
say("run at : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("files  : ", length(prose_files))
if (!file.exists(manuscript)) say("note   : manuscript not found, its checks are skipped")
if (!file.exists(response)) say("note   : response letter not found, its checks are skipped")
say("")

# ------------------------------------------------- Part A: paths named ----

say("--- Part A: every path named in prose exists ---")
path_pat <- "(?:dev|tests|inst|R|src|man|vignettes|tools|scripts)/[A-Za-z0-9_./-]+"
n_paths <- 0L
for (f in prose_files) {
  txt <- read_prose(f)
  hits <- unique(regmatches(txt, gregexpr(path_pat, txt, perl = TRUE))[[1]])
  for (h in hits) {
    p <- sub("[.,;:)]+$", "", h)
    n_paths <- n_paths + 1L
    if (!file.exists(p) && !dir.exists(p)) {
      fail(basename(f), " names ", p, ", which does not exist")
    }
  }
}
say("  ", n_paths, " path references checked")
say("")

# --------------------------------------------- Part B: functions named ----

say("--- Part B: every function named in prose is defined here, in a dependency or in base R ---")
desc0 <- read.dcf("DESCRIPTION")
nmsp <- readLines("NAMESPACE", warn = FALSE)
exported <- gsub('[")]', "", sub(".*export\\(", "", nmsp[grepl("^export\\(", nmsp)]))
exported <- trimws(exported)
exported <- exported[nzchar(exported)]

internal <- unlist(lapply(list.files("R", "\\.R$", full.names = TRUE), function(f) {
  l <- readLines(f, warn = FALSE)
  regmatches(l, regexpr("^\\.?[A-Za-z_][A-Za-z0-9._]*(?=\\s*<-\\s*function)", l, perl = TRUE))
}))
defined <- unique(c(exported, internal))

# The prose legitimately names functions of the packages this one depends on,
# for instance pbivnorm() in NEWS.md, and those namespaces are not attached in a
# plain session, so their exports have to be collected explicitly.
dep_fields <- c("Depends", "Imports", "Suggests")
dep_pkgs <- unlist(lapply(dep_fields, function(k) {
  if (!(k %in% colnames(desc0))) return(character(0))
  v <- trimws(strsplit(desc0[1, k], ",")[[1]])
  sub("\\s*\\(.*", "", v)
}))
dep_pkgs <- setdiff(unique(dep_pkgs[nzchar(dep_pkgs)]), c("R", ""))
dep_unavailable <- character(0)
dep_exports <- unlist(lapply(dep_pkgs, function(pkg) {
  e <- try(getNamespaceExports(pkg), silent = TRUE)
  if (inherits(e, "try-error")) {
    dep_unavailable <<- c(dep_unavailable, pkg)
    character(0)
  } else {
    e
  }
}))

is_known <- function(nm) {
  if (nm %in% defined) return(TRUE)
  if (nm %in% dep_exports) return(TRUE)
  for (env in search()) {
    if (exists(nm, where = env, inherits = FALSE)) return(TRUE)
  }
  FALSE
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
say("  ", n_names, " function references checked against ",
    length(defined), " names defined here and ",
    length(dep_exports), " exported by ", length(dep_pkgs) - length(dep_unavailable),
    " dependencies")
if (length(dep_unavailable)) {
  say("  note  not installed, so their exports could not be checked: ",
      paste(dep_unavailable, collapse = ", "))
}
say("")

# --------------------------------------------------- Part C: counts -------

say("--- Part C: counts stated in prose against the thing they count ---")

WORD_NUM <- c(one = 1, two = 2, three = 3, four = 4, five = 5, six = 6,
              seven = 7, eight = 8, nine = 9, ten = 10, eleven = 11, twelve = 12)

# Returns the integer a file states in the first match of pattern, or NA.
stated_in <- function(path, pattern, words = FALSE) {
  if (!file.exists(path)) return(NA_integer_)
  txt <- read_prose(path)
  m <- regmatches(txt, regexpr(pattern, txt, perl = TRUE))
  if (length(m) == 0L) return(NA_integer_)
  if (words) {
    w <- regmatches(m, regexpr(paste(names(WORD_NUM), collapse = "|"), m))
    if (length(w) == 0L) return(NA_integer_)
    return(as.integer(WORD_NUM[[w]]))
  }
  d <- regmatches(m, regexpr("[0-9]+", m))
  if (length(d) == 0L) return(NA_integer_)
  as.integer(d)
}

check_count <- function(label, stated, actual) {
  if (length(stated) != 1L || is.na(stated)) {
    say("  skip  ", label, " (not stated in the text)")
  } else if (stated != actual) {
    fail(label, ": the text says ", stated, ", the source has ", actual)
  } else {
    say("  ok    ", label, " = ", actual)
  }
}

news <- read_prose("NEWS.md")
sec_111 <- sub("(?s)# twoCoprimary 1\\.1\\.0.*", "", news, perl = TRUE)
bugs_block <- sub("(?s)## Performance.*", "", sec_111, perl = TRUE)
bug_bullets <- length(regmatches(bugs_block, gregexpr("\n\\* ", bugs_block))[[1]])

# "non-base packages" in the article means the Imports field minus the packages
# that ship with R itself, so stats and its siblings do not count.
BASE_PKGS <- c("base", "compiler", "datasets", "graphics", "grDevices", "grid",
               "methods", "parallel", "splines", "stats", "stats4", "tcltk",
               "tools", "utils")
imports <- trimws(strsplit(desc0[1, "Imports"], ",")[[1]])
imports <- sub("\\s*\\(.*", "", imports)
n_imports_nonbase <- sum(!(imports %in% BASE_PKGS))

check_count("defects fixed in 1.1.1, cran-comments vs NEWS bullets",
            stated_in("cran-comments.md", "fixes \\d+ defects"), bug_bullets)
check_count("exported functions, manuscript vs NAMESPACE",
            stated_in(manuscript, "\\d+ exported functions"), length(exported))
check_count("non-base packages, manuscript vs DESCRIPTION Imports less base R",
            stated_in(manuscript, "[a-z]+ non-base packages", words = TRUE),
            n_imports_nonbase)
check_count("vignettes, manuscript vs vignettes/",
            stated_in(manuscript, "\\d+ vignettes"),
            length(list.files("vignettes", "\\.Rmd$")))

# Cell counts come from the reproduction logs rather than from a directory.
for (csv in list.files(out_dir, "^reproduce_.*\\.csv$", full.names = TRUE)) {
  d <- try(utils::read.csv(csv, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || !all(c("table", "difference") %in% names(d))) {
    say("  skip  ", basename(csv), " (unexpected layout)")
    next
  }
  has_verdict <- "verdict" %in% names(d)
  for (tab in unique(d$table)) {
    k <- d$table == tab
    raw <- sum(k & d$difference != 0)
    if (has_verdict) {
      beyond <- sum(k & toupper(d$verdict) != "MATCH")
      say("  info  ", tab, ": ", sum(k), " cells, ", raw,
          " differ from the published value, ", beyond, " beyond tolerance")
    } else {
      say("  info  ", tab, ": ", sum(k), " cells, ", raw, " differ")
    }
  }
}
say("")

# --------------------------------- Part D: sentences shared by two files ----

say("--- Part D: sentences that appear in more than one file ---")
say("  A correction applied to one copy and not the other is how the same error")
say("  survived in the manuscript and a vignette. Confirm the copies still agree.")
say("  Bibliography entries and roxygen argument boilerplate are left out; they")
say("  are shared by construction and carry no claim.")
say("")

sentences_of <- function(f, drop_rd_fields = FALSE) {
  txt <- strip_code(read_prose(f))
  if (drop_rd_fields && grepl("\\.Rd$", f)) {
    txt <- gsub("\\\\arguments\\{(?:[^{}]|\\{[^{}]*\\})*\\}", " ", txt, perl = TRUE)
    txt <- gsub("\\\\value\\{(?:[^{}]|\\{[^{}]*\\})*\\}", " ", txt, perl = TRUE)
  }
  txt <- gsub("\\s+", " ", txt)
  s <- trimws(unlist(strsplit(txt, "(?<=[.]) (?=[A-Z])", perl = TRUE)))
  s[nchar(s) >= 60]
}
norm <- function(s) gsub("[^a-z0-9 ]", "", tolower(s))

# A shared sentence only matters if it makes a claim. A bibliography entry that
# five vignettes cite, and the boilerplate that roxygen repeats across argument
# lists, are shared by construction and would bury the ones worth reading.
is_boilerplate <- function(s) {
  grepl("\\\\item\\{|\\\\itemize\\{", s) |
    grepl("[0-9]+\\([0-9]+\\), *[0-9]+ *[-\u2013] *[0-9]+", s) |
    grepl("\\bdoi:|\\\\doi\\{", s)
}

pieces <- lapply(prose_files, function(f) {
  s <- sentences_of(f, drop_rd_fields = TRUE)
  s <- s[!is_boilerplate(s)]
  if (length(s) == 0L) return(NULL)
  data.frame(file = basename(f), sentence = s, key = norm(s), stringsAsFactors = FALSE)
})
all_s <- do.call(rbind, pieces[!vapply(pieces, is.null, logical(1))])

if (is.null(all_s) || nrow(all_s) == 0L) {
  say("  no sentences collected")
} else {
  dup_keys <- unique(all_s$key[duplicated(all_s$key)])
  shown <- 0L
  for (k in dup_keys) {
    grp <- all_s[all_s$key == k, ]
    if (length(unique(grp$file)) < 2L) next
    shown <- shown + 1L
    say("  * ", paste(unique(grp$file), collapse = " + "))
    say("    ", substr(grp$sentence[1], 1, 200))
  }
  if (shown == 0L) say("  none")
}
say("")

# -------------------- Part E: sentences to be read against their source ----

say("--- Part E: sentences carrying a number or a quantifier ---")
say("  Read each against the table or the log that produces the value. The list")
say("  is generated, so it cannot omit a sentence; the judgement stays human.")
say("  In Rd files the \\arguments and \\value fields are left out, since they")
say("  restate the function signature, which R CMD check already verifies.")
say("")

quant <- paste0("\\b(all|every|always|never|none|only|exactly|identical|unchanged|",
                "smallest|largest|each|reproduces?|agrees?|at most|at least)\\b")
numeric_claim <- "(?<![A-Za-z0-9._])[0-9]+(\\.[0-9]+)?(?![A-Za-z0-9._])"

for (f in prose_files) {
  s <- sentences_of(f, drop_rd_fields = TRUE)
  if (length(s) == 0L) next
  keep <- grepl(quant, s, perl = TRUE, ignore.case = TRUE) | grepl(numeric_claim, s, perl = TRUE)
  s <- s[keep]
  if (length(s) == 0L) next
  say("  ", basename(f), "  (", length(s), " sentences)")
  for (x in s) say("      - ", substr(x, 1, 200))
  say("")
}

# -------------------------------------------- Part F: release readiness ----

say("--- Part F: the three things that are settled last ---")
say("  Each of these is decided at the final knit and is therefore the class of")
say("  item that gets discovered at the last minute. The first is checked here;")
say("  the other two are printed every run so that they cannot be forgotten.")
say("")

rpkg_file <- file.path(repo_root, "04_Manuscript", "_Rpackages.txt")
if (file.exists(manuscript) && file.exists(rpkg_file)) {
  named <- unique(regmatches(
    read_prose(manuscript),
    gregexpr("(?<=CRANpkg\\{)[A-Za-z0-9.]+(?=\\})", read_prose(manuscript), perl = TRUE)
  )[[1]])
  listed <- trimws(readLines(rpkg_file, warn = FALSE))
  listed <- listed[nzchar(listed)]
  missing <- setdiff(named, listed)
  if (length(missing)) {
    fail("_Rpackages.txt omits ", paste(missing, collapse = ", "),
         ", which the article names with \\CRANpkg{}")
  } else {
    say("  all ", length(named), " packages the article names are listed in _Rpackages.txt")
  }
} else {
  say("  note  manuscript or _Rpackages.txt not found, package listing not checked")
}

if (file.exists(manuscript)) {
  m <- read_prose(manuscript)
  d <- regmatches(m, regexpr('(?m)^date:\\s*"[^"]+', m, perl = TRUE))
  d <- sub('^date:\\s*"', "", d)
  v <- regmatches(m, regexpr("(?<=\\(version )[0-9.]+", m, perl = TRUE))
  say("  manuscript date field    : ", if (length(d)) d else "not found",
      "   (set this to the resubmission date at the final knit)")
  say("  version string in the text: ", if (length(v)) v else "not found",
      "   (change to 1.1.1 only after CRAN accepts it)")
}
say("")

say("==============================================================================")
say("Summary")
say("==============================================================================")
say("  automatic failures in Parts A to C : ", n_fail)
say("  Parts D and E are lists to read, not verdicts.")
say("  log written to ", log_path)
close(con)

if (n_fail > 0L) stop("audit_prose_claims: ", n_fail, " failure(s); see ", log_path)
invisible(NULL)
