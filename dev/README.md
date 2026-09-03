# `dev/` — verification scripts

Nothing in this directory is part of the package. `dev` is listed in
`.Rbuildignore`, so none of it is shipped to CRAN, and `dev/out/` is listed in
`.gitignore`, so the results are not committed either. Every script is run from
the package root, against the **installed** package, as

```r
source("dev/<name>.R")
```

and writes its results into `dev/out/` as a log and one or more CSV files, so
that a number quoted anywhere else can be traced back to a file rather than to a
memory of a console session.

The prefix says what kind of script it is.

| prefix | what it does |
|---|---|
| `audit_` | sweeps the whole package against its own contracts |
| `reproduce_` | reproduces the published tables of one source article |
| `verify_` | checks one specific repair against an independent source |
| `measure_` | measures running time |
| `article_` | produces the numbers the R Journal article quotes |
| `release_` | the checks to run immediately before a submission |

There is one script outside that scheme. `audit_prose_claims.R` checks the
English rather than the code: the paths, function names and counts that the
vignettes, `NEWS.md`, `README.md`, `cran-comments.md`, the manuscript, the
response letter and the motivating letter state; it lists the sentences shared
between two files and the sentences carrying a number or a quantifier; and its
last part checks what is settled at the submission, that every package the
article names with `\CRANpkg{}` appears in `_Rpackages.txt` and that the
article's version string and the response letter's mention of 1.1.1 have been
updated together. It exists because the code was swept exhaustively while the
prose never was, which is how a claim corrected in the manuscript stayed wrong
in a vignette.

## Order to run before a release

1. `fuzz_all_functions.R` — does anything error, warn or return nonsense
   anywhere in the parameter space?
2. `audit_all_functions.R` — does every function still honour its contract?
3. the five `reproduce_*.R` scripts — does the package still return what the
   articles it implements printed?
4. the `verify_*.R` scripts for anything that changed in this release
5. `audit_prose_claims.R` — does every statement the prose makes about paths,
   names and counts still hold, and has every shared sentence been corrected in
   every copy? Run it after the `reproduce_*` scripts, whose CSV files it reads.
   Its last part fails on purpose while the response letter names version 1.1.1
   and the article still reports 1.1.0; that failure clears when CRAN accepts
   1.1.1 and the article's version string is updated with it.
6. `release_spelling.R` and `release_checklist.R` last, because Parts D and E of
   the checklist inspect what the earlier scripts left in `dev/out/`

## The scripts

### `audit_all_functions.R`
Exhaustive contract audit, about 470 checks in thirteen parts: documented
arguments against formals and documented return values against the actual
columns; the structure of every shape of object the package returns; the
probability identities that any correct power function satisfies; the
minimality of every sample size function over every test and allocation; the
unified interfaces and `design_table()` against the functions they dispatch to;
the shape, size and ordering of the rejection regions; the bivariate binomial
against its own marginals and moments; the correlation bounds against
independent constructions; argument validation; every `plot()` type against
every shape of object; `print()`; the reproducibility of the Monte Carlo paths;
and a running time budget.
Writes `audit_all_functions.log`, `.csv` and `_timing.csv`.

### `fuzz_all_functions.R`
Exhaustive sweep of every exported function over its whole documented domain,
about twenty-five thousand calls with the grids made dense at the edges. Every
call must return without an unexpected error, without a warning, and with a
value that satisfies the contract of its type; calls listed as invalid must be
refused; and each call runs under an elapsed time limit, so a loop that does
not terminate is reported rather than freezing the session. This is the file
that replaces reading the sources for edge cases, after six defects were found
by reading that `audit_all_functions.R` and the reproduce scripts had both
passed over. Writes `fuzz_all_functions.log` and `.csv`.

### `reproduce_Sozu_et_al_2011.R`
Table 1 of Sozu, Sugimoto and Hamasaki (2011), two continuous endpoints, and
its two single-endpoint columns. Exact reproduction expected.

### `reproduce_Sozu_et_al_2010.R`
Table III of Sozu, Sugimoto and Hamasaki (2010), two binary endpoints under the
four asymptotic methods, with the achieved power at every reproduced size. The
Fisher column of the same table is compared with the exact calculation for the
settings small enough to compute. Tolerance two subjects, the original having
been computed in SAS.
Also writes `reproduce_Sozu_et_al_2010_power.csv`.

### `reproduce_Sozu_et_al_2012.R`
Table 2 (the PREMIER illustration), Table 1 and Supporting Information Table S5
of Sozu, Sugimoto and Hamasaki (2012), mixed continuous and binary endpoints.
The Fisher column of Table S5 is handled by `verify_mixed_fisher_option.R`
instead, because a Monte Carlo power is better checked at the published sample
size than by re-running the search.

### `reproduce_Homma_and_Yoshida_2024.R`
Tables 1 and 2 of Homma and Yoshida (2024), mixed count and continuous
endpoints, both allocations, together with the correlation between the two test
statistics, which is recovered by inverting the bivariate normal probability.
Exact reproduction expected.

### `reproduce_Homma_and_Yoshida_2025.R`
Tables 3 and 4 of Homma and Yoshida (2025), two binary endpoints by exact
methods. Exact reproduction expected, with one known exception: the entry of
Table 4 at Z-pool, alpha = 0.05, r = 2 and rho = 0.3 was produced by a
rejection region that mishandled tied outcomes, and the corrected value is
checked separately. This is the slowest script; restrict `ALPHAS` to 0.025 for
a shorter run.

### `verify_two_continuous_family.R`
Eight-part verification of the two continuous endpoint family, including the
Wishart construction for unknown variance and a simulation check of the power.

### `verify_exact_binary_ties.R`
The corrected tie handling in `rr1Binary()`, checked cell by cell against the
rejection regions of the **Exact** package at three significance levels.

### `verify_tie_impact_on_manuscript.R`
Whether the tie correction moves any sample size the R Journal article reports,
and which convention agrees with independent software.

### `verify_exact_binary_n_grid.R`
The `n_grid` argument: what a finer nuisance parameter grid does to the
p-values, the rejection regions, the sample sizes and the running time.

### `verify_bibinom_cpp_kernel.R`
The compiled `dbibinom_g` kernel against the R reference implementation
retained in the package, over agreement, identities, sample sizes and timing.

### `verify_mixed_fisher_option.R`
The repaired Fisher option for mixed continuous and binary endpoints, including
Table S5 of Sozu et al. (2012), where the power at the published sample size is
compared with the published empirical power.

### `verify_1_1_1_fixes.R`
Every repair made in version 1.1.1, plus a recomputation of every validation
grid the R Journal article reports, so that the article's numbers can be diffed
against the current PDF.

### `measure_countcont_speed.R`
Where the time goes in `ss2MixedCountContinuous()`. Section 1 measures the
current total; sections 3 and 4 deliberately call the exported power function,
which still validates the correlation bounds on every call, in order to show
what the search would cost without the bounds being hoisted out of it.

### `measure_exact_binary_speed.R`
Accuracy and timing of the matrix product form of the exact co-primary power.

### `article_usage_examples.R`
The `design_table()` and `plot()` output that the Usage examples section of the
R Journal article shows.

### `article_correlation_examples.R`
The correlation estimation examples added to the article at the request of
Reviewer 4.

### `article_HY2025_tie_disclosure.R`
Which published results of Homma and Yoshida (2025) the tie correction changes,
computed three ways: as printed, from the fixed package, and from a
reimplementation of the pre-fix convention.

### `release_spelling.R`
Triage of `spelling::spell_check_package()` into words already in
`inst/WORDLIST`, words to add, and words to correct, so that none has to be
judged by eye.

### `release_checklist.R`
The final consistency check before a CRAN submission: that the attached build
is newer than every source file, the worked examples printed in `README.md`, the
citations in the help pages against the reference list, the version against the
`NEWS.md` heading, this README against the contents of `dev/`, and every log in
`dev/out/` against the installed version. The first part stops the script,
because a failed or forgotten reinstall makes every script in `dev/` report on
code that is no longer on disk and the version number cannot detect it, a
rebuild during a release carrying the same version. The last two parts exist so
that neither the index nor the results can go quietly out of date.

### `reproduce_helpers.R`
Shared logging and comparison machinery for the five `reproduce_*.R` scripts.
Not run on its own.

## `dev/out/`

Everything in `dev/out/` is regenerable and none of it is committed. Delete the
whole directory and rerun whenever the package is rebuilt: a result file written
by an earlier version is the one way a number quoted in the article or the
response letter can be silently stale. Part E of `release_checklist.R` checks
that every log there names the installed version, and no script reads another
script's output, so deleting the directory loses nothing.
