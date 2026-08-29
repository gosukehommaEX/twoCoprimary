## R CMD check results

0 errors | 0 warnings | 0 notes

## Submission

This is a minor release. It fixes three defects, improves the performance of the
exact methods for binary endpoints, and adds one argument.

* `rr1Binary()` now assigns the same p-value to outcomes that share the same
  value of the ordering statistic in the two exact unconditional tests. The
  resulting rejection regions agree with those implied by the `Exact` package.

* `power2Continuous()` with `known_var = FALSE` now draws the Wishart matrix
  from the correlation matrix of the standardized endpoints, so power is
  invariant to a common rescaling of the effects and the standard deviations.

* `power2MixedContinuousBinary()` with `Test = "Fisher"` now returns the
  simulated result and uses a group specific threshold for the latent binary
  variable. It reproduces Table S5 of the Supporting Information of
  Sozu et al. (2012).

* The package now contains compiled code. `Rcpp` was added to `Imports` and
  `LinkingTo`, and `src/twoCoprimary.cpp` holds a single function that evaluates
  the conditional probability of the bivariate binomial distribution.
  `NeedsCompilation` therefore changes from no to yes.

* `n_grid` was added to the functions for the exact binary methods, exposing the
  number of nuisance parameter grid points that was previously fixed at 100
  internally. The default reproduces the results of version 1.0.0.

## Test environments

* local Windows 11 install, R 4.6.0
* win-builder: R-devel (2026-08-27 r90452), R-release (4.6.1)
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  windows-latest (R-release), macos-latest (R-release)

## Downstream dependencies

There are currently no downstream dependencies for this package.
