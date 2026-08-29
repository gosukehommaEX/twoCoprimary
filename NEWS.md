# twoCoprimary 1.1.0

## Bug fixes

* `rr1Binary()` now gives every outcome in a tie group the same p-value for the
  two exact unconditional tests, `"Z-pool"` and `"Boschloo"`. The tail event of
  such a test is the set of outcomes at least as extreme as the observed one, so
  outcomes sharing the same value of the ordering statistic must share a tail
  probability. A plain cumulative sum gave the earlier members of a group a
  smaller p-value, and which member came first depended on the order that
  `order()` happened to return. The rejection region could therefore be larger
  than the standard test's, and the result was not reproducible. The corrected
  rejection regions agree with those implied by the `Exact` package at
  `alpha` = 0.01, 0.025 and 0.05 for both tests.

  This changes one entry of Table 4 of Homma and Yoshida (2025): for the
  Z-pooled test with `alpha` = 0.05, `r` = 2 and `rho` = 0.3 the required total
  sample size is 147 rather than the published 144. The remaining 47 entries of
  that table, and every design behind its Figures 1 to 3 but one, are unchanged.

* `power2Continuous()` with `known_var = FALSE` now draws the Wishart matrix
  from the correlation matrix of the standardized endpoints rather than from the
  variance-covariance matrix. The test statistic is already divided by `sd1` and
  `sd2`, so a factor of the standard deviation was left in the critical value
  and canceled only when `sd1 = sd2 = 1`. Power is now invariant to a common
  rescaling of the effects and the standard deviations, as it must be, and
  agrees with a direct simulation of the trial.

* `power2MixedContinuousBinary()` with `Test = "Fisher"` now returns the
  simulated result. The block computing the asymptotic result was not inside an
  `else`, so it ran for every value of `Test` and overwrote what the Fisher
  branch had computed; a request for `"Fisher"` returned the `"ASc"` value. The
  latent binary variable is also now centered in both groups and dichotomized at
  a group specific threshold, following Sozu et al. (2012). Both groups had
  previously been dichotomized at a single threshold, which pinned the control
  response probability at 0.5 whatever `p2` was and reversed the direction of
  the effect. The corrected implementation reproduces all 24 Fisher entries of
  Table S5 in the Supporting Information of Sozu et al. (2012) to within Monte
  Carlo error.

## Performance

* The co-primary power in `power2BinaryExact()` is now evaluated as two matrix
  products of order `n + 1` instead of forming two matrices of order `K`, the
  number of outcomes in the rejection region. Since `K` grows like `n ^ 2` the
  cost fell from order `n ^ 4` to order `n ^ 3` in both time and memory. At
  `n1 = n2 = 200` the previous expression required about 4 GB and could not be
  evaluated on a typical machine.

* The conditional probability of the bivariate binomial distribution used by
  `dbibinom()` is now computed in C++. The powers and binomial coefficients that
  the sum needs are tabulated once, which reduces the number of calls to `pow`
  and `choose` from order `N ^ 3` to order `N ^ 2`. Building one probability
  mass matrix at `N = 200` takes about 0.07 seconds rather than 2.1 seconds.

* Together these make a sample size search for two binary endpoints with exact
  methods around four times faster at moderate sample sizes, and possible at all
  at large ones.

## New features

* `rr1Binary()`, `power2BinaryExact()`, `ss2BinaryExact()`,
  `twoCoprimary2BinaryExact()` and `design_table()` gain an `n_grid` argument
  giving the number of points at which the null tail probability is maximized
  over the nuisance parameter in the two exact unconditional tests. The value
  was previously fixed at 100 internally, which remains the default, so existing
  results are reproduced exactly. Values below 10 are rejected, since a coarse
  grid can miss the maximum and return an anti-conservative p-value.

## Documentation

* `power2Continuous()`, `power2MixedContinuousBinary()` and
  `ss2MixedContinuousBinary()` now state that their Monte Carlo branches draw
  random numbers and that a seed should be set for reproducible results.

* `rr1Binary()` documents how tied outcomes are handled.

## Internal

* `Rcpp` is a new dependency, in `Imports` and `LinkingTo`.

* Test coverage was extended with invariance checks that apply across all five
  endpoint type combinations: the co-primary power factorizes into the product
  of the marginal powers at zero correlation, power depends on a continuous
  endpoint only through its standardized effect, and the returned sample size is
  the smallest one reaching the target. The exact unconditional tests are also
  checked against the `Exact` package where it is installed.

* A spelling check was added (`tests/spelling.R` with `inst/WORDLIST`), so the
  British and American spelling inconsistencies corrected in this version cannot
  reappear unnoticed.

# twoCoprimary 1.0.0

## Initial Release

This is the first release of twoCoprimary, providing comprehensive tools for sample size and power calculation in clinical trials with two co-primary endpoints.

### Features

* **Two Continuous Endpoints**
  - `ss2Continuous()`, `power2Continuous()`, `twoCoprimary2Continuous()`
  - Based on Sozu et al. (2011)
  - Supports known and unknown variance cases

* **Two Binary Endpoints (Asymptotic Methods)**
  - `ss2BinaryApprox()`, `power2BinaryApprox()`, `twoCoprimary2BinaryApprox()`
  - Based on Sozu et al. (2010)
  - Four test methods: AN, ANc, AS, ASc

* **Two Binary Endpoints (Exact Methods)**
  - `ss2BinaryExact()`, `power2BinaryExact()`, `twoCoprimary2BinaryExact()`
  - Based on Homma and Yoshida (2025)
  - Five exact tests: Chisq, Fisher, Fisher-midP, Z-pool, Boschloo

* **Mixed Continuous and Binary Endpoints**
  - `ss2MixedContinuousBinary()`, `power2MixedContinuousBinary()`, `twoCoprimary2MixedContinuousBinary()`
  - Based on Sozu et al. (2012)
  - Supports biserial correlation structure

* **Mixed Count and Continuous Endpoints**
  - `ss2MixedCountContinuous()`, `power2MixedCountContinuous()`, `twoCoprimary2MixedCountContinuous()`
  - Based on Homma and Yoshida (2024)
  - Handles overdispersed count data with negative binomial distribution

### Utility Functions

* `corrbound2Binary()` - Calculate valid correlation bounds for binary endpoints
* `corrbound2MixedCountContinuous()` - Calculate valid correlation bounds for count and continuous endpoints
* `design_table()` - Create comprehensive design comparison tables
* `plot.twoCoprimary()` - Visualize sample size vs correlation relationships

### Documentation

* Six comprehensive vignettes covering all methodologies
* Complete function documentation with examples
* Validation against published results

### Testing

* Comprehensive test suite with testthat
* Tests for all major functions
* Validation against published tables and results
