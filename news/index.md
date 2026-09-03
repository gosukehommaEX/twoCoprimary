# Changelog

## twoCoprimary 1.1.1

### Bug fixes

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) no longer
  fails on a sample size object when `type = "power_curve"` is
  requested. The reference-line block selected its branch on the
  presence of an `n2` column, which both object shapes have, so a sample
  size object entered the branch written for power objects and read a
  `powerCoprimary` column that is not there. The branches now select on
  `powerCoprimary`, which distinguishes the two shapes, and the branch
  already written for sample size objects is reachable for the first
  time.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) no longer
  fails on a power object when `type = "sample_size_rho"` is requested.
  A power object carries neither the allocation ratio `r`, the type II
  error rate `beta` nor the total `N`. These are now derived from the
  columns it does carry: the realized allocation `n1 / n2`, the achieved
  power as the target, and the sum of the two group sizes.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) no longer
  fails on results from
  [`ss2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2BinaryExact.md)
  and
  [`power2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryExact.md).
  The plotting helpers passed the object’s `Test` value into the
  approximate binary functions, which do not accept the exact test
  names, so the default call failed for every exact test. The helpers
  now route exact test names to
  [`power2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryExact.md)
  and
  [`ss2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2BinaryExact.md).

- [`power2BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryApprox.md)
  now validates `Test`. An unrecognized value previously fell through
  every branch and the function failed with an internal error about a
  missing object instead of a message naming the valid values.

- [`ss2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2MixedCountContinuous.md)
  now requires `r1 < r2` and `mu1 < mu2`. Treatment benefit on both
  endpoints is a lower value, so with the effects reversed the power
  decreases as the sample size grows, no sample size attains the target,
  and the sequential search did not terminate.

- [`power2Continuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2Continuous.md)
  with `known_var = FALSE` now returns zero power when `n1 + n2 - 2 < 1`
  instead of failing inside
  [`rWishart()`](https://rdrr.io/r/stats/rWishart.html). The variance
  cannot be estimated from fewer than three patients, so the t-test does
  not exist. `ss2Continuous(known_var = FALSE)` previously failed for
  standardized effect sizes above about 4 for this reason.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a result
  from
  [`ss1Continuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1Continuous.md),
  [`ss1Count()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1Count.md)
  or
  [`ss1BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1BinaryApprox.md)
  now explains that the co-primary plot methods need a two-endpoint
  result, rather than reporting that the endpoint type could not be
  determined.

- [`ss1BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1BinaryApprox.md)
  returns the smallest sample size at which
  [`power2BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryApprox.md)
  attains the target power, for all four asymptotic methods. The
  closed-form expression is kept only as the starting value of the
  sequential search, so the size returned and the power reported by the
  package can no longer disagree. Four faults are removed by the change.
  The arcsine multiplier was `(1 + kappa) / kappa`, the form written for
  `kappa = n1 / n2`, while the function defines `kappa = n2 / n1`, so at
  `r = 2` the function returned roughly twice the required sample size.
  The `"ASc"` continuity correction moved the two groups apart rather
  than toward each other, so it reduced the sample size below the
  uncorrected value. The corrected methods `"ANc"` and `"ASc"` solved
  their fixed point by an iteration that stopped once two successive
  values differed by one, which can return a size that does not satisfy
  its own correction and can cycle without terminating, for instance at
  `p1` = 0.95, `p2` = 0.05 and `r` = 5. And the `"ASc"` variance was
  `1 / (4 n)` rather than the corrected-proportion form of Sozu et
  al. (2012), Supporting Information Table 1, which the power function
  uses; the two disagreed for about half of the designs of that table.
  The twelve published `"ANc"` and `"ASc"` single-endpoint sizes of
  Table S5 are now reproduced exactly.

- `ss1BinaryApprox(Test = "Fisher")` steps back down after its upward
  search, so the size returned is the smallest one attaining the target
  power rather than the first one reached from the starting value.

- `power2BinaryApprox(Test = "ASc")` returns zero power where the
  continuity correction carries a probability out of `(0, 1)`, which
  happens at the smallest sample sizes. The corrected variance then
  vanished or was negative, and the reported power was one half or
  `NaN`.

- [`power2BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryApprox.md)
  evaluates the co-primary power at the Frechet bounds of its two
  marginals when the correlation between the test statistics reaches
  plus or minus one, rather than calling the bivariate normal
  distribution function at a correlation it does not accept. The power
  is the smaller marginal at plus one and the Bonferroni bound at minus
  one.

- `power2MixedContinuousBinary(Test = "ASc")` gains the same guard on
  the corrected proportions, and its `"Fisher"` branch returns zero
  power when `n1 + n2 - 2 < 1` rather than dividing by zero degrees of
  freedom.

- [`power2Continuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2Continuous.md)
  validates its arguments. It previously accepted a negative standard
  deviation, a group size of zero and a non-integer group size, and
  returned a number for each.

- [`dbibinom()`](https://gosukehommaex.github.io/twoCoprimary/reference/dbibinom.md)
  returns the comonotone limit at the upper Prentice bound when the two
  marginal probabilities are equal. The dependence parameter of the
  bivariate binomial distribution diverges there and the returned masses
  were `NaN`.

- [`ss1Count()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1Count.md)
  requires `r1 < r2`. With the two rates reversed the power decreases as
  the sample size grows, so no sample size attains the target.

- `rr1Binary(Test = "Z-pool")` and `rr1Binary(Test = "Boschloo")` build
  the rejection region when a group has a single subject. Only one
  distinct value of the ordering statistic is then ordered, and the grid
  of null probabilities was simplified from a matrix to a vector, so the
  function failed with an error about an incorrect number of dimensions.

- `plot(type = "effect_contour")` grids standardized effect sizes, which
  is what its axis labels state. It previously used a grid of mean
  differences and labeled them as standardized, so the contours were
  misplaced for any standard deviation other than one. The marked point
  is now on the same scale.

- `plot(type = "power_curve")` no longer inverts its default window for
  a design of fewer than twenty per group. The lower end of the window
  was a fixed floor of ten, above the upper end for a small design.

- [`design_table()`](https://gosukehommaex.github.io/twoCoprimary/reference/design_table.md)
  now forwards `nMC` to the mixed continuous-binary functions. The
  argument was documented and honored for continuous endpoints but
  silently ignored on that path, which ran at the default of 10000 in
  the function it calls.

- The sequential search no longer adds one back to `n2` when the
  downward search stops at `n2 = 1` with the power still at or above
  target.

- [`design_table()`](https://gosukehommaex.github.io/twoCoprimary/reference/design_table.md)
  uses [`seq_len()`](https://rdrr.io/r/base/seq.html) over the parameter
  grid, so a zero-row grid no longer iterates.

- [`power2BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryApprox.md)
  now validates `n1`, `n2`, the four probabilities and `alpha`, as
  [`power2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryExact.md)
  already did. A group size of zero or a probability outside the unit
  interval previously returned a number rather than a message.

- [`corrbound2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/corrbound2MixedCountContinuous.md)
  no longer depends on `mu` or `sd`. The bounds are a correlation, which
  is invariant to the location and the scale of the continuous endpoint,
  but they were computed by a quadrature over the whole real line in the
  untransformed variable. That quadrature returns zero when the mean
  lies far from the origin relative to the standard deviation: at `mu`
  of -50 and `sd` of 0.5 both bounds collapsed to zero, and
  [`power2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2MixedCountContinuous.md)
  and
  [`ss2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2MixedCountContinuous.md)
  then refused every correlation, zero included. The integrals are taken
  in the standardized variable instead.

- `power2Continuous(known_var = FALSE)` returns a power for a design of
  three patients in total.
  [`rWishart()`](https://rdrr.io/r/stats/rWishart.html) requires as many
  degrees of freedom as the dimension of its scale matrix, so
  `n1 + n2 - 2 = 1` failed inside it. The Wishart matrix is the outer
  product of a single normal draw there, and only its diagonal enters
  the calculation, so that case is drawn directly.

- The bivariate normal distribution function is no longer evaluated at
  arguments it cannot handle. A standardized effect of several hundred,
  which a small standard deviation produces, sends an argument far
  enough into a tail that `pbivnorm()` returns `NaN`, and the reported
  co-primary power was not a number. The arguments are held at eight
  standard deviations, beyond which the normal distribution function is
  zero or one to within 1e-15.

- The co-primary power is held inside the interval that its two marginal
  powers allow, in every power function that evaluates a bivariate
  normal distribution function. Rounding could otherwise return a value
  a few times 1e-19 below zero, and in
  `power2Continuous(known_var = FALSE)`, where the marginals are exact
  and the joint probability is simulated, a small `nMC` could put the
  simulated value above the smaller marginal.

- [`ss1Continuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1Continuous.md),
  [`ss1Count()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1Count.md)
  and
  [`ss1BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1BinaryApprox.md)
  return at least one subject per group. Their closed forms give zero
  when the target power does not exceed the size of the test, at `alpha`
  of 0.1 and `beta` of 0.9 for instance, and
  [`ss1BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1BinaryApprox.md)
  then evaluated its power at a group of no patients.

The two-endpoint sample size functions are unaffected by the
[`ss1BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1BinaryApprox.md)
corrections, because they use it only as the starting value for a
sequential search that converges to the minimum sample size from any
starting point.

### Performance

- [`ss2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2MixedCountContinuous.md)
  evaluates the correlation bounds once instead of twice at every step
  of its sequential search. The bounds are a quadrature over the support
  of the negative binomial distribution and depend only on the
  parameters held fixed by the search. A search at the default settings
  takes about 0.2 seconds rather than about 2.7.

### Documentation

- [`twoCoprimary2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/twoCoprimary2BinaryExact.md)
  documented `Test = "Z-pooled"`, which the code rejects, and omitted
  `"Fisher-midP"`, which it accepts.
- [`design_table()`](https://gosukehommaex.github.io/twoCoprimary/reference/design_table.md)
  now documents the five exact test methods it dispatches on.
- [`power2MixedContinuousBinary()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2MixedContinuousBinary.md)
  documents the `nMC` column it returns.
- [`power2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2MixedCountContinuous.md)
  writes the test statistics and the correlation between them in the
  indexing used by the package, group 1 for treatment and group 2 for
  control, rather than the indexing of the source article, and states
  that the event rates `r1` and `r2` are distinct from the allocation
  ratio `r`. The variance component `V_a` is now defined where it is
  used.
- [`ss2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2MixedCountContinuous.md)
  writes the mean count as `lambda_j = r_j * t`.
- [`ss1BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1BinaryApprox.md)
  documents the requirement that `p1` exceed `p2`, and no longer
  describes a binomial calculation as hypergeometric. Its arcsine
  formula matches the implementation.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) documents the
  endpoint-specific default for `rho_range` and the cost of `n_points`
  for exact and Monte Carlo based objects.
- [`ss1Count()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss1Count.md)
  documents that `r1` must be less than `r2`.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) documents
  that the `"effect_contour"` axes are standardized effect sizes.
- The vignettes carry the notation of the articles they follow. The
  variance of the log rate ratio in `mixed-count-continuous` is written
  with the event rates rather than the mean counts, which is the form of
  equation 8 of Homma and Yoshida (2024) and the form the code uses; its
  Case B design parameters name the treatment and control rates in the
  order the code passes them. The non-centrality parameter in
  `two-continuous-endpoints` is written with the standardized effect
  size. `mixed-continuous-binary` states the allocation ratio `r` of the
  package alongside the `kappa` of Sozu et al. (2012), uses the
  upper-tail quantile in the power formula and its two critical values,
  closes the parenthesis in the latent normal distribution, and no
  longer writes the standardized statistic of the continuous endpoint in
  terms of proportions. `overview` separates the patient-level
  correlation the user supplies from the correlation between test
  statistics that enters the power formula. `two-binary-endpoints-exact`
  writes the null and alternative hypotheses with the subscripts used by
  the other vignettes, and denotes the multinomial cell counts by `N`
  rather than by `Z`, which is the test statistic elsewhere.
- The vignettes state what their reproductions of the published tables
  actually show. Two footnotes attributed the differences to the
  bivariate normal distribution function of SAS differing from that of
  R, which is not the reason. In Sozu et al. (2010) Table III the twelve
  values of the block printed at a first probability of 0.87 are the
  only ones that do not reproduce, and recomputing that block at 0.868
  reproduces all twelve exactly. In Table 5 of the Supporting
  Information of Sozu et al. (2012) four of the forty-eight values
  differ by one subject, and for each of the four a standardized effect
  size inside the rounding window of the printed three decimal places
  returns the published size. Three further statements are corrected:
  the number of ties at forty subjects per group in
  `two-binary-endpoints-exact`, the grid sizes at which the rejection
  regions were compared in the same vignette, and how many of the
  published Fisher sample sizes the test suite checks.

## twoCoprimary 1.1.0

CRAN release: 2026-08-29

### Bug fixes

- [`rr1Binary()`](https://gosukehommaex.github.io/twoCoprimary/reference/rr1Binary.md)
  now gives every outcome in a tie group the same p-value for the two
  exact unconditional tests, `"Z-pool"` and `"Boschloo"`. The tail event
  of such a test is the set of outcomes at least as extreme as the
  observed one, so outcomes sharing the same value of the ordering
  statistic must share a tail probability. A plain cumulative sum gave
  the earlier members of a group a smaller p-value, and which member
  came first depended on the order that
  [`order()`](https://rdrr.io/r/base/order.html) happened to return. The
  rejection region could therefore be larger than the standard test’s,
  and the result was not reproducible. The corrected rejection regions
  agree with those implied by the `Exact` package at `alpha` = 0.01,
  0.025 and 0.05 for both tests.

  This changes one entry of Table 4 of Homma and Yoshida (2025): for the
  Z-pooled test with `alpha` = 0.05, `r` = 2 and `rho` = 0.3 the
  required total sample size is 147 rather than the published 144. The
  remaining 95 entries of that table, and every design behind its
  Figures 1 to 3 but one, are unchanged.

- [`power2Continuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2Continuous.md)
  with `known_var = FALSE` now draws the Wishart matrix from the
  correlation matrix of the standardized endpoints rather than from the
  variance-covariance matrix. The test statistic is already divided by
  `sd1` and `sd2`, so a factor of the standard deviation was left in the
  critical value and canceled only when `sd1 = sd2 = 1`. Power is now
  invariant to a common rescaling of the effects and the standard
  deviations, as it must be, and agrees with a direct simulation of the
  trial.

- [`power2MixedContinuousBinary()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2MixedContinuousBinary.md)
  with `Test = "Fisher"` now returns the simulated result. The block
  computing the asymptotic result was not inside an `else`, so it ran
  for every value of `Test` and overwrote what the Fisher branch had
  computed; a request for `"Fisher"` returned the `"ASc"` value. The
  latent binary variable is also now centered in both groups and
  dichotomized at a group specific threshold, following Sozu et
  al. (2012). Both groups had previously been dichotomized at a single
  threshold, which pinned the control response probability at 0.5
  whatever `p2` was and reversed the direction of the effect. The
  corrected implementation reproduces all 24 Fisher entries of Table S5
  in the Supporting Information of Sozu et al. (2012) to within Monte
  Carlo error.

### Performance

- The co-primary power in
  [`power2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryExact.md)
  is now evaluated as two matrix products of order `n + 1` instead of
  forming two matrices of order `K`, the number of outcomes in the
  rejection region. Since `K` grows like `n ^ 2` the cost fell from
  order `n ^ 4` to order `n ^ 3` in both time and memory. At
  `n1 = n2 = 200` the previous expression required about 4 GB and could
  not be evaluated on a typical machine.

- The conditional probability of the bivariate binomial distribution
  used by
  [`dbibinom()`](https://gosukehommaex.github.io/twoCoprimary/reference/dbibinom.md)
  is now computed in C++. The powers and binomial coefficients that the
  sum needs are tabulated once, which reduces the number of calls to
  `pow` and `choose` from order `N ^ 3` to order `N ^ 2`. Building one
  probability mass matrix at `N = 200` takes about 0.07 seconds rather
  than 2.1 seconds.

- Together these make a sample size search for two binary endpoints with
  exact methods around four times faster at moderate sample sizes, and
  possible at all at large ones.

### New features

- [`rr1Binary()`](https://gosukehommaex.github.io/twoCoprimary/reference/rr1Binary.md),
  [`power2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryExact.md),
  [`ss2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2BinaryExact.md),
  [`twoCoprimary2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/twoCoprimary2BinaryExact.md)
  and
  [`design_table()`](https://gosukehommaex.github.io/twoCoprimary/reference/design_table.md)
  gain an `n_grid` argument giving the number of points at which the
  null tail probability is maximized over the nuisance parameter in the
  two exact unconditional tests. The value was previously fixed at 100
  internally, which remains the default, so existing results are
  reproduced exactly. Values below 10 are rejected, since a coarse grid
  can miss the maximum and return an anti-conservative p-value.

### Documentation

- [`power2Continuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2Continuous.md),
  [`power2MixedContinuousBinary()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2MixedContinuousBinary.md)
  and
  [`ss2MixedContinuousBinary()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2MixedContinuousBinary.md)
  now state that their Monte Carlo branches draw random numbers and that
  a seed should be set for reproducible results.

- [`rr1Binary()`](https://gosukehommaex.github.io/twoCoprimary/reference/rr1Binary.md)
  documents how tied outcomes are handled.

### Internal

- `Rcpp` is a new dependency, in `Imports` and `LinkingTo`.

- Test coverage was extended with invariance checks that apply across
  all five endpoint type combinations: the co-primary power factorizes
  into the product of the marginal powers at zero correlation, power
  depends on a continuous endpoint only through its standardized effect,
  and the returned sample size is the smallest one reaching the target.
  The exact unconditional tests are also checked against the `Exact`
  package where it is installed.

- A spelling check was added (`tests/spelling.R` with `inst/WORDLIST`),
  so the British and American spelling inconsistencies corrected in this
  version cannot reappear unnoticed.

## twoCoprimary 1.0.0

CRAN release: 2025-11-21

### Initial Release

This is the first release of twoCoprimary, providing comprehensive tools
for sample size and power calculation in clinical trials with two
co-primary endpoints.

#### Features

- **Two Continuous Endpoints**
  - [`ss2Continuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2Continuous.md),
    [`power2Continuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2Continuous.md),
    [`twoCoprimary2Continuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/twoCoprimary2Continuous.md)
  - Based on Sozu et al. (2011)
  - Supports known and unknown variance cases
- **Two Binary Endpoints (Asymptotic Methods)**
  - [`ss2BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2BinaryApprox.md),
    [`power2BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryApprox.md),
    [`twoCoprimary2BinaryApprox()`](https://gosukehommaex.github.io/twoCoprimary/reference/twoCoprimary2BinaryApprox.md)
  - Based on Sozu et al. (2010)
  - Four test methods: AN, ANc, AS, ASc
- **Two Binary Endpoints (Exact Methods)**
  - [`ss2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2BinaryExact.md),
    [`power2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryExact.md),
    [`twoCoprimary2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/twoCoprimary2BinaryExact.md)
  - Based on Homma and Yoshida (2025)
  - Five exact tests: Chisq, Fisher, Fisher-midP, Z-pool, Boschloo
- **Mixed Continuous and Binary Endpoints**
  - [`ss2MixedContinuousBinary()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2MixedContinuousBinary.md),
    [`power2MixedContinuousBinary()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2MixedContinuousBinary.md),
    [`twoCoprimary2MixedContinuousBinary()`](https://gosukehommaex.github.io/twoCoprimary/reference/twoCoprimary2MixedContinuousBinary.md)
  - Based on Sozu et al. (2012)
  - Supports biserial correlation structure
- **Mixed Count and Continuous Endpoints**
  - [`ss2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2MixedCountContinuous.md),
    [`power2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2MixedCountContinuous.md),
    [`twoCoprimary2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/twoCoprimary2MixedCountContinuous.md)
  - Based on Homma and Yoshida (2024)
  - Handles overdispersed count data with negative binomial distribution

#### Utility Functions

- [`corrbound2Binary()`](https://gosukehommaex.github.io/twoCoprimary/reference/corrbound2Binary.md) -
  Calculate valid correlation bounds for binary endpoints
- [`corrbound2MixedCountContinuous()`](https://gosukehommaex.github.io/twoCoprimary/reference/corrbound2MixedCountContinuous.md) -
  Calculate valid correlation bounds for count and continuous endpoints
- [`design_table()`](https://gosukehommaex.github.io/twoCoprimary/reference/design_table.md) -
  Create comprehensive design comparison tables
- [`plot.twoCoprimary()`](https://gosukehommaex.github.io/twoCoprimary/reference/plot.twoCoprimary.md) -
  Visualize sample size vs correlation relationships

#### Documentation

- Six comprehensive vignettes covering all methodologies
- Complete function documentation with examples
- Validation against published results

#### Testing

- Comprehensive test suite with testthat
- Tests for all major functions
- Validation against published tables and results
