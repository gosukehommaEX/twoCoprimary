# Sample Size Calculation for Two Co-Primary Endpoints: One Count and One Continuous

Determines the sample size for a two-arm superiority trial with two
co-primary endpoints where one is a count (negative binomial) and one is
continuous (normal), to achieve a specified power at a given
significance level.

## Usage

``` r
ss2MixedCountContinuous(
  r1,
  r2,
  nu,
  t,
  mu1,
  mu2,
  sd,
  r,
  rho1,
  rho2,
  alpha,
  beta
)
```

## Arguments

- r1:

  Mean count rate in group 1 for the count endpoint

- r2:

  Mean count rate in group 2 for the count endpoint

- nu:

  Dispersion parameter for the negative binomial distribution (nu \> 0).
  Smaller values indicate greater overdispersion

- t:

  Follow-up time period

- mu1:

  Mean for group 1 for the continuous endpoint

- mu2:

  Mean for group 2 for the continuous endpoint

- sd:

  Common standard deviation for the continuous endpoint

- r:

  Allocation ratio n1/n2 where n1 is sample size for group 1

- rho1:

  Correlation between count and continuous endpoints in group 1

- rho2:

  Correlation between count and continuous endpoints in group 2

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- beta:

  Type II error rate (typically 0.1 or 0.2). Power = 1 - beta

## Value

A data frame with the following columns:

- r1, r2:

  Count rates

- nu:

  Dispersion parameter

- t:

  Follow-up time

- mu1, mu2:

  Means for continuous endpoint

- sd:

  Standard deviation for continuous endpoint

- r:

  Allocation ratio

- rho1, rho2:

  Correlations

- alpha:

  One-sided significance level

- beta:

  Type II error rate

- n1:

  Required sample size for group 1

- n2:

  Required sample size for group 2

- N:

  Total sample size (n1 + n2)

## Details

This function implements the sample size calculation for mixed
count-continuous co-primary endpoints following the methodology in Homma
and Yoshida (2024).

The sequential search algorithm (Homma and Yoshida 2025, Algorithm 1) is
used:

**Step 1:** Initialize with sample sizes from single endpoint formulas.

**Step 2:** Use sequential search:

- Calculate power at initial sample size

- If power \>= target: decrease n2 until power \< target, then add 1
  back

- If power \< target: increase n2 until power \>= target

**Step 3:** Return final sample sizes.

**Negative Binomial Distribution:** The count endpoint follows a
negative binomial distribution NB(lambda, nu) where:

- lambda = r \* t is the mean count

- nu is the dispersion parameter

- Variance = lambda + lambda^2 / nu

**Correlation:** The correlations rho1 and rho2 must satisfy feasibility
constraints that depend on the parameters. Use
[`corrbound2MixedCountContinuous`](https://gosukehommaEX.github.io/twoCoprimary/reference/corrbound2MixedCountContinuous.md)
to check valid correlation bounds.

## References

Homma, G., & Yoshida, T. (2024). Sample size calculation for count and
continuous multiple co-primary endpoints. *Pharmaceutical Statistics*,
23(3), 372-388.

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(1), 1-19.

## Examples

``` r
# Sample size calculation for count and continuous endpoints
ss2MixedCountContinuous(
  r1 = 1.0,
  r2 = 1.25,
  nu = 0.8,
  t = 1,
  mu1 = -50,
  mu2 = 0,
  sd = 250,
  r = 1,
  rho1 = 0.4,
  rho2 = 0.4,
  alpha = 0.025,
  beta = 0.2
)
#> 
#> Sample size calculation for mixed count and continuous co-primary endpoints
#> 
#>              n1 = 711
#>              n2 = 711
#>               N = 1422
#>              sd = 250
#>            rate = 1, 1.25
#>              nu = 0.8
#>               t = 1
#>              mu = -50, 0
#>             rho = 0.4, 0.4
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#> 

# With different dispersion parameter (more overdispersion)
ss2MixedCountContinuous(
  r1 = 1.0,
  r2 = 1.25,
  nu = 0.5,
  t = 1,
  mu1 = -50,
  mu2 = 0,
  sd = 250,
  r = 1,
  rho1 = 0.4,
  rho2 = 0.4,
  alpha = 0.025,
  beta = 0.2
)
#> 
#> Sample size calculation for mixed count and continuous co-primary endpoints
#> 
#>              n1 = 924
#>              n2 = 924
#>               N = 1848
#>              sd = 250
#>            rate = 1, 1.25
#>              nu = 0.5
#>               t = 1
#>              mu = -50, 0
#>             rho = 0.4, 0.4
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#> 
```
