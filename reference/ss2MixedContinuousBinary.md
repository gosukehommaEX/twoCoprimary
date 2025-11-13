# Sample Size Calculation for Two Co-Primary Endpoints: One Continuous and One Binary

Determines the sample size for a two-arm superiority trial with two
co-primary endpoints where one is continuous and one is binary, to
achieve a specified power at a given significance level.

## Usage

``` r
ss2MixedContinuousBinary(
  delta,
  sd,
  p1,
  p2,
  rho,
  r,
  alpha,
  beta,
  Test,
  nMC = 10000
)
```

## Arguments

- delta:

  Mean difference for the continuous endpoint (group 1 - group 2)

- sd:

  Common standard deviation for the continuous endpoint

- p1:

  Probability of response in group 1 for the binary endpoint (0 \< p1 \<
  1)

- p2:

  Probability of response in group 2 for the binary endpoint (0 \< p2 \<
  1)

- rho:

  Biserial correlation between the latent continuous variable underlying
  the binary endpoint and the observed continuous endpoint

- r:

  Allocation ratio n1/n2 where n1 is sample size for group 1

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- beta:

  Type II error rate (typically 0.1 or 0.2). Power = 1 - beta

- Test:

  Statistical testing method for the binary endpoint. One of:

  - `"AN"`: Asymptotic normal method without continuity correction

  - `"ANc"`: Asymptotic normal method with continuity correction

  - `"AS"`: Arcsine method without continuity correction

  - `"ASc"`: Arcsine method with continuity correction

  - `"Fisher"`: Fisher's exact test (uses sequential search)

- nMC:

  Number of Monte Carlo replications when Test = "Fisher" (default:
  10000)

## Value

A data frame with the following columns:

- delta:

  Mean difference for continuous endpoint

- sd:

  Standard deviation for continuous endpoint

- p1:

  Response probability in group 1 for binary endpoint

- p2:

  Response probability in group 2 for binary endpoint

- rho:

  Biserial correlation

- r:

  Allocation ratio

- alpha:

  One-sided significance level

- beta:

  Type II error rate

- Test:

  Testing method used for binary endpoint

- nMC:

  Number of Monte Carlo replications (NA if Test != "Fisher")

- n1:

  Required sample size for group 1

- n2:

  Required sample size for group 2

- N:

  Total sample size (n1 + n2)

## Details

This function implements the sample size calculation for mixed
continuous-binary co-primary endpoints following the methodology in Sozu
et al. (2012).

The sequential search algorithm (Homma and Yoshida 2025, Algorithm 1) is
used for all testing methods:

**Step 1:** Initialize with sample sizes from single endpoint formulas.

**Step 2:** Use sequential search:

- Calculate power at initial sample size

- If power \>= target: decrease n2 until power \< target, then add 1
  back

- If power \< target: increase n2 until power \>= target

**Step 3:** Return final sample sizes.

**Biserial Correlation:** The biserial correlation rho represents the
correlation between the latent continuous variable underlying the binary
endpoint and the observed continuous endpoint. This is not the same as
the point-biserial correlation observed in the data.

## References

Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination
in clinical trials with multiple co-primary endpoints including mixed
continuous and binary variables. *Biometrical Journal*, 54(5), 716-729.

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(1), 1-19.

## Examples

``` r
# Sample size calculation using asymptotic normal method
ss2MixedContinuousBinary(
  delta = 0.5,
  sd = 1,
  p1 = 0.6,
  p2 = 0.4,
  rho = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.1,
  Test = 'AN'
)
#> 
#> Sample size calculation for mixed continuous and binary co-primary endpoints
#> 
#>              n1 = 135
#>              n2 = 135
#>               N = 270
#>           delta = 0.5
#>              sd = 1
#>               p = 0.6, 0.4
#>             rho = 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#>            Test = AN
#> 

# With continuity correction
ss2MixedContinuousBinary(
  delta = 0.5,
  sd = 1,
  p1 = 0.6,
  p2 = 0.4,
  rho = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.1,
  Test = 'ANc'
)
#> 
#> Sample size calculation for mixed continuous and binary co-primary endpoints
#> 
#>              n1 = 143
#>              n2 = 143
#>               N = 286
#>           delta = 0.5
#>              sd = 1
#>               p = 0.6, 0.4
#>             rho = 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#>            Test = ANc
#> 

# \donttest{
# Fisher's exact test (computationally intensive)
ss2MixedContinuousBinary(
  delta = 0.5,
  sd = 1,
  p1 = 0.6,
  p2 = 0.4,
  rho = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.1,
  Test = 'Fisher',
  nMC = 5000
)
#> 
#> Sample size calculation for mixed continuous and binary co-primary endpoints
#> 
#>              n1 = 143
#>              n2 = 143
#>               N = 286
#>           delta = 0.5
#>              sd = 1
#>               p = 0.6, 0.4
#>             rho = 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#>            Test = Fisher
#>             nMC = 5000
#> 
# }
```
