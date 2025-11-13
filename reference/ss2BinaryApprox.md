# Sample Size Calculation for Two Co-Primary Binary Endpoints (Asymptotic Approximation)

Calculates the required sample size for a two-arm superiority trial with
two co-primary binary endpoints using asymptotic normal approximation or
arcsine transformation, as described in Sozu et al. (2010).

## Usage

``` r
ss2BinaryApprox(p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test)
```

## Arguments

- p11:

  True probability of responders in group 1 for the first outcome (0 \<
  p11 \< 1)

- p12:

  True probability of responders in group 1 for the second outcome (0 \<
  p12 \< 1)

- p21:

  True probability of responders in group 2 for the first outcome (0 \<
  p21 \< 1)

- p22:

  True probability of responders in group 2 for the second outcome (0 \<
  p22 \< 1)

- rho1:

  Correlation between the two outcomes for group 1

- rho2:

  Correlation between the two outcomes for group 2

- r:

  Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r
  \> 0)

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- beta:

  Target type II error rate (typically 0.1 or 0.2)

- Test:

  Statistical testing method. One of:

  - `"AN"`: Asymptotic normal method without continuity correction

  - `"ANc"`: Asymptotic normal method with continuity correction

  - `"AS"`: Arcsine method without continuity correction

  - `"ASc"`: Arcsine method with continuity correction

## Value

A data frame with the following columns:

- p11, p12, p21, p22:

  Response probabilities

- rho1, rho2:

  Correlations

- r:

  Allocation ratio

- alpha:

  One-sided significance level

- beta:

  Type II error rate

- Test:

  Testing method used

- n1:

  Required sample size for group 1

- n2:

  Required sample size for group 2

- N:

  Total sample size (n1 + n2)

## Details

This function uses a sequential search algorithm (Homma and Yoshida
2025, Algorithm 1) to find the minimum sample size:

**Step 1:** Initialize with sample sizes from single endpoint formulas.

**Step 2:** Use sequential search:

- Calculate power at initial sample size

- If power \>= target: decrease n2 until power \< target, then add 1
  back

- If power \< target: increase n2 until power \>= target

**Step 3:** Return final sample sizes.

The asymptotic normal (AN) and arcsine (AS) methods use normal
approximation with or without continuity correction. For small sample
sizes or extreme probabilities, consider using exact methods via
[`ss2BinaryExact`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2BinaryExact.md).

## References

Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination
in clinical trials with multiple co-primary binary endpoints.
*Statistics in Medicine*, 29(21), 2169-2179.

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(1), 1-19.

## Examples

``` r
# Sample size calculation using asymptotic normal method
ss2BinaryApprox(
  p11 = 0.5,
  p12 = 0.4,
  p21 = 0.3,
  p22 = 0.2,
  rho1 = 0.5,
  rho2 = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  Test = 'AN'
)
#> 
#> Sample size calculation for two binary co-primary endpoints
#> 
#>              n1 = 109
#>              n2 = 109
#>               N = 218
#>     p (group 1) = 0.5, 0.4
#>     p (group 2) = 0.3, 0.2
#>             rho = 0.5, 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = AN
#> 

# With continuity correction
ss2BinaryApprox(
  p11 = 0.5,
  p12 = 0.4,
  p21 = 0.3,
  p22 = 0.2,
  rho1 = 0.5,
  rho2 = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  Test = 'ANc'
)
#> 
#> Sample size calculation for two binary co-primary endpoints
#> 
#>              n1 = 119
#>              n2 = 119
#>               N = 238
#>     p (group 1) = 0.5, 0.4
#>     p (group 2) = 0.3, 0.2
#>             rho = 0.5, 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = ANc
#> 

# Using arcsine transformation
ss2BinaryApprox(
  p11 = 0.5,
  p12 = 0.4,
  p21 = 0.3,
  p22 = 0.2,
  rho1 = 0.5,
  rho2 = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  Test = 'AS'
)
#> 
#> Sample size calculation for two binary co-primary endpoints
#> 
#>              n1 = 109
#>              n2 = 109
#>               N = 218
#>     p (group 1) = 0.5, 0.4
#>     p (group 2) = 0.3, 0.2
#>             rho = 0.5, 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = AS
#> 
```
