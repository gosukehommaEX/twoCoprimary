# Exact Sample Size Calculation for Two Co-Primary Binary Endpoints

Calculates the required sample size for a two-arm superiority trial with
two co-primary binary endpoints using exact methods, as described in
Homma and Yoshida (2025).

## Usage

``` r
ss2BinaryExact(p11, p12, p21, p22, rho1, rho2, r, alpha, beta, Test)
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

  - `"Chisq"`: One-sided Pearson chi-squared test

  - `"Fisher"`: Fisher exact test

  - `"Fisher-midP"`: Fisher mid-p test

  - `"Z-pool"`: Z-pooled exact unconditional test

  - `"Boschloo"`: Boschloo exact unconditional test

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

This function uses a sequential search algorithm to find the minimum
sample size that achieves the target power:

**Step 1:** Initialize with sample size from approximate method (AN).
This provides a good starting point for the exact calculation.

**Step 2:** Use sequential search algorithm (Homma and Yoshida 2025,
Algorithm 1):

- Calculate power at initial sample size

- If power \>= target: decrease n2 until power \< target, then add 1
  back

- If power \< target: increase n2 until power \>= target

**Step 3:** Return final sample sizes.

Note: Due to the saw-tooth nature of exact power (power does not
increase monotonically with sample size), this sequential search ensures
the minimum sample size that achieves the target power.

## References

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(1), 1-19.

## Examples

``` r
# Quick example with Chi-squared test (faster)
ss2BinaryExact(
  p11 = 0.6,
  p12 = 0.5,
  p21 = 0.4,
  p22 = 0.3,
  rho1 = 0.3,
  rho2 = 0.3,
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  Test = "Chisq"
)
#> 
#> Sample size calculation for two binary co-primary endpoints
#> 
#>              n1 = 123
#>              n2 = 123
#>               N = 246
#>     p (group 1) = 0.6, 0.5
#>     p (group 2) = 0.4, 0.3
#>             rho = 0.3, 0.3
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = Chisq
#> 

# \donttest{
# More computationally intensive example with Fisher test
ss2BinaryExact(
  p11 = 0.5,
  p12 = 0.4,
  p21 = 0.3,
  p22 = 0.2,
  rho1 = 0.5,
  rho2 = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  Test = "Fisher"
)
#> 
#> Sample size calculation for two binary co-primary endpoints
#> 
#>              n1 = 117
#>              n2 = 117
#>               N = 234
#>     p (group 1) = 0.5, 0.4
#>     p (group 2) = 0.3, 0.2
#>             rho = 0.5, 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = Fisher
#> 
# }
```
