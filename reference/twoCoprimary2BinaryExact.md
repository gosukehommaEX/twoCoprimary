# Unified Interface for Two Co-Primary Binary Endpoints (Exact Methods)

This function provides a unified interface for both power calculation
and sample size determination for two co-primary binary endpoints using
exact inference methods.

## Usage

``` r
twoCoprimary2BinaryExact(
  n1 = NULL,
  n2 = NULL,
  p11,
  p12,
  p21,
  p22,
  rho1,
  rho2,
  power = NULL,
  r = NULL,
  alpha = 0.025,
  Test = "Fisher"
)
```

## Arguments

- n1:

  Sample size for group 1 (treatment group). If NULL, will be
  calculated.

- n2:

  Sample size for group 2 (control group). If NULL, will be calculated.

- p11:

  True response probability for endpoint 1 in group 1

- p12:

  True response probability for endpoint 2 in group 1

- p21:

  True response probability for endpoint 1 in group 2

- p22:

  True response probability for endpoint 2 in group 2

- rho1:

  Correlation between endpoints 1 and 2 in group 1

- rho2:

  Correlation between endpoints 1 and 2 in group 2

- power:

  Target power (1 - beta). If NULL, will be calculated.

- r:

  Allocation ratio (n1/n2). Required when calculating sample size.

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- Test:

  Test method: "Fisher" (Fisher's exact test), "Chisq" (Chi-squared
  test), "Z-pooled" (Z-pooled exact unconditional test), or "Boschloo"
  (Boschloo's exact unconditional test)

## Value

An object of class "twoCoprimary" containing either:

- Power calculation results (when n1 and n2 are specified)

- Sample size calculation results (when power and r are specified)

## Details

This function serves as a unified interface similar to
[`power.prop.test()`](https://rdrr.io/r/stats/power.prop.test.html). The
function determines the operation mode based on which parameters are
NULL.

Exactly one of {(n1, n2), (power, r)} must be NULL.

Note: Exact methods are computationally intensive and may take
considerable time, especially for large sample sizes.

## Examples

``` r
# \donttest{
# Calculate power given sample sizes
twoCoprimary2BinaryExact(
  n1 = 50, n2 = 50,
  p11 = 0.5, p12 = 0.4,
  p21 = 0.3, p22 = 0.2,
  rho1 = 0.5, rho2 = 0.5,
  alpha = 0.025, Test = "Fisher"
)
#> 
#> Power calculation for two binary co-primary endpoints
#> 
#>              n1 = 50
#>              n2 = 50
#>     p (group 1) = 0.5, 0.4
#>     p (group 2) = 0.3, 0.2
#>             rho = 0.5, 0.5
#>           alpha = 0.025
#>            Test = Fisher
#>          power1 = 0.46345
#>          power2 = 0.515232
#>  powerCoprimary = 0.321793
#> 

# Calculate sample size given target power
twoCoprimary2BinaryExact(
  p11 = 0.5, p12 = 0.4,
  p21 = 0.3, p22 = 0.2,
  rho1 = 0.5, rho2 = 0.5,
  power = 0.8, r = 1,
  alpha = 0.025, Test = "Chisq"
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
#>            Test = Chisq
#> 
# }
```
