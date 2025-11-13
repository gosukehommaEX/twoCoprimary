# Unified Interface for Mixed Continuous and Binary Co-Primary Endpoints

This function provides a unified interface for both power calculation
and sample size determination for trials with one continuous and one
binary co-primary endpoint.

## Usage

``` r
twoCoprimary2MixedContinuousBinary(
  n1 = NULL,
  n2 = NULL,
  delta,
  sd,
  p1,
  p2,
  rho,
  power = NULL,
  r = NULL,
  alpha = 0.025,
  Test = "AN",
  nMC = 10000
)
```

## Arguments

- n1:

  Sample size for group 1 (treatment group). If NULL, will be
  calculated.

- n2:

  Sample size for group 2 (control group). If NULL, will be calculated.

- delta:

  Mean difference for the continuous endpoint

- sd:

  Common standard deviation for the continuous endpoint

- p1:

  True response probability for the binary endpoint in group 1

- p2:

  True response probability for the binary endpoint in group 2

- rho:

  Biserial correlation between the continuous endpoint and the latent
  continuous variable underlying the binary endpoint

- power:

  Target power (1 - beta). If NULL, will be calculated.

- r:

  Allocation ratio (n1/n2). Required when calculating sample size.

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- Test:

  Test method for the binary endpoint: "AN" (asymptotic normal), "ANc"
  (with continuity correction), "AS" (arcsine), "ASc" (arcsine with
  continuity correction), or "Fisher" (Fisher's exact test)

- nMC:

  Number of Monte Carlo simulations when Test = "Fisher". Default is
  10000.

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

The biserial correlation rho represents the correlation between the
observed continuous endpoint and the latent continuous variable
underlying the binary endpoint.

## Examples

``` r
# Calculate power given sample sizes
twoCoprimary2MixedContinuousBinary(
  n1 = 100, n2 = 100,
  delta = 0.5, sd = 1,
  p1 = 0.6, p2 = 0.4,
  rho = 0.5,
  alpha = 0.025, Test = "AN"
)
#> 
#> Power calculation for mixed continuous and binary co-primary endpoints
#> 
#>              n1 = 100
#>              n2 = 100
#>           delta = 0.5
#>              sd = 1
#>               p = 0.6, 0.4
#>             rho = 0.5
#>           alpha = 0.025
#>            Test = AN
#>       powerCont = 0.942438
#>        powerBin = 0.812291
#>  powerCoprimary = 0.781111
#> 

# Calculate sample size given target power
twoCoprimary2MixedContinuousBinary(
  delta = 0.5, sd = 1,
  p1 = 0.6, p2 = 0.4,
  rho = 0.5,
  power = 0.8, r = 1,
  alpha = 0.025, Test = "AN"
)
#> 
#> Sample size calculation for mixed continuous and binary co-primary endpoints
#> 
#>              n1 = 105
#>              n2 = 105
#>               N = 210
#>           delta = 0.5
#>              sd = 1
#>               p = 0.6, 0.4
#>             rho = 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = AN
#> 
```
