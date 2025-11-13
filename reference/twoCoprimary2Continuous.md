# Unified Interface for Two Co-Primary Continuous Endpoints

This function provides a unified interface for both power calculation
and sample size determination for two co-primary continuous endpoints.
Depending on which parameters are provided (sample sizes or power), the
function automatically determines whether to calculate power or sample
size.

## Usage

``` r
twoCoprimary2Continuous(
  n1 = NULL,
  n2 = NULL,
  delta1,
  delta2,
  sd1,
  sd2,
  rho,
  power = NULL,
  r = NULL,
  alpha = 0.025,
  known_var = TRUE,
  nMC = 10000
)
```

## Arguments

- n1:

  Sample size for group 1 (treatment group). If NULL, will be
  calculated.

- n2:

  Sample size for group 2 (control group). If NULL, will be calculated.

- delta1:

  Mean difference for the first endpoint

- delta2:

  Mean difference for the second endpoint

- sd1:

  Common standard deviation for the first endpoint

- sd2:

  Common standard deviation for the second endpoint

- rho:

  Common correlation between the two outcomes

- power:

  Target power (1 - beta). If NULL, will be calculated.

- r:

  Allocation ratio (n1/n2). Required when calculating sample size.

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- known_var:

  Logical indicating whether variance is known (TRUE) or unknown
  (FALSE). Default is TRUE.

- nMC:

  Number of Monte Carlo simulations when known_var = FALSE. Default is
  10000.

## Value

An object of class "twoCoprimary" containing either:

- Power calculation results (when n1 and n2 are specified)

- Sample size calculation results (when power and r are specified)

## Details

This function serves as a unified interface similar to
[`power.prop.test()`](https://rdrr.io/r/stats/power.prop.test.html). The
function determines the operation mode based on which parameters are
NULL:

- If n1 and n2 are provided and power is NULL: calculates power

- If power and r are provided and n1/n2 are NULL: calculates sample size

Exactly one of {(n1, n2), (power, r)} must be NULL.

## Examples

``` r
# Calculate power given sample sizes
twoCoprimary2Continuous(
  n1 = 100, n2 = 100,
  delta1 = 0.5, delta2 = 0.5,
  sd1 = 1, sd2 = 1,
  rho = 0.3, alpha = 0.025,
  known_var = TRUE
)
#> 
#> Power calculation for two continuous co-primary endpoints
#> 
#>              n1 = 100
#>              n2 = 100
#>           delta = 0.5, 0.5
#>              sd = 1, 1
#>             rho = 0.3
#>           alpha = 0.025
#>       known_var = TRUE
#>          power1 = 0.942438
#>          power2 = 0.942438
#>  powerCoprimary = 0.893807
#> 

# Calculate sample size given target power
twoCoprimary2Continuous(
  delta1 = 0.5, delta2 = 0.5,
  sd1 = 1, sd2 = 1,
  rho = 0.3, power = 0.8,
  r = 1, alpha = 0.025,
  known_var = TRUE
)
#> 
#> Sample size calculation for two continuous co-primary endpoints
#> 
#>              n1 = 81
#>              n2 = 81
#>               N = 162
#>           delta = 0.5, 0.5
#>              sd = 1, 1
#>             rho = 0.3
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>       known_var = TRUE
#> 
```
