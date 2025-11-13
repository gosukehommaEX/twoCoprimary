# Unified Interface for Mixed Count and Continuous Co-Primary Endpoints

This function provides a unified interface for both power calculation
and sample size determination for trials with one count endpoint
(modeled by negative binomial distribution) and one continuous endpoint.

## Usage

``` r
twoCoprimary2MixedCountContinuous(
  n1 = NULL,
  n2 = NULL,
  r1,
  r2,
  nu,
  t,
  mu1,
  mu2,
  sd,
  rho1,
  rho2,
  power = NULL,
  r = NULL,
  alpha = 0.025
)
```

## Arguments

- n1:

  Sample size for group 1 (treatment group). If NULL, will be
  calculated.

- n2:

  Sample size for group 2 (control group). If NULL, will be calculated.

- r1:

  Event rate per unit time for the count endpoint in group 1

- r2:

  Event rate per unit time for the count endpoint in group 2

- nu:

  Dispersion parameter for the negative binomial distribution (nu \> 0).
  When nu approaches infinity, the distribution converges to Poisson.

- t:

  Follow-up period (time unit)

- mu1:

  Mean for the continuous endpoint in group 1

- mu2:

  Mean for the continuous endpoint in group 2

- sd:

  Common standard deviation for the continuous endpoint

- rho1:

  Correlation between count and continuous endpoints in group 1

- rho2:

  Correlation between count and continuous endpoints in group 2

- power:

  Target power (1 - beta). If NULL, will be calculated.

- r:

  Allocation ratio (n1/n2). Required when calculating sample size.

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

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

The count endpoint is modeled using a negative binomial distribution to
account for overdispersion. The dispersion parameter nu controls the
variance: Var = lambda + lambda^2/nu.

## Examples

``` r
# Calculate power given sample sizes
twoCoprimary2MixedCountContinuous(
  n1 = 300, n2 = 300,
  r1 = 1.0, r2 = 1.25,
  nu = 0.8, t = 1,
  mu1 = -50, mu2 = 0, sd = 250,
  rho1 = 0.5, rho2 = 0.5,
  alpha = 0.025
)
#> 
#> Power calculation for mixed count and continuous co-primary endpoints
#> 
#>              n1 = 300
#>              n2 = 300
#>              sd = 250
#>            rate = 1, 1.25
#>              nu = 0.8
#>               t = 1
#>              mu = -50, 0
#>             rho = 0.5, 0.5
#>           alpha = 0.025
#>       powerCont = 0.687765
#>      powerCount = 0.461715
#>  powerCoprimary = 0.389192
#> 

# Calculate sample size given target power
twoCoprimary2MixedCountContinuous(
  r1 = 1.0, r2 = 1.25,
  nu = 0.8, t = 1,
  mu1 = -50, mu2 = 0, sd = 250,
  rho1 = 0.5, rho2 = 0.5,
  power = 0.8, r = 1,
  alpha = 0.025
)
#> 
#> Sample size calculation for mixed count and continuous co-primary endpoints
#> 
#>              n1 = 705
#>              n2 = 705
#>               N = 1410
#>              sd = 250
#>            rate = 1, 1.25
#>              nu = 0.8
#>               t = 1
#>              mu = -50, 0
#>             rho = 0.5, 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#> 
```
