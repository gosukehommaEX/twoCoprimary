# Power Calculation for Two Co-Primary Endpoints (Count and Continuous)

Calculates the power for a two-arm superiority trial with one
overdispersed count co-primary endpoint and one continuous co-primary
endpoint, as described in Homma and Yoshida (2024).

## Usage

``` r
power2MixedCountContinuous(
  n1,
  n2,
  r1,
  r2,
  nu,
  t,
  mu1,
  mu2,
  sd,
  rho1,
  rho2,
  alpha
)
```

## Arguments

- n1:

  Sample size for group 1 (test group)

- n2:

  Sample size for group 2 (control group)

- r1:

  Mean rate (events per unit time) for the treatment group (count
  endpoint)

- r2:

  Mean rate (events per unit time) for the control group (count
  endpoint)

- nu:

  Common dispersion parameter for the negative binomial distribution (nu
  \> 0)

- t:

  Common follow-up time period

- mu1:

  Mean for group 1 (continuous endpoint)

- mu2:

  Mean for group 2 (continuous endpoint)

- sd:

  Common standard deviation for the continuous endpoint

- rho1:

  Correlation between count and continuous outcomes for treatment group

- rho2:

  Correlation between count and continuous outcomes for control group

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

## Value

A data frame with the following columns:

- n1:

  Sample size for group 1

- n2:

  Sample size for group 2

- r1:

  Mean rate in group 1 for count endpoint

- r2:

  Mean rate in group 2 for count endpoint

- nu:

  Dispersion parameter

- t:

  Follow-up time

- mu1:

  Mean in group 1 for continuous endpoint

- mu2:

  Mean in group 2 for continuous endpoint

- sd:

  Standard deviation for continuous endpoint

- rho1:

  Correlation for group 1

- rho2:

  Correlation for group 2

- alpha:

  One-sided significance level

- powerCount:

  Power for the count endpoint alone

- powerCont:

  Power for the continuous endpoint alone

- powerCoprimary:

  Power for both co-primary endpoints

## Details

The test statistics are (equation 7 in Homma and Yoshida 2024): \$\$Z_1
= \frac{\hat{\beta}\_1}{\sqrt{Var(\hat{\beta}\_1)}}, \quad Z_2 =
\frac{\hat{\delta}}{\sigma\sqrt{(1+\kappa)/(\kappa n_0)}}\$\$

The joint distribution of (Z1, Z2) follows an asymptotic bivariate
normal distribution with correlation gamma (equation 11): \$\$\gamma =
\sum\_{j=0,1} \frac{n_0 \rho_j \sqrt{1+\lambda_j/\nu}} {n_j
\sqrt{\lambda_j V_a} \sqrt{(1+\kappa)/\kappa}}\$\$

where \\\lambda_j = r_j \times t\\.

The correlation bounds are automatically checked using
[`corrbound2MixedCountContinuous`](https://gosukehommaEX.github.io/twoCoprimary/reference/corrbound2MixedCountContinuous.md).

## References

Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical
trials with two co-primary endpoints including overdispersed count and
continuous outcomes. *Pharmaceutical Statistics*, 23(1), 46-59.

## Examples

``` r
# Power calculation with moderate correlation
power2MixedCountContinuous(
  n1 = 300,
  n2 = 300,
  r1 = 1.0,
  r2 = 1.25,
  nu = 0.8,
  t = 1,
  mu1 = -50,
  mu2 = 0,
  sd = 250,
  rho1 = 0.5,
  rho2 = 0.5,
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

# Power calculation with no correlation
power2MixedCountContinuous(
  n1 = 350,
  n2 = 350,
  r1 = 1.0,
  r2 = 1.5,
  nu = 1,
  t = 1,
  mu1 = -40,
  mu2 = 0,
  sd = 200,
  rho1 = 0,
  rho2 = 0,
  alpha = 0.025
)
#> 
#> Power calculation for mixed count and continuous co-primary endpoints
#> 
#>              n1 = 350
#>              n2 = 350
#>              sd = 200
#>            rate = 1, 1.5
#>              nu = 1
#>               t = 1
#>              mu = -40, 0
#>             rho = 0, 0
#>           alpha = 0.025
#>       powerCont = 0.753576
#>      powerCount = 0.977329
#>  powerCoprimary = 0.736492
#> 

# Unbalanced design
power2MixedCountContinuous(
  n1 = 400,
  n2 = 200,
  r1 = 1,
  r2 = 1.25,
  nu = 1,
  t = 1,
  mu1 = -50,
  mu2 = 0,
  sd = 250,
  rho1 = 0.6,
  rho2 = 0.6,
  alpha = 0.025
)
#> 
#> Power calculation for mixed count and continuous co-primary endpoints
#> 
#>              n1 = 400
#>              n2 = 200
#>              sd = 250
#>            rate = 1, 1.25
#>              nu = 1
#>               t = 1
#>              mu = -50, 0
#>             rho = 0.6, 0.6
#>           alpha = 0.025
#>       powerCont = 0.636619
#>      powerCount = 0.470483
#>  powerCoprimary = 0.393625
#> 
```
