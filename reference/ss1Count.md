# Sample Size Calculation for a Single Count Endpoint (Negative Binomial)

Calculates the required sample size for a two-arm superiority trial with
a single overdispersed count endpoint following a negative binomial
distribution, as described in Homma and Yoshida (2024).

## Usage

``` r
ss1Count(r1, r2, nu, t, r, alpha, beta)
```

## Arguments

- r1:

  Mean rate (events per unit time) for the treatment group

- r2:

  Mean rate (events per unit time) for the control group

- nu:

  Common dispersion parameter for the negative binomial distribution (nu
  \> 0)

- t:

  Common follow-up time period

- r:

  Allocation ratio (treatment:control = r:1, where r \> 0)

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- beta:

  Target type II error rate (typically 0.1 or 0.2)

## Value

A data frame with the following columns:

- r1:

  Mean rate for treatment group

- r2:

  Mean rate for control group

- nu:

  Dispersion parameter

- t:

  Follow-up time

- r:

  Allocation ratio

- alpha:

  One-sided significance level

- beta:

  Type II error rate

- n1:

  Required sample size for treatment group

- n2:

  Required sample size for control group

- N:

  Total sample size (n2 + n1)

## Details

The test statistic for the negative binomial rate ratio is: \$\$Z_1 =
\frac{\hat{\beta}\_1}{\sqrt{Var(\hat{\beta}\_1)}}\$\$ where
\\\hat{\beta}\_1 = \log(\bar{Y}\_1) - \log(\bar{Y}\_2)\\ and the
variance is: \$\$Var(\hat{\beta}\_1) =
\frac{1}{n_2}\left\[\frac{1}{t}\left(\frac{1}{r_2} + \frac{1}{r \cdot
r_1}\right) + \frac{1+r}{\nu \cdot r}\right\]\$\$

This is equation (8) in Homma and Yoshida (2024).

## References

Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical
trials with two co-primary endpoints including overdispersed count and
continuous outcomes. *Pharmaceutical Statistics*, 23(1), 46-59.

## Examples

``` r
# Sample size for count endpoint with nu = 0.8
ss1Count(r1 = 1.0, r2 = 1.25, nu = 0.8, t = 1, r = 1,
         alpha = 0.025, beta = 0.1)
#> 
#> Sample size calculation for single count endpoint
#> 
#>              n1 = 908
#>              n2 = 908
#>               N = 1816
#>            rate = 1, 1.25
#>              nu = 0.8
#>               t = 1
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#> 

# Unbalanced design with 2:1 allocation
ss1Count(r1 = 1.0, r2 = 1.5, nu = 1.0, t = 1, r = 2,
         alpha = 0.025, beta = 0.2)
#> 
#> Sample size calculation for single count endpoint
#> 
#>              n1 = 256
#>              n2 = 128
#>               N = 384
#>            rate = 1, 1.5
#>              nu = 1
#>               t = 1
#>      allocation = 2
#>           alpha = 0.025
#>            beta = 0.2
#> 

# Higher dispersion
ss1Count(r1 = 1.5, r2 = 2.0, nu = 3.0, t = 1, r = 1,
         alpha = 0.025, beta = 0.1)
#> 
#> Sample size calculation for single count endpoint
#> 
#>              n1 = 233
#>              n2 = 233
#>               N = 466
#>            rate = 1.5, 2
#>              nu = 3
#>               t = 1
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#> 
```
