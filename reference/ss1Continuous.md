# Sample Size Calculation for a Single Continuous Endpoint

Calculates the required sample size for a two-arm superiority trial with
a single continuous endpoint using the standard formula for normally
distributed outcomes.

## Usage

``` r
ss1Continuous(delta, sd, r, alpha, beta)
```

## Arguments

- delta:

  Mean difference between treatment groups (treatment effect)

- sd:

  Common standard deviation for the continuous endpoint

- r:

  Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r
  \> 0)

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- beta:

  Target type II error rate (typically 0.1 or 0.2)

## Value

A data frame with the following columns:

- delta:

  Mean difference (treatment effect)

- sd:

  Common standard deviation

- r:

  Allocation ratio

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

The required sample size for group 2 is calculated using the standard
formula: \$\$n_2 = \left\lceil \frac{(1 + 1/r) \sigma^2 (z\_\alpha +
z\_\beta)^2}{\delta^2} \right\rceil\$\$ where \\z\_\alpha\\ and
\\z\_\beta\\ are the quantiles of the standard normal distribution
corresponding to the one-sided significance level \\\alpha\\ and type II
error rate \\\beta\\, respectively. The sample size for group 1 is \\n_1
= \lceil r \times n_2 \rceil\\.

## Examples

``` r
# Balanced design with 1:1 allocation
ss1Continuous(delta = 0.4, sd = 1, r = 1, alpha = 0.025, beta = 0.1)
#> 
#> Sample size calculation for single continuous endpoint
#> 
#>              n1 = 132
#>              n2 = 132
#>               N = 264
#>           delta = 0.4
#>              sd = 1
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#> 

# Unbalanced design with 2:1 allocation
ss1Continuous(delta = 0.5, sd = 1.2, r = 2, alpha = 0.025, beta = 0.2)
#> 
#> Sample size calculation for single continuous endpoint
#> 
#>              n1 = 136
#>              n2 = 68
#>               N = 204
#>           delta = 0.5
#>              sd = 1.2
#>      allocation = 2
#>           alpha = 0.025
#>            beta = 0.2
#> 

# Large treatment effect
ss1Continuous(delta = 0.8, sd = 1, r = 1, alpha = 0.025, beta = 0.1)
#> 
#> Sample size calculation for single continuous endpoint
#> 
#>              n1 = 33
#>              n2 = 33
#>               N = 66
#>           delta = 0.8
#>              sd = 1
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#> 
```
