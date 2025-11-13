# Sample Size Calculation for Two Co-Primary Continuous Endpoints (Approximate)

Calculates the required sample size for a two-arm superiority trial with
two co-primary continuous endpoints using sequential search algorithm.

## Usage

``` r
ss2Continuous(
  delta1,
  delta2,
  sd1,
  sd2,
  rho,
  r,
  alpha,
  beta,
  known_var = TRUE,
  nMC = 10000
)
```

## Arguments

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

- r:

  Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r
  \> 0)

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- beta:

  Target type II error rate (typically 0.1 or 0.2)

- known_var:

  Logical value indicating whether variance is known (TRUE) or unknown
  (FALSE). If TRUE, power is calculated analytically; if FALSE, Monte
  Carlo simulation is used

- nMC:

  Number of Monte Carlo simulations when known_var = FALSE (default is
  10000)

## Value

A data frame with the following columns:

- delta1, delta2:

  Mean differences

- sd1, sd2:

  Standard deviations

- rho:

  Correlation

- r:

  Allocation ratio

- alpha:

  One-sided significance level

- beta:

  Type II error rate

- known_var:

  Variance assumption

- nMC:

  Number of Monte Carlo simulations (NA if known_var = TRUE)

- n1:

  Required sample size for group 1

- n2:

  Required sample size for group 2

- N:

  Total sample size (n1 + n2)

## Details

This function uses a sequential search algorithm (Homma and Yoshida
2025, Algorithm 1) for both known and unknown variance cases:

**Step 1:** Initialize with sample sizes from single endpoint formulas.

**Step 2:** Use sequential search:

- Calculate power at initial sample size

- If power \>= target: decrease n2 until power \< target, then add 1
  back

- If power \< target: increase n2 until power \>= target

**Step 3:** Return final sample sizes.

For known variance, the standardized test statistics are: \$\$Z_k =
\frac{\delta_k}{\sigma_k \sqrt{1/n_1 + 1/n_2}}\$\$ For unknown variance,
t-statistics with \\\nu = n_1 + n_2 - 2\\ degrees of freedom are used,
and power is calculated using Monte Carlo simulation following Sozu et
al. (2011).

## References

Sozu, T., Sugimoto, T., & Hamasaki, T. (2011). Sample size determination
in superiority clinical trials with multiple co-primary correlated
endpoints. *Journal of Biopharmaceutical Statistics*, 21(4), 650-668.

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(1), 1-19.

## Examples

``` r
# Sample size calculation with known variance
ss2Continuous(
  delta1 = 0.2,
  delta2 = 0.2,
  sd1 = 1,
  sd2 = 1,
  rho = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.1,
  known_var = TRUE
)
#> 
#> Sample size calculation for two continuous co-primary endpoints
#> 
#>              n1 = 626
#>              n2 = 626
#>               N = 1252
#>           delta = 0.2, 0.2
#>              sd = 1, 1
#>             rho = 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#>       known_var = TRUE
#> 

# Sample size calculation with unequal allocation
ss2Continuous(
  delta1 = 0.3,
  delta2 = 0.25,
  sd1 = 1,
  sd2 = 1,
  rho = 0.3,
  r = 2,
  alpha = 0.025,
  beta = 0.2,
  known_var = TRUE
)
#> 
#> Sample size calculation for two continuous co-primary endpoints
#> 
#>              n1 = 418
#>              n2 = 209
#>               N = 627
#>           delta = 0.3, 0.25
#>              sd = 1, 1
#>             rho = 0.3
#>      allocation = 2
#>           alpha = 0.025
#>            beta = 0.2
#>       known_var = TRUE
#> 

# \donttest{
# Sample size calculation with unknown variance (uses sequential search)
ss2Continuous(
  delta1 = 0.5,
  delta2 = 0.4,
  sd1 = 1,
  sd2 = 1,
  rho = 0.4,
  r = 1,
  alpha = 0.025,
  beta = 0.1,
  known_var = FALSE,
  nMC = 10000
)
#> 
#> Sample size calculation for two continuous co-primary endpoints
#> 
#>              n1 = 138
#>              n2 = 138
#>               N = 276
#>           delta = 0.5, 0.4
#>              sd = 1, 1
#>             rho = 0.4
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#>       known_var = FALSE
#>             nMC = 10000
#> 
# }
```
