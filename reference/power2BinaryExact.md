# Exact Power Calculation for Two Co-Primary Binary Endpoints

Calculates the exact power for a two-arm superiority trial with two
co-primary binary endpoints using the bivariate binomial distribution,
as described in Homma and Yoshida (2025).

## Usage

``` r
power2BinaryExact(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test)
```

## Arguments

- n1:

  Sample size for group 1 (test group)

- n2:

  Sample size for group 2 (control group)

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

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- Test:

  Statistical testing method. One of:

  - `"Chisq"`: One-sided Pearson chi-squared test

  - `"Fisher"`: Fisher exact test

  - `"Fisher-midP"`: Fisher mid-p test

  - `"Z-pool"`: Z-pooled exact unconditional test

  - `"Boschloo"`: Boschloo exact unconditional test

## Value

A data frame with the following columns:

- n1:

  Sample size for group 1

- n2:

  Sample size for group 2

- p11, p12, p21, p22:

  Response probabilities

- rho1, rho2:

  Correlations

- alpha:

  One-sided significance level

- Test:

  Testing method used

- power1:

  Power for the first endpoint alone

- power2:

  Power for the second endpoint alone

- powerCoprimary:

  Exact power for both co-primary endpoints

## Details

This function calculates exact power using equation (9) in Homma and
Yoshida (2025): \$\$power_A(\theta) =
\sum\_{(a\_{1,1},a\_{2,1})\in\mathcal{A}\_1}
\sum\_{(a\_{1,2},a\_{2,2})\in\mathcal{A}\_2} f(a\_{1,1}\|N_1,p\_{1,1})
\times f(a\_{2,1}\|N_2,p\_{2,1}) \times
g(a\_{1,2}\|a\_{1,1},N_1,p\_{1,1},p\_{1,2},\gamma_1) \times
g(a\_{2,2}\|a\_{2,1},N_2,p\_{2,1},p\_{2,2},\gamma_2)\$\$

where \\\mathcal{A}\_k\\ is the rejection region for endpoint k, and
\\(Y\_{j,1}, Y\_{j,2}) \sim BiBin(N_j, p\_{j,1}, p\_{j,2}, \gamma_j)\\
follows the bivariate binomial distribution.

The correlation bounds are automatically checked using
[`corrbound2Binary`](https://gosukehommaEX.github.io/twoCoprimary/reference/corrbound2Binary.md).

## References

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(1), 1-19.

## Examples

``` r
# Exact power calculation using Boschloo test
power2BinaryExact(
  n1 = 100,
  n2 = 50,
  p11 = 0.5,
  p12 = 0.4,
  p21 = 0.3,
  p22 = 0.2,
  rho1 = 0.7,
  rho2 = 0.7,
  alpha = 0.025,
  Test = 'Boschloo'
)
#> 
#> Power calculation for two binary co-primary endpoints
#> 
#>              n1 = 100
#>              n2 = 50
#>     p (group 1) = 0.5, 0.4
#>     p (group 2) = 0.3, 0.2
#>             rho = 0.7, 0.7
#>           alpha = 0.025
#>            Test = Boschloo
#>          power1 = 0.651316
#>          power2 = 0.70332
#>  powerCoprimary = 0.563055
#> 

# Exact power with Fisher exact test
power2BinaryExact(
  n1 = 80,
  n2 = 80,
  p11 = 0.6,
  p12 = 0.5,
  p21 = 0.4,
  p22 = 0.3,
  rho1 = 0.5,
  rho2 = 0.5,
  alpha = 0.025,
  Test = 'Fisher'
)
#> 
#> Power calculation for two binary co-primary endpoints
#> 
#>              n1 = 80
#>              n2 = 80
#>     p (group 1) = 0.6, 0.5
#>     p (group 2) = 0.4, 0.3
#>             rho = 0.5, 0.5
#>           alpha = 0.025
#>            Test = Fisher
#>          power1 = 0.658351
#>          power2 = 0.677271
#>  powerCoprimary = 0.517687
#> 

# \donttest{
# Larger sample sizes (computationally intensive)
power2BinaryExact(
  n1 = 200,
  n2 = 100,
  p11 = 0.5,
  p12 = 0.4,
  p21 = 0.3,
  p22 = 0.2,
  rho1 = 0.6,
  rho2 = 0.6,
  alpha = 0.025,
  Test = 'Chisq'
)
#> 
#> Power calculation for two binary co-primary endpoints
#> 
#>              n1 = 200
#>              n2 = 100
#>     p (group 1) = 0.5, 0.4
#>     p (group 2) = 0.3, 0.2
#>             rho = 0.6, 0.6
#>           alpha = 0.025
#>            Test = Chisq
#>          power1 = 0.9219
#>          power2 = 0.949665
#>  powerCoprimary = 0.892527
#> 
# }
```
