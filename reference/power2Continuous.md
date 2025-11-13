# Power Calculation for Two Co-Primary Continuous Endpoints

Calculates the power for a two-arm superiority trial with two co-primary
continuous endpoints, as described in Sozu et al. (2011).

## Usage

``` r
power2Continuous(
  n1,
  n2,
  delta1,
  delta2,
  sd1,
  sd2,
  rho,
  alpha,
  known_var = TRUE,
  nMC = 10000
)
```

## Arguments

- n1:

  Sample size for group 1 (test group)

- n2:

  Sample size for group 2 (control group)

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

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- known_var:

  Logical value indicating whether variance is known (TRUE) or unknown
  (FALSE). If TRUE, power is calculated analytically; otherwise, Monte
  Carlo simulation is used for unknown variance

- nMC:

  Number of Monte Carlo simulations when known_var = FALSE (default is
  10000)

## Value

A data frame with the following columns:

- n1:

  Sample size for group 1

- n2:

  Sample size for group 2

- delta1:

  Mean difference for endpoint 1

- delta2:

  Mean difference for endpoint 2

- sd1:

  Standard deviation for endpoint 1

- sd2:

  Standard deviation for endpoint 2

- rho:

  Correlation between endpoints

- alpha:

  One-sided significance level

- known_var:

  Variance assumption

- nMC:

  Number of Monte Carlo simulations (NA if known_var = TRUE)

- power1:

  Power for the first endpoint alone

- power2:

  Power for the second endpoint alone

- powerCoprimary:

  Power for both co-primary endpoints

## Details

For known variance, the power is calculated using the bivariate normal
distribution as described in Sozu et al. (2011). The test statistics
are: \$\$Z_k = \frac{\delta_k}{\sigma_k \sqrt{1/n_1 + 1/n_2}}\$\$ for k
= 1, 2. The co-primary power is: \$\$1 - \beta =
\Phi_2\left(-z\_{1-\alpha} + Z_1, -z\_{1-\alpha} + Z_2 \mid
\rho\right)\$\$ where \\\Phi_2\\ is the cumulative distribution function
of the bivariate standard normal distribution.

For unknown variance, Monte Carlo simulation is used with
Wishart-distributed variance-covariance matrices to account for variance
estimation uncertainty, following equation (6) in Sozu et al. (2011):
\$\$\text{Power} = E_W\left\[\Phi_2(-c_1^\*\sqrt{w\_{11}},
-c_2^\*\sqrt{w\_{22}} \| \rho)\right\]\$\$ where \\c_k^\* =
t\_{\alpha,\nu}\sqrt{\frac{1}{\nu}} - \frac{Z_k}{\sqrt{w\_{kk}}}\\ and
\\W\\ follows a Wishart distribution with \\\nu = n_1 + n_2 - 2\\
degrees of freedom.

## References

Sozu, T., Sugimoto, T., & Hamasaki, T. (2011). Sample size determination
in superiority clinical trials with multiple co-primary correlated
endpoints. *Journal of Biopharmaceutical Statistics*, 21(4), 650-668.

## Examples

``` r
# Example parameters for comparison across methods
n1_ex <- 100
n2_ex <- 100
delta1_ex <- 0.5
delta2_ex <- 0.5
sd1_ex <- 1
sd2_ex <- 1
rho_ex <- 0.3
alpha_ex <- 0.025

# Power calculation with known variance
power2Continuous(
  n1 = n1_ex,
  n2 = n2_ex,
  delta1 = delta1_ex,
  delta2 = delta2_ex,
  sd1 = sd1_ex,
  sd2 = sd2_ex,
  rho = rho_ex,
  alpha = alpha_ex,
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

# \donttest{
# Power calculation with unknown variance (Monte Carlo)
power2Continuous(
  n1 = n1_ex,
  n2 = n2_ex,
  delta1 = delta1_ex,
  delta2 = delta2_ex,
  sd1 = sd1_ex,
  sd2 = sd2_ex,
  rho = rho_ex,
  alpha = alpha_ex,
  known_var = FALSE,
  nMC = 10000
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
#>       known_var = FALSE
#>             nMC = 10000
#>          power1 = 0.940427
#>          power2 = 0.940427
#>  powerCoprimary = 0.890195
#> 
# }
```
