# Create Design Comparison Table for Two Co-Primary Endpoints

Generates a comprehensive table comparing sample sizes or power across
different parameter combinations and correlation values. This function
is useful for sensitivity analyses and exploring how design parameters
affect statistical properties.

## Usage

``` r
design_table(
  param_grid,
  rho_values = c(0, 0.3, 0.5, 0.8),
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  endpoint_type = c("continuous", "binary", "mixed_cont_binary", "mixed_count_cont"),
  Test = "AN",
  known_var = TRUE,
  nMC = 1000,
  output_var = NULL
)
```

## Arguments

- param_grid:

  A data.frame containing parameter combinations. Required columns
  depend on endpoint_type and calculation_mode:

  - For continuous endpoints (sample size): delta1, delta2, sd1, sd2

  - For continuous endpoints (power): n1, n2, delta1, delta2, sd1, sd2

  - For binary endpoints (sample size): p11, p12, p21, p22

  - For binary endpoints (power): n1, n2, p11, p12, p21, p22

  - For mixed continuous-binary (sample size): delta, sd, p1, p2

  - For mixed continuous-binary (power): n1, n2, delta, sd, p1, p2

  - For mixed count-continuous (sample size): r1, r2, nu, t, mu1, mu2,
    sd

  - For mixed count-continuous (power): n1, n2, r1, r2, nu, t, mu1, mu2,
    sd

- rho_values:

  Numeric vector of correlation values to evaluate. Default is c(0, 0.3,
  0.5, 0.8).

- r:

  Allocation ratio (n1/n2). Required for sample size calculation.
  Default is 1.

- alpha:

  One-sided significance level. Default is 0.025.

- beta:

  Type II error rate (1 - power). Required for sample size calculation.
  Default is 0.2 (power = 0.8).

- endpoint_type:

  Character string specifying endpoint type: "continuous", "binary",
  "mixed_cont_binary", or "mixed_count_cont".

- Test:

  Test method for binary endpoints: "AN" (asymptotic normal), "ANc"
  (with continuity correction), "AS" (arcsine), or "ASc". Default is
  "AN". Only used for binary and mixed_cont_binary endpoints.

- known_var:

  Logical indicating whether variance is known for continuous endpoints.
  Default is TRUE.

- nMC:

  Number of Monte Carlo simulations for certain calculations. Default is
  1000.

- output_var:

  Character string specifying which variable to output in the result
  columns: "N" (total sample size, default for sample size calculation)
  or "powerCoprimary" (co-primary power, default for power calculation).

## Value

A data.frame of class "twoCoprimary_table" with:

- Parameter columns (from param_grid)

- Result columns for each correlation value (rho_0.0, rho_0.3, etc.)

## Details

This function performs systematic calculations across all combinations
of parameters specified in param_grid and correlation values in
rho_values.

The calculation mode (sample size vs power) is automatically determined:

- If param_grid contains n1 and n2: calculates power

- Otherwise: calculates sample size (requires r, alpha, beta)

For binary endpoints with two correlations (rho1, rho2), both are set to
the same value from rho_values for each calculation.

The output format follows the style of Sozu et al. (2011), with
parameters displayed in the leftmost columns and results for each
correlation in subsequent columns.

## References

Sozu, T., Kanou, T., Hamada, C., & Yoshimura, I. (2011). Power and
sample size calculations in clinical trials with multiple primary
variables. Japanese Journal of Biometrics, 27, 83-96.

## Examples

``` r
# Sample size calculation for continuous endpoints
param_grid <- expand.grid(
  delta1 = c(0.3, 0.5),
  delta2 = c(0.1, 0.2, 0.3),
  sd1 = c(1.0, 1.5),
  sd2 = c(1.0, 1.5)
)

result <- design_table(
  param_grid = param_grid,
  rho_values = c(0, 0.3, 0.5, 0.8),
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  endpoint_type = "continuous"
)
print(result)
#> 
#> Design Comparison Table for Two Co-Primary Endpoints
#> ======================================================
#> 
#>  delta1 delta2 sd1 sd2 rho_0.0 rho_0.3 rho_0.5 rho_0.8
#>     0.3    0.1 1.0 1.0    3140    3140    3140    3140
#>     0.5    0.1 1.0 1.0    3140    3140    3140    3140
#>     0.3    0.2 1.0 1.0     804     798     794     786
#>     0.5    0.2 1.0 1.0     786     786     786     786
#>     0.3    0.3 1.0 1.0     460     448     436     408
#>     0.5    0.3 1.0 1.0     352     352     350     350
#>     0.3    0.1 1.5 1.0    3142    3140    3140    3140
#>     0.5    0.1 1.5 1.0    3140    3140    3140    3140
#>     0.3    0.2 1.5 1.0    1032    1006     980     916
#>     0.5    0.2 1.5 1.0     792     790     788     786
#>     0.3    0.3 1.5 1.0     804     798     794     786
#>     0.5    0.3 1.5 1.0     418     408     398     374
#>     0.3    0.1 1.0 1.5    7064    7064    7064    7064
#>     0.5    0.1 1.0 1.5    7064    7064    7064    7064
#>     0.3    0.2 1.0 1.5    1768    1768    1766    1766
#>     0.5    0.2 1.0 1.5    1766    1766    1766    1766
#>     0.3    0.3 1.0 1.5     804     798     794     786
#>     0.5    0.3 1.0 1.5     786     786     786     786
#>     0.3    0.1 1.5 1.5    7064    7064    7064    7064
#>     0.5    0.1 1.5 1.5    7064    7064    7064    7064
#>     0.3    0.2 1.5 1.5    1808    1794    1784    1768
#>     0.5    0.2 1.5 1.5    1766    1766    1766    1766
#>     0.3    0.3 1.5 1.5    1032    1006     980     916
#>     0.5    0.3 1.5 1.5     792     790     788     786

# Power calculation for continuous endpoints
param_grid_power <- expand.grid(
  n1 = c(50, 100),
  n2 = c(50, 100),
  delta1 = 0.5,
  delta2 = 0.5,
  sd1 = 1.0,
  sd2 = 1.0
)

result_power <- design_table(
  param_grid = param_grid_power,
  rho_values = c(0, 0.3, 0.5, 0.8),
  alpha = 0.025,
  endpoint_type = "continuous"
)
print(result_power)
#> 
#> Design Comparison Table for Two Co-Primary Endpoints
#> ======================================================
#> 
#>   n1  n2 delta1 delta2 sd1 sd2   rho_0.0   rho_0.3   rho_0.5   rho_0.8
#>   50  50    0.5    0.5   1   1 0.4976088 0.5352018 0.5634879 0.6173565
#>  100  50    0.5    0.5   1   1 0.6772986 0.7002326 0.7190959 0.7573394
#>   50 100    0.5    0.5   1   1 0.6772986 0.7002326 0.7190959 0.7573394
#>  100 100    0.5    0.5   1   1 0.8881885 0.8938066 0.8997323 0.9141060

# Binary endpoints
param_grid_binary <- expand.grid(
  p11 = c(0.6, 0.7),
  p12 = c(0.4, 0.5),
  p21 = c(0.4, 0.5),
  p22 = c(0.2, 0.3)
)

result_binary <- design_table(
  param_grid = param_grid_binary,
  rho_values = c(0.3, 0.5, 0.7),
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  endpoint_type = "binary",
  Test = "AN"
)
print(result_binary)
#> 
#> Design Comparison Table for Two Co-Primary Endpoints
#> ======================================================
#> 
#>  p11 p12 p21 p22 rho_0.3 rho_0.5 rho_0.7
#>  0.6 0.4 0.4 0.2     230     224      NA
#>  0.7 0.4 0.4 0.2     168     166      NA
#>  0.6 0.5 0.4 0.2     196     196      NA
#>  0.7 0.5 0.4 0.2     104     100      NA
#>  0.6 0.4 0.5 0.2     776     776      NA
#>  0.7 0.4 0.5 0.2     224     218      NA
#>  0.6 0.5 0.5 0.2     776     776      NA
#>  0.7 0.5 0.5 0.2     188     188      NA
#>  0.6 0.4 0.4 0.3     714     712      NA
#>  0.7 0.4 0.4 0.3     712     712      NA
#>  0.6 0.5 0.4 0.3     244     238     228
#>  0.7 0.5 0.4 0.3     190     188      NA
#>  0.6 0.4 0.5 0.3     952     928      NA
#>  0.7 0.4 0.5 0.3     714     712      NA
#>  0.6 0.5 0.5 0.3     776     776      NA
#>  0.7 0.5 0.5 0.3     238     232      NA
```
