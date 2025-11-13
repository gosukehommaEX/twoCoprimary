# Power Calculation for Two Co-Primary Mixed Endpoints

Calculates the power for a two-arm superiority trial with two co-primary
endpoints where one is continuous and one is binary, as described in
Sozu et al. (2012).

## Usage

``` r
power2MixedContinuousBinary(
  n1,
  n2,
  delta,
  sd,
  p1,
  p2,
  rho,
  alpha,
  Test,
  nMC = 10000
)
```

## Arguments

- n1:

  Sample size for group 1 (test group)

- n2:

  Sample size for group 2 (control group)

- delta:

  Mean difference for the continuous endpoint (group 1 - group 2)

- sd:

  Common standard deviation for the continuous endpoint

- p1:

  Probability of response in group 1 for the binary endpoint (0 \< p1 \<
  1)

- p2:

  Probability of response in group 2 for the binary endpoint (0 \< p2 \<
  1)

- rho:

  Biserial correlation between the latent continuous variable underlying
  the binary endpoint and the observed continuous endpoint

- alpha:

  One-sided significance level (typically 0.025 or 0.05)

- Test:

  Statistical testing method for the binary endpoint. One of:

  - `"AN"`: Asymptotic normal method without continuity correction

  - `"ANc"`: Asymptotic normal method with continuity correction

  - `"AS"`: Arcsine method without continuity correction

  - `"ASc"`: Arcsine method with continuity correction

  - `"Fisher"`: Fisher's exact test (Monte Carlo simulation required)

- nMC:

  Number of Monte Carlo replications when Test = "Fisher" (default:
  10000)

## Value

A data frame with the following columns:

- n1:

  Sample size for group 1

- n2:

  Sample size for group 2

- delta:

  Mean difference for continuous endpoint

- sd:

  Standard deviation for continuous endpoint

- p1:

  Response probability in group 1 for binary endpoint

- p2:

  Response probability in group 2 for binary endpoint

- rho:

  Biserial correlation

- alpha:

  One-sided significance level

- Test:

  Testing method used for binary endpoint

- powerCont:

  Power for the continuous endpoint alone

- powerBin:

  Power for the binary endpoint alone

- powerCoprimary:

  Power for both co-primary endpoints

## Details

This function implements the power calculation for mixed endpoints (one
continuous and one binary) as described in Sozu et al. (2012). The
method assumes that the binary variable is derived from a latent
continuous variable via dichotomization at a threshold point.

**For Fisher's exact test**, Monte Carlo simulation is used because
exact calculation is computationally intensive. The continuous endpoint
is analyzed using t-test, and the binary endpoint uses Fisher's exact
test.

**For asymptotic methods (AN, ANc, AS, ASc)**, analytical formulas are
used based on bivariate normal approximation. The correlation between
test statistics depends on the biserial correlation rho and the specific
testing method.

**Biserial Correlation:** The biserial correlation rho represents the
correlation between the latent continuous variable underlying the binary
endpoint and the observed continuous endpoint. This is not the same as
the point-biserial correlation observed in the data.

## References

Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination
in clinical trials with multiple co-primary endpoints including mixed
continuous and binary variables. *Biometrical Journal*, 54(5), 716-729.

## Examples

``` r
# Power calculation using asymptotic normal method
power2MixedContinuousBinary(
  n1 = 100,
  n2 = 100,
  delta = 0.5,
  sd = 1,
  p1 = 0.6,
  p2 = 0.4,
  rho = 0.5,
  alpha = 0.025,
  Test = 'AN'
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

# \donttest{
# Power calculation with Fisher's exact test (computationally intensive)
power2MixedContinuousBinary(
  n1 = 50,
  n2 = 50,
  delta = 0.5,
  sd = 1,
  p1 = 0.6,
  p2 = 0.4,
  rho = 0.5,
  alpha = 0.025,
  Test = 'Fisher',
  nMC = 5000
)
#> 
#> Power calculation for mixed continuous and binary co-primary endpoints
#> 
#>              n1 = 50
#>              n2 = 50
#>           delta = 0.5
#>              sd = 1
#>               p = 0.6, 0.4
#>             rho = 0.5
#>           alpha = 0.025
#>            Test = Fisher
#>       powerCont = 0.705414
#>        powerBin = 0.440109
#>  powerCoprimary = 0.364131
#> 
# }
```
