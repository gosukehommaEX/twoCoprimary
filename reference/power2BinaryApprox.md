# Power Calculation for Two Co-Primary Binary Endpoints (Approximate)

Calculates the power for a two-arm superiority trial with two co-primary
binary endpoints using various asymptotic normal approximation methods,
as described in Sozu et al. (2010).

## Usage

``` r
power2BinaryApprox(n1, n2, p11, p12, p21, p22, rho1, rho2, alpha, Test)
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

  - `"AN"`: Asymptotic normal method without continuity correction

  - `"ANc"`: Asymptotic normal method with continuity correction

  - `"AS"`: Arcsine method without continuity correction

  - `"ASc"`: Arcsine method with continuity correction

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

  Power for both co-primary endpoints

## Details

This function implements four approximate power calculation methods:

**Asymptotic Normal (AN):** Uses the standard normal approximation
without continuity correction (equations 3-4 in Sozu et al. 2010).

**Asymptotic Normal with Continuity Correction (ANc):** Includes Yates's
continuity correction (equation 5 in Sozu et al. 2010).

**Arcsine (AS):** Uses arcsine transformation without continuity
correction (equation 6 in Sozu et al. 2010).

**Arcsine with Continuity Correction (ASc):** Arcsine method with
continuity correction (equation 7 in Sozu et al. 2010).

The correlation between test statistics for the two endpoints depends on
the method:

**For AN and ANc methods:** \$\$\rho\_{nml} = \frac{\sum\_{j=1}^{2}
\rho_j \sqrt{\nu\_{j1}\nu\_{j2}}/n_j} {se_1 \times se_2}\$\$ where
\\\nu\_{jk} = p\_{jk}(1-p\_{jk})\\.

**For AS method:** \$\$\rho\_{arc} = \frac{n_2 \rho_1 + n_1
\rho_2}{n_1 + n_2}\$\$ This is the weighted average of the correlations
from both groups.

**For ASc method:** \$\$\rho\_{arc,c} = \frac{1}{se_1 \times se_2}
\left(\frac{\rho_1
\sqrt{\nu\_{11}\nu\_{12}}}{4n_1\sqrt{\nu\_{11,c}\nu\_{12,c}}} +
\frac{\rho_2
\sqrt{\nu\_{21}\nu\_{22}}}{4n_2\sqrt{\nu\_{21,c}\nu\_{22,c}}}\right)\$\$
where \\\nu\_{jk,c} = (p\_{jk} + c_j)(1 - p\_{jk} - c_j)\\, \\c_1 =
-1/(2n_1)\\, and \\c_2 = 1/(2n_2)\\.

The correlation bounds are automatically checked using
[`corrbound2Binary`](https://gosukehommaEX.github.io/twoCoprimary/reference/corrbound2Binary.md).

## References

Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination
in clinical trials with multiple co-primary binary endpoints.
*Statistics in Medicine*, 29(21), 2169-2179.

## Examples

``` r
# Power calculation using asymptotic normal method
power2BinaryApprox(
  n1 = 200,
  n2 = 100,
  p11 = 0.5,
  p12 = 0.4,
  p21 = 0.3,
  p22 = 0.2,
  rho1 = 0.7,
  rho2 = 0.7,
  alpha = 0.025,
  Test = 'AN'
)
#> 
#> Power calculation for two binary co-primary endpoints
#> 
#>              n1 = 200
#>              n2 = 100
#>     p (group 1) = 0.5, 0.4
#>     p (group 2) = 0.3, 0.2
#>             rho = 0.7, 0.7
#>           alpha = 0.025
#>            Test = AN
#>          power1 = 0.91929
#>          power2 = 0.949617
#>  powerCoprimary = 0.894946
#> 

# Power calculation with continuity correction
power2BinaryApprox(
  n1 = 200,
  n2 = 100,
  p11 = 0.5,
  p12 = 0.4,
  p21 = 0.3,
  p22 = 0.2,
  rho1 = 0.7,
  rho2 = 0.7,
  alpha = 0.025,
  Test = 'ANc'
)
#> 
#> Power calculation for two binary co-primary endpoints
#> 
#>              n1 = 200
#>              n2 = 100
#>     p (group 1) = 0.5, 0.4
#>     p (group 2) = 0.3, 0.2
#>             rho = 0.7, 0.7
#>           alpha = 0.025
#>            Test = ANc
#>          power1 = 0.898088
#>          power2 = 0.933117
#>  powerCoprimary = 0.867311
#> 

# Power calculation using arcsine method
power2BinaryApprox(
  n1 = 150,
  n2 = 150,
  p11 = 0.6,
  p12 = 0.5,
  p21 = 0.4,
  p22 = 0.3,
  rho1 = 0.5,
  rho2 = 0.5,
  alpha = 0.025,
  Test = 'AS'
)
#> 
#> Power calculation for two binary co-primary endpoints
#> 
#>              n1 = 150
#>              n2 = 150
#>     p (group 1) = 0.6, 0.5
#>     p (group 2) = 0.4, 0.3
#>             rho = 0.5, 0.5
#>           alpha = 0.025
#>            Test = AS
#>          power1 = 0.936701
#>          power2 = 0.945629
#>  powerCoprimary = 0.897574
#> 
```
