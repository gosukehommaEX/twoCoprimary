# Sample Size Calculation for a Single Binary Endpoint

Calculates the required sample size for a two-arm superiority trial with
a single binary endpoint using various statistical testing methods.

## Usage

``` r
ss1BinaryApprox(p1, p2, r, alpha, beta, Test = "AN")
```

## Arguments

- p1:

  True probability of responders in group 1 (0 \< p1 \< 1)

- p2:

  True probability of responders in group 2 (0 \< p2 \< 1)

- r:

  Allocation ratio of group 1 to group 2 (group 1:group 2 = r:1, where r
  \> 0)

- alpha:

  One-sided significance level (typically 0.025)

- beta:

  Target type II error rate (typically 0.1 or 0.2)

- Test:

  Statistical testing method. One of:

  - `"AN"`: Asymptotic normal method without continuity correction
    (default)

  - `"ANc"`: Asymptotic normal method with continuity correction

  - `"AS"`: Arcsine transformation without continuity correction

  - `"ASc"`: Arcsine transformation with continuity correction

  - `"Fisher"`: Fisher's exact test with iterative sample size
    determination

## Value

A data frame with the following columns:

- p1:

  Probability of responders in group 1

- p2:

  Probability of responders in group 2

- r:

  Allocation ratio

- alpha:

  One-sided significance level

- beta:

  Type II error rate

- Test:

  Testing method used

- n1:

  Required sample size for group 1

- n2:

  Required sample size for group 2

- N:

  Total sample size (n1 + n2)

## Details

This function implements sample size calculations for single binary
endpoint trials using five different methods.

**Important:** This function is designed for a **single binary
endpoint**. For co-primary endpoints, use
[`ss2BinaryApprox`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2BinaryApprox.md)
(for approximate methods) or
[`ss2BinaryExact`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2BinaryExact.md)
(for exact methods).

**Notation:**

- r = n1/n2: allocation ratio (group 1 to group 2)

- kappa = 1/r = n2/n1: inverse allocation ratio

- p1, p2: response probabilities

- theta1 = 1 - p1, theta2 = 1 - p2: non-response probabilities

- delta = p1 - p2: treatment effect

**AN (Asymptotic Normal) Method:** Uses the standard normal
approximation with pooled variance under H0: \$\$n_2 = \left\lceil
\frac{(1 + \kappa)}{(\pi_1 - \pi_2)^2} \left(z\_{1-\alpha}
\sqrt{\bar{\pi}(1-\bar{\pi})} + z\_{1-\beta} \sqrt{\kappa\pi_1\theta_1 +
\pi_2\theta_2}\right)^2 / \kappa \right\rceil\$\$ where \\\bar{\pi} =
(r\pi_1 + \pi_2)/(1 + r)\\ is the pooled proportion.

**ANc Method:** Adds continuity correction to the AN method. Uses
iterative calculation because the correction term depends on sample
size. Converges when the difference between successive iterations is
less than or equal to 1.

**AS (Arcsine) Method:** Uses the variance-stabilizing arcsine
transformation: \$\$n_2 = \left\lceil \frac{(z\_{1-\alpha} +
z\_{1-\beta})^2}{4(\sin^{-1}\sqrt{\pi_1} - \sin^{-1}\sqrt{\pi_2})^2}
\times \frac{1 + \kappa}{\kappa} \right\rceil\$\$

**ASc Method:** Applies continuity correction to the arcsine method.
Uses iterative procedure with convergence criterion.

**Fisher Method:** Fisher's exact test does not have a closed-form
sample size formula. This method:

1.  Starts with the AN method's sample size as initial value

2.  Incrementally increases n2 by 1

3.  Calculates exact power using hypergeometric distribution

4.  Stops when power is greater than or equal to 1 - beta

Note: Due to the saw-tooth nature of exact power (power does not
increase monotonically with sample size), a sequential search approach
is used. The incremental approach ensures the minimum sample size that
achieves the target power.

## References

Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination
in clinical trials with multiple co-primary binary endpoints.
*Statistics in Medicine*, 29(21), 2169-2179.

## Examples

``` r
# Balanced design with 1:1 allocation (AN method)
ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = 1, alpha = 0.025, beta = 0.1, Test = "AN")
#> 
#> Sample size calculation for single binary endpoint
#> 
#>              n1 = 130
#>              n2 = 130
#>               N = 260
#>               p = 0.6, 0.4
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#>            Test = AN
#> 

# Unbalanced design with 2:1 allocation (ANc method)
ss1BinaryApprox(p1 = 0.5, p2 = 0.3, r = 2, alpha = 0.025, beta = 0.2, Test = "ANc")
#> 
#> Sample size calculation for single binary endpoint
#> 
#>              n1 = 156
#>              n2 = 78
#>               N = 234
#>               p = 0.5, 0.3
#>      allocation = 2
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = ANc
#> 

# Arcsine transformation method
ss1BinaryApprox(p1 = 0.55, p2 = 0.35, r = 1, alpha = 0.025, beta = 0.1, Test = "AS")
#> 
#> Sample size calculation for single binary endpoint
#> 
#>              n1 = 129
#>              n2 = 129
#>               N = 258
#>               p = 0.55, 0.35
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#>            Test = AS
#> 

# Arcsine with continuity correction
ss1BinaryApprox(p1 = 0.65, p2 = 0.45, r = 1, alpha = 0.025, beta = 0.1, Test = "ASc")
#> 
#> Sample size calculation for single binary endpoint
#> 
#>              n1 = 121
#>              n2 = 121
#>               N = 242
#>               p = 0.65, 0.45
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.1
#>            Test = ASc
#> 

# Fisher's exact test
ss1BinaryApprox(p1 = 0.6, p2 = 0.4, r = 2, alpha = 0.025, beta = 0.1, Test = "Fisher")
#> 
#> Sample size calculation for single binary endpoint
#> 
#>              n1 = 206
#>              n2 = 103
#>               N = 309
#>               p = 0.6, 0.4
#>      allocation = 2
#>           alpha = 0.025
#>            beta = 0.1
#>            Test = Fisher
#> 
```
