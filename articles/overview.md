# Overview of Two Co-Primary Endpoints Analysis

## Introduction

The `twoCoprimary` package provides comprehensive tools for sample size
calculation and power analysis in clinical trials with two co-primary
endpoints. This package implements methodologies from multiple
publications and supports various endpoint types.

``` r
library(twoCoprimary)
```

## What are Co-Primary Endpoints?

In clinical trials, **co-primary endpoints** are multiple primary
endpoints that must **all** show statistically significant treatment
effects for the trial to be considered successful. This is in contrast
to multiple primary endpoints where demonstrating an effect on **any
one** endpoint is sufficient.

### Statistical Properties

- **No multiplicity adjustment needed** for Type I error control
- The overall Type I error rate is controlled at level $\alpha$ without
  Bonferroni correction
- The overall power is the joint probability: Power =
  $\Pr\left( {\text{Reject}\mspace{6mu}}\text{H}_{01}{\mspace{6mu}\text{and Reject}\mspace{6mu}}\text{H}_{02} \right)$
- Accounting for **correlation** between endpoints can improve
  efficiency

### Hypotheses Structure

For two co-primary endpoints, we test:

**Null hypothesis**: $\text{H}_{0} = \text{H}_{01} \cup \text{H}_{02}$
(at least one null hypothesis is true)

**Alternative hypothesis**:
$\text{H}_{1} = \text{H}_{11} \cap \text{H}_{12}$ (both alternative
hypotheses are true)

We reject $\text{H}_{0}$ only if **both** $\text{H}_{01}$ and
$\text{H}_{02}$ are rejected at level $\alpha$.

## Statistical Framework

### Intersection-Union Test (IUT)

The co-primary endpoint framework is based on the intersection-union
test principle. Let $Z_{1}$ and $Z_{2}$ be the test statistics for
endpoints 1 and 2, respectively.

**Decision rule**: Reject $\text{H}_{0}$ if and only if:
$$Z_{1} > z_{1 - \alpha}{\mspace{6mu}\text{and}\mspace{6mu}}Z_{2} > z_{1 - \alpha}$$

where $z_{1 - \alpha}$ is the $(1 - \alpha)$-th quantile of the standard
normal distribution.

### Type I Error Control

The overall Type I error rate is:

$$\alpha_{\text{overall}} = \Pr\left( {\text{Reject}\mspace{6mu}}\text{H}_{0} \mid \text{H}_{0}{\mspace{6mu}\text{true}} \right)$$

Under the intersection-union test, this is automatically controlled at
level $\alpha$ without adjustment:

$$\alpha_{\text{overall}} \leq \alpha$$

This is because rejecting when **both** statistics exceed the threshold
is more conservative than rejecting when **either** statistic exceeds
the threshold.

### Overall Power

Under the alternative hypothesis $\text{H}_{1}$, the overall power is:

$$1 - \beta = \Pr\left( Z_{1} > z_{1 - \alpha}{\mspace{6mu}\text{and}\mspace{6mu}}Z_{2} > z_{1 - \alpha} \mid \text{H}_{1} \right)$$

When $\left( Z_{1},Z_{2} \right)$ follow a bivariate normal distribution
with correlation $\rho$:

$$1 - \beta = \Phi_{2}\left( - z_{1 - \alpha} + \omega_{1}, - z_{1 - \alpha} + \omega_{2} \mid \rho \right)$$

where:

- $\Phi_{2}( \cdot , \cdot \mid \rho)$ is the bivariate normal
  cumulative distribution function with correlation $\rho$
- $\omega_{1}$ and $\omega_{2}$ are the non-centrality parameters under
  $\text{H}_{1}$

### Impact of Correlation

The correlation $\rho$ between test statistics affects the overall
power:

- **Positive correlation** ($\rho > 0$): Increases power and reduces
  required sample size
- **Zero correlation** ($\rho = 0$): Test statistics are independent
- **Negative correlation** ($\rho < 0$): Decreases power and increases
  required sample size

**Key insight**: Accounting for positive correlation between endpoints
can lead to substantial sample size reductions (typically 5-15% for
$\rho = 0.5$-0.8) compared to assuming independence.

## Supported Endpoint Types

The package supports five combinations of co-primary endpoints:

### 1. Two Continuous Endpoints

**Use case**: Trials measuring two continuous outcomes (e.g., systolic
and diastolic blood pressure)

**Statistical model**: Both endpoints follow normal distributions

**Key functions**:

- [`ss2Continuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2Continuous.md):
  Sample size calculation
- [`power2Continuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2Continuous.md):
  Power calculation
- [`twoCoprimary2Continuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2Continuous.md):
  Unified interface for both

**Reference**: Sozu et al. (2011)

``` r
# Example: Two continuous endpoints with correlation rho = 0.5
ss2Continuous(
  delta1 = 0.5,  # Standardized effect size for endpoint 1
  delta2 = 0.5,  # Standardized effect size for endpoint 2
  sd1 = 1,       # Standard deviation for endpoint 1
  sd2 = 1,       # Standard deviation for endpoint 2
  rho = 0.5,     # Correlation between endpoints
  r = 1,         # Balanced allocation
  alpha = 0.025,
  beta = 0.2,
  known_var = TRUE
)
#> 
#> Sample size calculation for two continuous co-primary endpoints
#> 
#>              n1 = 79
#>              n2 = 79
#>               N = 158
#>           delta = 0.5, 0.5
#>              sd = 1, 1
#>             rho = 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>       known_var = TRUE
```

### 2. Two Binary Endpoints (Asymptotic Approximation)

**Use case**: Trials with two binary outcomes using normal approximation
(large sample)

**Statistical model**: Binary endpoints with asymptotic normal
approximation

**Key functions**:

- [`ss2BinaryApprox()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2BinaryApprox.md):
  Sample size calculation
- [`power2BinaryApprox()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2BinaryApprox.md):
  Power calculation
- [`twoCoprimary2BinaryApprox()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2BinaryApprox.md):
  Unified interface for both

**Reference**: Sozu et al. (2010)

**Supported test methods**:

- AN: Asymptotic normal test without continuity correction
- ANc: Asymptotic normal test with continuity correction
- AS: Arcsine transformation without continuity correction
- ASc: Arcsine transformation with continuity correction

**When to use**: Large sample sizes (typically $N > 200$) and
probabilities not too extreme ($0.1 < p < 0.9$)

``` r
# Example: Two binary endpoints
ss2BinaryApprox(
  p11 = 0.7, p12 = 0.6,  # Endpoint 1 and 2 in treatment group
  p21 = 0.4, p22 = 0.3,  # Endpoint 1 and 2 in control group
  rho1 = 0.5,            # Correlation in treatment group
  rho2 = 0.5,            # Correlation in control group
  r = 1,                 # Balanced allocation
  alpha = 0.025,
  beta = 0.2,
  Test = "AN"
)
#> 
#> Sample size calculation for two binary co-primary endpoints
#> 
#>              n1 = 52
#>              n2 = 52
#>               N = 104
#>     p (group 1) = 0.7, 0.6
#>     p (group 2) = 0.4, 0.3
#>             rho = 0.5, 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = AN
```

### 3. Two Binary Endpoints (Exact Methods)

**Use case**: Small to medium sample sizes requiring exact inference

**Statistical model**: Binary endpoints with exact tests

**Key functions**:

- [`ss2BinaryExact()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2BinaryExact.md):
  Sample size calculation using exact tests
- [`power2BinaryExact()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2BinaryExact.md):
  Exact power calculation
- [`twoCoprimary2BinaryExact()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2BinaryExact.md):
  Unified interface for both

**Reference**: Homma and Yoshida (2025)

**Supported tests**:

- Chisq: Chi-squared test
- Fisher: Fisher’s exact test (conditional test)
- Fisher-midP: Fisher’s mid-p test
- Z-pool: Z-pooled exact unconditional test
- Boschloo: Boschloo’s exact unconditional test

**When to use**: Small/medium samples ($N < 200$), extreme probabilities
($p < 0.1$ or $p > 0.9$), or when strict Type I error control is
required

``` r
# Example: Exact methods for small samples
ss2BinaryExact(
  p11 = 0.7, p12 = 0.6,
  p21 = 0.4, p22 = 0.3,
  rho1 = 0.5, rho2 = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  Test = "Fisher"  # or "Chisq", "Fisher-midP", "Z-pool", "Boschloo"
)
#> 
#> Sample size calculation for two binary co-primary endpoints
#> 
#>              n1 = 59
#>              n2 = 59
#>               N = 118
#>     p (group 1) = 0.7, 0.6
#>     p (group 2) = 0.4, 0.3
#>             rho = 0.5, 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = Fisher
```

### 4. Mixed Continuous and Binary Endpoints

**Use case**: Trials with one continuous and one binary outcome

**Statistical model**: Normal distribution for continuous endpoint,
Bernoulli for binary endpoint, with biserial correlation

**Key functions**:

- [`ss2MixedContinuousBinary()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2MixedContinuousBinary.md):
  Sample size calculation
- [`power2MixedContinuousBinary()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2MixedContinuousBinary.md):
  Power calculation
- [`twoCoprimary2MixedContinuousBinary()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2MixedContinuousBinary.md):
  Unified interface for both

**Reference**: Sozu et al. (2012)

**Supported test methods for binary endpoint**:

- AN: Asymptotic normal test without continuity correction
- ANc: Asymptotic normal test with continuity correction
- AS: Arcsine transformation without continuity correction
- ASc: Arcsine transformation with continuity correction
- Fisher: Fisher’s exact test (simulation-based)

**Correlation structure**: Uses biserial correlation between observed
continuous variable and latent continuous variable underlying the binary
outcome

``` r
# Example: Continuous + Binary endpoints
ss2MixedContinuousBinary(
  delta = 0.5,           # Effect size for continuous endpoint
  sd = 1,                # Standard deviation
  p1 = 0.7,              # Success probability in treatment group
  p2 = 0.4,              # Success probability in control group
  rho = 0.5,             # Biserial correlation
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  Test = "AN"
)
#> 
#> Sample size calculation for mixed continuous and binary co-primary endpoints
#> 
#>              n1 = 68
#>              n2 = 68
#>               N = 136
#>           delta = 0.5
#>              sd = 1
#>               p = 0.7, 0.4
#>             rho = 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = AN
```

### 5. Mixed Count and Continuous Endpoints

**Use case**: Trials with overdispersed count data (e.g., exacerbations)
and continuous outcomes (e.g., lung function)

**Statistical model**: Negative binomial distribution for count
endpoint, normal distribution for continuous endpoint

**Key functions**:

- [`ss2MixedCountContinuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2MixedCountContinuous.md):
  Sample size calculation
- [`power2MixedCountContinuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2MixedCountContinuous.md):
  Power calculation
- [`twoCoprimary2MixedCountContinuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2MixedCountContinuous.md):
  Unified interface for both
- [`corrbound2MixedCountContinuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/corrbound2MixedCountContinuous.md):
  Calculate valid correlation bounds

**Reference**: Homma and Yoshida (2024)

**Special considerations**:

- The negative binomial distribution accommodates overdispersion
  (variance $>$ mean) common in count data
- **Treatment effects must be in the negative direction for both
  endpoints**: Lower event rate (count endpoint) and lower/better
  continuous values indicate treatment benefit. For example, reduction
  in exacerbation rate (`r1 < r2`) and improvement in lung function
  (e.g., `mu1 < mu2` when lower is better, or larger negative change
  from baseline)

``` r
# Example: Count (exacerbations) + Continuous (FEV1)
ss2MixedCountContinuous(
  r1 = 1.0, r2 = 1.25,   # Count rates (events per unit time)
  nu = 0.8,              # Dispersion parameter
  t = 1,                 # Follow-up time
  mu1 = -50, mu2 = 0,    # Continuous means
  sd = 250,              # Standard deviation
  rho1 = 0.5, rho2 = 0.5, # Correlations
  r = 1,
  alpha = 0.025,
  beta = 0.2
)
#> 
#> Sample size calculation for mixed count and continuous co-primary endpoints
#> 
#>              n1 = 705
#>              n2 = 705
#>               N = 1410
#>              sd = 250
#>            rate = 1, 1.25
#>              nu = 0.8
#>               t = 1
#>              mu = -50, 0
#>             rho = 0.5, 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
```

## Impact of Correlation

A key advantage of accounting for correlation between co-primary
endpoints is the potential for **sample size reduction**.

The table below illustrates this for two continuous endpoints:

``` r
# Sample size at different correlation levels
correlations <- c(0, 0.3, 0.5, 0.8)
results <- sapply(correlations, function(rho) {
  ss2Continuous(
    delta1 = 0.5, delta2 = 0.5,
    sd1 = 1, sd2 = 1,
    rho = rho, r = 1,
    alpha = 0.025, beta = 0.2,
    known_var = TRUE
  )$N
})

data.frame(
  Correlation = correlations,
  Total_N = results,
  Reduction = paste0(round((1 - results/results[1]) * 100, 1), "%")
)
#>   Correlation Total_N Reduction
#> 1         0.0     166        0%
#> 2         0.3     162      2.4%
#> 3         0.5     158      4.8%
#> 4         0.8     148     10.8%
```

As correlation increases, the required sample size decreases. At
$\rho = 0.8$, approximately 11% reduction in sample size can be achieved
compared to $\rho = 0$.

### Why Does Correlation Matter?

The correlation between endpoints affects the joint distribution of test
statistics. When endpoints are positively correlated:

1.  **Test statistics tend to move together**: If $Z_{1}$ is large,
    $Z_{2}$ is also likely to be large
2.  **Higher probability of rejecting both nulls**:
    $\Pr\left( Z_{1} > c,Z_{2} > c \right)$ increases with $\rho$
3.  **Sample size reduction**: Fewer subjects needed to achieve target
    power

Mathematically, for bivariate normal $\left( Z_{1},Z_{2} \right)$ with
correlation $\rho$:

$$\Pr\left( Z_{1} > c,Z_{2} > c \mid \rho \right) > \Pr\left( Z_{1} > c,Z_{2} > c \mid \rho = 0 \right)$$

when $\rho > 0$ and both endpoints have positive treatment effects.

## Choosing the Right Method

| Endpoint Types          | Sample Size       | Method     | Key Considerations                           |
|-------------------------|-------------------|------------|----------------------------------------------|
| Both continuous         | Any               | Asymptotic | Simple, well-established                     |
| Both binary             | Large ($N > 200$) | Asymptotic | Fast computation, probabilities moderate     |
| Both binary             | Small/Medium      | Exact      | Better Type I error control, exact inference |
| 1 Continuous + 1 Binary | Any               | Asymptotic | Handles mixed types, biserial correlation    |
| 1 Count + 1 Continuous  | Any               | Asymptotic | Accounts for overdispersion                  |

### Decision Guidelines

**For binary endpoints**:

- Use **asymptotic methods** when: $N > 200$, $0.1 < p < 0.9$,
  computational efficiency important
- Use **exact methods** when: $N < 200$, extreme probabilities
  ($p < 0.1$ or $p > 0.9$), regulatory requirements for exact tests

**For mixed endpoint types**:

- Ensure correlation structure is appropriate (e.g., biserial for
  continuous-binary)
- Consider clinical plausibility of correlation magnitude
- Use conservative estimates if correlation is uncertain

## Sample Size Calculation Approach

All methods in this package follow a similar computational approach:

1.  **Specify design parameters**: Effect sizes, probabilities, standard
    deviations, etc.
2.  **Specify correlation**: Between endpoints
3.  **Specify error rates**: Type I error $\alpha$ (typically 0.025 for
    one-sided) and Type II error $\beta$ (typically 0.2 for 80% power)
4.  **Calculate sample size**: Using iterative algorithms
5.  **Verify power**: Confirm that calculated sample size achieves
    target power

## Detailed Vignettes

For detailed methodology, examples, and validation against published
results, please see:

1.  [`vignette("two-continuous-endpoints")`](https://gosukehommaEX.github.io/twoCoprimary/articles/two-continuous-endpoints.md):
    Two continuous endpoints
2.  [`vignette("two-binary-endpoints-approx")`](https://gosukehommaEX.github.io/twoCoprimary/articles/two-binary-endpoints-approx.md):
    Two binary endpoints (asymptotic)
3.  [`vignette("two-binary-endpoints-exact")`](https://gosukehommaEX.github.io/twoCoprimary/articles/two-binary-endpoints-exact.md):
    Two binary endpoints (exact)
4.  [`vignette("mixed-continuous-binary")`](https://gosukehommaEX.github.io/twoCoprimary/articles/mixed-continuous-binary.md):
    Mixed continuous and binary
5.  [`vignette("mixed-count-continuous")`](https://gosukehommaEX.github.io/twoCoprimary/articles/mixed-count-continuous.md):
    Mixed count and continuous

## References

Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical
trials with two co-primary endpoints including overdispersed count and
continuous outcomes. *Pharmaceutical Statistics*, 23(1), 46-59.

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(1), 1-19.

Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination
in clinical trials with multiple co-primary binary endpoints.
*Statistics in Medicine*, 29(21), 2169-2179.

Sozu, T., Sugimoto, T., & Hamasaki, T. (2011). Sample size determination
in superiority clinical trials with multiple co-primary correlated
endpoints. *Journal of Biopharmaceutical Statistics*, 21(4), 650-668.

Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination
in clinical trials with multiple co-primary endpoints including mixed
continuous and binary variables. *Biometrical Journal*, 54(5), 716-729.
