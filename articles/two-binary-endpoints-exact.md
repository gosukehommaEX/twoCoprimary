# Two Binary Co-Primary Endpoints (Exact Methods)

## Overview

This vignette demonstrates exact sample size calculation and power
analysis for clinical trials with two co-primary binary endpoints. The
methodology is based on Homma and Yoshida (2025), which provides exact
inference methods using the bivariate binomial distribution.

``` r

library(twoCoprimary)
library(dplyr)
library(tidyr)
library(knitr)
```

## Background

### When to Use Exact Methods

Exact methods are recommended when:

- **Small to medium sample sizes** ($`N < 200`$)
- **Extreme probabilities** ($`p < 0.10`$ or $`p > 0.90`$)
- **Strict Type I error control** is required
- **Regulatory requirements** for exact inference

Asymptotic methods may not maintain the nominal Type I error rate in
these situations.

### Advantages of Exact Methods

1.  **Accurate Type I error control**: Exact tests guarantee
    $`\alpha \leq`$ nominal level
2.  **Better small-sample performance**: No reliance on asymptotic
    approximations
3.  **Valid for extreme probabilities**: No restrictions on $`p`$ values
4.  **Regulatory acceptance**: Often preferred by regulatory agencies

### Disadvantages

1.  **Computational intensity**: Requires enumeration of possible
    outcomes
2.  **Conservatism**: Discrete nature can lead to conservatism
3.  **Implementation complexity**: More complex than asymptotic methods

## Statistical Framework

### Model and Assumptions

Consider a two-arm parallel-group superiority trial comparing treatment
(group 1) with control (group 2). Let $`n_{1}`$ and $`n_{2}`$ denote the
sample sizes in groups 1 and 2, respectively.

For patient $`i`$ in group $`j`$ ($`j = 1`$: treatment, $`j = 2`$:
control), we observe two binary outcomes:

**Endpoint $`k`$** ($`k = 1, 2`$):
``` math
X_{i,j,k} \in \{0, 1\}
```

where $`X_{i,j,k} = 1`$ if patient $`i`$ in group $`j`$ is a responder
for endpoint $`k`$, and 0 otherwise.

**True response probabilities**:
``` math
p_{j,k} = \text{P}(X_{i,j,k} = 1)
```

where $`0 < p_{j,k} < 1`$ for each $`j`$ and $`k`$.

### Joint Distribution of Binary Outcomes

The paired binary outcomes $`(X_{i,j,1}, X_{i,j,2})`$ for patient $`i`$
in group $`j`$ follow a multinomial distribution with four possible
outcomes:

**Per-trial probabilities**:

- $`p_{j}^{(1,1)} = \phi_{j}`$: Both endpoints successful
- $`p_{j}^{(1,0)} = p_{j,1} - \phi_{j}`$: Only endpoint 1 successful
- $`p_{j}^{(0,1)} = p_{j,2} - \phi_{j}`$: Only endpoint 2 successful
- $`p_{j}^{(0,0)} = 1 - p_{j,1} - p_{j,2} + \phi_{j}`$: Both endpoints
  unsuccessful

where $`\phi_{j} = \text{P}(X_{i,j,1} = 1, X_{i,j,2} = 1)`$.

Let $`Z_{j}^{(\ell,m)}`$ denote the random variable representing the
number of times $`\{(X_{i,j,1}, X_{i,j,2}) : i = 1, \ldots, n_{j}\}`$
takes the value $`(\ell, m)`$ for $`\ell, m \in \{0, 1\}`$. Then:

``` math
(Z_{j}^{(0,0)}, Z_{j}^{(1,0)}, Z_{j}^{(0,1)}, Z_{j}^{(1,1)}) \sim \text{Multinomial}(n_{j}; p_{j}^{(0,0)}, p_{j}^{(1,0)}, p_{j}^{(0,1)}, p_{j}^{(1,1)})
```

### Number of Responders

Let $`Y_{j,k} = \sum_{i=1}^{n_{j}} X_{i,j,k}`$ represent the number of
responders in group $`j`$ for endpoint $`k`$. Then:

- $`Y_{j,1} = Z_{j}^{(1,1)} + Z_{j}^{(1,0)}`$
- $`Y_{j,2} = Z_{j}^{(1,1)} + Z_{j}^{(0,1)}`$

### Bivariate Binomial Distribution

Following Homma and Yoshida (2025), the joint distribution of
$`(Y_{j,1}, Y_{j,2})`$ can be expressed as a **bivariate binomial
distribution**:

``` math
(Y_{j,1}, Y_{j,2}) \sim \text{BiBin}(n_{j}, p_{j,1}, p_{j,2}, \gamma_{j})
```

where $`\gamma_{j}`$ is a dependence parameter related to the
correlation $`\rho_{j}`$ between $`X_{i,j,1}`$ and $`X_{i,j,2}`$.

**Probability mass function** (Equation 3 in Homma and Yoshida, 2025):

``` math
\text{P}(Y_{j,1} = y_{j,1}, Y_{j,2} = y_{j,2} \mid n_{j}, p_{j,1}, p_{j,2}, \gamma_{j}) = f(y_{j,1} \mid n_{j}, p_{j,1}) \times g(y_{j,2} \mid y_{j,1}, n_{j}, p_{j,1}, p_{j,2}, \gamma_{j})
```
For more details, please see Homma and Yoshida (2025).

### Correlation Structure

The **correlation** $`\rho_{j}`$ between $`X_{i,j,1}`$ and $`X_{i,j,2}`$
is:

``` math
\rho_{j} = \text{Cor}(X_{i,j,1}, X_{i,j,2}) = \frac{\phi_{j} - p_{j,1} p_{j,2}}{\sqrt{p_{j,1}(1 - p_{j,1}) p_{j,2}(1 - p_{j,2})}}
```

The dependence parameter $`\gamma_{j}`$ is related to $`\rho_{j}`$
through (Equation 4 in Homma and Yoshida, 2025):

``` math
\gamma_{j} = \gamma(\rho_{j}, p_{j,1}, p_{j,2}) = \rho_{j} \sqrt{\frac{p_{j,2}(1 - p_{j,2})}{p_{j,1}(1 - p_{j,1})}} \left(1 - \rho_{j} \sqrt{\frac{p_{j,2}(1 - p_{j,2})}{p_{j,1}(1 - p_{j,1})}}\right)^{-1}
```

**Important property**: The correlation between $`Y_{j,1}`$ and
$`Y_{j,2}`$ equals $`\rho_{j}`$, the same as the correlation between
$`X_{i,j,1}`$ and $`X_{i,j,2}`$.

**Marginal distributions**:
``` math
Y_{j,k} \sim \text{Bin}(n_{j}, p_{j,k})
```

**Correlation bounds**: Due to $`0 < p_{j,k} < 1`$, the correlation
$`\rho_{j}`$ is bounded:

``` math
\rho_{j} \in [L(p_{j,1}, p_{j,2}), U(p_{j,1}, p_{j,2})] \subseteq [-1, 1]
```

where:

``` math
L(p_{j,1}, p_{j,2}) = \max\left\{-\sqrt{\frac{p_{j,1} p_{j,2}}{(1 - p_{j,1})(1 - p_{j,2})}}, -\sqrt{\frac{(1 - p_{j,1})(1 - p_{j,2})}{p_{j,1} p_{j,2}}}\right\}
```

``` math
U(p_{j,1}, p_{j,2}) = \min\left\{\sqrt{\frac{p_{j,1}(1 - p_{j,2})}{p_{j,2}(1 - p_{j,1})}}, \sqrt{\frac{p_{j,2}(1 - p_{j,1})}{p_{j,1}(1 - p_{j,2})}}\right\}
```

**Special cases**:

- If $`p_{j,1} = p_{j,2}`$, then $`U(p_{j,1}, p_{j,2}) = 1`$
- If $`p_{j,1} + p_{j,2} = 1`$, then $`L(p_{j,1}, p_{j,2}) = -1`$

## Hypothesis Testing

### Superiority Hypotheses

Since higher values of both endpoints indicate treatment benefit, we
test:

**For endpoint 1**:
``` math
\text{H}_{0}^{(1)}: p_{1,1} \leq p_{2,1} \text{ vs. } \text{H}_{1}^{(1)}: p_{1,1} > p_{2,1}
```

**For endpoint 2**:
``` math
\text{H}_{0}^{(2)}: p_{1,2} \leq p_{2,2} \text{ vs. } \text{H}_{1}^{(2)}: p_{1,2} > p_{2,2}
```

### Co-Primary Endpoints (Intersection-Union Test)

The trial succeeds only if superiority is demonstrated for **both**
endpoints simultaneously:

**Null hypothesis**:
$`\text{H}_{0} = \text{H}_{0}^{(1)} \cup \text{H}_{0}^{(2)}`$ (at least
one null is true)

**Alternative hypothesis**:
$`\text{H}_{1} = \text{H}_{1}^{(1)} \cap \text{H}_{1}^{(2)}`$ (both
alternatives are true)

**Decision rule**: Reject $`\text{H}_{0}`$ at level $`\alpha`$ if and
only if **both** $`\text{H}_{0}^{(1)}`$ and $`\text{H}_{0}^{(2)}`$ are
rejected at level $`\alpha`$ without multiplicity adjustment.

## Statistical Tests

Homma and Yoshida (2025) consider five exact test methods:

### Method 1: One-sided Pearson Chi-squared Test (Chisq)

For endpoint $`k`$, the test statistic is:

``` math
Z(y_{1,k}, y_{2,k}) = \frac{\hat{p}_{1,k} - \hat{p}_{2,k}}{\sqrt{\hat{p}_{k}(1 - \hat{p}_{k})\left(\frac{1}{n_{1}} + \frac{1}{n_{2}}\right)}}
```

where:

- $`\hat{p}_{j,k} = y_{j,k} / n_{j}`$ is the sample proportion
- $`\hat{p}_{k} = \frac{n_{1} \hat{p}_{1,k} + n_{2} \hat{p}_{2,k}}{n_{1} + n_{2}}`$
  is the pooled proportion

Reject $`\text{H}_{0}^{(k)}`$ if $`Z(y_{1,k}, y_{2,k}) > z_{1-\alpha}`$,
where $`z_{1-\alpha}`$ is the $`(1-\alpha)`$-quantile of the standard
normal distribution.

### Method 2: Fisher’s Exact Test (Fisher)

**Conditional test**: Conditions on the total number of successes
$`y_{1,k} + y_{2,k}`$.

Under $`\text{H}_{0}^{(k)}`$, $`Y_{1,k}`$ follows a hypergeometric
distribution given $`Y_{1,k} + Y_{2,k} = y_{k}`$.

**One-sided p-value**:

``` math
p_{k}^{\text{Fisher}} = \sum_{y=y_{1,k}}^{\min(n_{1}, y_{k})} \frac{\binom{n_{1}}{y} \binom{n_{2}}{y_{k} - y}}{\binom{n_{1} + n_{2}}{y_{k}}}
```

Reject $`\text{H}_{0}^{(k)}`$ if $`p_{k}^{\text{Fisher}} < \alpha`$.

### Method 3: Fisher’s Mid-P Test (Fisher-midP)

Reduces conservatism by adding half the probability of the observed
outcome:

``` math
p_{k}^{\text{mid-p}} = p_{k}^{\text{Fisher}} - \frac{1}{2} \times \frac{\binom{n_{1}}{y_{1,k}} \binom{n_{2}}{y_{k} - y_{1,k}}}{\binom{n_{1} + n_{2}}{y_{k}}}
```
Note: The `twoCoprimary` package can implement the Fisher’s Mid-P Test,
but Homma and Yoshida (2025) has not investigated this test.

### Method 4: Z-pooled Exact Unconditional Test (Z-pool)

**Unconditional test**: the $`p`$-value is the null probability of the
outcomes at least as extreme as the observed one, maximized over the
nuisance parameter $`p_{k}`$, the common success probability under
$`\text{H}_{0}`$.

The ordering statistic is the pooled $`Z`$ statistic. Writing
$`T(y_{1}, y_{2})`$ for that statistic and $`t`$ for its observed value,

``` math
p_{k} = \max_{0 \leq \theta \leq 1} \; \Pr\left( T(Y_{1}, Y_{2}) \geq t \mid \theta \right),
```

where $`Y_{1}`$ and $`Y_{2}`$ are independent binomial variables with
common probability $`\theta`$.

### Method 5: Boschloo’s Exact Unconditional Test (Boschloo)

Same construction with Fisher’s exact $`p`$-value as the ordering
statistic, so the tail event is the set of outcomes whose Fisher
$`p`$-value is no larger than the observed one.

**Most powerful** of the exact unconditional tests, but the most
demanding to compute.

### Tied outcomes

The tail event is defined by the inequality $`T \geq t`$, so every
outcome sharing the observed value of the ordering statistic belongs to
it. Outcomes with the same value of $`T`$ therefore receive the same
$`p`$-value. Accumulating the null probabilities along an ordering of
the outcomes without grouping the ties would give the earlier members of
a tie group a smaller $`p`$-value than the later ones, and which member
comes first would depend on an arbitrary sort order. Ties are common: at
$`n_{1} = n_{2} = 40`$ roughly half of the outcomes with a positive
$`Z`$ statistic share their value with another outcome.

### The nuisance parameter grid

The maximization over $`\theta`$ is carried out on a finite grid of
equally spaced values on $`[0, 1]`$, not analytically. The number of
grid points is the `n_grid` argument of
[`rr1Binary()`](https://gosukehommaex.github.io/twoCoprimary/reference/rr1Binary.md),
[`power2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/power2BinaryExact.md),
[`ss2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/ss2BinaryExact.md)
and
[`twoCoprimary2BinaryExact()`](https://gosukehommaex.github.io/twoCoprimary/reference/twoCoprimary2BinaryExact.md),
and its default is 100.

A coarse grid can only understate the maximum, so it can only make the
$`p`$-value smaller and the rejection region larger. In practice the
default is ample: for $`n_{1} = n_{2} = 20`$ and $`40`$ at
$`\alpha = 0.025`$, the rejection regions of both unconditional tests
are identical for every grid size from 25 to 4000 points. Cost grows
roughly in proportion to `n_grid`, so a finer grid is inexpensive to try
when a design sits close to a decision boundary.

``` r

# The default reproduces a much finer grid
identical(
  rr1Binary(n1 = 30, n2 = 30, alpha = 0.025, Test = "Boschloo"),
  rr1Binary(n1 = 30, n2 = 30, alpha = 0.025, Test = "Boschloo", n_grid = 1000)
)
#> [1] TRUE
```

## Exact Power Calculation

### Power Formula

The exact power for test method $`A`$ is (Equation 9 in Homma and
Yoshida, 2025):

``` math
\text{power}_{A}(\boldsymbol{\theta}) = \text{P}\left[\bigcap_{k=1}^{2} \{p_{A}(y_{1,k}, y_{2,k}) < \alpha\} \mid \text{H}_{1}\right]
```

``` math
= \sum_{(a_{1,1}, a_{2,1}) \in \mathcal{A}_{1}} \sum_{(a_{1,2}, a_{2,2}) \in \mathcal{A}_{2}} f(a_{1,1} \mid n_{1}, p_{1,1}) \times f(a_{2,1} \mid n_{2}, p_{2,1}) \times g(a_{1,2} \mid a_{1,1}, n_{1}, p_{1,1}, p_{1,2}, \gamma_{1}) \times g(a_{2,2} \mid a_{2,1}, n_{2}, p_{2,1}, p_{2,2}, \gamma_{2})
```

where:

- $`\boldsymbol{\theta} = (p_{1,1}, p_{2,1}, p_{1,2}, p_{2,2}, n_{1}, n_{2}, \gamma_{1}, \gamma_{2})`$
  is the parameter vector
- $`\mathcal{A}_{k}`$ is the rejection region for endpoint $`k`$
- $`\mathcal{A}_{k} = \{(y_{1,k}, y_{2,k}) : p_{A}(y_{1,k}, y_{2,k}) < \alpha\}`$

### Sample Size Calculation

The required sample size $`n_{2}`$ to achieve target power $`1 - \beta`$
is (Equation 10 in Homma and Yoshida, 2025):

``` math
n_{2} = \arg\min_{n_{2} \in \mathbb{Z}} \{\text{power}_{A}(\boldsymbol{\theta}) \geq 1 - \beta\}
```

This cannot be expressed as a closed-form formula due to:

1.  Discreteness of binary outcomes
2.  Non-monotonic “saw-tooth” power curve

**Algorithm**: Sequential search starting from asymptotic normal
approximation (AN method) as initial value.

## Replicating Homma and Yoshida (2025) Table 4

Table 4 from Homma and Yoshida (2025) shows sample sizes for various
correlations using the Chisq, Fisher, Z-pool, and Boschloo. Note that
the following sample code compute only scenario for $`\alpha=0.025`$.

The notation used in the function is: `p11` = $`p_{1,1}`$, `p12` =
$`p_{1,2}`$, `p21` = $`p_{2,1}`$, `p22` = $`p_{2,2}`$, where the first
subscript denotes the group (1 = treatment, 2 = control) and the second
subscript denotes the endpoint (1 or 2).

``` r

# Recreate Homma and Yoshida (2025) Table 4
library(dplyr)
library(tidyr)
library(readr)

param_grid_bin_exact_ss <- tibble(
  p11 = 0.54, 
  p12 = 0.54,
  p21 = 0.25,
  p22 = 0.25
)

result_bin_exact_ss <- do.call(
  bind_rows,
  lapply(c("Chisq", "Fisher", "Z-pool", "Boschloo"), function(test) {
    do.call(
      bind_rows,
      lapply(1:2, function(r) {
        design_table(
          param_grid = param_grid_bin_exact_ss,
          rho_values = c(0, 0.3, 0.5, 0.8),
          r = r,
          alpha = 0.025,
          beta = 0.1,
          endpoint_type = "binary",
          Test = test
        ) %>% 
          mutate(alpha = 0.025, r = r, Test = test)
      })
    )
  })
) %>% 
  pivot_longer(
    cols = starts_with("rho_"),
    names_to = "rho",
    values_to = "N",
    names_transform = list(rho = parse_number)
  ) %>% 
  select(r, rho, Test, N) %>% 
  pivot_wider(names_from = Test,  values_from = N) %>% 
  as.data.frame()

kable(result_bin_exact_ss,
      caption = "Table 4: Total Sample Size (N) for Two Co-Primary Binary Endpoints (α = 0.025, 1-β = 0.90)^a,b^",
      digits = 1,
      col.names = c("r", "ρ", "Chisq", "Fisher", "Z-pool", "Boschloo"))
```

|   r |   ρ | Chisq | Fisher | Z-pool | Boschloo |
|----:|----:|------:|-------:|-------:|---------:|
|   1 | 0.0 |   142 |    152 |    144 |      144 |
|   1 | 0.3 |   142 |    150 |    142 |      142 |
|   1 | 0.5 |   140 |    150 |    140 |      140 |
|   1 | 0.8 |   128 |    144 |    134 |      134 |
|   2 | 0.0 |   162 |    174 |    180 |      162 |
|   2 | 0.3 |   159 |    174 |    180 |      159 |
|   2 | 0.5 |   156 |    171 |    177 |      156 |
|   2 | 0.8 |   147 |    159 |    168 |      150 |

Table 4: Total Sample Size (N) for Two Co-Primary Binary Endpoints (α =
0.025, 1-β = 0.90)^(a,b) {.table}

^(a) Chisq denotes the one-sided Pearson chi-squared test. Fisher stands
for Fisher’s exact test. Z-pool represents the Z-pooled exact
unconditional test. Boschloo signifies Boschloo’s exact unconditional
test.

^(b) The required sample sizes were obtained by assuming that
$`p_{1,1} = p_{1,2} = 0.54`$ and $`p_{2,1} = p_{2,2} = 0.25`$.

## Practical Examples

### Example 1: Basic Exact Power Calculation

``` r

# Calculate exact power using Fisher's exact test
result_fisher <- power2BinaryExact(
  n1 = 50,
  n2 = 50,
  p11 = 0.70, p12 = 0.65,
  p21 = 0.50, p22 = 0.45,
  rho1 = 0.5, rho2 = 0.5,
  alpha = 0.025,
  Test = "Fisher"
)

print(result_fisher)
#> 
#> Power calculation for two binary co-primary endpoints
#> 
#>              n1 = 50
#>              n2 = 50
#>     p (group 1) = 0.7, 0.65
#>     p (group 2) = 0.5, 0.45
#>             rho = 0.5, 0.5
#>           alpha = 0.025
#>            Test = Fisher
#>          power1 = 0.46345
#>          power2 = 0.46196
#>  powerCoprimary = 0.297231
```

**Interpretation**:

- `power1`: Power for endpoint 1 alone
- `power2`: Power for endpoint 2 alone
- `powerCoprimary`: Exact power for both co-primary endpoints

### Example 2: Sample Size Calculation

``` r

# Calculate required sample size using Boschloo's test
result_ss <- ss2BinaryExact(
  p11 = 0.70, p12 = 0.65,
  p21 = 0.50, p22 = 0.45,
  rho1 = 0.5, rho2 = 0.5,
  r = 1,
  alpha = 0.025,
  beta = 0.2,
  Test = "Boschloo"
)

print(result_ss)
#> 
#> Sample size calculation for two binary co-primary endpoints
#> 
#>              n1 = 120
#>              n2 = 120
#>               N = 240
#>     p (group 1) = 0.7, 0.65
#>     p (group 2) = 0.5, 0.45
#>             rho = 0.5, 0.5
#>      allocation = 1
#>           alpha = 0.025
#>            beta = 0.2
#>            Test = Boschloo
```

### Example 3: Comparison of Test Methods

``` r

# Compare different exact test methods
test_methods <- c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")

comparison <- lapply(test_methods, function(test) {
    result <- ss2BinaryExact(
        p11 = 0.50, p12 = 0.40,
        p21 = 0.20, p22 = 0.10,
        rho1 = 0.7, rho2 = 0.6,
        r = 1,
        alpha = 0.025,
        beta = 0.2,
        Test = test
    )
    data.frame(
        Test = test,
        n2 = result$n2,
        N = result$N
    )
})

comparison_table <- bind_rows(comparison)

kable(comparison_table,
      caption = "Sample Size Comparison Across Test Methods",
      col.names = c("Test Method", "n per group", "N total"))
```

| Test Method | n per group | N total |
|:------------|------------:|--------:|
| Chisq       |          42 |      84 |
| Fisher      |          49 |      98 |
| Fisher-midP |          43 |      86 |
| Z-pool      |          43 |      86 |
| Boschloo    |          43 |      86 |

Sample Size Comparison Across Test Methods {.table}

## Impact of Correlation

### Example 4: Correlation Effect

``` r

# Calculate sample size for different correlation values
rho_values <- c(0, 0.3, 0.5, 0.8)

correlation_effect <- lapply(rho_values, function(rho) {
    result <- ss2BinaryExact(
        p11 = 0.70, p12 = 0.60,
        p21 = 0.40, p22 = 0.30,
        rho1 = rho, rho2 = rho,
        r = 1,
        alpha = 0.025,
        beta = 0.2,
        Test = "Fisher"
    )
    data.frame(
        rho = rho,
        n2 = result$n2,
        N = result$N
    )
})

rho_table <- bind_rows(correlation_effect)

kable(rho_table,
      caption = "Impact of Correlation on Sample Size (Fisher's Test)",
      col.names = c("ρ", "n per group", "N total"))
```

|   ρ | n per group | N total |
|----:|------------:|--------:|
| 0.0 |          61 |     122 |
| 0.3 |          60 |     120 |
| 0.5 |          59 |     118 |
| 0.8 |          56 |     112 |

Impact of Correlation on Sample Size (Fisher’s Test) {.table}

**Key finding**: Higher positive correlation reduces required sample
size.

## Comparison: Exact vs Asymptotic

### Example 5: Exact vs AN Method

``` r

# Exact method (Chisq)
exact_result <- ss2BinaryExact(
    p11 = 0.60, p12 = 0.40,
    p21 = 0.30, p22 = 0.10,
    rho1 = 0.5, rho2 = 0.5,
    r = 1,
    alpha = 0.025,
    beta = 0.1,
    Test = "Chisq"
)

# Asymptotic method (AN)
asymp_result <- ss2BinaryApprox(
    p11 = 0.60, p12 = 0.40,
    p21 = 0.30, p22 = 0.10,
    rho1 = 0.5, rho2 = 0.5,
    r = 1,
    alpha = 0.025,
    beta = 0.1,
    Test = "AN"
)

comparison_exact_asymp <- data.frame(
  Method = c("Exact (Chisq)", "Asymptotic (AN)"),
  n_per_group = c(exact_result$n2, asymp_result$n2),
  N_total = c(exact_result$N, asymp_result$N),
  Difference = c(0, asymp_result$N - exact_result$N)
)

kable(comparison_exact_asymp,
      caption = "Comparison: Exact vs Asymptotic Methods",
      col.names = c("Method", "n per group", "N total",
                    "Difference in N total"))
```

| Method          | n per group | N total | Difference in N total |
|:----------------|------------:|--------:|----------------------:|
| Exact (Chisq)   |          59 |     118 |                     0 |
| Asymptotic (AN) |          60 |     120 |                     2 |

Comparison: Exact vs Asymptotic Methods {.table}

## Practical Recommendations

### Test Method Selection

1.  **Fisher’s exact test**:
    - Most widely used and accepted
    - Conservative but guarantees Type I error control
    - Recommended for regulatory submissions
2.  **Boschloo’s test**:
    - Most powerful among exact tests
    - Best choice when computational resources permit
    - Recommended for final analysis
3.  **Chi-squared test**:
    - Less conservative than Fisher
    - May be anti-conservative for small samples
    - Use with caution for $`N < 200`$
4.  **Z-pooled and Fisher-midP**:
    - Intermediate between Fisher and chi-squared
    - Reduce conservatism while maintaining validity

### When to Use Each Method

**Sample size guidelines**:

1.  **$`N < 100`$**: Always use exact methods

2.  **$`100 \leq N < 200`$**: Exact methods preferred, especially if:

    - Extreme probabilities ($`p < 0.1`$ or $`p > 0.9`$)
    - Strict Type I error control required

3.  **$`N \geq 200`$ and $`0.1 < p < 0.9`$**: Asymptotic methods
    acceptable

### Correlation Estimation

- Use pilot data or historical information
- Be conservative if uncertain (use $`\rho = 0`$)
- Consider sensitivity analysis across plausible range

### Allocation Ratio

- Balanced design ($`r = 1`$) generally most efficient
- Unbalanced designs may be justified by:
  - Limited control group availability
  - Ethical considerations
  - Cost constraints

## Computational Considerations

Two quantities dominate the cost of an exact sample size search. The
rejection region has to be built once for each candidate sample size,
and the bivariate binomial probability mass function has to be evaluated
over the whole outcome grid.

### Software Implementation

The `twoCoprimary` package uses:

- the bivariate binomial distribution (`dbibinom`), whose conditional
  probability is evaluated in compiled code, with the powers and
  binomial coefficients that the sum needs tabulated once so that the
  cost grows like $`N^{2}`$ rather than $`N^{3}`$;
- the rejection region (`rr1Binary`), built by vectorized operations
  over the outcome grid;
- a matrix product for the co-primary power. Summing over the outcomes
  of the second group first turns the double sum over the rejection
  region into two matrix products of order $`n + 1`$. Evaluating that
  double sum directly would cost $`O(n^{4})`$ in both time and memory,
  which is prohibitive beyond moderate sample sizes.

A complete sample size search with Boschloo’s test at $`N`$ around 370
takes on the order of a second on a current desktop machine.

## References

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(11), 2183-2201.
