# Rejection Region for Two-Arm Trials with a Single Binary Endpoint

Calculates the rejection region for two-arm trials with a single binary
endpoint using various exact statistical tests, as described in Homma
and Yoshida (2025).

## Usage

``` r
rr1Binary(n1, n2, alpha, Test)
```

## Arguments

- n1:

  Sample size for group 1 (test group)

- n2:

  Sample size for group 2 (control group)

- alpha:

  One-sided significance level (typically 0.025)

- Test:

  Type of statistical test. One of:

  - `"Chisq"`: One-sided Pearson chi-squared test

  - `"Fisher"`: Fisher exact test

  - `"Fisher-midP"`: Fisher mid-p test

  - `"Z-pool"`: Z-pooled exact unconditional test

  - `"Boschloo"`: Boschloo exact unconditional test

## Value

A logical matrix of dimensions (n1+1) x (n2+1), where TRUE indicates
rejection of the null hypothesis. Rows correspond to the number of
responders in group 1 (0 to n1), and columns correspond to the number of
responders in group 2 (0 to n2).

## Details

This function computes the rejection region for five different one-sided
tests:

**Chi-squared test:** Uses the asymptotic normal approximation of the
chi-squared statistic.

**Fisher exact test:** Uses the hypergeometric distribution to calculate
exact p-values conditional on the total number of successes.

**Fisher mid-p test:** Modification of Fisher's exact test that adds
half the probability of the observed outcome to reduce conservatism.

**Z-pooled test:** Exact unconditional test that maximizes p-values over
all possible values of the nuisance parameter (common success
probability under H0).

**Boschloo test:** Exact unconditional test similar to Z-pooled but
based on Fisher's exact p-values, maximizing over the nuisance
parameter.

## References

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(1), 1-19.

## Examples

``` r
# Simple example with small sample sizes
n1 <- 5
n2 <- 5
alpha <- 0.025
RR <- rr1Binary(n1, n2, alpha, Test = 'Chisq')
print(dim(RR))  # Should be (6, 6)
#> [1] 6 6

# Fisher exact test
RR_fisher <- rr1Binary(n1 = 10, n2 = 10, alpha = 0.025, Test = 'Fisher')

# \donttest{
# More computationally intensive: Boschloo test
n1 <- 20
n2 <- 10
alpha <- 0.025
RR <- rr1Binary(n1, n2, alpha, Test = 'Boschloo')
print(RR)
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
#>  [1,] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [2,] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [3,] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [4,] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [5,] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [6,] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [7,] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [8,]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [9,]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [10,]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [11,]  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [12,]  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [13,]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [14,]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [15,]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [16,]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [17,]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
#> [18,]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
#> [19,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE
#> [20,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE
#> [21,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE
# }
```
