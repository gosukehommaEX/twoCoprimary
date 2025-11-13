# Probability Mass Function of Bivariate Binomial Distribution

Calculates the probability mass function of the bivariate binomial
distribution for given parameters, as described in Homma and Yoshida
(2025).

## Usage

``` r
dbibinom(N, y1, y2, p1, p2, rho)
```

## Arguments

- N:

  Sample size (number of trials)

- y1:

  Observed value(s) of the first random variable (0 to N)

- y2:

  Observed value(s) of the second random variable (0 to N)

- p1:

  True probability of responders for the first outcome (0 \< p1 \< 1)

- p2:

  True probability of responders for the second outcome (0 \< p2 \< 1)

- rho:

  Correlation coefficient between the two binary outcomes

## Value

Probability mass function value(s) for the bivariate binomial
distribution. If y1 and y2 are vectors, returns a vector of
probabilities.

## Details

The bivariate binomial distribution BiBin(N, p1, p2, gamma) has
probability mass function given by equation (3) in Homma and Yoshida
(2025): \$\$P(Y_1 = y_1, Y_2 = y_2) = f(y_1\|N, p_1) \times g(y_2\|y_1,
N, p_1, p_2, \gamma)\$\$ where \$\$g(y_2\|y_1, N, p_1, p_2, \gamma) =
\frac{1}{(1+\gamma)^N} \sum\_{m \in \mathcal{M}} \binom{y_1}{m}
\binom{N-y_1}{y_2-m} (\xi+\gamma)^m (1-\xi)^{y_1-m} \xi^{y_2-m}
(1-\xi+\gamma)^{N-y_1-(y_2-m)}\$\$ with \\\xi = p_2 + \gamma(p_2 -
p_1)\\ and \\\mathcal{M} = \\m : m = \max(0, y_2-(N-y_1)), \ldots,
\min(y_1, y_2)\\\\.

## References

Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical
trials with two co-primary binary endpoints. *Statistical Methods in
Medical Research*, 34(1), 1-19.

## Examples

``` r
# Calculate single probability mass
dbibinom(N = 100, y1 = 30, y2 = 50, p1 = 0.3, p2 = 0.5, rho = 0.5)
#> [1] 0.007981836

# Verify that probabilities sum to 1
N <- 20
p1 <- 0.3
p2 <- 0.5
rho <- 0.5
sum(outer(0:N, 0:N, function(x, y) dbibinom(N, x, y, p1, p2, rho)))
#> [1] 1
```
