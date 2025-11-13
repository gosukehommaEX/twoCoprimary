# Calculate Correlation Bounds Between Count and Continuous Outcomes

Computes the lower and upper bounds of the correlation coefficient
between an overdispersed count outcome (negative binomial) and a
continuous outcome (normal), as described in Homma and Yoshida (2024).

## Usage

``` r
corrbound2MixedCountContinuous(lambda, nu, mu, sd)
```

## Arguments

- lambda:

  Mean parameter for the negative binomial distribution (lambda \> 0)

- nu:

  Dispersion parameter for the negative binomial distribution (nu \> 0)

- mu:

  Mean for the continuous outcome

- sd:

  Standard deviation for the continuous outcome (sd \> 0)

## Value

A named numeric vector with two elements:

- L_bound:

  Lower bound of the correlation

- U_bound:

  Upper bound of the correlation

## Details

The correlation bounds are calculated using the Frechet-Hoeffding bounds
for copulas, as described in Trivedi and Zimmer (2007). The negative
binomial distribution has mean lambda and variance: \$\$Var(Y_1) =
\lambda + \frac{\lambda^2}{\nu}\$\$

The variance of the negative binomial distribution is: Var(Y1) =
lambda + lambda^2/nu

## References

Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical
trials with two co-primary endpoints including overdispersed count and
continuous outcomes. *Pharmaceutical Statistics*, 23(1), 46-59.

Trivedi, P. K., & Zimmer, D. M. (2007). Copula modeling: an introduction
for practitioners. *Foundations and Trends in Econometrics*, 1(1),
1-111.

## Examples

``` r
# Calculate correlation bounds for NB(1.25, 0.8) and N(0, 250)
corrbound2MixedCountContinuous(lambda = 1.25, nu = 0.8, mu = 0, sd = 250)
#>    L_bound    U_bound 
#> -0.8457747  0.8457747 

# Higher dispersion parameter
corrbound2MixedCountContinuous(lambda = 2.0, nu = 2.0, mu = 50, sd = 200)
#>    L_bound    U_bound 
#> -0.9209669  0.9209671 

# Different follow-up time
corrbound2MixedCountContinuous(lambda = 1.0 * 2, nu = 1.0, mu = 0, sd = 300)
#>    L_bound    U_bound 
#> -0.8812028  0.8812028 
```
