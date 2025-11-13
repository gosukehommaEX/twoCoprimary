# Calculate Correlation Bounds Between Two Binary Outcomes

Computes the lower and upper bounds of the correlation coefficient
between two binary outcomes based on their marginal probabilities, as
described in Prentice (1988).

## Usage

``` r
corrbound2Binary(p1, p2)
```

## Arguments

- p1:

  True probability of responders for the first outcome (0 \< p1 \< 1)

- p2:

  True probability of responders for the second outcome (0 \< p2 \< 1)

## Value

A named numeric vector with two elements:

- L_bound:

  Lower bound of the correlation

- U_bound:

  Upper bound of the correlation

## Details

For two binary outcomes with marginal probabilities p1 and p2, the
correlation coefficient rho is bounded by: \$\$\rho \in \[L(p_1, p_2),
U(p_1, p_2)\]\$\$ where \$\$L(p_1, p_2) = \max\left\\-\sqrt{\frac{p_1
p_2}{(1-p_1)(1-p_2)}}, -\sqrt{\frac{(1-p_1)(1-p_2)}{p_1
p_2}}\right\\\$\$ \$\$U(p_1, p_2) =
\min\left\\\sqrt{\frac{p_1(1-p_2)}{p_2(1-p_1)}},
\sqrt{\frac{p_2(1-p_1)}{p_1(1-p_2)}}\right\\\$\$

## References

Prentice, R. L. (1988). Correlated binary regression with covariates
specific to each binary observation. *Biometrics*, 44(4), 1033-1048.

## Examples

``` r
# Calculate correlation bounds for two binary outcomes
corrbound2Binary(p1 = 0.3, p2 = 0.5)
#>    L_bound    U_bound 
#> -0.6546537  0.6546537 

# When probabilities are equal, upper bound is 1
corrbound2Binary(p1 = 0.4, p2 = 0.4)
#>    L_bound    U_bound 
#> -0.6666667  1.0000000 

# When p1 + p2 = 1, lower bound is -1
corrbound2Binary(p1 = 0.3, p2 = 0.7)
#>    L_bound    U_bound 
#> -1.0000000  0.4285714 
```
