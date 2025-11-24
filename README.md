# twoCoprimary <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/twoCoprimary/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/twoCoprimary/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/twoCoprimary)](https://CRAN.R-project.org/package=twoCoprimary)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/twoCoprimary)](https://CRAN.R-project.org/package=twoCoprimary)
[![Monthly downloads](https://cranlogs.r-pkg.org/badges/twoCoprimary)](https://CRAN.R-project.org/package=twoCoprimary)
<!-- badges: end -->

## Overview

`twoCoprimary` provides comprehensive tools for sample size and power calculation in clinical trials with two co-primary endpoints. In co-primary endpoint trials, treatment success requires demonstrating statistically significant effects on **all** primary endpoints simultaneously. This package implements state-of-the-art methodologies that properly account for correlation between endpoints, leading to more efficient trial designs.

📖 **[Documentation Website](https://gosukehommaEX.github.io/twoCoprimary/)**

## Key Features

The package supports five combinations of co-primary endpoints:

- **Two continuous endpoints** - Normal distribution with correlation (Sozu et al., 2011)
- **Two binary endpoints (asymptotic)** - Large sample approximation methods (Sozu et al., 2010)
- **Two binary endpoints (exact)** - Exact inference for small to medium samples (Homma & Yoshida, 2025)
- **Mixed continuous and binary** - Biserial correlation structure (Sozu et al., 2012)
- **Mixed count and continuous** - Negative binomial for overdispersed counts (Homma & Yoshida, 2024)

All methods provide:
- ✅ Sample size calculation given target power
- ✅ Power calculation given sample size
- ✅ Proper Type I error control without multiplicity adjustment
- ✅ Accounting for correlation between endpoints
- ✅ Support for unbalanced allocation ratios

## Installation

Install from CRAN:

``` r
install.packages("twoCoprimary")
```

Or install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("gosukehommaEX/twoCoprimary")
```

## Quick Start

### Example 1: Two Continuous Endpoints

Calculate sample size for a trial with two continuous co-primary endpoints:

``` r
library(twoCoprimary)

# Sample size calculation
result <- ss2Continuous(
  delta1 = 0.5,      # Standardized effect size for endpoint 1
  delta2 = 0.4,      # Standardized effect size for endpoint 2
  rho = 0.3,         # Correlation between endpoints
  alpha = 0.025,     # One-sided significance level
  power = 0.80,      # Target power
  r = 1              # Allocation ratio (1:1)
)

print(result)
#> 
#> Sample size calculation for two continuous co-primary endpoints
#> 
#> Parameters:
#>   Effect sizes: delta1 = 0.5, delta2 = 0.4
#>   Correlation: rho = 0.3
#>   Significance level: alpha = 0.025 (one-sided)
#>   Target power: 1 - beta = 0.8
#>   Allocation ratio: r = 1
#> 
#> Results:
#>   Sample size per group: n1 = n2 = 130
#>   Total sample size: N = 260
```

### Example 2: Two Binary Endpoints (Exact Method)

For small to medium sample sizes, use exact methods:

``` r
# Sample size with exact inference
result_exact <- ss2BinaryExact(
  p11 = 0.30, p12 = 0.15,    # Response rates for endpoint 1
  p21 = 0.50, p22 = 0.30,    # Response rates for endpoint 2
  rho1 = 0.3, rho2 = 0.3,    # Within-group correlations
  alpha = 0.025,             # One-sided significance level
  power = 0.80,              # Target power
  r = 1,                     # Allocation ratio
  test_method = "Fisher"     # Exact test method
)

print(result_exact)
```

### Example 3: Mixed Count and Continuous Endpoints

For COPD/asthma trials with exacerbation count and lung function:

``` r
# Exacerbation rates (events per year)
r1 <- 0.80  # Treatment group
r2 <- 1.25  # Control group

# Sample size calculation
result_mixed <- ss2MixedCountContinuous(
  r1 = r1, r2 = r2,          # Event rates
  nu = 1.0,                  # Dispersion parameter
  t = 1,                     # Follow-up period (years)
  mu1 = 250, mu2 = 200,      # Mean FEV1 (mL)
  sd = 300,                  # Common SD
  rho1 = 0.5, rho2 = 0.5,    # Correlations
  alpha = 0.025,
  power = 0.80,
  r = 1
)

print(result_mixed)
```

## Documentation

Comprehensive vignettes are available:

- [Overview of co-primary endpoints](https://gosukehommaEX.github.io/twoCoprimary/articles/overview.html)
- [Two continuous endpoints](https://gosukehommaEX.github.io/twoCoprimary/articles/two-continuous-endpoints.html)
- [Two binary endpoints (asymptotic)](https://gosukehommaEX.github.io/twoCoprimary/articles/two-binary-endpoints-approx.html)
- [Two binary endpoints (exact)](https://gosukehommaEX.github.io/twoCoprimary/articles/two-binary-endpoints-exact.html)
- [Mixed continuous and binary](https://gosukehommaEX.github.io/twoCoprimary/articles/mixed-continuous-binary.html)
- [Mixed count and continuous](https://gosukehommaEX.github.io/twoCoprimary/articles/mixed-count-continuous.html)

## References

1. Homma, G., & Yoshida, T. (2024). Sample size calculation for clinical trials with co‐primary outcomes: Negative binomial and continuous outcomes. *Pharmaceutical Statistics*, 23(3), 368-392. https://doi.org/10.1002/pst.2337

2. Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical trials with two co-primary binary endpoints. *Statistical Methods in Medical Research*, 34(1). https://doi.org/10.1177/09622802251368697

3. Sozu, T., Sugimoto, T., Hamasaki, T., & Evans, S. R. (2010). Sample size determination in superiority clinical trials with multiple co-primary correlated endpoints. *Statistics in Medicine*, 29(21), 2219-2227. https://doi.org/10.1002/sim.3972

4. Sozu, T., Sugimoto, T., & Hamasaki, T. (2011). Sample size determination in superiority clinical trials with multiple co-primary correlated endpoints. *Journal of Biopharmaceutical Statistics*, 21(4), 650-668. https://doi.org/10.1080/10543406.2011.551329

5. Sozu, T., Sugimoto, T., Hamasaki, T., & Evans, S. R. (2012). Sample size determination in clinical trials with multiple co-primary binary endpoints including mixed binary and continuous endpoints. *Biometrical Journal*, 54(5), 716-729. https://doi.org/10.1002/bimj.201100221



## Citation

```r
citation("twoCoprimary")
```

## Getting Help

- For bug reports and feature requests, please use the [GitHub issue tracker](https://github.com/gosukehommaEX/twoCoprimary/issues)
- For questions about usage, see the [documentation website](https://gosukehommaEX.github.io/twoCoprimary/)

## License

MIT © Gosuke Homma
