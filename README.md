# twoCoprimary

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/twoCoprimary/workflows/R-CMD-check/badge.svg)](https://github.com/gosukehommaEX/twoCoprimary/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/twoCoprimary)](https://CRAN.R-project.org/package=twoCoprimary)
<!-- badges: end -->

## Overview

`twoCoprimary` provides comprehensive tools for sample size and power calculation in clinical trials with two co-primary endpoints. In co-primary endpoint trials, treatment success requires demonstrating statistically significant effects on **all** primary endpoints simultaneously. This package implements state-of-the-art methodologies that properly account for correlation between endpoints, leading to more efficient trial designs.

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

Install the development version from GitHub:

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
  delta2 = 0.5,      # Standardized effect size for endpoint 2
  sd1 = 1,           # Standard deviation for endpoint 1
  sd2 = 1,           # Standard deviation for endpoint 2
  rho = 0.5,         # Correlation between endpoints
  r = 1,             # Balanced allocation (1:1)
  alpha = 0.025,     # One-sided significance level
  beta = 0.2,        # Type II error (power = 0.8)
  known_var = TRUE
)

print(result)
#> 
#>  Sample size calculation for two continuous co-primary endpoints
#> 
#>      delta1 = 0.5
#>      delta2 = 0.5
#>         sd1 = 1
#>         sd2 = 1
#>         rho = 0.5
#>           r = 1
#>       alpha = 0.025
#>        beta = 0.2
#>   known_var = TRUE
#>          n1 = 85
#>          n2 = 85
#>           N = 170
```

### Example 2: Two Binary Endpoints (Exact Methods)

Calculate sample size using exact inference for small sample sizes:

``` r
# Exact sample size calculation
result <- ss2BinaryExact(
  p11 = 0.7,         # Success probability for endpoint 1 in treatment
  p12 = 0.6,         # Success probability for endpoint 2 in treatment
  p21 = 0.4,         # Success probability for endpoint 1 in control
  p22 = 0.3,         # Success probability for endpoint 2 in control
  rho1 = 0.5,        # Correlation in treatment group
  rho2 = 0.5,        # Correlation in control group
  r = 1,             # Balanced allocation
  alpha = 0.025,
  beta = 0.2,
  Test = "Fisher"    # Fisher's exact test
)

print(result)
#> 
#>  Sample size calculation for two binary co-primary endpoints (exact)
#> 
#>         p11 = 0.7
#>         p12 = 0.6
#>         p21 = 0.4
#>         p22 = 0.3
#>        rho1 = 0.5
#>        rho2 = 0.5
#>           r = 1
#>       alpha = 0.025
#>        beta = 0.2
#>        Test = Fisher
#>          n1 = 56
#>          n2 = 56
#>           N = 112
```

## Documentation

For detailed methodology, examples, and validation against published results, see the comprehensive vignettes:

- `vignette("overview")` - Package overview and statistical framework
- `vignette("two-continuous-endpoints")` - Continuous endpoints methodology
- `vignette("two-binary-endpoints-approx")` - Binary endpoints (asymptotic methods)
- `vignette("two-binary-endpoints-exact")` - Binary endpoints (exact methods)
- `vignette("mixed-continuous-binary")` - Mixed continuous and binary endpoints
- `vignette("mixed-count-continuous")` - Mixed count and continuous endpoints

## Why Account for Correlation?

When endpoints are positively correlated, accounting for this correlation can substantially reduce required sample sizes (typically 5-15% reduction for ρ = 0.5-0.8) compared to assuming independence, while maintaining proper Type I error control.

## References

**Methodology papers:**

- Homma, G., & Yoshida, T. (2025). Exact power and sample size in clinical trials with two co-primary binary endpoints. *Statistical Methods in Medical Research*, 34(1), 1-19. https://doi.org/10.1177/09622802251368697

- Homma, G., & Yoshida, T. (2024). Sample size calculation in clinical trials with two co-primary endpoints including overdispersed count and continuous outcomes. *Pharmaceutical Statistics*, 23(1), 46-59. https://doi.org/10.1002/pst.2337

- Sozu, T., Sugimoto, T., & Hamasaki, T. (2012). Sample size determination in clinical trials with multiple co-primary endpoints including mixed continuous and binary variables. *Biometrical Journal*, 54(5), 716-729. https://doi.org/10.1002/bimj.201100221

- Sozu, T., Sugimoto, T., & Hamasaki, T. (2011). Sample size determination in superiority clinical trials with multiple co-primary correlated endpoints. *Journal of Biopharmaceutical Statistics*, 21(4), 650-668. https://doi.org/10.1080/10543406.2011.551329

- Sozu, T., Sugimoto, T., & Hamasaki, T. (2010). Sample size determination in clinical trials with multiple co-primary binary endpoints. *Statistics in Medicine*, 29(21), 2169-2179. https://doi.org/10.1002/sim.3972

## License

MIT License. See [LICENSE](LICENSE) file for details.

## Author

**Gosuke Homma**  
Email: my.name.is.gosuke@gmail.com  
GitHub: [@gosukehommaEX](https://github.com/gosukehommaEX)

## Citation

If you use this package in your research, please cite the relevant methodology papers listed above.

---

For questions, bug reports, or feature requests, please open an issue on [GitHub](https://github.com/gosukehommaEX/twoCoprimary/issues).
