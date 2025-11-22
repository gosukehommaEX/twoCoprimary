# Changelog

## twoCoprimary 1.0.0

CRAN release: 2025-11-21

### Initial Release

This is the first release of twoCoprimary, providing comprehensive tools
for sample size and power calculation in clinical trials with two
co-primary endpoints.

#### Features

- **Two Continuous Endpoints**
  - [`ss2Continuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2Continuous.md),
    [`power2Continuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2Continuous.md),
    [`twoCoprimary2Continuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2Continuous.md)
  - Based on Sozu et al. (2011)
  - Supports known and unknown variance cases
- **Two Binary Endpoints (Asymptotic Methods)**
  - [`ss2BinaryApprox()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2BinaryApprox.md),
    [`power2BinaryApprox()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2BinaryApprox.md),
    [`twoCoprimary2BinaryApprox()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2BinaryApprox.md)
  - Based on Sozu et al. (2010)
  - Four test methods: AN, ANc, AS, ASc
- **Two Binary Endpoints (Exact Methods)**
  - [`ss2BinaryExact()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2BinaryExact.md),
    [`power2BinaryExact()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2BinaryExact.md),
    [`twoCoprimary2BinaryExact()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2BinaryExact.md)
  - Based on Homma and Yoshida (2025)
  - Five exact tests: Chisq, Fisher, Fisher-midP, Z-pool, Boschloo
- **Mixed Continuous and Binary Endpoints**
  - [`ss2MixedContinuousBinary()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2MixedContinuousBinary.md),
    [`power2MixedContinuousBinary()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2MixedContinuousBinary.md),
    [`twoCoprimary2MixedContinuousBinary()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2MixedContinuousBinary.md)
  - Based on Sozu et al. (2012)
  - Supports biserial correlation structure
- **Mixed Count and Continuous Endpoints**
  - [`ss2MixedCountContinuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/ss2MixedCountContinuous.md),
    [`power2MixedCountContinuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/power2MixedCountContinuous.md),
    [`twoCoprimary2MixedCountContinuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/twoCoprimary2MixedCountContinuous.md)
  - Based on Homma and Yoshida (2024)
  - Handles overdispersed count data with negative binomial distribution

#### Utility Functions

- [`corrbound2Binary()`](https://gosukehommaEX.github.io/twoCoprimary/reference/corrbound2Binary.md) -
  Calculate valid correlation bounds for binary endpoints
- [`corrbound2MixedCountContinuous()`](https://gosukehommaEX.github.io/twoCoprimary/reference/corrbound2MixedCountContinuous.md) -
  Calculate valid correlation bounds for count and continuous endpoints
- [`design_table()`](https://gosukehommaEX.github.io/twoCoprimary/reference/design_table.md) -
  Create comprehensive design comparison tables
- [`plot.twoCoprimary()`](https://gosukehommaEX.github.io/twoCoprimary/reference/plot.twoCoprimary.md) -
  Visualize sample size vs correlation relationships

#### Documentation

- Six comprehensive vignettes covering all methodologies
- Complete function documentation with examples
- Validation against published results

#### Testing

- Comprehensive test suite with testthat
- Tests for all major functions
- Validation against published tables and results
