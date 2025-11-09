# twoCoprimary 1.0.0

## Initial Release

This is the first release of twoCoprimary, providing comprehensive tools for sample size and power calculation in clinical trials with two co-primary endpoints.

### Features

* **Two Continuous Endpoints**
  - `ss2Continuous()`, `power2Continuous()`, `twoCoprimary2Continuous()`
  - Based on Sozu et al. (2011)
  - Supports known and unknown variance cases

* **Two Binary Endpoints (Asymptotic Methods)**
  - `ss2BinaryApprox()`, `power2BinaryApprox()`, `twoCoprimary2BinaryApprox()`
  - Based on Sozu et al. (2010)
  - Four test methods: AN, ANc, AS, ASc

* **Two Binary Endpoints (Exact Methods)**
  - `ss2BinaryExact()`, `power2BinaryExact()`, `twoCoprimary2BinaryExact()`
  - Based on Homma and Yoshida (2025)
  - Five exact tests: Chisq, Fisher, Fisher-midP, Z-pool, Boschloo

* **Mixed Continuous and Binary Endpoints**
  - `ss2MixedContinuousBinary()`, `power2MixedContinuousBinary()`, `twoCoprimary2MixedContinuousBinary()`
  - Based on Sozu et al. (2012)
  - Supports biserial correlation structure

* **Mixed Count and Continuous Endpoints**
  - `ss2MixedCountContinuous()`, `power2MixedCountContinuous()`, `twoCoprimary2MixedCountContinuous()`
  - Based on Homma and Yoshida (2024)
  - Handles overdispersed count data with negative binomial distribution

### Utility Functions

* `corrbound2Binary()` - Calculate valid correlation bounds for binary endpoints
* `corrbound2MixedCountContinuous()` - Calculate valid correlation bounds for count and continuous endpoints
* `design_table()` - Create comprehensive design comparison tables
* `plot.twoCoprimary()` - Visualize sample size vs correlation relationships

### Documentation

* Six comprehensive vignettes covering all methodologies
* Complete function documentation with examples
* Validation against published results

### Testing

* Comprehensive test suite with testthat
* Tests for all major functions
* Validation against published tables and results
