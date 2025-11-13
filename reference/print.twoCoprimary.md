# Print Method for twoCoprimary Objects

Provides a clean, formatted display of sample size and power calculation
results from the twoCoprimary package.

## Usage

``` r
# S3 method for class 'twoCoprimary'
print(x, ...)
```

## Arguments

- x:

  An object of class "twoCoprimary"

- ...:

  Additional arguments (currently unused)

## Value

Invisibly returns the original object

## Details

This print method provides a formatted output that displays key
parameters and results in an easy-to-read format. The specific format
adapts to the type of calculation (sample size vs power) and the type of
endpoints involved.
