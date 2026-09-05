## R CMD check results

0 errors | 0 warnings | 0 notes on every environment listed below.

## Submission

This is a patch release. It fixes 27 defects. Twenty-two were found by a
systematic audit of the package and by reproducing the tables of the five
articles the methods come from. The remaining five were found by an exhaustive
sweep of the argument space of every exported function, written after version
1.1.0 was published and kept in the package's public repository rather than in
the sources submitted here.

Version 1.1.0 reached CRAN on 2026-08-29, so this submission comes sooner after
the previous one than the convention on update frequency allows. We would not
normally submit this soon. Three classes of defect led us to judge that an early
fix was warranted, and we followed the guidance to submit a corrected version
and to explain the reason here.

* Several documented calls to `plot()` failed with an error rather than
  producing a figure. A sample size object with `type = "power_curve"`, a power
  object with `type = "sample_size_rho"` and any object returned by the exact
  binary functions each entered a branch written for a different object shape.

* Two functions returned incorrect numbers. `ss1BinaryApprox()` used the arcsine
  multiplier written for the reciprocal of its own allocation ratio, so at an
  allocation of two to one it returned about twice the required sample size, and
  its `"ASc"` continuity correction moved the two groups apart rather than
  toward each other. `corrbound2MixedCountContinuous()` computed the correlation
  bounds by a quadrature in the untransformed continuous endpoint, so the bounds
  collapsed to zero when the mean lay far from the origin relative to the
  standard deviation, and the two functions that use them then rejected every
  admissible correlation.

* `rr1Binary()` failed with an error about an incorrect number of dimensions
  when either group held a single subject, and `power2Continuous()` accepted a
  negative standard deviation, a group size of zero and a non-integer group size
  and returned a number for each.

The four asymptotic methods of `ss1BinaryApprox()` are now a sequential search
against `power2BinaryApprox()`, so the sample size returned and the power
reported by the package refer to the same expression. The twelve published
single-endpoint sizes of Table S5 of the Supporting Information of Sozu et al.
(2012) are reproduced exactly. The two-endpoint sample size functions use
`ss1BinaryApprox()` only as the starting value of a search that converges to the
minimum sample size from any starting point, so their results are unchanged.

`ss2MixedCountContinuous()` evaluates the correlation bounds once instead of at
every step of its search, which takes a call at the default settings from about
2.7 seconds to about 0.2.

`NeedsCompilation` remains yes, as in 1.1.0. There are no user visible interface
changes and no new dependencies.

## Reverse dependencies

There are currently no downstream dependencies for this package.

## Test environments

* local Windows 11 x64 install, R 4.6.0
* win-builder: R-devel (2026-09-04 r90492), R-release (4.6.1),
  R-oldrelease (4.5.3)
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  windows-latest (R-release), macos-latest (R-release)
