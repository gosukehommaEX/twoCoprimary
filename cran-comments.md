## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Resubmission

This is a resubmission. In this version I have:

* Fixed URL field formatting in DESCRIPTION: Ensured proper comma-space 
  separation with no line breaks
  
* Wrapped all computationally intensive examples in ss2BinaryExact(), 
  power2BinaryExact(), and twoCoprimary2BinaryExact() with \donttest{} 
  to reduce check time

The spelling NOTEs for "Homma", "Sozu", "Yoshida", "et", and "al" are author 
names from published references and are spelled correctly.

## Test environments

* local Windows install, R 4.5.0
* GitHub Actions (ubuntu-latest, windows-latest, macos-latest): R-devel, 
  R-release, R-oldrel-1

## Downstream dependencies

There are currently no downstream dependencies for this package.
