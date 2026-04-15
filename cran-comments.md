## Update notes for betaARMA 1.1.0
This is a minor release that fixes a mathematical bug in the score vector derivation for non-logit link functions and refactors the underlying parameter architecture for improved stability. We have also exported several internal functions for advanced users.

## Test environments
* local Ubuntu 24.04, R 4.4.2
* win-builder (devel and release)
* macOS (release) via devtools

## R CMD check results
0 errors | 0 warnings | 2 notes

* NOTE 1: Found the following (possibly) invalid URLs: https://www.scimagojr.com/... (Status: 403)
  * Explanation: This is a false positive. The URL points to the SCImago Journal Rank page for the journal TEST (Springer). The link is valid and works perfectly in a standard web browser, but SCImago's server blocks automated ping requests (like CRAN's urlchecker) with a 403 Forbidden error due to their bot protection. This same link is used and approved in our other CRAN packages.

* NOTE 2: unable to verify current time
  * Explanation: Local Ubuntu environment issue connecting to the time server during the check.

## Reverse dependencies
There are currently no downstream dependencies for this package.