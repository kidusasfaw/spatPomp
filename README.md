# `spatPomp` has moved to [spatPomp-org/spatPomp](https://github.com/spatPomp-org/spatPomp/)

[![Development Release](https://img.shields.io/github/release/spatPomp-org/spatPomp.svg)](https://github.com/spatPomp-org/spatPomp/releases/latest)
[![CRAN Status](https://www.r-pkg.org/badges/version/spatPomp?color=blue)](https://cran.r-project.org/package=spatPomp)
[![Last CRAN release date](https://www.r-pkg.org/badges/last-release/spatPomp?color=blue)](https://cran.r-project.org/package=spatPomp)

[![R-CMD-check](https://github.com/spatPomp-org/spatPomp/actions/workflows/r-cmd-check.yml/badge.svg)](https://github.com/spatPomp-org/spatPomp/actions/workflows/r-cmd-check.yml)
[![test-coverage](https://github.com/spatPomp-org/spatPomp/actions/workflows/test-coverage.yml/badge.svg)](https://github.com/spatPomp-org/spatPomp/actions/workflows/test-coverage.yml)
[![codecov](https://codecov.io/github/spatPomp-org/spatPomp/graph/badge.svg?token=O97GJYUGNH)](https://codecov.io/github/spatPomp-org/spatPomp)

![CRAN mirror monthly downloads](https://cranlogs.r-pkg.org/badges/last-month/spatPomp?color=yellow)
![CRAN mirror total downloads](https://cranlogs.r-pkg.org/badges/grand-total/spatPomp?color=yellow)


## What is this package?
The [`spatPomp` package](https://spatPomp-org.github.io/spatPomp/) provides facilities for inference on panel data using spatiotemporal partially-observed Markov process (SpatPOMP) models.
To do so, it relies on and extends a number of facilities that the [`pomp` package](https://kingaa.github.io/pomp/) provides for inference on time series data using partially-observed Markov process (POMP) models.

## Why use `spatPomp`?
The [`spatPomp` package](https://spatPomp-org.github.io/spatPomp/) concerns models consisting of a collection of interacting units.
The methods in `spatPomp` may be applicable whether or not these units correspond to spatial locations.

## Installing `spatPomp`
To install the package from this GitHub source, type the following command into your console (assumes you have the `devtools` package installed):
`devtools::install_github('spatPomp-org/spatPomp')`.
The CRAN version of `spatPomp` can be installed by `install.packages("spatPomp")`.

## Documentation
Details on the motivation and use of this package can be found in Asfaw et al. (2021) (https://arxiv.org/abs/2101.01157).
Further information can be found on the [`spatPomp` website](https://spatPomp-org.github.io/spatPomp/).
