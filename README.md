# The `spatPomp` package

[![R-CMD-check](https://github.com/kidusasfaw/spatPomp/actions/workflows/r-cmd-check.yml/badge.svg)](https://github.com/kidusasfaw/spatPomp/actions/workflows/r-cmd-check.yml)[![test-coverage](https://github.com/kidusasfaw/spatPomp/actions/workflows/test-coverage.yml/badge.svg)](https://github.com/kidusasfaw/spatPomp/actions/workflows/test-coverage.yml)

## What is this package?
The `spatPomp` package provides facilities for inference on panel data using spatiotemporal partially-observed Markov process (SpatPOMP) models.
To do so, it relies on and extends a number of facilities that the `pomp` package provides for inference on time series data using partially-observed Markov process (POMP) models.

## Why use `spatPomp`?
The `spatPomp` package concerns models consisting of a collection of interacting units.
The methods in `spatPomp` may be applicable whether or not these units correspond to spatial locations.

## Installing `spatPomp`
To install the package from this GitHub source, type the following command into your console (assumes you have the `devtools` package installed):
`devtools::install_github('kidusasfaw/spatPomp')`

## Documentation
Details on the motivation and use of this package can be found in Asfaw et al. (2021) (https://arxiv.org/abs/2101.01157).

