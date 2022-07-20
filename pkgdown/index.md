
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

# smouseGPS <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: end -->

The goal of smouseGPS is to smooth and interpolate noisy GPS data. It
enables to deal with outliers by indentifying them thanks to a
collection of b-splines. This package is inspired by the work of Jeffrey
J. Early and Adam M. Sykulski
(<https://jeffreyearly.com/smoothing-and-interpolating-noisy-gps-data/>).
It is a translation from Matlab language to R language.

## Installation

You can install the development version of smouseGPS like so:

``` r
# install.packages("remotes")
remotes::install_github("alexmerot/smouseGPS")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(smouseGPS)
## basic example code
```