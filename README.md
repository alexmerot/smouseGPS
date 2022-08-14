
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

# smouseGPS <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: end -->

<div class="alert alert-warning"
style="font-variant: small-caps;font-size: 2em;text-align: center;margin-top: 10%;">

**Work in Progress**

</div>

The goal of the R package is to smooth and interpolate noisy GPS data.
It enables to deal with outliers by indentifying them thanks to a
collection of b-splines. This package is based on the work of Jeffrey J.
Early and Adam M. Sykulski
(<https://jeffreyearly.com/smoothing-and-interpolating-noisy-gps-data/>).
It is a translation from Matlab language to R language. The
documentation is on my github page
[alexmerot.github.io/smouseGPS/](https://alexmerot.github.io/smouseGPS/).

## Installation

You can install the development version of smouseGPS like so:

``` r
# install.packages("remotes")
remotes::install_github("alexmerot/smouseGPS")
```
