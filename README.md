
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

# smouseGPS <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: end -->

<div class="alert alert-warning"
style="font-variant: small-caps; font-weight: bold; font-size: 2em; color: #a56404; border-radius: 33px; text-align: center;content-align: center; margin: 20% 2 2 2; padding: 2 0 0 0;"
role="alert">

<svg xmlns="http://www.w3.org/2000/svg" width="10%" heigth="10%" viewBox="0 0 576 512">
<!--! Font Awesome Pro 6.1.2 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license (Commercial License) Copyright 2022 Fonticons, Inc. -->
<path fill="#a56404" d="M272 95.93c26.5 0 47.99-21.47 47.99-47.97S298.5 0 272 0C245.5 0 224 21.47 224 47.97S245.5 95.93 272 95.93zM209.7 357.3c-25.75-17.25-52.25-33.24-79.5-48.11L58.62 270.2L1.246 471.1c-4.875 16.1 4.1 34.74 22 39.62s34.63-4.998 39.5-21.99l36.63-128.1l60.63 40.37v78.86c0 17.62 14.38 31.99 32 31.99s32-14.37 32-31.99l.0022-95.93C224 373.2 218.6 363.2 209.7 357.3zM311.1 416c-13.88 0-25.95 8.863-30.33 21.86l-24.75 74.07h319.9l-101.9-206.3c-11.38-22.49-43.1-23.63-56.1-2.01l-31.89 54.21l-65.26-35.64l-24-121.2C288.1 161.3 263.2 127.7 227.1 109.7c-1-.4999-2.125-.625-3.125-1.125c-2.25-1.125-4.752-1.1-7.252-2.625C201.5 99.85 185.2 95.98 168.7 95.98H95.1c-9.25 0-18.05 4.061-24.18 10.93l-55.95 63.92c-.75 .9998-1.5 2.124-2.25 3.249c-8.875 13.1-3 32.87 11.63 40.74l336.6 184.3l-9.837 16.87H311.1zM105.9 204.1l-23.5-12.87l28.13-32.12h34.38L105.9 204.1zM199.5 256.1l34.9-41.28l13.5 67.61L199.5 256.1z"/>
</svg>

Work in Progress

</div>

The goal of the R package `smouseGPS` is to smooth and interpolate noisy
GPS data. It enables to deal with outliers by indentifying them thanks
to a collection of b-splines. This package is based on the work of
Jeffrey J. Early and Adam M. Sykulski
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
