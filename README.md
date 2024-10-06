
# BayesianBombCarbon
Bayesian methods including penalized B-splines for fitting bomb radiocarbon reference series and assessing aging bias

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![GPL3
license](https://img.shields.io/badge/License-GPL3-blue.svg)](https://github.com/zsiders/BayesianBombCarbon/blob/main/LICENSE)
<!-- badges: end -->

## About

The `BayesianBombCarbon` R package provides Bayesian models in `stan` to estimate ∆14C reference series, estimate aging bias of samples with unvalidated formation/birth years given a reference series, and assess if the aging of samples are valid via hypothesis testing. Currently, reference series are fit using a Bayesian penalized B-spline with a random walk prior to control the degree of smoothing, first presented in [Chamberlin et al. 2023](https://doi.org/10.1038/s41598-023-34680-0). 


## Installing the package

The development version of `BayesianBombCarbon` can be installed from GitHub as shown below. Please note that the package is still in early development and may be subject to breaking changes.

``` r
remotes::install_github("zsiders/BayesianBombCarbon", dependencies = TRUE)
```

`BayesianBombCarbon` also requires CmdStan to be installed on your system. This can be done using the `install_cmdstan()` function from `cmdstanr`. If you experience any problems installing CmdStan, see the [cmdstanr vignette](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) for help.

``` r
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan(cores = 2) # use more cores to speed up
```

The stan models used by BayesianBombCarbon need to be compiled for your device. This is only necessary once - after installing or updating the package - and can be done using the `bombcarbon_compile()` function.

``` r
BayesianBombCarbon::bombcarbon_compile()
```

If the models are not successfully compiled, please ensure that
`cmdstan` is properly set up and try updating it to a newer version
using `cmdstanr::install_cmdstan()`. If the problem persists, please run
`BayesianBombCarbon::bombcarbon_compile(verbose = TRUE)` and post the output in a new issue on GitHub, along with your `cmdstanr::cmdstan_version()`.


