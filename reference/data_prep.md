# Prepares data for Bayesian penalized B-spline

Forces a standardized format for inputting bomb radidocarbon data into a
list object for estimating with 'cmdstanr'

## Usage

``` r
data_prep(
  df_ref,
  df_unk,
  model = "bspline",
  likelihood = "normal",
  coupled_function = "exponential",
  ll_wt = 5,
  pred.by = list(min.by = 1940, max.by = 2020, inc.by = 1),
  adj_prior = NULL,
  bspline.control = list(knot.min = 10, knot.adj = 4, fixed.knot = NULL, spline.degree =
    3, pad.spline = 0.01)
)
```

## Arguments

- df_ref:

  a data.frame with two named columns for the reference series data; 1)
  BY - vector of *known* formation (birth) years of the reference
  series; 2) C14 - vector of ∆14C values of the reference series

- df_unk:

  a data.frame with two named columns for samples with unknown true
  birth year; 1) BY - vector of estimated formation (birth) years; 2)
  C14 - vector of ∆14C values of the samples

- model:

  [character](https://rdrr.io/r/base/character.html) of
  'coupled_function','bspline'; defaults to bspline

- likelihood:

  [character](https://rdrr.io/r/base/character.html) of
  'normal','student_t'; defaults to normal

- coupled_function:

  [character](https://rdrr.io/r/base/character.html) of
  'exponential','linear'; defaults to exponential

- ll_wt:

  weighting value of the reference series relative to the sample values
  for the integrated method, default is 5 (tested to be sufficient
  unless there is many samples).

- pred.by:

  a vector of formation years or a named list object for predicting the
  reference series ∆14C at. The named list needs: min.by = start year of
  prediction sequence, max.by = end year of prediction sequence, inc.by
  = year increment of the prediction sequence.

- adj_prior:

  a vector of two values specifying a mean and sd of a normal
  distribution for the bias adjustment prior; only used with df_unk is
  provided otherwise ignored. Default is NULL which attempts to specify
  a weakly informative normal prior.

- bspline.control:

  a named list of B-spline control parameters (see Details)

## Value

A named list object containing:

- **flag** character string indicating model type

- **data** a named list object that matches the DATA section in the STAN
  model

## Details

The weakly informative prior of **adj_prior** is set with a mean of zero
and a standard deviation equal to the difference in range of test sample
birth years divided by 4.

- **knot.min** minimum number of knots, default is 10.

- **knot.adj** divisor of the number of observations to set the number
  of knots, set a default of 4 with number of knots =
  nrow(df_ref)/knot.adj. Increasing this value decreases the number of
  knots (more smoothing in the spline).

- **fixed.knot** (optional) set number of knots that overrides number of
  knots = nrow(df_ref)/knot.adj. Must be a integer followed by L to be
  considered.

- **spline.degree** degree of polynomial spline, set at a default of 3.
  Must be 0 or greater.

- **pad.spline** amount of years to pad spline knot locations, default
  is 0.01

## Examples

``` r
#default BY_pred, reference only
df <- data_prep(sim_ref)
#> [1] "B-spline model indicated, adding additional data values for spline control"
#> No values provided in df_unk, skipping validation and estimating reference series only

#default BY_pred, integrated model
df <- data_prep(sim_ref, sim_unk)
#> [1] "B-spline model indicated, adding additional data values for spline control"

#custom BY_pred
df <- data_prep(sim_ref, sim_unk, pred.by = c(1956, 1959, 1988, 1991))
#> [1] "B-spline model indicated, adding additional data values for spline control"

#custom BY_pred with sequence function
df <- data_prep(sim_ref, sim_unk, pred.by = list(min.by=1900, max.by=2020, inc.by=0.5))
#> [1] "B-spline model indicated, adding additional data values for spline control"
```
