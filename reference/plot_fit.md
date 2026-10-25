# Plotting Bayesian B-spline model

function for plotting Bayesian B-spline model for estimating reference
series and/or assessing aging bias

## Usage

``` r
plot_fit(
  df,
  ext,
  probs = c(0.05, 0.5, 0.95),
  post.den = TRUE,
  legend = TRUE,
  min.BY,
  max.BY,
  inset.control = list(width = 0.003, height = 0.125, sig_low = 0.3, adj_low = 0.05,
    xbump = 0, lab_adj = 0.18),
  plot_bias = TRUE,
  legend.lab = c("Ref.", "Obs.")
)
```

## Arguments

- df:

  data list object from `data_prep` and used to estimate the STAN model

- ext:

  posterior draws extraction list from `extract_draws` function

- probs:

  vector of three probabilities for lower CI, centrality, upper CI

- post.den:

  logical flag to plot posterior densities of certain parameters

- legend:

  logical flag to plot legend

- inset.control:

  a list with values of width (width of the inset posterior in
  proportions of xlim), height (height of the inset posterior in
  proportions of ylim), sig_low (lower position of sigma posterior in
  proportions of ylim), adj_low (lower position of BY_adj posterior in
  proportions of ylim), xbump (translation of posteriors along x-axis in
  proportions of xlim).

- plot_bias:

  logical flag whether to plot the formation year adjustment (only for
  models with test samples)

- legend.lab:

  labels for the reference and test observations

## Value

figure only, invisible

## Examples

``` r
#REFERENCE-ONLY
df <- data_prep(sim_ref)
#> [1] "B-spline model indicated, adding additional data values for spline control"
#> No values provided in df_unk, skipping validation and estimating reference series only

#estimate model
fit <- est_model(df, show_messages = FALSE, show_exceptions = FALSE)
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.

#extract draws
draws <- extract_draws(fit)
#> Error: object 'fit' not found

#plot fit
plot_fit(df, draws)
#> Error: object 'draws' not found

#INTEGRATED
if (FALSE) { # \dontrun{
df_int <- data_prep(sim_ref, sim_unk)
fit_int <- est_model(df_int)
draws_int <- extract_draws(fit_int)
plot_fit(df_int, draws_int)
} # }
```
