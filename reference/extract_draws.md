# Function for extracting Bayesian B-spline posteriors

wrapper function for extracting posterior samples from the Bayesian
B-spline model

## Usage

``` r
extract_draws(
  fit,
  df,
  pars = c("adj", "sigma_ref", "sigma_obs", "tau", "nu", "C14_pred")
)
```

## Arguments

- fit:

  cmdstan fit object of a Bayesian B-spline model

- df:

  data prepared using the data_prep function used to fit the model

- pars:

  list of parameters to extract. defaults to
  c('adj','sigma_ref','sigma_obs','tau','C14_pred')

## Value

a list containing draws_matrix of posterior draws for parameters of
interest; parameters return depend on the model fit

## Examples

``` r
#reference-only
df <- data_prep(sim_ref)
#> [1] "B-spline model indicated, adding additional data values for spline control"
#> No values provided in df_unk, skipping validation and estimating reference series only

#estimate model
fit <- est_model(df, show_messages = FALSE, show_exceptions = FALSE)
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.

#extract draws
draws <- extract_draws(fit)
#> Error: object 'fit' not found
```
