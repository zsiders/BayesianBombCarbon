# Wrapper for model estimation based on data_prep flag

this is a wrapper function around compiling and sampling from the
Bayesian bomb radiocarbon model sets. Flags generated in the data_prep
function are used to compile the appropriate model specific to your
machine's C++ compiler and then sample from the model.

## Usage

``` r
est_model(data, save_dir, ...)
```

## Arguments

- data:

  data prepared using the data_prep function

- save_dir:

  directory to save cmdstan outputs

- ...:

  arguments to 'cmdstanr' sample function

## Value

A list containing the STAN model and samples if fitted successfully.

- **model**: the CmdStanModel object used to sample from

- **fitted**: the returned samples from calling `sample` on the
  CmdStanModel object

## Examples

``` r
#REFERENCE-ONLY
df <- data_prep(sim_ref)
#> [1] "B-spline model indicated, adding additional data values for spline control"
#> No values provided in df_unk, skipping validation and estimating reference series only

#estimate model
fit <- est_model(df, show_messages = FALSE, show_exceptions = FALSE)
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.

#estimate model with args to cmdstanr sample function
if (FALSE) { # \dontrun{
fit <- est_model(df, iter_warmup = 3000, iter_sampling = 500, parallel_chains = 4)
} # }

#INTEGRATED MODEL
df_int <- data_prep(sim_ref, sim_unk)
#> [1] "B-spline model indicated, adding additional data values for spline control"

#estimate model
if (FALSE) { # \dontrun{
fit_int <- est_model(df_int, parallel_chains = 4, iter_warmup=3000, iter_sampling = 250)
} # }

#FIXED NUMBER OF KNOTS
if (FALSE) { # \dontrun{
df <- data_prep(sim_ref, fixed.knot = 12L)
fit <- est_model(df, parallel_chains = 4, iter_warmup=1000, iter_sampling = 250)
} # }
```
