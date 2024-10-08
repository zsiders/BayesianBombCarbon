% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_model.R
\name{est_model}
\alias{est_model}
\title{Wrapper for model estimation based on data_prep flag}
\usage{
est_model(data, save_dir, ...)
}
\arguments{
\item{data}{data prepared using the data_prep function}

\item{save_dir}{directory to save cmdstan outputs}

\item{...}{arguments to 'cmdstanr' sample function}
}
\value{
A list containing the STAN model and samples if fitted successfully.
\itemize{
\item \strong{model}: the CmdStanModel object used to sample from
\item \strong{fitted}: the returned samples from calling \code{sample} on the CmdStanModel object
}
}
\description{
this is a wrapper function around compiling and sampling from the Bayesian bomb radiocarbon model sets. Flags generated in the data_prep function are used to compile the appropriate model specific to your machine's C++ compiler and then sample from the model.
}
\examples{
#REFERENCE-ONLY
df <- data_prep(sim_ref)

#estimate model
fit <- est_model(df, show_messages = FALSE, show_exceptions = FALSE)

#estimate model with args to cmdstanr sample function
\dontrun{
fit <- est_model(df, iter_warmup = 3000, iter_sampling = 500, parallel_chains = 4)
}

#INTEGRATED MODEL
df_int <- data_prep(sim_ref, sim_unk)

#estimate model
\dontrun{
fit_int <- est_model(df_int, parallel_chains = 4, iter_warmup=3000, iter_sampling = 250)
}

#FIXED NUMBER OF KNOTS
\dontrun{
df <- data_prep(sim_ref, fixed.knot = 12L)
fit <- est_model(df, parallel_chains = 4, iter_warmup=1000, iter_sampling = 250)
}
}
