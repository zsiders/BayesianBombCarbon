#' Wrapper for model estimation based on data_prep flag
#' 
#' @description this is a wrapper function around compiling and sampling from the Bayesian bomb radiocarbon model sets. Flags generated in the data_prep function are used to compile the appropriate model specific to your machine's C++ compiler and then sample from the model.
#' 
#' @param data data prepared using the data_prep function
#' @param save_dir directory to save cmdstan outputs
#' @param ... arguments to 'cmdstanr' sample function
#' 
#' @return A list containing the STAN model and samples if fitted successfully. 
#' \itemize{
#' 	\item \strong{model}: the CmdStanModel object used to sample from
#' \item \strong{fitted}: the returned samples from calling \code{sample} on the CmdStanModel object
#' }
#' 
#' @export
#' 
#' @examples
#' #REFERENCE-ONLY
#' df <- data_prep(sim_ref)
#' 
#' #estimate model
#' fit <- est_model(df, show_messages = FALSE, show_exceptions = FALSE)
#' 
#' #estimate model with args to cmdstanr sample function
#' \dontrun{
#' fit <- est_model(df, iter_warmup = 3000, iter_sampling = 500, parallel_chains = 4)
#' }
#' 
#' #INTEGRATED MODEL
#' df_int <- data_prep(sim_ref, sim_unk)
#' 
#' #estimate model
#' \dontrun{
#' fit_int <- est_model(df_int, parallel_chains = 4, iter_warmup=3000, iter_sampling = 250)
#' }
#' 
#' #FIXED NUMBER OF KNOTS
#' \dontrun{
#' df <- data_prep(sim_ref, fixed.knot = 12L)
#' fit <- est_model(df, parallel_chains = 4, iter_warmup=1000, iter_sampling = 250)
#' }

est_model <- function(data, save_dir, ...){
	if(data$flag == 'reference-only'){
		mod <- get_stan_model('ref-only')
	}else if(data$flag == 'integrated'){
		mod <- get_stan_model('integrated')
	}
	fitting_successful <- FALSE
	bspline_fit <- list()
	bspline_fit$model <- mod$load_model[[1]]()
	fit <- tryCatch(
		{
		  fit <- bspline_fit$model$sample(data = data$data,...)
		  if(length(fit$warnings) == 0) {
		    fitting_successful <- TRUE
		  }else{
		    cat("\n")
		    cli::cli_warn(
		      paste(
		        "There was an error while fitting the model.",
		        "Only the model input is returned."
		      )
		    )
		    fit <- list(
		      errors = unlist(lapply(
		        fit$warnings, function(x) stringr::str_remove(x$message, "\n")
		      )),
		      sampler_output = fit$value$output()
		    )
		  }
		  fit
		},
		error = function(err) {
		  cat("\n")
		  cli::cli_warn(c(
		    paste(
		      "There was an error while fitting the model.",
		      "Only the model input is returned."
		    ),
		    err$message
		  ))
		  return(list(errors = err, sampler_output = NULL))
		}
	)
	if(fitting_successful){
		bspline_fit$fitted <- fit
	}
	if(!missing(save_dir) & fitting_successful){
		bspline_fit$save_output_files(dir=save_dir)
	}
	return(bspline_fit)
}
	