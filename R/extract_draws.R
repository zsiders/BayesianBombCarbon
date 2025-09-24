#' Function for extracting Bayesian B-spline posteriors
#' 
#' @description wrapper function for extracting posterior samples from the Bayesian B-spline model
#' @param fit cmdstan fit object of a Bayesian B-spline model
#' @param df data prepared using the data_prep function used to fit the model
#' @param pars list of parameters to extract. defaults to c('adj','sigma_ref','sigma_obs','tau','C14_pred')
#' 
#' @return a list containing draws_matrix of posterior draws for parameters of interest; parameters return depend on the model fit
 
#' @export
#' 
#' @examples
#' #reference-only
#' df <- data_prep(sim_ref)
#' 
#' #estimate model
#' fit <- est_model(df, show_messages = FALSE, show_exceptions = FALSE)
#' 
#' #extract draws
#' draws <- extract_draws(fit)
#' 
extract_draws <- function(fit, df, pars = c('adj','sigma_ref','sigma_obs','tau','C14_pred')){
	fit_vars <- unique(gsub("\\[[[:digit:]]+\\]|\\[[[:digit:]]+\\,[[:digit:]]+\\]","",fit$fitted$metadata()$variables))
	pars <- pars[pars %in% fit_vars]
	bspline_ext <- extract(fit = fit$fitted,
                           pars = pars)
	if('adj' %in% pars & !missing(df)){
		bspline_ext$BY_adj <- matrix(bspline_ext$adj,
                                nrow = nrow(bspline_ext$adj),
                                ncol = df$data$Nobs) + 
						  matrix(df$data$BY_obs,
						         nrow = nrow(bspline_ext$adj),
						         ncol = df$data$Nobs,byrow=TRUE)
	}
	
	return(bspline_ext)
}	
extract <- function(fit, pars){
	ext <- lapply(pars, function(x)fit$draws(variables=x, format='draws_matrix'))
	names(ext) <- pars
	ext
}