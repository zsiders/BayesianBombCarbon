#' Prepares data for Bayesian penalized B-spline
#' 
#' @description Forces a standardized format for inputting bomb radidocarbon data into a list object for estimating with 'cmdstanr'
#' 
#' @param df_ref a data.frame with two named columns for the reference series data; 1) BY - vector of *known* formation (birth) years of the reference series; 2) C14 - vector of ∆14C values of the reference series
#' @param df_unk a data.frame with two named columns for samples with unknown true birth year; 1) BY - vector of estimated formation (birth) years; 2) C14 - vector of ∆14C values of the samples
#' @param knot.min minimum number of knots, default is 10.
#' @param knot.adj divisor of the number of observations to set the number of knots, set a default of 4 with number of knots = nrow(df_ref)/knot.adj. Increasing this value decreases the number of knots (more smoothing in the spline).
#' @param fixed.knot (optional) set number of knots that overrides number of knots = nrow(df_ref)/knot.adj. Must be a integer followed by L to be considered.
#' @param spline.degree degree of polynomial spline, set at a default of 3. Must be 0 or greater.
#' @param pad.spline amount of years to pad spline knot locations, default is 0.01
#' @param ll_wt weighting value of the reference series relative to the sample values for the integrated method, default is 5 (tested to be sufficient unless there is many samples).
#' @param pred.by a vector of formation years or a named list object for predicting the reference series ∆14C at. The named list needs: min.by = start year of prediction sequence, max.by = end year of prediction sequence, inc.by = year increment of the prediction sequence.
#' 
#' @return A named list object containing:
#' \itemize{
#' \item \strong{flag} character string indicating model type
#' \item \strong{data} a named list object that matches the DATA section in the STAN model
#' }
#' 
#' @export
#' 
#' @examples
#' #default BY_pred, reference only
#' df <- data_prep(sim_ref)
#' 
#' #default BY_pred, integrated model
#' df <- data_prep(sim_ref, sim_unk)
#' 
#' #custom BY_pred
#' df <- data_prep(sim_ref, sim_unk, pred.by = c(1956, 1959, 1988, 1991))
#' 
#' #custom BY_pred with sequence function
#' df <- data_prep(sim_ref, sim_unk, pred.by = list(min.by=1900, max.by=2020, inc.by=0.5))

data_prep <- function(df_ref, df_unk, knot.min = 10, knot.adj = 4, fixed.knot, spline.degree = 3, pad.spline=0.01, ll_wt = 5, pred.by=list(min.by=1940, max.by=2020, inc.by=1)){
	if(any(!c('BY','C14') %in% colnames(df_ref))) stop('Missing BY or C14 named column in df_ref, check column names')
	if(any(is.na(df_ref[,c('BY','C14')]))){
		df_ref <- na.omit(df_ref[,c('BY','C14')])
		warning('df_ref contains NAs, using na.omit to remove')
	}
	v <- list(flag = NULL, data = NULL)
	v$data <- list(BY_mu = mean(df_ref$BY)) #mean birth year
	if(is.list(pred.by)){
		pred_seq <- seq(pred.by$min.by,pred.by$max.by,pred.by$inc.by) #predicted sequence
	}else if(is.vector(pred.by)){
		pred_seq <- pred.by
	}else{
		warning('pred.by format is invalid; Setting reference series prediction years to range of sampled formation years and increment of one year')
		pred_seq <- seq(floor(min(df_ref$BY)),ceiling(max(df_ref$BY)),by=1)
	}
	
	v$data <- c(v$data,list(Nref = nrow(df_ref), #number of observations
		     BY_ref = df_ref$BY-v$data$BY_mu, #mean-centered formation year
		     C14_ref = df_ref$C14, #∆14C reference values
		     Np = length(pred_seq), #length of prediction series
		     BY_pred = pred_seq-v$data$BY_mu)) #mean-centered prediction series
	v$data$Nk <- pmax(knot.min,floor(v$data$Nref/knot.adj)) #can adjust knots with this
	if(!missing(fixed.knot)){
		if(is.integer(fixed.knot)){
			v$data$Nk <- fixed.knot
		}else if(!is.integer(fixed.knot)){
			warning('Tried to supply non-integer fixed.knot value; reverting to knot.adj calculation for number of knots')
		}
	}
	dBY <- diff(range(v$data$BY_ref)) #difference in range of sampled formation years
	if(pad.spline > 0.25){
		warning('User input of pad.spline > 0.25, exceeding recommended tolerance')
	}
	#start and end points of knot locations
	st_point <- dBY * -pad.spline + min(v$data$BY_ref)
	end_point <- dBY * pad.spline + max(v$data$BY_ref)
	#knot locations
	v$data$knots <- unname(quantile(c(st_point,v$data$BY_ref,end_point),
								probs=seq(from=0, to=1,
								          length.out = v$data$Nk)
								))
	v$data$spline_degree <- pmax(0,3)

	#validation samples
	if(!missing(df_unk)){
		if(any(!c('BY','C14') %in% colnames(df_unk))){
			warning('Missing BY or C14 named column in df_unk, check column names, skipping unknown FY/BY samples')
		}else{
			if(any(is.na(df_unk[,c('BY','C14')]))){
				df_unk <- na.omit(df_unk[,c('BY','C14')])
				warning('df_unk contains NAs, using na.omit to remove')
			}
			v$data$Nobs <- nrow(df_unk) #number of observations
			v$data$BY_obs <- df_unk$BY - v$data$BY_mu #estimated birth year
			v$data$C14_obs <- df_unk$C14 #measured ∆14C
			v$data$wt <- ll_wt #weight
			v$flag <- 'integrated'
		}
	}else{
		message('No values provided in df_unk, skipping validation and estimating reference series only')
		v$flag <- 'reference-only'
	}

	return(v)
}