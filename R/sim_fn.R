#' @examples
#' exp.decay <- sim_fn(seq(1940,2020, length.out=125), mode='exp',plot=T)
#' ln.decay <- sim_fn(seq(1940,2020, length.out=125),mode='linear', gamma=3, plot=T)
#' wig.decay <- sim_fn(seq(1940,2020, length.out=125),mode='wiggle', gamma=3, A=6, damp=0.01, omega=0.4, plot=T)

sim_fn <- function(BY,mode='exp',rand_BY=FALSE, BY_n=50, lambda=-50, kappa=250, mu=1959, gamma=0.018, sigmaN=4, sigmaE=5, sigmaO=10, A = 1, damp = 0.01, omega=0.2, plot=FALSE){
	#lambda is pre-bomb delta 14C 
	#kappa is asymptote delta 14C (logistic model)
	#mu is peak year
	#gamma is exponential decay rate
	#sigmaN is gaussian SD of pulse
	#sigmaE is observation error in reference
	#sigmaO is observation error in samples
	if(rand_BY){
		BY <- sort(rnorm(BY_n,seq(min(BY),max(BY),length.out=BY_n),2))
	}
	if(mode=='exp'){
		cdf_mu <- mu+(sigmaN^2)*gamma
		mean <- lambda + 
				kappa * exp((mu*gamma) + (sigmaN^2*gamma^2)/2) *
				exp(-gamma*BY) * pnorm(BY, cdf_mu, sigmaN)
		rand <- rnorm(length(BY),mean,sigmaE)

		df <- data.frame(BY=BY,
		                 mean = mean,
		                 rng = rand)
	}else if(mode=='linear'){
		cdf_mu <- mu
		mean <- lambda + 
				kappa*pnorm(BY, cdf_mu, sigmaN) + 
				pnorm(BY, cdf_mu, sigmaN)*(gamma*(mu-BY))
		rand <- rnorm(length(BY),mean,sigmaE)

		df <- data.frame(BY=BY,
		                 mean = mean,
		                 rng = rand)
	}else if(mode=='wiggle'){
		cdf_mu <- mu
		mean <- lambda + 
				kappa*pnorm(BY, cdf_mu, sigmaN) + 
				pnorm(BY, cdf_mu, sigmaN)*(gamma*(mu-BY)) + 
				(A * exp(damp*(mu-BY))*
				 (cos(omega*(mu-BY))+sin(omega*(mu-BY)))) *
				pnorm(BY, cdf_mu, sigmaN)
		rand <- rnorm(length(BY),mean,sigmaE)
		df <- data.frame(BY=BY,
		                 mean = mean,
		                 rng = rand)
	}
	
	if(plot){
		plot(df$BY, df$rng, pch=16,
		     xlab = "Year",
		     ylab = expression(paste(Delta^14,"C")))
		df.ord <- df[order(df$BY),]
		lines(df.ord$BY, df.ord$mean, lwd=2)
	}
	return(df)
}