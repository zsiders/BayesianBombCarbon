data{
	int<lower=1> Nref;             // number of data points
	int<lower=1> Nobs;             // number of data points
	int<lower=1> Np;			// num of predictions
	int<lower=1> wt; //wt of ref LL
	array[Nref] real C14_ref; //C14 from reference points
	array[Nref] real BY_ref; //Birth year from reference points
	array[Nobs] real C14_obs; //C14 from observer points
	array[Nobs] real BY_obs; //Birth year from observed points
	//for predictions
	array[Np] real BY_pred; //Birth year to predict to
	vector[2] lambda_prior;
	vector[2] mu_prior;
	vector[2] kappa_prior;
	vector[2] gamma_prior;
	vector[2] sigmaN_prior;
	vector[2] sigma_ref_prior;
	vector[2] sigma_obs_prior;
	vector[2] adj_prior;
}
parameters{
	real lambda; //pre-bomb delta 14C
	real<lower=0> kappa; //asymptote of delta 14C
	real mu; //peak year of delta 14C
	real<lower=0> gamma; //exponential decay rate post-peak decline
	real<lower=0> sigmaN; //Gaussian SD of pulse
	real<lower=0> sigma_ref; //variability of the reference
	real<lower=0> sigma_obs; //variability of the observations

	real adj; //birth year adjustment
}
transformed parameters{
	vector[Nref] C14_ref_hat; //predicted reference
	vector[Nobs] C14_obs_hat;  //predicted observations
	array[Nobs] real BY_adj; //repeated adjustment
	BY_adj = to_array_1d(to_vector(BY_obs) + rep_vector(adj,Nobs));
	for(i in 1:Nref){
		C14_ref_hat[i] = lambda+kappa *normal_cdf(BY_ref[i] | mu, sigmaN) + (gamma*(mu-BY_ref[i]))*normal_cdf(BY_ref[i] | mu, sigmaN);
	}
	for(i in 1:Nobs){
		C14_obs_hat[i] = lambda+kappa *normal_cdf(BY_adj[i] | mu, sigmaN) + (gamma*(mu-BY_adj[i]))*normal_cdf(BY_adj[i] | mu, sigmaN);
	}
}
model{
	//priors
		lambda ~ normal(lambda_prior[1],lambda_prior[2]); 
		mu ~ normal(mu_prior[1],mu_prior[2]);
		kappa ~ normal(kappa_prior[1],kappa_prior[2]);
		gamma ~ normal(gamma_prior[1],gamma_prior[2]);
		sigmaN ~ normal(sigmaN_prior[1],sigmaN_prior[2]);
		sigma_ref ~ normal(sigma_ref_prior[1],sigma_ref_prior[2]);
		sigma_obs ~ normal(sigma_obs_prior[1],sigma_obs_prior[2]);
		adj ~ normal(adj_prior[1], adj_prior[2]);
	//model
		target += normal_lpdf(C14_ref|C14_ref_hat, sigma_ref)*wt; //reference
		target += normal_lpdf(C14_obs|C14_obs_hat, sigma_obs); //observations
}
generated quantities{
	vector[Nref] log_lik_ref;
	vector[Nobs] log_lik_obs;
	vector[Np] C14_pred; //predicted C14 including variablity
	vector[Np] C14_pred_hat; //predicted mean C14

	//likelihood
	for(ii in 1:Nref){
		log_lik_ref[ii] = normal_lpdf(C14_ref[ii]|C14_ref_hat[ii],sigma_ref);
	}
	for(ii in 1:Nobs){
		log_lik_obs[ii] = normal_lpdf(C14_obs[ii]|C14_obs_hat[ii],sigma_obs);
	}

	for(i in 1:Np){
		C14_pred_hat[i] = lambda+kappa *normal_cdf(BY_pred[i] | mu, sigmaN) + (gamma*(mu-BY_pred[i]))*normal_cdf(BY_pred[i] | mu, sigmaN);
		C14_pred[i] = normal_rng(C14_pred_hat[i],	sigma_ref);
	}
}
