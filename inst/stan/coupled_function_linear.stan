data{
	int<lower=1> Nref;             // number of data points
	int<lower=1> Np;			// num of predictions
	array[Nref] real C14_ref;
	array[Nref]real BY_ref;
	array[Np] real BY_pred;
	vector[2] lambda_prior;
	vector[2] mu_prior;
	vector[2] kappa_prior;
	vector[2] gamma_prior;
	vector[2] sigmaN_prior;
	vector[2] sigma_ref_prior;
}
parameters{
	real lambda; //pre-bomb delta 14C
	real<lower=0> kappa; //asymptote of delta 14C
	real mu; //peak year of delta 14C
	real<lower=0> gamma; //exponential decay rate post-peak decline
	real<lower=0> sigmaN; //Gaussian SD of pulse
	real<lower=0> sigma_ref; //observation error
}
transformed parameters{
	vector[Nref] C14_ref_hat;

	for(i in 1:Nref){
		C14_ref_hat[i] = lambda+kappa * normal_cdf(BY_ref[i] | mu, sigmaN) +(gamma*(mu-BY_ref[i]))* normal_cdf(BY_ref[i] | mu, sigmaN);
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
	//model
		C14_ref ~ normal(C14_ref_hat, sigma_ref);
}
generated quantities{
	vector[Nref] log_lik_ref;
	vector[Np] C14_pred;
	vector[Np] C14_pred_hat;

	//likelihood
	for(ii in 1:Nref){
		log_lik_ref[ii] = normal_lpdf(C14_ref[ii]|C14_ref_hat[ii],sigma_ref);
	}

	//predictions
	for(ii in 1:Np){
		C14_pred_hat[ii] = lambda+kappa * normal_cdf(BY_pred[ii] | mu, sigmaN) +(gamma*(mu-BY_pred[ii]))* normal_cdf(BY_pred[ii] | mu, sigmaN);
		C14_pred[ii] = normal_rng(C14_pred_hat[ii], sigma_ref);
	}
}
