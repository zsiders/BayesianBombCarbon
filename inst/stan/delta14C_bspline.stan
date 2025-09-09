functions {
	vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order);
	vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order) {
	// INPUTS:
	//    t:          the points at which the b_spline is calculated
	//    ext_knots:  the set of extended knots
	//    ind:        the index of the b_spline
	//    order:      the order of the b-spline
	vector[size(t)] b_spline;
	vector[size(t)] w1 = rep_vector(0, size(t));
	vector[size(t)] w2 = rep_vector(0, size(t));
	if(order==1){
		for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
		b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]); 
	}else{
		if (ext_knots[ind] != ext_knots[ind+order-1])
		w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / 
		     (ext_knots[ind+order-1] - ext_knots[ind]);
		if (ext_knots[ind+1] != ext_knots[ind+order])
		w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) / 
		         (ext_knots[ind+order] - ext_knots[ind+1]);
		// Calculating the B-spline recursively as linear interpolation of two lower-order splines 
		b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) + 
		         w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
	}
	return b_spline;
	}
}
data {
	int<lower=1> Nref;             // number of data points
	int<lower=1> Nk;            // num of knots
	int<lower=1> Np;			// num of predictions
	vector[Nk] knots;  // the sequence of knots
	int spline_degree;        // the degree of spline (is equal to order - 1)
	array[Nref] real C14_ref;
	array[Nref]real BY_ref;
	array[Np] real BY_pred;
}
transformed data {
	int num_basis = Nk + spline_degree - 1; // total number of B-splines
	matrix[num_basis, Nref] B;  // matrix of B-splines
	vector[spline_degree + Nk] ext_knots_temp;
	vector[2*spline_degree + Nk] ext_knots; // set of extended knots
	ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
	ext_knots = append_row(ext_knots_temp, rep_vector(knots[Nk], spline_degree));
	for (ii in 1:num_basis){
		B[ii,:] = to_row_vector(build_b_spline(BY_ref, to_array_1d(ext_knots), ii, spline_degree + 1));
	}
	
	B[Nk + spline_degree - 1, Nref] = 1; 
}
parameters {
	row_vector[num_basis] a_raw; 
	real a0;  // intercept
	real<lower=0> sigma_ref; 
	real<lower=0> tau;   
}
transformed parameters {
	row_vector[num_basis] a;
	vector[Nref] C14_ref_hat;
	a[1] = a_raw[1];
	for (i in 2:num_basis){
		a[i] = a[i-1] + a_raw[i]*tau; 
	}
	C14_ref_hat = a0*to_vector(BY_ref) + to_vector(a*B);
}
model {
	// Priors
	a_raw ~ std_normal();
	a0 ~ std_normal();
	tau ~ std_normal();
	sigma_ref ~ std_normal();

	//Likelihood
	C14_ref ~ normal(C14_ref_hat, sigma_ref);
}
generated quantities{
	vector[Nref] log_lik_ref;
	vector[Np] C14_pred;
	vector[Np] C14_pred_hat;
	matrix[num_basis, Np] B_pred;  // matrix of B-splines

	//likelihood
	for(ii in 1:Nref){
		log_lik_ref[ii] = normal_lpdf(C14_ref[ii]|C14_ref_hat[ii],sigma_ref);
	}

	//predictions
	for (ii in 1:num_basis){
		B_pred[ii,:] = to_row_vector(build_b_spline(BY_pred, to_array_1d(ext_knots), ii, spline_degree + 1));
	}
	B_pred[num_basis, Np] = 1;
	C14_pred_hat = a0*to_vector(BY_pred) + to_vector(a*B_pred);
	for(ii in 1:Np){
		C14_pred[ii] = normal_rng(C14_pred_hat[ii], sigma_ref);
	}
}
