# Bayesian Negative Binomial regression for single-cell RNA-seq

data {
	int<lower=0> N;					# number of outcomes
	int<lower=0> K;					# number of predictors
	matrix<lower=0>[N,K] x;	# predictor matrix 
	int y[N];						    # outcomes
}

parameters {
	vector<lower=1>[K] beta;	# coefficients
	real<lower=0.001> r;			# overdispersion
}

model {	
	vector<lower=0.001>[N] mu;
	vector<lower=1.001>[N] rv;

	# priors
	r ~ cauchy(0, 1);
	beta ~ pareto(1, 1.5);

	# vectorize the scale
	for (n in 1:N) {
		rv[n] <- square(r + 1) - 1;
	}

	# regression
	mu <- x * (beta - 1) + 0.001;
	y ~ neg_binomial(mu ./ rv, 1 / rv[1]);
}
