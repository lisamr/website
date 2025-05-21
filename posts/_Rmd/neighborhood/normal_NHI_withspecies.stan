//gaussian model with global intercept and kernel covariate that does incorporate species
data{
  int<lower=1> N; //number of plants
  int<lower=1> S; //number of species
  vector[N] y;//response variable (size of plant)
  int obs; //number of non-NA pairwise distances
  int<lower=1, upper=S> spID_vector[obs]; //species identity of plants
  vector[obs] dist_vector; //pairwise distances^2
  int n_nonNA[N]; //number of neighbors for each plant
  int pos[N]; //position demarcating where each plant starts in dist_vector 
}
parameters{
  real a0; 
  real<lower=0> sigma; 
  vector<lower=0>[S] eta; 
  real<lower=0> rho; 
}
transformed parameters{ 
  vector[N] mu; 
  vector[N] NHI; 
  vector<lower=0>[S] eta_sq; 
  //NHI kernel
  for(s in 1:S){ //only real numbers, not vectors can be squared or used with power function
    eta_sq[s] = pow(eta[s], 2); 
  }
  for(i in 1:N){
    NHI[i] = sum(eta_sq[segment(spID_vector, pos[i], n_nonNA[i])] .* exp(-segment(dist_vector, pos[i], n_nonNA[i]) /
      (2 * rho^2))); // vector multiplication uses .* while division with real numbers (i.e. rho) just needs /
  }
  //main model
  mu = a0 - NHI; 
}
model{
  //priors
  a0 ~ normal(0, 5);
  eta ~ normal(0, 3);
  rho ~ normal(0, 5);
  sigma ~ exponential(1);
  
  //likelihood
  y ~ normal(mu, sigma);
}
