//gaussian model with global intercept and kernel covariate that doesn't incorporate species
data{
  int<lower=1> N; //number of plants
  vector[N] y;//response variable (size of plant)
  int obs; //number of non-NA pairwise distances
  vector[obs] dist_vector; //pairwise distances^2
  int n_nonNA[N]; //number of neighbors for each plant
  int pos[N]; //position demarcating where each plant starts in dist_vector 
}
parameters{
  real a0; 
  real<lower=0> sigma;
  real<lower=0> eta;
  real<lower=0> rho;
}
transformed parameters{ 
  vector[N] mu;
  vector[N] NHI;
  //NHI kernel
  for(i in 1:N){
    NHI[i] = sum(eta^2 * exp(
      -segment(dist_vector, pos[i], n_nonNA[i]) /
      (2 * rho^2)));
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
