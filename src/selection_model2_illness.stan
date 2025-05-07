// Define data
data {
  int N_miss;     // number of missing (-99) values in f
  int index_miss[N_miss]; // indices of missing (-99) values
  
  int N_subset;          // Number of elements in the subset
  int subset_indices[N_subset]; // Indices of elements in the subset
  
  int F;          // Number of households
  int L;          // Number of locations
  int N;          // Number of rows in data
  int K;          // Number of predictors
  int<lower = 1, upper = F> hh[N];      // household identifiers
  int<lower = 1, upper = L> loc[N];     // location identifiers
  
  matrix[N, K] X; // Matrix of factual predictors
  matrix[N, K] X_cntr; // Matrix of counterfactual predictors
  int R[N]; // binary missingness indicator vector
  vector[N] f;    // Outcome variable: -99 represents missing value
  
  //EMPIRICAL BAYES
  vector[2] phi_prior_mean;
  matrix[2,2] phi_prior_var;
  real kappa;
}

// Define parameters
parameters {
  // DATA MODEL
  vector[F] alpha_raw;          // Random intercept -noncentered param
  vector[L] tau_raw;            // random intercept for location
  vector[N_miss] f_imp;         // vector of imputed/missing f values
  vector[K] beta;               // Coefficients on X
  //vector<lower = 0>[K] sigma_beta;
  //real<lower = 0> tau_beta;
  
  vector<lower = 0>[N] u;       // "Inefficiency" term
  real<lower = 0> sigma_v;      // Standard deviation on noise term
  
  // RANDOM GROUP MODEL
  real<lower = 0> sigma_a;      // household standard deviation
  real mu_a;                    // intercept 
  real<lower = 0> sigma_l;      // location standard deviation 

  // ASYMMETRIC TERM MODEL
  real gamma1;                    // Coefficient in scale of asymmetric 
  vector[7] gamma2;                    // Coefficients in scale of asymmetric 
  real<lower = 0, upper = 1> rho; // Probability of asymmetric inclusion
  
  // MISSINGNESS MODEL
  vector[2] phi;                      // pg 305 Ma, Chen

  //vector[K] phi2; //taking out for now...
  
}

transformed parameters{
  vector[F] alpha;                // Random intercept -noncentered param
  vector[L] tau; // random intercept for location
  vector<lower=0>[N] lambda;
  vector[N] Data; // will contain observed, imputed quantities
  
  //replace missing values (-99) with imputed/parameters
  Data = f;
  Data[index_miss] = f_imp;   
  // noncentered alpha
  alpha = mu_a + sigma_a * alpha_raw;
  tau = sigma_l*tau_raw; 
  lambda =      exp(gamma1 + X[,1:7]*gamma2);
  
}

// Define model
model {
  
  alpha_raw ~ normal(0,1);
  tau_raw ~ normal(0,1);
  mu_a ~ normal(0,2.5);
  
  sigma_v ~ inv_gamma(20, .1); 
  sigma_l ~ inv_gamma(20, .1); 
  sigma_a ~ inv_gamma(20, .1); 

  //sigma_v ~ normal(0,1);
  //sigma_l ~ normal(0,1);
  
  //tried horseshoe prior - shrinks most coefficients pretty rigorously - given our low dimension and goal of quantifying effects, this is probably not that useful
  // go back to standard ridge prior
  //issue: assumes independence among coefficients
  //this is a problem if we believe in heredity principals
  // e.g.)Hierarchical Shrinkage Priors for Regression Models by Griffin, Brown
  beta ~ normal(0, 1);
 // sigma_beta ~ cauchy(0,1);
  //tau_beta ~ cauchy(0,1);
  
  gamma1 ~ normal(0,1);
  gamma2 ~ normal(0,1);
  rho ~ beta(2,4);
  //EMPIRICAL BAYES based prior
  phi ~ multi_normal(phi_prior_mean,kappa*phi_prior_var);
      u ~ exponential(1 ./ lambda);
      R ~ bernoulli_logit(phi[1] + phi[2]*Data);
  // Likelihood
  for (i in 1:N) {
    real mu = alpha[hh[i]] + tau[loc[i]] + X[i] * beta;
    target += log_mix(
      rho,
      normal_lpdf(Data[i] | mu - u[i], sigma_v),
      normal_lpdf(Data[i] | mu, sigma_v)
    );
  }
}

// Generate predicted capabilities
generated quantities {
    vector[N_subset] c;
    vector[N_subset] c_cntr;
    vector[N_subset] f_pred;
    vector[N_subset] f_pred_cntr;
    real u_cntr[N_subset];


    for (n in 1:N_subset) {
        int i = subset_indices[n];
        real alpha_tau = alpha[hh[i]]+tau[loc[i]];
        c[n] = alpha_tau + X[i] * beta; // factual potential
        c_cntr[n] = alpha_tau + X_cntr[i] * beta; // counterfactual potential
        u_cntr[n] = exponential_rng(1. / exp(gamma1 + X_cntr[i, 1:7] * gamma2)); // counterfactual deviation

        f_pred[n] = normal_rng(c[n] - bernoulli_rng(rho) * u[i], sigma_v);
        f_pred_cntr[n] = normal_rng(c_cntr[n] - bernoulli_rng(rho) * u_cntr[n], sigma_v);
    }
}
  