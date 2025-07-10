// Define input data
data {
  int N_miss;                      // number of missing (-99) values in f
  int index_miss[N_miss];          // indices of missing (-99) values
  int F;                           // Number of households
  int L;                           // Number of locations
  int N;                           // Number of rows in data
  int K;                           // Number of predictors
  int<lower = 1, upper = F> hh[N]; // household identifiers
  int<lower = 1, upper = L> loc[N];// location identifiers
  matrix[N, K] X;                  // Matrix of factual predictors
  matrix[N, K] X_cntr;             // Matrix of counterfactual predictors
  int R[N];                        // binary missingness indicator vector
  vector[N] f;                     // Outcome variable: -99 represents missing value
  
  //EMPIRICAL BAYES
  vector[2] phi_prior_mean;
  matrix[2,2] phi_prior_var;
  real kappa;                //scaling term on empirical Bayes prior
}

// Define parameters
parameters {
  // DATA MODEL
  vector[F] alpha_raw;          // Random intercept -noncentered param
  vector[L] tau_raw;            // random intercept for location
  vector[N_miss] f_imp;         // vector of imputed/missing f values
  vector[K] beta;               // Coefficients on X
  vector<lower = 0>[N] u;       // "Inefficiency" term
  real<lower = 0> sigma_v;      // Standard deviation on noise term
  
  // RANDOM GROUP MODEL
  real<lower = 0> sigma_a;      // household standard deviation
  real mu_a;                    // intercept 
  real<lower = 0> sigma_l;      // location standard deviation 

  // ASYMMETRIC TERM MODEL
  real gamma1;                    // Coefficient in scale of asymmetric 
  vector[7] gamma2;               // Coefficients in scale of asymmetric 
  real<lower = 0, upper = 1> rho; // Probability of asymmetric inclusion
  
  // MISSINGNESS MODEL
  vector[2] phi;                    
  
}

transformed parameters{
  vector[F] alpha;                // random intercept -noncentered param
  vector[L] tau;                  // random intercept for location
  vector<lower=0>[N] lambda;
  vector[N] Data;                 // will contain observed and imputed quantities
  vector[N] log_lik;              // for calculation of LOO-based ELPD
  vector[N] alphataubeta;
  
  //replace missing values (-99) with imputed/parameters
  Data = f;
  Data[index_miss] = f_imp;
  
  // random household, location effects
  alpha = mu_a + sigma_a * alpha_raw;
  tau = sigma_l*tau_raw; 
  
  // exponential scale term
  lambda = exp(gamma1 + X[,1:7]*gamma2);
  
  //data model
  for (i in 1:N){
  alphataubeta[i] = alpha[hh[i]] +tau[loc[i]] + X[i]*beta;
  log_lik[i] = log_mix(rho,
            normal_lpdf(Data[i] |  alphataubeta[i] - u[i], sigma_v),
            normal_lpdf(Data[i] |  alphataubeta[i],        sigma_v));
  } 
}

// Define model
model {
  
  alpha_raw ~ normal(0,1);
  tau_raw ~ normal(0,1);
  mu_a ~ normal(0,5);
  
  sigma_a ~ normal(0,2.5);
  sigma_v ~ normal(0,2.5);
  sigma_l ~ normal(0,2.5);
  
  beta ~ normal(0, 2.5);
  
  gamma1 ~ normal(0,1);
  gamma2 ~ normal(0,1);
  rho ~ beta(1,1);
  
  //EMPIRICAL BAYES based prior
  phi ~ multi_normal(phi_prior_mean, kappa*phi_prior_var);
  
  //distribution on deviation from potential
  u ~ exponential(1 ./ lambda);
  
  //missingness model
  R ~ bernoulli_logit(phi[1] + phi[2]*Data);
  
  //log likelihood calculation
   target += sum(log_lik);

}

// Generate predicted capabilities
generated quantities {
  vector[N] c;
  vector[N] c_cntr;
  real u_cntr[N];
  real f_pred[N];
  real f_pred_cntr[N];
  real alphatau[N];

    for (i in 1:N){
    alphatau[i] = alpha[hh[i]] + tau[loc[i]];
    c[i] =      alphatau[i] + X[i]*beta; //factual potential
    c_cntr[i] = alphatau[i] + X_cntr[i]*beta; //counterfactual potential
    u_cntr[i] = exponential_rng(1 ./ exp(gamma1 + X_cntr[i, 1:7]*gamma2)); //counterfactual deviation from potential
    f_pred[i] =      normal_rng(c[i]-      bernoulli_rng(rho)*u[i],      sigma_v); //factual functioning
    f_pred_cntr[i] = normal_rng(c_cntr[i]- bernoulli_rng(rho)*u_cntr[i], sigma_v); //counterfactual functioning
    }
  
  }
  