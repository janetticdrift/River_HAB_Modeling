//This model include all abiotic and biotic factors

data {
  int uniqueID; //Total number of weeks down the years
  int firstdays[uniqueID]; //Days to skip modeling, first day of the year
  int Npredictors; //Number of predictor variables (includes species int + env effect)
  int<lower=0, upper=1> is_obs[uniqueID];
  
  int<lower=0> Toxins[uniqueID]; //Vector of raw toxin concentrations
  matrix[uniqueID, Npredictors] X2; //Design matrix of all predictors
}

parameters {
  
  real<lower= 0> sigma_p; //var w/ process model
  
  vector[uniqueID] tox;   //latent anatoxin state on log scale
  
  real Beta0;            
  real Beta1;            
  real Beta2;            
  real Beta3;            
  real Beta4;            
  
  // real BetaGreen;        
  // real BetaMicro;        
  real BetaAna;          
  // real BetaNFix;         
  
  real Ntheta;
  real Ptheta; 
  real Atheta;
  real Dtheta; 
  real Ttheta; 
  real Ctheta; 
  real Rtheta; 
}

transformed parameters {
  
 vector[Npredictors] beta;
  
  beta[1] = Beta0; //Intercept
  beta[2] = Beta1;
  beta[3] = Beta2;
  beta[4] = Beta3;
  beta[5] = Beta4;

  beta[6] = Ntheta;
  beta[7] = Ptheta;
  beta[8] = Atheta;
  beta[9] = Dtheta;
  beta[10] = Ttheta;
  beta[11] = Ctheta;
  beta[12] = Rtheta;
  
  vector[1] beta_lag; //vector [4] before
  
  // beta_lag[1] = BetaGreen;
  // beta_lag[2] = BetaMicro;
  beta_lag[1] = BetaAna;
  // beta_lag[4] = BetaNFix;
}

model {
  
  // Priors
  
  sigma_p ~ normal(0,0.3);
  
  Beta0 ~ normal(0,0.3);
  Beta1 ~ normal(0,0.3);
  Beta2 ~ normal(0,0.3);
  Beta3 ~ normal(0,0.3);
  Beta4 ~ normal(0,0.3);

  // BetaGreen ~ normal(0,0.3);
  // BetaMicro ~ normal(0,0.3);
  BetaAna ~ normal(0,0.3);
  // BetaNFix ~ normal(0,0.3);  
  
  Ntheta ~ normal(0,0.3);
  Ptheta ~ normal(0,0.3);
  Atheta ~ normal(0,0.3);
  Dtheta ~ normal(0,0.3);
  Ttheta ~ normal(0,0.3);
  Ctheta ~ normal(0,0.3);
  Rtheta ~ normal(0,0.3);
  
  //Process model
  
  // Initial states for timesteps the lag skips
  tox[1] ~ normal(0, sigma_p);
  tox[2] ~ normal(0, sigma_p);
  
  for(t in 3:uniqueID){
    if(firstdays[t] == 1){
      tox[t] ~ normal(0, sigma_p); 
      continue;
    }
    
    tox[t] ~ normal(X2[t-1,]*beta + X2[t-2, 4]*beta_lag, //used to be 2:5 for X2[t-2, 2:5] 
    sigma_p);
  }
  
  //Observation model
  
  for(t in 1:uniqueID){
    if(is_obs[t] == 1){ //If the timestep is one where we have observed field data
      Toxins[t] ~ poisson_log(tox[t]);
    }
  }
}

generated quantities {
  
  vector[uniqueID] tox_raw; //Backtransform logged values
  vector[uniqueID] log_lik; //Calculate log-likelihood
  
  for (t in 1:uniqueID) {
    
    tox_raw[t] = exp(tox[t]);
    
    if (is_obs[t] == 1) {
      log_lik[t] = poisson_log_lpmf(Toxins[t] | tox[t]); //We can only calculate likelihood when there's observed data to compare against
    } else {
      log_lik[t] = 0;
    }
  }
}
