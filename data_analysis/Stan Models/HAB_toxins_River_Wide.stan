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
  
  real<lower= 0> sigma_p;       //Standard deviation w/ process model
  real<lower=0, upper=1> phi;   //Probability of excess zeros, constrained between 0 and 1
  
  vector[uniqueID] tox;   //latent anatoxin state on log scale
  
  real Beta0;            
  real Beta1;            
  real Beta2;            
  real Beta3;            
  real Beta4;            
  
  real BetaAna;          
  
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
  // beta[13] = DINtheta;
  
  vector[1] beta_lag; 
  
  beta_lag[1] = BetaAna;
}

model {
  
  // Priors
  
  sigma_p ~ normal(0,0.3);
  phi ~ beta(1,1);        // A beta prior naturally restricts to between 0 and 1
  
  Beta0 ~ normal(0,0.3);
  Beta1 ~ normal(0,0.3);
  Beta2 ~ normal(0,0.3);
  Beta3 ~ normal(0,0.3);
  Beta4 ~ normal(0,0.3);

  BetaAna ~ normal(0,0.3);
  
  Ntheta ~ normal(0,0.3);
  Ptheta ~ normal(0,0.3);
  Atheta ~ normal(0,0.3);
  Dtheta ~ normal(0,0.3);
  Ttheta ~ normal(0,0.3);
  Ctheta ~ normal(0,0.3);
  Rtheta ~ normal(0,0.3);
  
  phi ~ beta(1,1);        // Weak prior
  
  //Process model
  
  // Initial states for timesteps the lag skips
  tox[1] ~ normal(0, sigma_p);
  tox[2] ~ normal(0, sigma_p);
  
  for(t in 3:uniqueID){
    if(firstdays[t] == 1){
      tox[t] ~ normal(0, sigma_p); 
      continue;
    }
    
    tox[t] ~ normal(X2[t-1,]*beta + X2[t-2, 4]*beta_lag, 
    sigma_p);
  }
  
  //Observation model (Zero-Inflated Poisson)
  
  for(t in 1:uniqueID){
    if(is_obs[t] == 1){    //If the timestep is one where we have observed data...
      if(Toxins[t] == 0){  //And if that timestep does NOT observe any toxins...
        target += log_sum_exp(bernoulli_lpmf(1 | phi),   //Parts of the data where there are simply just many zeros (structural)
                              bernoulli_lpmf(0 | phi) +
                              poisson_log_lpmf(Toxins[t] | tox[t]));  // Parts of the data where zeros could still occur in the middle of the process
                  //Sum both structural and process zeros on the log scale (log_sum_exp)
                              
      } else {
        
        target += bernoulli_lpmf(0 | phi) + poisson_log_lpmf(Toxins[t] | tox[t]);
        //target += is the longhand version of the ~ sampling notation
          //Needed here because there is no other left-hand side of the equation
      }
    }
  }
}

generated quantities {
  
  vector[uniqueID] tox_raw; //Backtransform logged toxin values
  vector[uniqueID] log_lik; //Calculate log-likelihood
  
  for (t in 1:uniqueID) {
    
    tox_raw[t] = exp(tox[t]);  //Exponentiate latent toxin values
    
    if (is_obs[t] == 1) {      //If the time series is a week we have observed data for...
      if(Toxins[t] == 0){      //And if there are no toxins at this week...
        log_lik[t] = log_sum_exp(bernoulli_lpmf(1 | phi),
                                 bernoulli_lpmf(0 | phi) + 
                                 poisson_log_lpmf(Toxins[t] | tox[t]));
      } else {
        
        //If observations are present but toxins are non-zero...
        log_lik[t] = bernoulli_lpmf(0 | phi) + 
                     poisson_log_lpmf(Toxins[t] | tox[t]);
      }
    } else {
      //If there were no observations at this time step...
      log_lik[t] = 0;
    }
  }
}
