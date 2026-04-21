//This model include all abiotic and biotic factors

data {
  int uniqueID; //Total number of weeks down the years
  int firstdays[uniqueID]; //Days to skip modeling, first day of the year
  int Npredictors; //Number of predictor variables (includes species int + env effect)
  int<lower=0, upper=1> is_obs[uniqueID];
  
  int<lower=0> Toxins[uniqueID]; //Vector of raw toxin concentrations
  matrix[uniqueID, Npredictors] X; //Design matrix of all predictors
}

parameters {
  
  real<lower= 0> sigma_p; //var w/ process model
  real<lower= 0> sigma_o; //var w/ observation model
  
  vector[uniqueID] tox_nc;  //estimated anatoxin state, non-centered and treated as on the log scale
  
  real Beta0;            // intercept
  real Beta1;            // species 1 effects
  real Beta2;            // species 2 effects
  real Beta3;            // species 3 effects
  
  real BetaAna;         // species 1 effects, lag t-2
  real BetaEpi;         // species 2 effects, lag t-2
  real BetaGeit;        // species 3 effects, lag t-2
  
  
  real Ntheta; //parameter for nitrate each species
  real Ptheta; //parameter for o phos each species
  real Atheta; //parameter for ammonium each species
  real Dtheta; //parameter for discharge each species
  real Ttheta; //parameter for temps each species
  real Ctheta; //parameter for conductivity each species
  real Rtheta; //parameter for shortwave radiation each species
}

transformed parameters {
  
  vector[uniqueID] tox;

  tox[1] = sigma_p * tox_nc[1];
  
  vector[Npredictors] beta;   // Beta vector: 1 intercept + 3 algae species + 7 env vars
  
  //Fill in empty beta vector:
  beta[1] = Beta0; //intercept
  beta[2] = Beta1; //Anabaena
  beta[3] = Beta2; //Epithemia
  beta[4] = Beta3; //Geitlerinema

  beta[5] = Ntheta;
  beta[6] = Ptheta;
  beta[7] = Atheta;
  beta[8] = Dtheta;
  beta[9] = Ttheta;
  beta[10] = Ctheta;
  beta[11] = Rtheta;
  
  vector[3] beta_lag; //Beta vector: for using a t-2 lag with the species
  
  beta_lag[1] = BetaAna; //Anabaena
  beta_lag[2] = BetaEpi; //Epithemia
  beta_lag[3] = BetaGeit; //Geitlerinema
  
  for(t in 3:uniqueID){
    if(firstdays[t]==1){
      tox[t] = tox_nc[t]; 
      continue;
    }
    tox[t] = X[t-1,]*beta + X[t-2, 2:4]*beta_lag + sigma_p*tox_nc[t];

  }

}

model {
  
  //priors
  sigma_p ~ normal(0,0.3); //process model var
  sigma_o ~ normal(0,0.3); //observation model var
  
  Beta0 ~ normal(0,0.3);    //Intercept
  BetaAna ~ normal(0,0.3);    //population coefficient
  BetaEpi ~ normal(0,0.3);    //population coefficient
  BetaGeit ~ normal(0,0.3);    //population coefficient
  
  Ntheta ~ normal(0,0.3);
  Ptheta ~ normal(0,0.3);
  Atheta ~ normal(0,0.3);
  Dtheta ~ normal(0,0.3);
  Ttheta ~ normal(0,0.3);
  Ctheta ~ normal(0,0.3);
  Rtheta ~ normal(0,0.3);
  
  // ----------------- Process model (NON-CENTERED) -----------------
    
    tox_nc ~ normal(0,1);
  //tox_nc is drawn here, then used to estimate tox in the transformed parameters block
  
  // ----------------- Observation model -----------------
    for(t in 1:uniqueID){
      if(is_obs[t] == 1){
        Toxins[t] ~ poisson_log(tox[t]);
      }
    }
}

generated quantities {
  vector[uniqueID] tox_raw;
  
  for (t in 1:uniqueID) {
    tox_raw[t] = exp(tox[t]);
  }
}
