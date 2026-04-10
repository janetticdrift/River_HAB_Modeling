//This model include all abiotic and biotic factors

data {
  int uniqueID; //Total number of weeks down the years
  int Nspecies; //Number of species in microscopy
  int firstdays[uniqueID]; //Days to skip modeling, first day of the year
  
  int<lower=0, upper=1> is_obs[uniqueID];
  int<lower=0> Toxins[uniqueID]; //Vector of raw toxin concentrations
  matrix[uniqueID, Nspecies] N; //microscopy abundances per week
  
  vector [uniqueID] nitrate; //Vector of nitrate levels, standardized
  vector [uniqueID] phos; //Vector of o phos levels, standardized
  vector [uniqueID] ammonium; //Vector of ammonium levels, standardized
  vector [uniqueID] discharge; //Vector of discharge levels, logged
  vector [uniqueID] temp; //Vector of temperatures, Celsius
  vector [uniqueID] cond; //Vector of conductivity, standardized
  vector [uniqueID] rad; //Vector of shortwave radiation, standardized
}

parameters {
  
  real<lower= 0> sigma_p; //var w/ process model
  real<lower= 0> sigma_o; //var w/ observation model
  
  vector[uniqueID] log_tox_nc;  //estimated anatoxin state, non-centered and treated as on the log scale
  
  real Beta0;            // intercept
  vector[Nspecies] Beta1;// species effects
  
  real Ntheta; //parameter for nitrate each species
  real Ptheta; //parameter for o phos each species
  real Atheta; //parameter for ammonium each species
  real Dtheta; //parameter for discharge each species
  real Ttheta; //parameter for temps each species
  real Ctheta; //parameter for conductivity each species
  real Rtheta; //parameter for shortwave radiation each species
}

transformed parameters {
  
  vector[uniqueID] log_tox;
  
  log_tox[1] = log_tox_nc[1];
  
  for(t in 2:uniqueID){
    if(firstdays[t]==1){
      log_tox[t] = log_tox_nc[t]; 
      continue;
    }
    log_tox[t] = Beta0 + dot_product(Beta1, N[t-1]) + 
      Ntheta*nitrate[t-1] + Ptheta*phos[t-1] + 
      Atheta*ammonium[t-1] + Dtheta*discharge[t-1] + 
      Ttheta*temp[t-1] + Ctheta*cond[t-1] + 
      Rtheta*rad[t-1] +
      sigma_p * log_tox_nc[t];
    
  }
  
}

model {
  
  //priors
  sigma_p ~ normal(0,0.3); //process model var
  sigma_o ~ normal(0,0.3); //observation model var
  
  Beta0 ~ normal(0,0.2);    //Intercept
  Beta1 ~ normal(0,0.2);    //population coefficient
  
  Ntheta ~ normal(0,1);
  Ptheta ~ normal(0,1);
  Atheta ~ normal(0,1);
  Dtheta ~ normal(0,1);
  Ttheta ~ normal(0,1);
  Ctheta ~ normal(0,1);
  Rtheta ~ normal(0,1);
  
  // ----------------- Process model (NON-CENTERED) -----------------
    
    log_tox_nc ~ normal(0,1);
  //tox_nc is drawn here, then used to estimate tox in the transformed parameters block
  
  // ----------------- Observation model -----------------
    for(t in 1:uniqueID){
      if(is_obs[t] == 1){
        Toxins[t] ~ poisson_log(log_tox[t]);
      }
    }
}

generated quantities {
  vector[uniqueID] tox_raw;
  
  for (t in 1:uniqueID) {
    tox_raw[t] = exp(log_tox[t]);
  }
}
