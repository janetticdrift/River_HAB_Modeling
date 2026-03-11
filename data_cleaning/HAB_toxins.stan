//This model include all abiotic and biotic factors

data {
  int uniqueID; //Total number of weeks down the years
  int Nspecies; //Number of species in microscopy
  vector[uniqueID] firstdays; //Days to skip modeling, first day of the year
  vector[uniqueID] Toxins; //Vector of known toxin concentrations
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
  
  vector[uniqueID] tox_nc;  //estimated anatoxin state, non-centered
  
  real Beta0;            // intercept
  vector[Nspecies] Beta1;// species effects
  real Beta_tox;         // Toxin effect
  
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

  tox[1] = tox_nc[1];

  for(t in 2:uniqueID){
    if(firstdays[t]==1){
      tox[t] = tox_nc[t]; 
      continue;
    }
    tox[t] = Beta0 + Beta_tox*tox[t-1] + dot_product(Beta1, N[t-1]) + 
                            Ntheta*nitrate[t-1] + Ptheta*phos[t-1] + 
                            Atheta*ammonium[t-1] + Dtheta*discharge[t-1] + 
                            Ttheta*temp[t-1] + Ctheta*cond[t-1] + 
                            Rtheta*rad[t-1] +
                            sigma_p * tox_nc[t];

  }

}

model {
	
  //priors
  sigma_p ~ inv_gamma(3,1); //process model var
  sigma_o ~ inv_gamma(3,1); //observation model var
  
  Beta0 ~ normal(0,1);    //Intercept
  Beta1 ~ normal(0,1);    //population coefficient
  Beta_tox ~ normal(0,1); //Toxin coefficient
  
  tox_nc ~ normal(0,1);
  
  Ntheta ~ normal(0,1);
  Ptheta ~ normal(0,1);
  Atheta ~ normal(0,1);
  Dtheta ~ normal(0,1);
  Ttheta ~ normal(0,1);
  Ctheta ~ normal(0,1);
  Rtheta ~ normal(0,1);


// ----------------- Observation model -----------------
    for(t in 1:uniqueID){
    if(Toxins[t] >= -3){ //If t is a week we actually have collected data for
      Toxins[t] ~ normal(tox[t], sigma_o);
        
      }
    }
}
