//This model include all abiotic and biotic factors

data {
  int uniqueID; //Total number of weeks down the years
  vector[uniqueID] firstdays; //Days to skip modeling, first day of the year
  vector[uniqueID] Toxins; //Vector of known toxin concentrations
  matrix[uniqueID, Nspecies] N; //microscopy abundances per week
  
  // vector [uniqueID] nitrate; //Vector of nitrate levels, standardized
  // vector [uniqueID] phos; //Vector of o phos levels, standardized
  // vector [uniqueID] ammonium; //Vector of ammonium levels, standardized
  // vector [uniqueID] discharge; //Vector of discharge levels, logged
  // vector [uniqueID] temp; //Vector of temperatures, Celsius
  // vector [uniqueID] cond; //Vector of conductivity, standardized
  // vector [uniqueID] rad; //Vector of shortwave radiation, standardized
}

parameters {
  
  real<lower= 0> sigma_p; //var w/ process model
  real<lower= 0> sigma_o; //var w/ observation model
  
  vector[uniqueID] tox; //estimated anatoxins
  
  real Beta0;
  vector[Nspecies] Beta1;
  vector[uniqueID] Toxtheta;
  
  // vector[Nspecies] Ntheta; //parameter for nitrate each week
  // vector[Nspecies] Ptheta; //parameter for o phos each week
  // vector[Nspecies] Atheta; //parameter for ammonium each week
  // vector[Nspecies] Dtheta; //parameter for discharge each week
  // vector[Nspecies] Ttheta; //parameter for temps each week
  // vector[Nspecies] Ctheta; //parameter for conductivity each week
  // vector[Nspecies] Rtheta; //parameter for shortwave radiation each week
}

model {
	
  //priors
  sigma_p ~ inv_gamma(3,1); //process model var
  sigma_o ~ inv_gamma(3,1); //observation model var
  
  Beta0 ~ normal(0,1);
  Beta1 ~ normal(0,1);
  Toxtheta ~ normal(0,1);
  
  // Ntheta ~ normal(0,1);
  // Ptheta ~ normal(0,1);
  // Atheta ~ normal(0,1);
  // Dtheta ~ normal(0,1);
  // Ttheta ~ normal(0,1);
  // Ctheta ~ normal(0,1);
  // Rtheta ~ normal(0,1);

  
  //Population models
  for(t in 2:uniqueID){
      if(firstdays[t]==1) continue; //continue ends current operation and returns to top of loop
       tox[t] ~ normal(Beta0 + Beta1*n + Toxtheta*tox[t-1], sigma_p);
}

    for(t in 1:uniqueID){
      if(Toxins[t] >= -3){ //if the week is a week we actually have sampled data for
        Toxins[t] ~ normal(tox[t], sigma_o); //for collected data
      }
    }
}