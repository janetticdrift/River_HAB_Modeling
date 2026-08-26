//This model include all abiotic and biotic factors

data {
  int uniqueID; //Total number of weeks down the years
  int Nspecies; //Total number of species
  vector[uniqueID] firstdays; //Days to skip modeling, first day of the year
  vector[uniqueID] is_obs;   //Is this a time step where toxin data was collected?
  int n_obs;  //How many observed time steps were there?
  
  matrix [uniqueID, Nspecies] N; //Percent cover at year per species
  
  vector [Nspecies] id; //Vector of 1s for ID matrix
  
  vector [uniqueID] nitrate; //Vector of nitrate levels, standardized
  vector [uniqueID] phos; //Vector of o phos levels, standardized
  vector [uniqueID] ammonium; //Vector of ammonium levels, standardized
  // vector [uniqueID] DIN; //Vector of DIN levels, standardized
  vector [uniqueID] discharge; //Vector of discharge levels, logged
  vector [uniqueID] temp; //Vector of temperatures, Celsius
  vector [uniqueID] cond; //Vector of conductivity, standardized
  vector [uniqueID] rad; //Vector of shortwave radiation, standardized
}



parameters {
  
  vector<lower= 0>[Nspecies] sigma_p; //var w/ process model
  vector<lower= 0>[Nspecies] sigma_o; //var w/ observation model

  vector<lower=0>[Nspecies] Alpha;
  
  //vector<lower=0,upper=1>[Nspecies] Beta_diag; //create diagonal vector
  //matrix[Nspecies, Nspecies] Beta_off; //create off diagonal matrix
  
  matrix<upper=99>[Nspecies, uniqueID] n; //percent cover each week at each reach
  
  vector[Nspecies] Ntheta; //parameter for nitrate each species
  vector[Nspecies] Ptheta; //parameter for o phos each species
  vector[Nspecies] Atheta; //parameter for ammonium each species
  // vector[Nspecies] DINtheta; //parameter for ammonium each species
  vector[Nspecies] Dtheta; //parameter for discharge each species
  vector[Nspecies] Ttheta; //parameter for temps each species
  vector[Nspecies] Ctheta; //parameter for conductivity each species
  vector[Nspecies] Rtheta; //parameter for shortwave radiation each species
}
transformed parameters{
  matrix[Nspecies, Nspecies] ID = diag_matrix(sigma_p);
  
   //matrix[Nspecies, Nspecies] Beta_d = diag_matrix(Beta_diag);
   
   //matrix[Nspecies, Nspecies] Beta;
   
   //for(i in 1:Nspecies){
     //for(j in 1:Nspecies){
      // Beta[i,j] = (Beta_d[i,j]==0) ? Beta_off[i,j] : Beta_d[i,j]; //if it's off diagonal, supply zero, otherwise keep diag
     //}
   //}
   
}

model {
	
  //priors
  
  sigma_p ~ inv_gamma(3,1); // process model var
  sigma_o ~ inv_gamma(3,1); // observation model var
  
  //gamma ~ normal(0,tauP); //random effect for site (later pop) //gamma[s]*tauP
  //omega ~ normal(0,tauT); //random effect for time //omega[t]*tauT if convergence issues


  Alpha ~ normal(0,1);
  
  //Beta_diag ~ normal(.5, .1) T[0,]; //T means Truncate, so bounded at zero
  //to_vector(Beta_off) ~ normal(0, .1);
  
  Ntheta ~ normal(0,1);
  Ptheta ~ normal(0,1);
  Atheta ~ normal(0,1);
  // DINtheta ~ normal(0,1);
  Dtheta ~ normal(0,1);
  Ttheta ~ normal(0,1);
  Ctheta ~ normal(0,1);
  Rtheta ~ normal(0,1);

  
  //Population models
  for(t in 2:uniqueID){
    //for(s in 1:Nspecies){
      
      if(firstdays[t]==1) continue; //continue ends current operation and returns to top of loop
       n[,t] ~ multi_normal(Alpha + Ntheta*nitrate[t-1] +
                            Ptheta*phos[t-1] + 
                            Atheta*ammonium[t-1] +
                            Dtheta*discharge[t-1] + Ttheta*temp[t-1] +
                            Ctheta*cond[t-1] + Rtheta*rad[t-1], ID);
                           
}
    for(t in 1:uniqueID){
      for(s in 1:Nspecies){
        
      if(N[t,s] >= -3){ //if the year is a year we actually have sampled data for
        N[t,s] ~ normal(n[s,t], sigma_o[s]); //for collected data, we apply poisson dist to use for estimating unknown weeks
          //N[t,r] ~ normal(exp(n[t,r]), sigma_o); 
      }
    }
  }
}

generated quantities {
  ///////////////Log-Likelihood///////////////

  matrix [n_obs-1, Nspecies] log_lik;   //Matrix storing log-likelihoods for n_obs, number of observed days
  matrix[Nspecies, uniqueID] mu;         //Matrix storing mu that is used to calculate LL
  int counter = 1;                   //Set starting counter to 1

  for (t in 2:uniqueID) {       //Loop through the entire timeseries...
    mu[,t] = Alpha +
    Ntheta*nitrate[t-1] + Ptheta*phos[t-1] + Atheta*ammonium[t-1] + 
    Dtheta*discharge[t-1] + Ttheta*temp[t-1] + 
    Ctheta*cond[t-1] + Rtheta*rad[t-1];
    
    if (is_obs[t] == 1) {       //... but only run the calculation when there's an observation
       for (s in 1:Nspecies){

        log_lik[counter,s] = normal_lpdf(N[t,s] | mu[s,t], sqrt((sigma_o[s]^2)+sigma_p[s]));
      
      }
      counter += 1;  //Move to the next observation day in the log-lik vector
    }
  }

}


