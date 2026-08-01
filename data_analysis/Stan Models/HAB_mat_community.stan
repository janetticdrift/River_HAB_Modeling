data {
  int uniqueID; //Total number of weeks down the years
  int Nspecies; //Total number of species
  int firstdays[uniqueID]; //Days to skip modeling, first day of the year
  
  matrix [uniqueID, Nspecies] N; //Proportion in assemblage at year per species
  
  vector [uniqueID] nitrate; //Vector of nitrate levels, standardized
  vector [uniqueID] phos; //Vector of o phos levels, standardized
  vector [uniqueID] ammonium; //Vector of ammonium levels, standardized
  vector [uniqueID] DIN; //Vector of DIN levels, standardized
  vector [uniqueID] discharge; //Vector of discharge levels, logged
  vector [uniqueID] temp; //Vector of temperatures, Celsius
  vector [uniqueID] cond; //Vector of conductivity, standardized
  vector [uniqueID] rad; //Vector of shortwave radiation, standardized
}



parameters {
  
  vector<lower=0>[Nspecies] sigma_p; //var w/ process model
  vector<lower=0>[Nspecies] sigma_o; //var w/ observation model

  vector[Nspecies] Alpha; //reparameterize constrained Alpha, unconstrain by transforming below
  vector[Nspecies] Beta_diag; //create diagonal vector for intra-interactions
  matrix[Nspecies, Nspecies] Beta_off; //create off diagonal matrix
  
  //Non-centered latent state, to be used in transformed parameters
  matrix[Nspecies, uniqueID] n_nc; //fill in with modeled data
  
  vector[Nspecies] Ntheta; //parameter for nitrate each week
  vector[Nspecies] Ptheta; //parameter for o phos each week
  vector[Nspecies] Atheta; //parameter for ammonium each week
  vector[Nspecies] DINtheta; //parameter for DIN each week 
  vector[Nspecies] Dtheta; //parameter for discharge each week
  vector[Nspecies] Ttheta; //parameter for temps each week
  vector[Nspecies] Ctheta; //parameter for conductivity each week
  vector[Nspecies] Rtheta; //parameter for shortwave radiation each week
  

}

transformed parameters{
  
   matrix[Nspecies, Nspecies] Beta_d = diag_matrix(Beta_diag);
   matrix[Nspecies, Nspecies] Beta;
   
   Beta = Beta_off;
   for (i in 1:Nspecies) {
     Beta[i,i] = Beta_diag[i];
   }
   
   //Reconstruct latent state matrix `n` from the non-centered innovations `n_nc`
   matrix[Nspecies, uniqueID] n;

  //t = 1: no previous state available. As a transformed param, every aspect of it must be 
  //explicity computed. n[,1] no longer has an implicit prior like it used to
  n[,1] = n_nc[,1];
    //Note the multiplication by sigma_p. Otherwise, n's distribution is N(0,1), instead
    //of N(0,sigma_p)

  //t >= 2
  for (t in 2:uniqueID) {
    if (firstdays[t] == 1){
      n[,t] = n_nc[,t];
     continue; //continue ends current operation and returns to top of loop
    }
    n[,t] = Alpha + Beta * n[,t-1] + 
            Ntheta*nitrate[t-1] + Ptheta*phos[t-1] + Atheta*ammonium[t-1] +
            DINtheta*DIN[t-1] +
            Dtheta*discharge[t-1] + Ttheta*temp[t-1] +
            Ctheta*cond[t-1] + Rtheta*rad[t-1] + 
            sigma_p .* n_nc[,t];
  }
  //multiplying a standard normal variable (n_nc) by sigma_p gives it variance sigma-squared
  //and adding it to the mean Alpha shifts the center of the distribution.
  //The first timesteps do not include Beta, because you do not want to inform the new year
  //with any of the previous year's information, and Beta was informed by the previous timestep
  
  //Instead of sampling the latent n state directly with its mean and variance, we sample 
  //a standardized version of it (n_nc) and then construct the latent state by scaling by
  //sigma and shifting the center by the added mean. This decouples sigma from n and reduces
  //divergences.
    
}

model {
	
  //priors
  
  sigma_p ~ inv_gamma(3,1); //process model var
  sigma_o ~ inv_gamma(3,1); //normal(2.5,1); //T[0,]; #observation model var, removed truncation bc log-scale

  Alpha ~ normal(0,1);
  
  Beta_diag ~ normal(0.5, 0.1);// T[0,]; //T means Truncate, so bounded at zero
  to_vector(Beta_off) ~ normal(0, 0.1); //input matrix reshaped to vector
  
  Ntheta ~ normal(0,1);
  Ptheta ~ normal(0,1);
  Atheta ~ normal(0,1);
  DINtheta ~ normal(0,1);
  Dtheta ~ normal(0,1);
  Ttheta ~ normal(0,1);
  Ctheta ~ normal(0,1);
  Rtheta ~ normal(0,1);
  
  //Population models
  
  // ----------------- Process model (NON-CENTERED) -----------------
  // For diagonal process covariance diag(sigma_p^2), the centered process:
  // n[,t] ~ multi_normal(Alpha + Beta * n[,t-1], diag_matrix(square(sigma_p)))
  // is equivalent to the NC form implemented here:
  to_vector(n_nc) ~ normal(0, 1);
    //n_nc is drawn here, then used to estimate n in the transformed parameters block
      
      
  // ----------------- Observation model -----------------
  for (t in 1:uniqueID) {
    for (s in 1:Nspecies) {
      if (N[t, s] > -99) { // sentinel check: if the week is a week we have sampled data 
        N[t, s] ~ normal(n[s, t], sigma_o[s]);   // observation model uses reconstructed n
      }
    }
  }

}
