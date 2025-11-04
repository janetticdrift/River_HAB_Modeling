data {
  int uniqueID; //Total number of weeks down the years
  int Nspecies; //Total number of species
  int firstdays[uniqueID]; //Days to skip modeling, first day of the year
  
  matrix [uniqueID, Nspecies] N; //Proportion in assemblage at year per species
}



parameters {
  
  vector<lower=0>[Nspecies] sigma_p; //var w/ process model
  vector<lower=0>[Nspecies] sigma_o; //var w/ observation model
  
  //vector[Nspecies] log_sigma_p;
  //vector[Nspecies] log_sigma_o;

  vector[Nspecies] Alpha; //reparameterize constrained Alpha, unconstrain by transforming below
  vector[Nspecies] Beta_diag; //create diagonal vector for intra-interactions
  matrix[Nspecies, Nspecies] Beta_off; //create off diagonal matrix
  
  //Non-centered latent state, to be used in transformed parameters
  matrix[Nspecies, uniqueID] n_nc; //fill in with modeled data
  

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
  n[,1] = Alpha + sigma_p .* n_nc[,1];

  //t >= 2
  for (t in 2:uniqueID) {
    if (firstdays[t] == 1){
      n[,t] = Alpha + sigma_p .* n_nc[,t];
     continue; //continue ends current operation and returns to top of loop
    }
    n[,t] = Alpha + Beta * n[,t-1] + sigma_p .* n_nc[,t];
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
  
  //log_sigma_p ~ normal(log(0.5), 0.1);  
  //log_sigma_o ~ normal(log(0.5), 0.1);
  
  sigma_p ~ inv_gamma(3,1); //process model var
  sigma_o ~ inv_gamma(3,1); //normal(2.5,1); //T[0,]; #observation model var, removed truncation bc log-scale

  Alpha ~ normal(0,1);
  
  Beta_diag ~ normal(0.5, 0.2);// T[0,]; //T means Truncate, so bounded at zero
  to_vector(Beta_off) ~ normal(0, 0.1); //input matrix reshaped to vector
  
  //Population models
  
  // ----------------- Process model (NON-CENTERED) -----------------
  // For diagonal process covariance diag(sigma_p^2), the centered process:
  // n[,t] ~ multi_normal(Alpha + Beta * n[,t-1], diag_matrix(square(sigma_p)))
  // is equivalent to the NC form implemented here:
  to_vector(n_nc) ~ normal(0, 1);
      
      
  // ----------------- Observation model -----------------
  for (t in 1:uniqueID) {
    for (s in 1:Nspecies) {
      if (N[t, s] > -98) { // sentinel check: if the week is a week we have sampled data 
        N[t, s] ~ normal(n[s, t], sigma_o[s]);   // observation model uses reconstructed n
      }
    }
  }

}
