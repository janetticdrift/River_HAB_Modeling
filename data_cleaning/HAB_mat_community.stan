data {
  int uniqueID; //Total number of weeks down the years
  int Nspecies; //Total number of species
  int Nreach; //Total number of reaches (3 each year, removed the added one in 2023)
  //vector[uniqueID] firstdays; //Days to skip modeling, first day of the year
  
  matrix [uniqueID, Nspecies] N [Nreach]; //Percent cover at reach and year per species
    //Declares as "dims found=(weeks,species,reach)
  
  //vector [Nspecies] id; //Vector of 1s for ID matrix
  
}



parameters {
  real<lower=0> gamma[Nreach]; //random effect reach
  real<lower=0> tauP; //var w/ gamma
  
  vector<lower=0>[Nspecies] sigma_p; //var w/ process model
  vector<lower=0>[Nspecies] sigma_o; //var w/ observation model

  vector<lower=0>[Nspecies] Alpha; //Is a beta per reach overkill? could turn this into a matrix
  vector<lower=0,upper=1>[Nspecies] Beta_diag; //create diagonal vector for intra-interac
  matrix[Nspecies, Nspecies] Beta_off; //create off diagonal matrix
  
  matrix<upper=99>[uniqueID, Nspecies] n[Nreach]; //fill in with modeled data

}

transformed parameters{
  matrix[Nspecies, Nspecies] ID = diag_matrix(sigma_p);
  
   matrix[Nspecies, Nspecies] Beta_d = diag_matrix(Beta_diag);
   
   matrix[Nspecies, Nspecies] Beta;
   
   for(i in 1:Nspecies){
     for(j in 1:Nspecies){
       Beta[i,j] = (Beta_d[i,j]==0) ? Beta_off[i,j] : Beta_d[i,j];
     }
   }
}

model {
	
  //priors
  
  sigma_p ~ inv_gamma(3,1); //process model var
  sigma_o ~ inv_gamma(3,1); //normal(2.5,1); //T[0,]; #observation model var, removed truncation bc log-scale
  
  tauP ~ inv_gamma(1,1); //reach random var
  gamma ~ normal(0,tauP); //random effect for reach //gamma[r]*tauP
  //omega ~ normal(0,tauT); //random effect for time //omega[t]*tauT if convergence issues


  Alpha ~ normal(0,1);
  
  Beta_diag ~ normal(.5, .2) T[0,]; //T means Truncate, so bounded at zero now
  to_vector(Beta_off) ~ normal(0, .2); //input matrix reshaped to vector
  
  //Population models
  
      //Process model
  for(r in 1:Nreach){
   // if(firstsample[r] <= -3) continue; way to skip to the first sampled date for reaches with missing data?
    for(t in 2:(uniqueID)){
    	
       n[r,,t] ~ multi_normal(Alpha + Beta*n[r,,t-1] + gamma[r], ID);
      
       
   }
}
    //Observation model
    for(r in 1:Nreach){
      for(t in 1:uniqueID){
        for(s in 1:Nspecies){
      
        if(N[r,t,s] >= -3){ //if the year is a year we actually have sampled data for
          N[r,t,s] ~ normal(n[r,t,s] + gamma[r], sigma_o[s]); //for collected data
                            //t and s may switch positions, pulling from all years stan
      }
    }  
  }
}
}
