data {
  int uniqueID; //Total number of weeks down the years
  int Nspecies; //Total number of species
  int firstdays[uniqueID]; //Days to skip modeling, first day of the year
  
  matrix [uniqueID, Nspecies] N; //Proportion in assemblage at year per species
}



parameters {
  //real<lower=0> gamma[Nreach]; //random effect reach
  //real<lower=0> tauP; //var w/ gamma
  
  vector<lower=0>[Nspecies] sigma_p; //var w/ process model
  vector<lower=0>[Nspecies] sigma_o; //var w/ observation model

  vector<lower=0>[Nspecies] Alpha; 
  vector<lower=0,upper=1>[Nspecies] Beta_diag; //create diagonal vector for intra-interactions
  matrix[Nspecies, Nspecies] Beta_off; //create off diagonal matrix
  
  matrix<upper=99>[Nspecies, uniqueID] n; //fill in with modeled data
  //vector<upper=99>[Nreach] n[Nspecies, uniqueID];

}

transformed parameters{
  matrix[Nspecies, Nspecies] ID = diag_matrix(sigma_p);
  
   matrix[Nspecies, Nspecies] Beta_d = diag_matrix(Beta_diag);
   
   matrix[Nspecies, Nspecies] Beta;
   
   Beta = Beta_off;
   for (i in 1:Nspecies) {
     Beta[i,i] = Beta_diag[i];
   }
   
   //for(i in 1:Nspecies){
   //  for(j in 1:Nspecies){
   //    Beta[i,j] = (Beta_d[i,j]==0) ? Beta_off[i,j] : Beta_d[i,j];
  //   }
 //  }
}

model {
	
  //priors
  
  sigma_p ~ gamma(2, 0.1);//inv_gamma(4,1); //process model var
  sigma_o ~ gamma(2, 0.1);//inv_gamma(4,1); //normal(2.5,1); //T[0,]; #observation model var, removed truncation bc log-scale
  
  //tauP ~ inv_gamma(1,1); //reach random var
  //gamma ~ normal(0,tauP); //random effect for reach //gamma[r]*tauP
  //omega ~ normal(0,tauT); //random effect for time //omega[t]*tauT if convergence issues


  Alpha ~ normal(1,1)T[0,];
  
  Beta_diag ~ normal(.5, .2) T[0,]; //T means Truncate, so bounded at zero now
  to_vector(Beta_off) ~ normal(0, .2); //input matrix reshaped to vector
  
  //Population models
  
      //Process model
  //for(r in 1:Nreach){
    for(t in 2:(uniqueID)){
   // if(firstsample[r] <= -3) continue; way to skip to the first sampled date for reaches with missing data?
    	if(firstdays[t]==1) continue; //continue ends current operation and returns to top of loop
       n[,t] ~ multi_normal(Alpha + Beta*n[,t-1], ID);
      
       
   }
//}
    //Observation model
    //for(r in 1:Nreach){
      for(t in 1:uniqueID){
        for(s in 1:Nspecies){
      
        if(N[t,s] >= -3){ //if the year is a year we actually have sampled data for
        N[t,s] ~ normal(n[s,t], sigma_o[s]); //for collected data, we apply poisson dist to use for estimating unknown weeks
      }
    }  
  }
//}
}
