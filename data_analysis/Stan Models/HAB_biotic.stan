//This model include only biotic factors

data {
  int uniqueID; //Total number of weeks down the years
  int Nspecies; //Total number of species
  vector[uniqueID] firstdays; //Days to skip modeling, first day of the year
  
  matrix [uniqueID, Nspecies] N; //Percent cover at year per species
  
  vector [Nspecies] id; //Vector of 1s for ID matrix
}



parameters {
  
  vector<lower= 0>[Nspecies] sigma_p; //var w/ process model
  vector<lower= 0>[Nspecies] sigma_o; //var w/ observation model

  vector<lower=0>[Nspecies] Alpha;
  
  vector<lower=0,upper=1>[Nspecies] Beta_diag; //create diagonal vector
  matrix[Nspecies, Nspecies] Beta_off; //create off diagonal matrix
  
  matrix<upper=99>[Nspecies, uniqueID] n; //percent cover each week at each reach
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
  //tauP ~ inv_gamma(1,1); site random var
  //tauT ~ inv_gamma(1,1); time random var
  
  sigma_p ~ inv_gamma(3,1); #process model var
  sigma_o ~ inv_gamma(3,1); //T[0,]; #observation model var
  
  //gamma ~ normal(0,tauP); //random effect for site (later pop) //gamma[s]*tauP
  //omega ~ normal(0,tauT); //random effect for time //omega[t]*tauT if convergence issues


  Alpha ~ normal(0,1);
  
  Beta_diag ~ normal(.5, .2) T[0,]; //T means Truncate, so bounded at zero now
  to_vector(Beta_off) ~ normal(0, .2);

  
  //Population models
  for(t in 2:uniqueID){
    //for(s in 1:Nspecies){
      
      if(firstdays[t]==1) continue;
       n[,t] ~ multi_normal(Alpha + Beta*n[,t-1], ID);
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

//Bare Biomass calculation

// generated quantities{
// vector[uniqueID] b; //100 minus everything else = bare
// 
//     for(t in 1:uniqueID){
//       b[t] = 100 - sum(exp(n[,t]));
//     }
// }
