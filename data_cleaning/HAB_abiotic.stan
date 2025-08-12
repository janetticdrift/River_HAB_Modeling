//This model includes only abiotic factors

data {
  int uniqueID; //Total number of weeks down the years
  int Nspecies; //Total number of species
  vector[uniqueID] firstdays; //Days to skip modeling, first day of the year
  
  matrix [uniqueID, Nspecies] N; //Percent cover at year per species
  
  vector [Nspecies] id; //Vector of 1s for ID matrix
  
  vector [uniqueID] nitrate; //Vector of nitrate levels, standardized
  vector [uniqueID] phos; //Vector of o phos levels, standardized
  vector [uniqueID] ammonium; //Vector of ammonium levels, standardized
  vector [uniqueID] discharge; //Vector of discharge levels, logged
  vector [uniqueID] temp; //Vector of temperatures, Celsius
  vector [uniqueID] cond; //Vector of conductivity, standardized
  vector [uniqueID] rad; //Vector of shortwave radiation, standardized
}



parameters {
  
  vector<lower= 0>[Nspecies] sigma_p; //var w/ process model
  vector<lower= 0>[Nspecies] sigma_o; //var w/ observation model

  vector<lower=0>[Nspecies] Alpha;
  
  vector<lower=0,upper=1>[Nspecies] Beta; //create intraspecific vector
  
  matrix<upper=99>[uniqueID, Nspecies] n; //percent cover each week at each reach
  
  vector[Nspecies] Ntheta; //parameter for nitrate each week
  vector[Nspecies] Ptheta; //parameter for o phos each week
  vector[Nspecies] Atheta; //parameter for ammonium each week
  vector[Nspecies] Dtheta; //parameter for discharge each week
  vector[Nspecies] Ttheta; //parameter for temps each week
  vector[Nspecies] Ctheta; //parameter for conductivity each week
  vector[Nspecies] Rtheta; //parameter for shortwave radiation each week
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
  
  Beta ~ normal(.5, .2) T[0,]; //T means Truncate, so bounded at zero now
  
  Ntheta ~ normal(0,1);
  Ptheta ~ normal(0,1);
  Atheta ~ normal(0,1);
  Dtheta ~ normal(0,1);
  Ttheta ~ normal(0,1);
  Ctheta ~ normal(0,1);
  Rtheta ~ normal(0,1);

  
  //Population models
  for(t in 2:uniqueID){
    for(s in 1:Nspecies){
      
      if(firstdays[t]==1) continue;
       n[t,s] ~ normal(Alpha[s] + Beta[s]*n[t-1, s] + Ntheta[s]*nitrate[t-1] +
                            Ptheta[s]*phos[t-1] + Atheta[s]*ammonium[t-1] +
                            Dtheta[s]*discharge[t-1] + Ttheta[s]*temp[t-1] +
                            Ctheta[s]*cond[t-1] + Rtheta[s]*rad[t-1], sigma_p[s]);
  }
}
    for(t in 1:uniqueID){
      for(s in 1:Nspecies){
        
      if(N[t,s] >= -3){ //if the year is a year we actually have sampled data for
        N[t,s] ~ normal(n[t,s], sigma_o[s]); //for collected data, we apply poisson dist to use for estimating unknown weeks
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

