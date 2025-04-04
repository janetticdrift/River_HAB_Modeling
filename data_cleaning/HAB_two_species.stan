data {
  int Nweeks; //Total number of weeks (some sites don't need as many weeks estimated)
  int Nspecies; //Total number of species
  
  matrix [Nweeks, Nspecies] N; //Percent cover observed
}



parameters {
  //real gamma[Nreach]; //random effect reach, ignore for now
  //real<lower= 0> tauP; //var w/ gamma
  
  //real omega[Nweeks]; //random effect for time, probably completely ignore
  //real<lower= 0> tauT; //var w/ time
  
  vector<lower= 0>[Nspecies] sigma_p; //var w/ process model
  vector<lower= 0>[Nspecies] sigma_o; //var w/ observation model

  vector<lower=0>[Nspecies] Beta0;
  vector<lower=0,upper=1>[Nspecies] Beta1;
  
  matrix[Nweeks, Nspecies] n; //percent cover each week
}

model {
	
  //priors
  //tauP ~ inv_gamma(1,1); site random var
  //tauT ~ inv_gamma(1,1); time random var
  
  sigma_p ~ inv_gamma(3,1); #process model var, this is extremely small and narrow
  sigma_o ~ inv_gamma(3,1); //T[0,]; #observation model var
  
  //gamma ~ normal(0,tauP); //random effect for site (later pop) //gamma[s]*tauP
  //omega ~ normal(0,tauT); //random effect for time //omega[t]*tauT if convergence issues


  Beta0 ~ normal(0,1);
  
  Beta1 ~ normal(.5, 1); //T[0,]; //T means Truncate, so bounded at zero now
  
  //Population models
  for(t in 2:Nweeks){
    for(s in 1:Nspecies){
    	
      n[t,s] ~ normal(Beta0[s] + Beta1[s]*n[t-1,s], sigma_p[s]); 
      
  }
}
    for(t in 1:Nweeks){
      for(s in 1:Nspecies){
        
      if(N[t,s] >= -3){ //if the year is a year we actually have sampled data for
        N[t,s] ~ normal(n[t,s], sigma_o[s]); //for collected data, we apply poisson dist to use for estimating unknown weeks
          //N[t,r] ~ normal(exp(n[t,r]), sigma_o); 
      }
      
    }

  }

 }

//Bare Biomass calculation

generated quantities{
vector[Nweeks] b; //100 minus everything else = bare
  
    for(t in 1:Nweeks){
      b[t] = 100 - sum(exp(n[t])); 
    }
}


