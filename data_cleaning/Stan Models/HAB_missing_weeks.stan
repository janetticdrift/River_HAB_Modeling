data {
  int Nweeks; //Total number of weeks (some sites don't need as many weeks estimated)
  int Nreach; //Total number of reaches (3 one year, 4 other years)
  
  matrix [Nweeks, Nreach] N; //Percent cover at year per reach
}



parameters {
  //real gamma[Nreach]; //random effect reach, ignore for now
  //real<lower= 0> tauP; //var w/ gamma
  
  //real omega[Nweeks]; //random effect for time, probably completely ignore
  //real<lower= 0> tauT; //var w/ time
  
  real<lower= 0> sigma_p; //var w/ process model
  real<lower= 0> sigma_o; //var w/ observation model

  real Beta0;
  //real Beta1;
  
  matrix[Nweeks, Nreach] n; //percent cover each week at each reach
}

model {
	
  //priors
  //tauP ~ inv_gamma(1,1); site random var
  //tauT ~ inv_gamma(1,1); time random var
  
  sigma_p ~ inv_gamma(1,1); #process model var
  sigma_o ~ normal(2.5,1); //T[0,]; #observation model var
  
  //gamma ~ normal(0,tauP); //random effect for site (later pop) //gamma[s]*tauP
  //omega ~ normal(0,tauT); //random effect for time //omega[t]*tauT if convergence issues


  Beta0 ~ normal(0,1);
  
  //Beta1 ~ normal(0,1)T[0,]; //T means Truncate, so bounded at zero now
  
  //Population models
  for(t in 2:(Nweeks)){
    for(r in 1:Nreach){
    	
      n[t,r] ~ normal(Beta0 + n[t-1,r], sigma_p); 
      //n[t,r] ~ normal(Beta0 + Beta1*n[t-1,r], sigma_p); with beta1
      

      if(N[t,r] >= 0){ //if the year is a year we actually have sampled data for... make sure it's greater than or EQUAL TO
        N[t,r] ~ normal(exp(n[t,r]), sigma_o); //for collected data, we apply poisson dist to use for estimating unknown weeks
          //consider that this may be a beta.
          //N[t,r] ~ normal(exp(n[t,r]), sigma_o); 
      }
      
    }
  }
  
}