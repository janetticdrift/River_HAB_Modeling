data {
  int Nyears; //Number of years we're estimating population size for
  int Nsites; //Number of sites 
  int Ncov; //Number of covars
  
  real areascale[Nsites];//area scaling parameter 
  
  real ski[Nsites]; //whether it's a ski area or not
  
  int N[Nyears, Nsites]; //Number of plants at year ndata //maybe matrix?? idk
  
  matrix[Ncov, Nsites] X [(Nyears)]; //Declares as "dims found=(14,2,10) "
  //matrix[(Nyears+1), Ncov] X [Nsites];   //Declares as "dims declared=(10,14,2)"
  //From manual "Continuing to add indices, m[1,2,3] is of type row_vector and denotes the third row of the matrix denoted by m[1,2] //This part is for incorporating env covariates
}


parameters {
  real gamma[Nsites]; //random effect for pop/site
  real<lower= 0> tauP; //var w/ gamma
  
  //real omega[Nyears]; //random effect for time
  //real<lower= 0> tauT; //var w/ time
  
  real<lower= 0> sigma; //var w/ proc
  
  real Beta0;
  real Alpha; //this was a categorical con for Sage specifically
  vector[Ncov] Beta;
  
  matrix[Nyears, Nsites] n; //number of plants each year at each sites, this is what JLD wants
}

model {
	
  //priors
  tauP ~ inv_gamma(1,1);
  //tauT ~ inv_gamma(1,1);
  
  sigma ~ inv_gamma(1,1);
  
  gamma ~ normal(0,tauP); //random effect for site (later pop) //gamma[s]*tauP
  //omega ~ normal(0,tauT); //random effect for time //omega[t]*tauT if convergence issues


  Beta0 ~ normal(0,1);
  Alpha ~ normal(0,1);
  Beta ~ normal(0,1); //setting them all to be the same I guess
  
  
  //the rest
  for(t in 2:(Nyears)){
    for(s in 1:Nsites){
    	
      //n[t,s] ~ normal(Beta0 + to_row_vector(X[t,,s])*Beta + n[t-1,s] + gamma[s] + omega[t],sigma); 
      //n[t,s] ~ normal(Beta0 + to_row_vector(X[t,,s])*Beta + n[t-1,s] ,sigma); //remove random effects to see if covariates show stronger signal
      //n[t,s] ~ normal(Beta0 + to_row_vector(X[t,,s])*Beta + n[t-1,s] + gamma[s],sigma);//incorporate back in site RE not time
      n[t,s] ~ normal(Beta0 + n[t-1,s] + gamma[s],sigma);//incorporate FE for ski area

      if(N[t,s] > 0){ //if the year is a year we actually have sampled data for.
        N[t,s] ~ poisson(exp(n[t,s]*areascale[s])); 
        //N[t,s] ~ poisson(exp(n[t,s])); 
      }
      
    }
  }
  
}

generated quantities{ //this part is model selection

  matrix[Nyears,Nsites] npreds; //new preds, based soley on process model, not restricted by data 
  
  for (t in 2:(Nyears)){
    for (s in 1:Nsites){
     npreds[t,s]= (Beta0 + Alpha*ski[s] + to_row_vector(X[t,,s])*Beta + n[t-1,s] + gamma[s]);
    }
  }
}

