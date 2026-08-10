//This model include all abiotic and biotic factors

data {
  int uniqueID;                             //Total number of weeks down the years
  int firstdays[uniqueID];                  //Days to skip modeling, first day of the year
  int Npredictors;                          //Number of predictor variables (includes species int + env effect)
  int<lower=0, upper=1> is_obs[uniqueID];   //Is this a time step where toxin data was collected?
  
  int<lower=0> Toxins[uniqueID];            //Vector of raw toxin concentrations
  matrix[uniqueID, Npredictors] X1;         //Design matrix of all predictors
}

parameters {
  
  real<lower= 0> sigma_p;       //Standard deviation w/ process model
  
  vector[uniqueID] tox;         //latent anatoxin state on log scale
  
  //Process model coefficients (How much toxin is produced?)
  real Beta0;            
  real Beta1;     //Anabaena      
  real Beta2;     //Epithemia       
  real Beta3;     //Geitlerinema       
  
  // real Ntheta; 
  // real Ptheta; 
  // real Atheta; 
  // real Dtheta; 
  // real Ttheta; 
  // real Ctheta; 
  // real Rtheta; 
  
  //Hurdle model coefficients (When do toxins start getting produced?)
  real Phi0;    //Probability intercept of timing of toxins initiating
  real PhiAna;  //Effect of t-2 Anabaena lag on timing of toxins initiating
}

transformed parameters {
  
 vector[Npredictors] beta;
  
  beta[1] = Beta0;
  beta[2] = Beta1;
  beta[3] = Beta2;
  beta[4] = Beta3;

  // beta[5] = Ntheta;
  // beta[6] = Ptheta;
  // beta[7] = Atheta;
  // beta[8] = Dtheta;
  // beta[9] = Ttheta;
  // beta[10] = Ctheta;
  // beta[11] = Rtheta;

}

model {
  
  // Priors
  
  sigma_p ~ normal(0,0.3);
  
  Beta0 ~ normal(0,1);
  Beta1 ~ normal(0,1);
  Beta2 ~ normal(0,1);
  Beta3 ~ normal(0,1);
  
  // Ntheta ~ normal(0,1);
  // Ptheta ~ normal(0,1);
  // Atheta ~ normal(0,1);
  // Dtheta ~ normal(0,1);
  // Ttheta ~ normal(0,1);
  // Ctheta ~ normal(0,1);
  // Rtheta ~ normal(0,1);
  
  //Priors for hurdle model
  Phi0 ~ normal(0,1);
  PhiAna ~ normal(0,1);
  
  //Process model
  
  // Initial states for timesteps the "t in 3" skips
  tox[1] ~ normal(0, sigma_p);
  tox[2] ~ normal(0, sigma_p);

  
  for(t in 3:uniqueID){
    if(firstdays[t] == 1){        //If this observation is the first day of a year, clear the previous year's history of effect
      tox[t] ~ normal(0, sigma_p); 
      continue;
    }
      //Latent toxin magnitudes are predicted by algal abundances t-1 and env drivers t-1
    tox[t] ~ normal(X1[t-1,]*beta, sigma_p);
  }
  
  //Observation model (Hurdle)
  
  for(t in 1:uniqueID){
    
    if(is_obs[t] == 1){    //If the timestep is one where we have observed data...
        
        real phi_t;        //phi_t is not a single estimated parameter, it is a calculated quantity that changes every time step. Should therefore not go in parameter block
        
        if(t <= 2 || firstdays[t] == 1 || firstdays[t-1] == 1){ //If the timestep is on the first day of the year, or day before current time step was the first day of the year...
          
          //Do not use the Anabaena t-2 lag
          phi_t = inv_logit(Phi0);
                  //inv_logit is the function that maps any number to a probability value between 0 and 1
          
        } else { //If the time step was not on the first or second day of the year, use the Anabaena t-2 lag
        
        phi_t =  inv_logit(Phi0 + PhiAna * X1[t-2,4]);  //This is the timestep-specific probability that toxins initiate (Intercept + lagged Anabaena)
         
        }                     
      if(Toxins[t] == 0){  //And if that observed timestep does NOT observe any toxins...
        
        //Absence: Hurdle not crossed
        target += bernoulli_lpmf(0 | phi_t); //Probability hurdle was not crossed (See 0 given the probability at timestep phi_t)
        //target += is the longhand version of the ~ sampling notation
          //Needed here because there is no other left-hand side of the equation
                              
      } else {
        
        //Presence: Hurdle crossed
        target += bernoulli_lpmf(1 | phi_t) +             //Probability that toxins occur given probability of timestep phi_t
                  poisson_lpmf(Toxins[t] | exp(tox[t])) -  //What is the observed toxin concentration then, given the estimated latent toxin concentration?
                  poisson_lccdf(0 | exp(tox[t]));          //Subtract the probability of the observed value being zero to get the prob of the value being non-zero
      }
    }
  }
}

generated quantities {
  
  vector[uniqueID] tox_raw; //Backtransform logged toxin values
  vector[uniqueID] log_lik; //Calculate log-likelihood
  
  for (t in 1:uniqueID) {
    
    tox_raw[t] = exp(tox[t]);  //Exponentiate latent toxin values
    
    if (is_obs[t] == 1) {      //If the time series is a week we have observed data for...
      
      real phi_t;              //Probability of hurdle being crossed
      
      if(t <= 2 || firstdays[t] == 1 || firstdays[t-1] == 1){  //Avoid indexing error by skipping first days
        phi_t = inv_logit(Phi0);
        
      } else {
      
      phi_t = inv_logit(Phi0 + PhiAna * X1[t-2,4]);
      
      }
      
      if(Toxins[t] == 0){      //And if there are no toxins at this week...
        log_lik[t] = bernoulli_lpmf(0 | phi_t);
        
      } else {
        
        //If observations are present but toxins are non-zero...
        log_lik[t] = bernoulli_lpmf(1 | phi_t) + 
                     poisson_log_lpmf(Toxins[t] | tox[t]) - 
                     poisson_lccdf(0 | exp(tox[t]));
      }
    } else {
      //If there were no observations at this time step at all...
      log_lik[t] = 0;
    }
  }
}
