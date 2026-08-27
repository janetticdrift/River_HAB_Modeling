##############
#FUNCTIONS
##############

#### Function for reading in and formatting NLDAS data
### Author: Jacob Schaperow 
## Data accessed from "https://disc.gsfc.nasa.gov/information/tools?title=Hydrology%20Data%20Rods"

get_NLDASv20_datarod <- function(start_date, end_date, lat, lon, var)
{
  base_url <- "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi"
  
  full_url <- paste0(base_url,
                     "?variable=NLDAS2:NLDAS_FORA0125_H_v2.0:", var,
                     "&startDate=", start_date, "T00&",
                     "endDate=", end_date, "T23&",
                     "location=GEOM:POINT(", lon, ",%20", lat, ")",
                     "&type=asc2"
  )
  print(full_url)
  
  x <- GET(full_url)
  
  # Parse response
  z <- content(x, as = "text")
  z1 <- strsplit(z, "\n")[[1]] # separate data by line
  
  # Get dates and values. First 13 lines are metadata.
  z1 <- z1[14:length(z1)]
  
  # Split each entry by the tab separator into datetime and value
  split_entries <- strsplit(z1, "\t")
  
  # Extract datetime and value columns
  thedatetime <- sapply(split_entries, function(x) x[1]) # First part is datetime
  thevalue <- sapply(split_entries, function(x) as.numeric(x[2])) # Second part is value (convert to numeric)
  
  # Create the data frame
  df <- data.frame(datetime = as.POSIXct(thedatetime, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"), value = thevalue)
  
  return(df)
}


##############################
#' Compute Box-Cox Transformation for a given lambda
##############################
#' @param x Vector of strictly positive data.
#' @param lambda Chosen lambda based on profile likelihood plot
#'
#' @return Box-Cox transformed x.
#' 
boxcox_transform <- function(x, lambda){
  
  if(sum(x <= 0) > 0){
    stop("Box-Cox Transform only appropriate for strictly positive data.
         Some observations <= 0.")
  }
  if(lambda == 0){
    return(log(x))
  } else{
    return((x^lambda - 1) / lambda)
  }
  
}

##############################
#####Calculate ammonium from ammonia
##############################
calculate_NH4 <- function(df) {
  df <- df %>%
    dplyr::mutate(
      pKa = 0.09018 + 2727.92 / (temp_C + 273.15),
      f = 1 / (10^(pKa - pH) + 1),
      ammonium_mg_N_L = round((1 - f) * ammonia, 5)
    )
}

##############################
####Calculate standard errors
##############################
calcSE <- function(x){
  x <- x[is.na(x)==F] 
  sd(x)/sqrt(length(x))
}

##############################
#####Calculate DIC
##############################
calcDIC <- function(fit, observed, positions){
  posteriors <- fit
  obs_data <- as.matrix(observed[, positions])
  is_obs <- observed$is_obs
          #####
          #Step 1
          #####
          #Calculate average of individual log-likelihoods, second part of pDIC
  log_lik <- posteriors$log_lik #Isolate log-likelihoods
  logLik_draws <- apply(log_lik,c(1,2),sum) #sum across species
  logLik_draws <- apply(logLik_draws,1,sum) #sum acros timesteps
  logLik_draws <- mean(logLik_draws) ## take the mean Log-likelihood
  
  mu <- posteriors$mu[,,-1]     #Isolate the mu term                            
  meanmu <- apply(mu,c(2,3),mean)    #Take average of the predictions
  SD <- apply(sqrt(posteriors$sigma_p + 
                     (posteriors$sigma_o^2)), 2, mean) #combining variance of process and observation models
          #####
          #Step 2
          #####
          #Average the posteriors of all parameters, recalculate the model, then calculate log-likelihood
  LogLikelihoodofMean <- matrix(NA,1,44)    #Create empty matrix
  
  #Fill in matrix with LLs using averaged posteriors
  for(t in 1:44){
    
    LogLikelihoodofMean[,t]<- sum(dnorm(N[t+1,], meanmu[,t], SD, log=T))
    
  }
  
  MLL <- sum(LogLikelihoodofMean[,which(is_obs[-1]==1)]) #Only add together the LLs that coincide with observed dates
  
  ##############################
  #Calculate DIC
  ##############################
  
  pDIC <- 2*MLL - logLik_draws
  dic <- -2*MLL + 2*pDIC
  
  print(dic)
}

##############################
#####Calculate ELPD
##############################
calcELPD <- function(fit) {
  log_lik <- extract_log_lik(fit, parameter_name = "log_lik")
  looall <- loo(log_lik)
  ICsummary <- loo_moment_match(x = fit, loo = looall)
  
  print(ICsummary)
}

##############################
#####Toxin Simulations for River-Wide TM Mats
##############################
#This function simulates then summarizes one year of data at a time
  #The model output only includes Anabaena(t-2) in the hurdle portion, not the entire model

simulate_toxin_year <- function(
    x,
    percentcover_predicted,
    year,
    start_timestep,
    time,
    week_offset) {
  
  #Pull first two iterations of latent state toxins to initalize the simulation
  tox_conc <- x[["tox_raw"]][,start_timestep:(start_timestep + 1)]
  
  #Extract model coefficients 
  Beta0 <- x[["Beta0"]]
  Beta1 <- x[["Beta1"]]   #Green Algae
  Beta2 <- x[["Beta2"]]   #Microcoleus
  Beta3 <- x[["Beta3"]]   #Anabaena
  Beta4 <- x[["Beta4"]]   #Other N Fixers
  
  sigma_p <- x[["sigma_p"]]   #Standard deviations
  
  Phi0 <- x[["Phi0"]]       #Hurdle initiation intercept
  PhiAna <- x[["PhiAna"]]   #Lagged Anabaena effect coefficient
  
  Ntheta <- x[["Ntheta"]]   #Nitrate
  Ptheta <- x[["Ptheta"]]   #Phosphate
  Atheta <- x[["Atheta"]]   #Ammonium
  Dtheta <- x[["Dtheta"]]   #Discharge
  Ttheta <- x[["Ttheta"]]   #Temperature
  Ctheta <- x[["Ctheta"]]   #Conductivity
  Rtheta <- x[["Rtheta"]]   #Radiation
  
  #Pull number of iterations in the model
  runs <- length(Beta1)
  
  #Index the latent states of the algal assemblage models
  percentcover_clean <- percentcover_latent[env_keep, -1]
  
  #Select the appropriate time period for this year
  year_starttime <- start_timestep:(start_timestep + time - 1)
  #Index the latent algal abundances by this period
  percentcover_use <- percentcover_clean[year_starttime,,drop = FALSE]
  #Index the environmental variables by this period
  nitrate_use <- nitrate_clean[year_starttime]
  phos_use <- phos_clean[year_starttime]
  amon_use <- amon_clean[year_starttime]
  discharge_use <- discharge_clean[year_starttime]
  temp_use <- temp_clean[year_starttime]
  cond_use <- cond_clean[year_starttime]
  rad_use <- rad_clean[year_starttime]
  
  #Build the design matrix for the current year
  X1 <- as.matrix(cbind(intercept = rep(1, time),
                        percentcover_use,  # Abundances are log-transformed
                        nitrate = nitrate_use,
                        phos = phos_use,
                        amon = amon_use,
                        discharge = discharge_use,
                        temp = temp_use,
                        cond = cond_use,
                        rad = rad_use))
  
  #Build the effect coefficient matrix
  
  beta_matrix <- cbind(Beta0, Beta1, Beta2, Beta3, Beta4, 
                       Ntheta, Ptheta, Atheta,
                       Dtheta, Ttheta, Ctheta, Rtheta)
  
  #Create an empty matrix to store toxin values in
  tox <- matrix(NA, nrow = runs, ncol = time)
  
  #####Run the Simulations#####
  
  for (z in 1:runs) {
    
    #Set initial toxin concentrations to fill the first two timesteps. As a 
    #non-autoregressive model, they do not actually inform the simulation.
    tox[z, 1] <- log(tox_conc[z, 1] + 1e-6)
    tox[z, 2] <- log(tox_conc[z, 2] + 1e-6)
    
    #Extract this iteration's coefficents
    beta <- beta_matrix[z, ]
    
    #Run simulation
    for (t in 3:time) {
      
      #Hurdle: Do toxins initiate?
      phi_t <- plogis(Phi0[z] + PhiAna[z] * X1[t - 2, 4])
      #plogis is an inverse logit function, returning a result as a 0 or 1 outcome
      
      #Random draw of if toxins initiate or not
      toxin_initiate <- rbinom(1, size = 1, prob = phi_t)
      
      #Process: How many toxins are created?
      if (toxin_initiate == 1) {
        
        #If initiation is yes, run the toxin simulation
        tox[z, t] <- rnorm(1, mean = beta %*% X1[t - 1, ], sd = sigma_p[z])
        
      } else {
        #If initiaion is no, set value near to zero
        tox[z, t] <- log(0.0001)
      }
    }
  }
  
  #Summarize the simulation outcomes, calculating median and upper/lower confidence intervals
  
  simsmedian <- as.data.frame(apply(exp(tox)/1000, 2, median)) %>% 
    dplyr::rename(toxins = 1) %>% 
    dplyr::mutate(time = 1:time)
  
  simslquant <- as.data.frame(apply(exp(tox)/1000, 2, quantile, probs = 0.025)) %>% 
    dplyr::rename(CIlower = 1) %>% 
    dplyr::mutate(time = 1:time)
  
  simsuquant <- as.data.frame(apply(exp(tox)/1000, 2, quantile, probs = 0.975)) %>% 
    dplyr::rename(CIupper = 1) %>% 
    dplyr::mutate(time = 1:time)
  
  #Merge summary calculations
  riversims <- dplyr::left_join(simsmedian, simslquant, by = "time") %>% 
    dplyr::left_join(simsuquant,by = "time") %>% 
    dplyr::mutate(real_week = time + week_offset, 
                  year = year) %>% #Only 2022's week indexing needs adding by 25 instead of 26
    dplyr::mutate(model_date = ceiling_date(ymd(paste(year,"01","01",sep = "-")) +
                                              (real_week - 1) * 7 - 1, "week", week_start = 7))
  
  #Create list of summarized predictions and unsummarized raw values to calculating
  #model fit indices with later
  return(list(predictions = riversims,
              raw = tox))
  
}