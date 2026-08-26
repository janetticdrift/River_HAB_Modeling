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
