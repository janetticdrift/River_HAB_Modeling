#Historical Predictions

#Pull out starting abundances 
x <- rstan::extract(fit.m4)
abundances <- as.data.frame(x[c("n","b")]) %>% 
  select(1:12) %>% #initiate model with first 3 weeks
  head(100) #Don't commit to all 15k iterations yet
alphas <- x[["Alpha"]]
betas <- as.array(x[["Beta"]]) %>% 
  head(100)
