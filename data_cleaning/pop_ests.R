#~~~~~~~~~~~~~~~
#Estimating population sizes
#SLE Fall 2023
#Top Level----
#~~~~~~~~~~~~~~~

#Functions----
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

log_lik<-function(data,pred,latent,sigma,iter, scale)
{
  keep<-which(data!=-99,arr.ind = T)
  log_lik_out<-array(NA,c(iter,dim(data)))
  log_lik_out_keep<-matrix(NA,iter,dim(keep)[1])
  
  for(i in 1:iter){
    for(t in 1:dim(data)[1]){
      for(s in 1:dim(data)[2]) {
        log_lik_out[i,t,s]<-dnorm(latent[i,t,s],pred[i,t,s],sigma[i],log=T)+dpois(data[t,s],exp(latent[i,t,s]*scale[s]),log=T)
        
      }}
    log_lik_out_keep[i,]<-c(log_lik_out[cbind(i,keep)])
  }
  
  return(log_lik_out_keep)
}

group_covars <- function(covars, prism, eo.names){
  lists <- vector("list", length(eo.names)) 
  
  for(i in 1:length(eo.names)){
    site <- prism[which(prism$EO == eo.names[i]),covars]
    lists[[i]] <- site
  }
  
  if(length(covars) >1){
    X <- array(c(unlist(lists)), dim = c(nrow(lists[[1]]),ncol(lists[[1]]),length(eo.names)))
    dimnames(X) <- list(c(2009:2024),covars,eo.names)
  }else{
    X <- array(c(unlist(lists)), dim = c(length(lists[[1]]),1,length(eo.names)))
    dimnames(X) <- list(c(2009:2024), covars, eo.names)
  }
  return(X)
}


#Packages----
library(rstan)
library(shinystan)
library(loo)


#Setup----
setwd("G:\\My Drive\\UNR MS 2022\\Population Dem of Draast")
#setwd("E:\\My Drive\\UNR MS 2022\\Population Dem of Draast") #On server

#Read in covariates
prism <- read.csv(".\\Covars\\prism_cleaned_24_v2.csv")
prism.og <- prism
  ##Winter 2009 means the winter 2008-2009 season. 

prism <- prism[which(prism$wnt.yr %in% 2009:2024),]

#Read in scale parameters
scale_params <- readRDS(".\\Scripts\\saved objects\\scale_param.RDS")

#Read in data
data <- read.csv(".\\Draba Data\\Draba_2010-2024_Plant_Data_cleaner.csv")
data.og <- data
data <- data[which(data$Use.Pop == 1),] 

data$Location[which(data$Location == "Relay Peak Road"|data$Location == "Tamarack")] <- "Mt Rose Wilderness"
data$Year <- substrRight(data$Survey.Date, 4)

summary(as.factor(data$Year))

#Calculate totals per SITE
totals.eo <- aggregate(Use.Pop ~ Year + EO, data = data, FUN = sum, na.rm = T)

totals.eo <- reshape(data =totals.eo, idvar = "Year", timevar = "EO", direction = "wide", sep = "_")
names(totals.eo) <- c("Year", sub(".+_(.+)","\\1",colnames(totals.eo[,2:11])))

totals.eo <- totals.eo[order(totals.eo$Year),]

#Silly to get indexing correct 
totals.eo <- rbind(matrix(NA, nrow = 7, ncol = 11, dimnames = list(c(2014, 2016:2021),colnames(totals.eo))), totals.eo)
totals.eo$Year[1:7] <- rownames(totals.eo)[1:7]
totals.eo <- totals.eo[order(totals.eo$Year),]

totals.eo[is.na(totals.eo)] <- -99 #NA's to -99 for if statment in model
totals.eo <- rbind(c(2009, rep(-99,10)),totals.eo) #Add a  year of NA's to use first year's of data
row.names(totals.eo) <- c(totals.eo$Year)

#Choose which sites to use
sites.keep <- c("1A", "1B", "1J", "2B", "2F", "3B", "Bonanza", "Bruce's")
totals.eo <- totals.eo[,c("Year", sites.keep)]

#Fix covariates to proper format
prism[,3:ncol(prism)] <- apply(prism[,3:ncol(prism)], 2, FUN = function(x) (x-mean(x))/sd(x)) #standardize covariates
cors <- (cor(prism[,3:ncol(prism)]))

cors[cors>0.70] <- NA

#Make a list with each site's covariates (to then be turned into an array)
eo.names <- colnames(totals.eo[,-1])

#H1 <- fit #Average conditions: Tmean.S, Tmean.W, PPT.S, SFD
#H2 <- fit #Maximum temperatures in Summer and Winter
#H4 <- fit #EWE through hot temperatures and DOS 
#H5 <- fit #Drought conditions (spei270d)

X.H1 <- group_covars(covars = c("TMEAN.W", "PPT.S", "SFD"), prism = prism, eo.names = eo.names)
X.H1.alt <- group_covars(covars = c("TMEAN.W", "PPT.S", "PPT.W"), prism = prism, eo.names = eo.names)
X.H2 <- group_covars(covars = c("TMAX.S", "TMAX.W"), prism = prism, eo.names = eo.names)
X.H4 <- group_covars(covars = c("TMAX.EH.S", "TMAX.EH.W", "SFD"), prism = prism, eo.names = eo.names)
X.H5 <- group_covars(covars = c("spei270d"), prism = prism, eo.names = eo.names)

###Gather data in proper forms for stan----

#Re-sort scaling parameter so ID is correct
scale_params <- scale_params[order(scale_params$EO),]
scale_params <- scale_params[scale_params$EO %in% eo.names,]

#Ski areas
ski <- data.frame(EO = unique(data$EO), yes = c(0,0,0,0,1,1,0,1,1,0)) #Ski area 0/1 for all pops
ski <- ski[ski$EO %in% eo.names,] #Subset to sites we're using

#Silly N's
Nyears <- length(2009:2024) #Number of years we're estimating for
totalyears <- 2009:2024
datayears <- as.numeric(totals.eo$Year)
Ndata <- length(datayears)
Nsites <- (ncol(totals.eo)-1)


#Model----

#Gather Stan Data
mod.data.h1 <- list("Nyears" = Nyears, "Nsites" = Nsites, "Ncov" = ncol(X.H1),
                 "N" = totals.eo[,-1], "X" = X.H1, "ski" = ski$yes, "areascale" = scale_params$alpha)
mod.data.h1.alt <- list("Nyears" = Nyears, "Nsites" = Nsites, "Ncov" = ncol(X.H1.alt),
                    "N" = totals.eo[,-1], "X" = X.H1.alt, "ski" = ski$yes, "areascale" = scale_params$alpha)
mod.data.h2 <- list("Nyears" = Nyears, "Nsites" = Nsites, "Ncov" = ncol(X.H2),
                    "N" = totals.eo[,-1], "X" = X.H2, "ski" = ski$yes, "areascale" = scale_params$alpha)
mod.data.h4 <- list("Nyears" = Nyears, "Nsites" = Nsites, "Ncov" = ncol(X.H4),
                    "N" = totals.eo[,-1], "X" = X.H4, "ski" = ski$yes, "areascale" = scale_params$alpha)
mod.data.h5 <- list("Nyears" = Nyears, "Nsites" = Nsites, "Ncov" = ncol(X.H5),
                    "N" = totals.eo[,-1], "X" = X.H5, "ski" = ski$yes, "areascale" = scale_params$alpha)

#Run Model
library(rstan)
setwd("G:\\My Drive\\UNR MS 2022\\Population Dem of Draast\\Scripts\\draast_pop")
#setwd("R:\\Users\\sageellis\\Documents\\GitHub\\draast_pop") #On server

options(mc.cores = parallel::detectCores())
fit.h1 <-  stan(file = c("model_pop_ests.stan"), data = mod.data.h1, chains = 3, iter = 10000, 
             warmup = 5000, refresh=10) #  control = list(max_treedepth = 12, adapt_delta = 0.95) 
fit.h1.alt <-  stan(file = c("model_pop_ests.stan"), data = mod.data.h1.alt, chains = 3, iter = 10000, 
                warmup = 5000, refresh=10) #  control = list(max_treedepth = 12, adapt_delta = 0.95) 

fit.h2 <-  stan(file = c("model_pop_ests.stan"), data = mod.data.h2, chains = 3, iter = 10000, 
                warmup = 5000, refresh=10)
fit.h4 <-  stan(file = c("model_pop_ests.stan"), data = mod.data.h4, chains = 3, iter = 10000, 
                warmup = 5000, refresh=10)
fit.h5 <-  stan(file = c("model_pop_ests.stan"), data = mod.data.h5, chains = 3, iter = 10000, 
                warmup = 5000, refresh=10)

#Interaction models
fit.h1int <- stan(file = c("model_pop_ests_interactions.stan"), data = mod.data.h1, chains = 3, iter = 10000, 
                  warmup = 5000, refresh=10)
fit.h2int <-  stan(file = c("model_pop_ests_interactions.stan"), data = mod.data.h2, chains = 3, iter = 10000, 
                warmup = 5000, refresh=10)
fit.h4int <-  stan(file = c("model_pop_ests_interactions.stan"), data = mod.data.h4, chains = 3, iter = 10000, 
                warmup = 5000, refresh=10)
fit.h5int <-  stan(file = c("model_pop_ests_interactions.stan"), data = mod.data.h5, chains = 3, iter = 10000, 
                warmup = 5000, refresh=10)

#Model Eval----

library(shinystan)
launch_shinystan(fit.h1)

#Store model fit in hypotheses
params.h1 <- rstan::extract(fit.h1, c('n', 'Alpha','Beta0', 'Beta', 'gamma', 'sigma', 'npreds'))
params.h1.alt <- rstan::extract(fit.h1.alt, c('n', 'Alpha','Beta0', 'Beta', 'gamma', 'sigma', 'npreds'))
params.h2 <- rstan::extract(fit.h2, c('n', 'Alpha','Beta0', 'Beta', 'gamma', 'sigma', 'npreds'))
params.h4 <- rstan::extract(fit.h4, c('n', 'Alpha','Beta0', 'Beta', 'gamma', 'sigma', 'npreds'))
params.h5 <- rstan::extract(fit.h5, c('n', 'Alpha','Beta0', 'Beta', 'gamma', 'sigma', 'npreds'))

params.h1int <- rstan::extract(fit.h1int, c('n', 'Alpha','Beta0', 'Beta1','Beta2', 'gamma', 'sigma', 'npreds'))
params.h2int <- rstan::extract(fit.h2int, c('n', 'Alpha','Beta0', 'Beta1','Beta2',  'gamma', 'sigma', 'npreds'))
params.h4int <- rstan::extract(fit.h4int, c('n', 'Alpha','Beta0', 'Beta1','Beta2',  'gamma', 'sigma', 'npreds'))
params.h5int <- rstan::extract(fit.h5int, c('n', 'Alpha','Beta0', 'Beta1','Beta2',  'gamma', 'sigma', 'npreds'))

#Model Selection----
#Calculate the log likelihood of each model
ll.h1 <- log_lik(data = mod.data.h1$N, pred = params.h1$npreds, latent = params.h1$n,sigma = params.h1$sigma, iter = length(params.h1$sigma),
        scale = mod.data.h1$areascale)
ll.h1.alt <- log_lik(data = mod.data.h1.alt$N, pred = params.h1.alt$npreds, latent = params.h1.alt$n,sigma = params.h1.alt$sigma, iter = length(params.h1.alt$sigma),
                     scale = mod.data.h1.alt$areascale)
ll.h2 <- log_lik(data = mod.data.h2$N, pred = params.h2$npreds, latent = params.h2$n,sigma = params.h2$sigma, iter = length(params.h2$sigma),
                 scale = mod.data.h2$areascale)
ll.h4 <- log_lik(data = mod.data.h4$N, pred = params.h4$npreds, latent = params.h4$n,sigma = params.h4$sigma, iter = length(params.h4$sigma),
                 scale = mod.data.h4$areascale)
ll.h5 <- log_lik(data = mod.data.h5$N, pred = params.h5$npreds, latent = params.h5$n,sigma = params.h5$sigma, iter = length(params.h5$sigma),
                 scale = mod.data.h5$areascale)

ll.h1int <- log_lik(data = mod.data.h1$N, pred = params.h1int$npreds, latent = params.h1int$n,sigma = params.h1int$sigma, iter = length(params.h1int$sigma),
                 scale = mod.data.h1$areascale)
ll.h2int <- log_lik(data = mod.data.h2$N, pred = params.h2int$npreds, latent = params.h2int$n,sigma = params.h2int$sigma, iter = length(params.h2int$sigma),
                 scale = mod.data.h2$areascale)
ll.h4int <- log_lik(data = mod.data.h4$N, pred = params.h4int$npreds, latent = params.h4int$n,sigma = params.h4int$sigma, iter = length(params.h4int$sigma),
                 scale = mod.data.h4$areascale)
ll.h5int <- log_lik(data = mod.data.h5$N, pred = params.h5int$npreds, latent = params.h5int$n,sigma = params.h5int$sigma, iter = length(params.h5int$sigma),
                 scale = mod.data.h5$areascale)


#wAIC
loo::waic(ll.h1)$estimates[3]
#loo::waic(ll.h1.alt)$estimates[3] #Tested alternative to SFD using winter precip
loo::waic(ll.h2)$estimates[3]
loo::waic(ll.h4)$estimates[3]
loo::waic(ll.h5)$estimates[3]

loo::waic(ll.h1int)$estimates[3]
loo::waic(ll.h2int)$estimates[3]
loo::waic(ll.h4int)$estimates[3]
loo::waic(ll.h5int)$estimates[3]


#Calculate Lambda and CI----
##Through latent values----
lambdas <- array(numeric(),c(Nyears,length(params.h1$sigma),Nsites)) #Storage array (t by i by s)??
lambda_stoch <- matrix(NA, nrow = length(params.h1$Alpha),Nsites )

i=1;s=1;t=2
for(s in 1:Nsites){
  for(i in 1:length(params.h1$sigma)){
    for(t in 2:(Nyears)){
      lambdas[t,i,s] <- params.h1$n[i,t,s]/params.h1$n[i,t-1,s]
    }
    lambda_stoch[i,s] <- exp(mean(log(lambdas[,i,s]), na.rm = T))
  }
}

lambdas <- lambdas[-1,,] #Get rid of NA's because no transition from 2008-2009 was calculated

median.lambdas <- apply(lambdas, c(1,3), median)
apply(lambdas[,,1],1,median) == median.lambdas[,1] #Check
lb.lambdas <- apply(lambdas, c(1,3), FUN = quantile, c(0.025))
ub.lambdas <- apply(lambdas, c(1,3), FUN = quantile, c(0.975))

#Scale n back up----
n.scaled <- array(numeric(),c(15000,Nyears,Nsites))
npreds.scaled <- array(numeric(),c(15000,Nyears,Nsites))

for(s in 1:length(eo.names)){
  n.scaled[,,s] <- params.h1$n[,,s]*scale_params$alpha[s]
  npreds.scaled[,,s] <- params.h1$npreds[,,s]*scale_params$alpha[s]
}

median.n <- apply(n.scaled, c(2,3), median) #;median(n.preds[,15,8])#Double check above is what I want - yes!
lb.n <- apply(n.scaled, c(2,3), FUN = quantile, c(0.025))
ub.n <- apply(n.scaled, c(2,3), FUN = quantile, c(0.975))

#Useful numbers for paper----

#Median & CI Lambdas
round(apply(lambda_stoch, 2, median),3)
round(apply(lambda_stoch,2,FUN = quantile, c(0.025)),3)
round(apply(lambda_stoch,2, FUN = quantile, c(0.975)),3)

#Median & CI Covariate effects
apply(params.h1$Beta,2,median)
apply(params.h1$Beta,2,FUN = quantile, c(0.025))
apply(params.h1$Beta,2, FUN = quantile, c(0.975))

#Ski effects
median(params.h1$Alpha)
quantile(params.h1$Alpha, c(0.025, 0.975))

#Save Obects for Figures----
save(params.h1, scale_params, median.lambdas, lb.lambdas, ub.lambdas, median.n, lb.n, ub.n,
       totals.eo, eo.names,datayears,half.x, file = ".\\Scripts\\saved objects\\figure_objects_popests_2.RData")

#Quick Figures----

par(mfrow = c(2,2))

for(s in 1:length(eo.names)){
  eo <- params.h1$n[,,s]*scale_params$alpha[s] #scale parameter again??
  colnames(eo) <- datayears
  name <- colnames(totals.eo[s+1])
  
  median.eo <- apply(eo, 2, median) #mean across iterations 
  UB.eo <- apply(eo, 2, quantile, prob = 0.975) #upper CI across iterations
  LB.eo <- apply(eo, 2, quantile, prob = 0.025) #lower CI across iterations
  
  plot(x = datayears, y = exp(median.eo), type = 'l', col = 'black', ylim = c(0,300), main = name,
       xlab = "Years", ylab = "Population Size")
  lines(x = datayears,exp(UB.eo), lty = 2, col = 'black')
  lines(x = datayears,exp(LB.eo), lty = 2, col = 'black')
  points(x = datayears, y = totals.eo[,(s+1)], pch = 16)
  
}

#Plot Lambdas
par(mfrow = c(2,2))

half.x <- datayears[1:(Nyears-1)]+0.5

for(s in 1:length(eo.names)){
  name <- colnames(totals.eo[s+1])
  plot(x = half.x, y = median.lambdas[,s], type = 'l', col = 'black', ylim = c(0.8,1.2), main = name,
       xlab = "Years", ylab = expression(paste("Population Growth Rate (" , lambda, ")")))
  lines(x = half.x, y = ub.lambdas[,s], lty = 2, col = 'black')
  lines(x = half.x, y = lb.lambdas[,s], lty = 2, col = 'black')
  abline(h = 1, col = "red")
}

#For ESA ----

##Lambda----
eo.names#remind myself of order of plots

half.x <- datayears[1:(Nyears-1)]+0.5

tiff("lambdas_S.tiff", width = 800)

par(mfrow = c(1,1), mar = c(6,5.5,3.9,1))

#Base plot
#Wilderness
plot(x = 2009:2024,ylim = c(0.79,1.23),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = '',
     cex.lab = 3, cex.axis = 2)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis =2)
abline(h = 1, lty = 2, lwd = 0.25)

points(x = half.x, y = median.lambdas[,1], type = 'l', col = '#A32952', lwd = 2)
lines(x = half.x, y = ub.lambdas[,1], lty = 2, col = '#A32952', lwd = 1)
lines(x = half.x, y = lb.lambdas[,1], lty = 2, col = '#A32952', lwd = 1)

points(x = half.x, y = median.lambdas[,2], type = 'l', col = 'palevioletred', lwd = 2)
lines(x = half.x, y = ub.lambdas[,2], lty = 2, col = 'palevioletred', lwd = 1)
lines(x = half.x, y = lb.lambdas[,2], lty = 2, col = 'palevioletred', lwd = 1)

points(x = half.x, y = median.lambdas[,3], type = 'l', col = '#EFBDCE', lwd = 2)
lines(x = half.x, y = ub.lambdas[,3], lty = 2, col = '#EFBDCE', lwd = 1)
lines(x = half.x, y = lb.lambdas[,3], lty = 2, col = '#EFBDCE', lwd = 1)

points(x = half.x, y = median.lambdas[,6], type = 'l', col = 'aquamarine4', lwd = 2)
lines(x = half.x, y = ub.lambdas[,6], lty = 2, col = 'aquamarine4', lwd = 1)
lines(x = half.x, y = lb.lambdas[,6], lty = 2, col = 'aquamarine4', lwd = 1)


dev.off()

tiff("lambdas_N.tiff", width = 800)
par(mfrow = c(1,1), mar = c(6,7,3.9,1))


#Ski areas
plot(x = 2009:2024,ylim = c(0.79,1.23),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = expression(paste("Rate of Growth ( ", lambda," )")),
     cex.lab = 3, cex.axis = 2)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis = 2)
abline(h = 1, lty = 2, lwd = 0.25)

points(x = half.x, y = median.lambdas[,4], type = 'l', col = 'coral', lwd = 2)
lines(x = half.x, y = ub.lambdas[,4], lty = 2, col = 'coral', lwd = 1)
lines(x = half.x, y = lb.lambdas[,4], lty = 2, col = 'coral', lwd = 1)

points(x = half.x, y = median.lambdas[,5], type = 'l', col = '#FFA585', lwd = 2)
lines(x = half.x, y = ub.lambdas[,5], lty = 2, col = '#FFA585', lwd = 1)
lines(x = half.x, y = lb.lambdas[,5], lty = 2, col = '#FFA585', lwd = 1)

points(x = half.x, y = median.lambdas[,7], type = 'l', col = 'lightskyblue3', lwd = 2)
lines(x = half.x, y = ub.lambdas[,7], lty = 2, col = 'lightskyblue3', lwd = 1)
lines(x = half.x, y = lb.lambdas[,7], lty = 2, col = 'lightskyblue3', lwd = 1)

points(x = half.x, y = median.lambdas[,8], type = 'l', col = '#5591B4', lwd = 2)
lines(x = half.x, y = ub.lambdas[,8], lty = 2, col = '#5591B4', lwd = 1)
lines(x = half.x, y = lb.lambdas[,8], lty = 2, col = '#5591B4', lwd = 1)

dev.off()

#4x4 Lambda Plot----
par(mfrow=c(2,2))
#par(mfrow=c(1,1))

plot(x = 2009:2024,ylim = c(0.79,1.23),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "", ylab = expression(paste("Rate of Growth ( ", lambda," )")),
     cex.lab = 1, cex.axis = 1)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis = 1)
abline(h = 1, lty = 2, lwd = 0.25)
points(x = half.x, y = median.lambdas[,4], type = 'l', col = 'coral', lwd = 2)
lines(x = half.x, y = ub.lambdas[,4], lty = 2, col = 'coral', lwd = 1)
lines(x = half.x, y = lb.lambdas[,4], lty = 2, col = 'coral', lwd = 1)

points(x = half.x, y = median.lambdas[,5], type = 'l', col = '#FFA585', lwd = 2)
lines(x = half.x, y = ub.lambdas[,5], lty = 2, col = '#FFA585', lwd = 1)
lines(x = half.x, y = lb.lambdas[,5], lty = 2, col = '#FFA585', lwd = 1)
legend(x = 2021, y = 1.27, bty = 'n', x.intersp = 0.25,
       legend = c("Population 1", "Population 2"), fill = c("coral", "#FFA585"))

plot(x = 2009:2024,ylim = c(0.79,1.23),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "", ylab = '',
     cex.lab = 1, cex.axis = 1)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis =1)
abline(h = 1, lty = 2, lwd = 0.25)
points(x = half.x, y = median.lambdas[,7], type = 'l', col = 'lightskyblue3', lwd = 2)
lines(x = half.x, y = ub.lambdas[,7], lty = 2, col = 'lightskyblue3', lwd = 1)
lines(x = half.x, y = lb.lambdas[,7], lty = 2, col = 'lightskyblue3', lwd = 1)

points(x = half.x, y = median.lambdas[,8], type = 'l', col = '#5591B4', lwd = 2)
lines(x = half.x, y = ub.lambdas[,8], lty = 2, col = '#5591B4', lwd = 1)
lines(x = half.x, y = lb.lambdas[,8], lty = 2, col = '#5591B4', lwd = 1)
legend(x = 2021, y = 1.27, bty = 'n',x.intersp = 0.25,
       legend = c("Population 1", "Population 2"), fill = c("#5591B4", "lightskyblue3"))

plot(x = 2009:2024,ylim = c(0.79,1.3),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = expression(paste("Rate of Growth ( ", lambda," )")),
     cex.lab = 1, cex.axis = 1)
abline(h = 1, lty = 2, lwd = 0.25)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis = 1)
points(x = half.x, y = median.lambdas[,1], type = 'l', col = '#A32952', lwd = 2)
lines(x = half.x, y = ub.lambdas[,1], lty = 2, col = '#A32952', lwd = 1)
lines(x = half.x, y = lb.lambdas[,1], lty = 2, col = '#A32952', lwd = 1)

points(x = half.x, y = median.lambdas[,2], type = 'l', col = 'palevioletred', lwd = 2)
lines(x = half.x, y = ub.lambdas[,2], lty = 2, col = 'palevioletred', lwd = 1)
lines(x = half.x, y = lb.lambdas[,2], lty = 2, col = 'palevioletred', lwd = 1)

points(x = half.x, y = median.lambdas[,3], type = 'l', col = '#EFBDCE', lwd = 2)
lines(x = half.x, y = ub.lambdas[,3], lty = 2, col = '#EFBDCE', lwd = 1)
lines(x = half.x, y = lb.lambdas[,3], lty = 2, col = '#EFBDCE', lwd = 1)
legend(x = 2018.5, y = 1.35, bty = 'n', ncol = 2, x.intersp = 0.25,text.width = 1.25,
       legend = c("Population 1", "Population 2", "Population 3"), 
       fill = c("#A32952", "palevioletred", "#EFBDCE"))

plot(x = 2009:2024,ylim = c(0.79,1.23),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = '',
     cex.lab = 1, cex.axis = 1)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis = 1)
abline(h = 1, lty = 2, lwd = 0.25)
points(x = half.x, y = median.lambdas[,6], type = 'l', col = 'aquamarine4', lwd = 2)
lines(x = half.x, y = ub.lambdas[,6], lty = 2, col = 'aquamarine4', lwd = 1)
lines(x = half.x, y = lb.lambdas[,6], lty = 2, col = 'aquamarine4', lwd = 1)
legend(x = 2021, y = 1.27, bty = 'n',x.intersp = 0.25,
       legend = c("Population 1"), fill = c("aquamarine4"))



##Plot N values----
eo.names 

#Base plot

tiff("nlat_1B.tiff", width = 800)
par(mfrow = c(1,1), mar = c(5.1,6,3.9,2.1))

plot(x = 2009:2024,ylim = c(0,400),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = 'Population Size',
     cex.lab = 3, cex.axis = 2)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis = 2)

points(x = datayears, y = exp(median.n[,2]), type = 'l', lwd = 2, col = 'palevioletred')
lines(x = datayears, y = exp(lb.n[,2]),lty = 2,lwd = 2, col = 'palevioletred')
lines(x = datayears, y = exp(ub.n[,2]), lty = 2, lwd = 2, col = 'palevioletred')
points(x = datayears, totals.eo[,2+1], col = 'palevioletred', pch = 19,cex = 1.5)

dev.off()


tiff("nlat_2B.tiff", width = 800)
par(mfrow = c(1,1), mar = c(5.1,6,3.9,2.1))

plot(x = 2009:2024,ylim = c(0,400),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = 'Population Size',
     cex.lab = 3, cex.axis = 2)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis = 2)

points(x = datayears, y = exp(median.n[,4]), type = 'l', lwd = 2, col = 'coral')
lines(x = datayears, y = exp(lb.n[,4]),lty = 2,lwd = 2, col = 'coral')
lines(x = datayears, y = exp(ub.n[,4]), lty = 2, lwd = 2, col = 'coral')
points(x = datayears, totals.eo[,4+1], col = 'coral', pch = 19,cex = 1.5)

dev.off()

tiff("nlat_3B.tiff", width = 800)
par(mfrow = c(1,1), mar = c(5.1,6,3.9,2.1))

plot(x = 2009:2024,ylim = c(0,400),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = 'Population Size',
     cex.lab = 3, cex.axis = 2)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis = 2)

points(x = datayears, y = exp(median.n[,6]), type = 'l', lwd = 2, col = 'aquamarine4')
lines(x = datayears, y = exp(lb.n[,6]),lty = 2,lwd = 2, col = 'aquamarine4')
lines(x = datayears, y = exp(ub.n[,6]), lty = 2, lwd = 2, col = 'aquamarine4')
points(x = datayears, totals.eo[,6+1], col = 'aquamarine4', pch = 19, cex = 1.5)
dev.off()

tiff("nlat_Bruce.tiff", width = 800)
par(mfrow = c(1,1), mar = c(5.1,6,3.9,2.1))
plot(x = 2009:2024,ylim = c(0,400),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = 'Population Size',
     cex.lab = 3, cex.axis = 2)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1, cex.axis = 2)

points(x = datayears, y = exp(median.n[,8]), type = 'l', lwd = 2, col = 'lightskyblue3')
lines(x = datayears, y = exp(lb.n[,8]),lty = 2,lwd = 2, col = 'lightskyblue3')
lines(x = datayears, y = exp(ub.n[,8]), lty = 2, lwd = 2, col = 'lightskyblue3')
points(x = datayears, totals.eo[,8+1], col = 'lightskyblue3', pch = 19,cex = 1.5)
dev.off()

##Alpha----
par(mar = c(5.1,4.2,4.1,2.1))

alpha.plot <- data.frame(vals = params.h1$Alpha) %>% 
                ggplot(.,aes(vals)) +
                geom_density(fill = "#91B1CD", color = "#91B1CD", alpha = 0.75)+
                geom_vline(aes(xintercept = quantile(params.h1$Alpha, 0.025)), col = "#91B1CD",
                           linetype = 'dashed', size = 1)+
                geom_vline(aes(xintercept = quantile(params.h1$Alpha, 0.975)), col = "#91B1CD",
                           linetype = 'dashed', size = 1)+
                scale_y_continuous(expand = c(0,0),
                                   limits = c(0,3)) +
                scale_x_continuous(limits = c(-0.7,0.6))+
                labs(x = "Alpha", y = "Density")+
                theme(text = element_text(size = 30),
                      plot.margin = margin(20,5.5,5.5,5.5, 'pt'),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))

alpha.plot

ggsave("alpha_plot.png", alpha.plot, dpi = 600)

##Plot violin effects----
covars 

cov_data <- data.frame("Average Winter Temperature" = params.h1$Beta[,1],
                       "Total Summer Precipitation" = params.h1$Beta[,2],
                       "Snow Free Day" = params.h1$Beta[,3],
                       check.names = F)

cov_data <- cov_data %>% pivot_longer(cols = colnames(cov_data),
                         names_to = "Covariate",
                         values_to = "Value")

data_summary <- function(x){
  x <- unlist(x)
  m <- median(x)
  ymin <- unname(quantile(x, 0.025))
  ymax <- unname(quantile(x, 0.975))
  return(c(y=m,ymin=ymin,ymax=ymax))
}


param.violin <- ggplot(cov_data, aes(y = Value, x =Covariate, fill=Covariate)) +
  geom_violin()+
  scale_fill_manual(values=c("#ed6a5a","#E1BE6A",  "#9bc1bc"), name = "Covariates")+
  stat_summary(fun = median, show.legend = F)+
  theme(text = element_text(size = 30),
        #legend.position = 'none',
        legend.key.spacing.y = unit(20, 'pt'),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

param.violin 

ggsave("covar_plot.png",param.violin, width = 10, height = 7,dpi = 800)

#Nice Figures----

eo.names#remind myself of order of plots

half.x <- datayears[1:(Nyears-1)]+0.5

##Combine Lambda figures ski vs not ski first----

#Base plot
plot(x = 2009:2024,ylim = c(0.85,1.15),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = expression(paste("Rate of Growth ( ", lambda," )")),
     cex.lab = 1.25)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1)
abline(h = 1, lty = 2, lwd = 0.25)

#Add ski plots
points(x = half.x, y = median.lambdas[,7], type = 'l', col = 'palevioletred')
lines(x = half.x, y = ub.lambdas[,7], lty = 2, col = 'palevioletred')
lines(x = half.x, y = lb.lambdas[,7], lty = 2, col = 'palevioletred')

points(x = half.x, y = median.lambdas[,8], type = 'l', col = 'palevioletred')
lines(x = half.x, y = ub.lambdas[,8], lty = 2, col = 'palevioletred')
lines(x = half.x, y = lb.lambdas[,8], lty = 2, col = 'palevioletred')

points(x = half.x, y = median.lambdas[,4], type = 'l', col = 'aquamarine4')
lines(x = half.x, y = ub.lambdas[,4], lty = 2, col = 'aquamarine4')
lines(x = half.x, y = lb.lambdas[,4], lty = 2, col = 'aquamarine4')

points(x = half.x, y = median.lambdas[,5], type = 'l', col = 'aquamarine4')
lines(x = half.x, y = ub.lambdas[,5], lty = 2, col = 'aquamarine4')
lines(x = half.x, y = lb.lambdas[,5], lty = 2, col = 'aquamarine4')

#Base plot
plot(x = 2009:2024,ylim = c(0.5,1.7),  xlim = c(2009,2024), xaxt = 'n',
     xlab = "Years", ylab = expression(paste("Rate of Growth ( ", lambda," )")),
     cex.lab = 1.25)
axis(1, at=2009:2024, labels=c(2009:2024), las = 1)
abline(h = 1, lty = 2, lwd = 0.25)

#Other plot
points(x = half.x, y = median.lambdas[,1], type = 'l', col = 'lightskyblue3')
lines(x = half.x, y = ub.lambdas[,1], lty = 2, col = 'lightskyblue3')
lines(x = half.x, y = lb.lambdas[,1], lty = 2, col = 'lightskyblue3')

points(x = half.x, y = median.lambdas[,2], type = 'l', col = 'lightskyblue3')
lines(x = half.x, y = ub.lambdas[,2], lty = 2, col = 'lightskyblue3')
lines(x = half.x, y = lb.lambdas[,2], lty = 2, col = 'lightskyblue3')

points(x = half.x, y = median.lambdas[,3], type = 'l', col = 'lightskyblue3')
lines(x = half.x, y = ub.lambdas[,3], lty = 2, col = 'lightskyblue3')
lines(x = half.x, y = lb.lambdas[,3], lty = 2, col = 'lightskyblue3')

points(x = half.x, y = median.lambdas[,6], type = 'l', col = 'coral')
lines(x = half.x, y = ub.lambdas[,6], lty = 2, col = 'coral')
lines(x = half.x, y = lb.lambdas[,6], lty = 2, col = 'coral')





