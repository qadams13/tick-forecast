#Load necessary libraries
library(readr)
library(rjags)
library(lubridate)
#library(rnoaa)
library(daymetr)
devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
library(ecoforecastR)

#Reading in data
ticks <- readr::read_csv(here::here(
  "./data/tick-targets.csv"
), guess_max = 1e6)

#NAing 2019 and deleting 2020
ticks$time <- as.POSIXct(ticks$time, format = "%m/%d/%Y")
ticks$Year <- format(ticks$time, format="%Y")

#We retain data no affected by COVID collection
ticks=subset(ticks,Year<2020)

#Data for cross-validation
ticks2=ticks[which(ticks$Year==2019),]
dim(ticks2)

#Eliminate data that we will use later for cross-validation
ticks[which(ticks$Year==2019), "amblyomma_americanum"]=NA
dim(ticks)
y = ticks$amblyomma_americanum
y = 
time = ticks$mmwrWeek
siteID = ticks$siteID

data <- list(y=log(y+1),n=length(y), nsite = length(unique(siteID)),## data
             x_ic=log(1000),tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

SpatialRandomWalk = "
model{

  #### Data Model
  for(t in 1:n){
    for(i in 1:nsite){
      y[i,t] ~ dpois(x[i,t])
    }
  }
  
  #### Process Model
  for(t in 2:n){
    for(i in 1:nsite){
      mu[i,t] <- x[i,t-1] + alpha * sum(adj[i,1:nsite]*x[1:nsite,t-1])
    }
    x[1:nsite,t] ~ dmnorm(mu[1:nsite,t],Omega_proc)
  }
  
  #### Priors
  for(i in 1:nsite){
    x[i,1] ~ dnorm(x_ic,tau_ic)
  }
  tau_obs ~ dgamma(a_obs,r_obs)
  Omega_proc ~ dwish(R,k)
  alpha ~ dbeta(1,20)
}
"

nchain = 3
init = list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),  ## initial guess on process precision
                    tau_obs=5/var(log(y.samp)))        ## initial guess on obs precision
}

j.model = jags.model (file = textConnection(SpatialRandomWalk),
                      data = data,
                      inits = init,
                      n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)