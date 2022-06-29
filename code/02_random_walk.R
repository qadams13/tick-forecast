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
time = ticks$mmwrWeek
siteID = ticks$siteID

data <- list(y=log(y+1),n=length(y), nsite = length(unique(siteID)),## data
             x_ic=log(1000),tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

# random walk 
RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:n){

    y[t] ~ dpois(x[t]) #t different y's (observations) and they are a function of x
  }
  
  
  #### Process Model
  for(t in 2:n){
    x[t]~pois(x[t-1])
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic) #initial condition for the state 
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),  ## initial guess on process precision
                    tau_obs=5/var(log(y.samp)))        ## initial guess on obs precision
}

j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)
plot(jags.out)

# larger samples
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)
summary(jags.out)

# plotting results
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags.out)         ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="tick density",xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), log='y', format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)
