library(rjags)
#library(rnoaa)
library(daymetr)
devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
library(ecoforecastR)

gflu = read.csv("data/gflu_data.txt",skip=11)
time = as.Date(gflu$Date)
y = gflu$Massachusetts
plot(time,y,type='l',ylab="Flu Index",lwd=2,log='y')

RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"

gflu = read.csv("data/gflu_data.txt",skip=11)
time = as.Date(gflu$Date)
y = gflu$Massachusetts
data <- list(y=log(y),n=length(y),      ## data
             x_ic=log(1000),tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),  ## initial guess on process precision
                    tau_obs=5/var(log(y.samp)))        ## initial guess on obs precision
}

j.model = jags.model(file = textConnection(RandomWalk),
                     data = data,
                     inits = init,
                     n.chains = 3)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)
plot(jags.out)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)

time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags.out)         ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Flu Index",log='y',xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)

hist(1/sqrt(out[,1]),main=colnames(out)[1])
hist(1/sqrt(out[,2]),main=colnames(out)[2])

plot(out[,1],out[,2],pch=".",xlab=colnames(out)[1],ylab=colnames(out)[2])
cor(out[,1:2])

