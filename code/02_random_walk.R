#Load necessary libraries
library(readr)
library(rjags)
library(lubridate)
#library(rnoaa)
library(daymetr)
#devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
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

# x = matrix(nrow = length(unique(ticks$mmwrWeek)),
#            ncol = length(unique(ticks$Year)))
# years_vec = unique(ticks$Year)
# week_vec = sort(unique(ticks$mmwrWeek))
# for(i in 1:length(years_vec)){
#   for(j in 1:length(week_vec)) {
#     x[j,i] = ticks[which(ticks$mmwrWeek == week_vec[j] & 
#                            ticks$amblyomma_americanum == years_vec[i]), 
#                    "amblyomma_americanum"]
#   }
# }

density = ticks$amblyomma_americanum
years_vec = as.integer(as.factor(ticks$Year))
site_id_vec = ticks$siteID
week_vec = as.integer(as.factor(ticks$mmwrWeek))
n_obs = length(density)
x_rows = sort(unique(week_vec))
x_cols = sort(unique(years_vec))

sort(unique(week_vec))

# create new df with NA's filled in 
length_new_df = length(unique(ticks$mmwrWeek)) * 
  length(unique(ticks$siteID)) * length(unique(ticks$Year))
new_df = data.frame(
  week = numeric(length = length_new_df),
  year = numeric(length = length_new_df),
  density = numeric(length = length_new_df),
  site = character(length = length_new_df)
)

first_week = numeric(length = length(unique(ticks$Year)))
last_week = numeric(length = length(unique(ticks$Year)))

j = 1
for(i in sort(unique(ticks$Year))) {
  first_week[j] = min(
    ticks[which(ticks$Year == i), "mmwrWeek"]
  )
  last_week[j] = max(
    ticks[which(ticks$Year == i), "mmwrWeek"]
  )
  j = j + 1
}

year_iter = 1
for(year in sort(unique(ticks$Year))) {
  first = first_week[year_iter]
  last = last_week[year_iter]
  for(week in first:last) {
    for(site in sort(unique(ticks$siteID))) {
      
      # check if value exists
      value = ifelse(nrow(ticks[which(
          ticks$mmwrWeek == week & 
            ticks$Year == year &
              ticks$siteID == site), ]) > 0, 
          ticks[which(
            ticks$mmwrWeek == week & 
              ticks$Year == year &
              ticks$siteID == site), "amblyomma_americanum"][[1]],
          NA)
      new_df[i, "week"] = as.numeric(week)
      new_df[i, "site"] = site
      new_df[i, "year"] = as.numeric(year)
      new_df[i, "density"] = as.numeric(value)
      i = i + 1
    }
  }
}

# test 
week = unique(ticks$mmwrWeek)[1]
site = unique(ticks$siteID)[1]
year = unique(ticks$Year)[1]


data <- list(y = new_df$density, 
             years_vec,
             site_id_vec,
             week_vec,
             n_obs,
             x_ic_alpha = 1.5, 
             x_beta=80, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

# look at distribution of the data
hist(ticks$amblyomma_americanum)

# random walk 
RandomWalk = "
model{

  #### Priors
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
  
  #### Process Model
  # loop for each year
  for(i in unique(years_vec)) {
        
        x[1,i] ~ dgamma(x_ic_alpha, x_beta) # give the first value as the prior
        alpha.y[i] ~ dnorm(0, tau_add)
        
    for(t in sort(unique(week_vec))[2:length(unique(week_vec))]) {
      t in sort(unique(week_vec))[2:length(unique(week_vec)
      log(x[t, i]) <- x[t-1, i] + alpha.y[i]
    }
  }
  
  #### Data Model
  
  for(i in 1:n_obs){
      y[i] ~ dpois() #t different y's (observations) and they are a function of x
  }
}
"

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),  ## initial guess on process precision
                    tau_obs=5/var(log(y.samp)))        ## initial guess on obs precision
}

j.model <- jags.model(file = textConnection(RandomWalk),
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
