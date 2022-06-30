covs <- read.csv("daymetChallengeSites.csv")

covs$Date <- as.Date(covs$Date)
#changing Date to time
names(covs)[names(covs) == "Date"] <- "time"

ticks <- readr::read_csv(
  "https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", 
  guess_max = 1e6)
#NAing 2019 and deleting 2020
ticks$time1 <- as.POSIXct(ticks$time, format = "%m/%d/%Y")
ticks$Year <- format(ticks$time1, format="%Y")
#We retain data no affected by COVID collection
ticks=subset(ticks,Year<2020)
#Data for cross-validation
ticks2=ticks[which(ticks$Year==2019),]

#Eliminate data that we will use later for cross-validation
ticks[which(ticks$Year==2019), "amblyomma_americanum"]=NA

ticks$time <- as.Date(ticks$time)

#merge with OG ticks data
tickscov <- merge(ticks, covs, by = c("time", "siteID"))