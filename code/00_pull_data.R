library(tidyverse)
library(here)

# pull data from website
data = readr::read_csv(
  "https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", 
  guess_max = 1e6)

# save to local
readr::write_csv(data, here::here("./data/tick-targets.csv"))
