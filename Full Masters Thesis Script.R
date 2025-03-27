library(rlang)
library(tidyverse)
library(lubridate)
library(tidybayes)
library(tidybayes.rethinking)
library(bayesplot)
library(devtools)
library(brms)
library(dplyr)
library(ggplot2)
library(rethinking)
library(ggpubr)
library(randomForest)
library(plyr)
library(mgcv)

#I first loaded the raw data files from the Utqiagvik ARM station on NOAA's website,
#then converted those into a data file that R could read.
#Data is already converted to .Rds file, accessible through the figshare link on Github.

read_minute_data <- function(file_name) {
  col_names <- c("site", "year", "month", "day", "hour", "minute", "wind_dir", "wind_spd", "wind_stead", "p", "t_2m", "t_10m", "t_top", "rh", "precip")
  
  col_types = "ciiiiiddddddddd"
  
  df <- read_table(file_name, 
                   col_names = col_names, col_types = col_types) %>%
    mutate(across(c(wind_dir, wind_spd, p, t_2m, t_10m, t_top), 
                  ~ifelse(.x <= -999, NA_real_, .x)), 
           wind_stead = ifelse(wind_stead <= -9, NA_real_, wind_stead), 
           rh = ifelse(rh <= -99, NA_real_, rh), 
           precip = ifelse(precip <= -99, NA_real_, precip)
    )
  invisible(df)
}


read_hour_data <- function(file_name) {
  col_names <- c("site", "year", "month", "day", "hour", "wind_dir", "wind_spd", "wind_stead", "p", "t_2m", "t_10m", "t_top", "rh", "precip")
  
  col_types = "ciiiiddddddddd"
  
  df <- read_table(file_name, 
                   col_names = col_names, col_types = col_types) %>%
    mutate(across(c(wind_dir, wind_spd, p, t_2m, t_10m, t_top), 
                  ~ifelse(.x <= -99, NA_real_, .x)), 
           wind_stead = ifelse(wind_stead <= -9, NA_real_, wind_stead), 
           rh = ifelse(rh <= -99, NA_real_, rh), 
           precip = ifelse(precip <= -99, NA_real_, precip)
    )
  invisible(df)
}

minute_data <- read_rds("data/minute_data.Rds")

updated_minute_summary <- minute_data %>% 
  mutate(wind_spd = ifelse(wind_spd < 0, NA, wind_spd)) %>%
  group_by(site, year, month, day, hour) %>% 
  summarize(wind_mean = mean(wind_spd), 
            gust_time = sum(wind_spd >= 10 & wind_spd >= (wind_mean + 0.5*wind_mean)), 
            wind_sd = sd(wind_spd), 
            wind_max = max(wind_spd), wind_min = min(wind_spd), 
            wind_med = median(wind_spd), nrow = n(), 
            na_ct = sum(is.na(wind_spd)),
            .groups = "drop")

updated_clean_minute_summary <- updated_minute_summary %>% filter(nrow == 60, na_ct == 0) %>% 
  mutate(date = make_datetime(year, month, day, hour))

full_dates <- updated_clean_minute_summary$date

hour_data <- read_rds("data/hour_data.Rds")

combined_data <- hour_data %>%
  left_join(updated_clean_minute_summary, by = c("site", "year",  "month", "day", "hour"))

combined_cc <- combined_data %>% drop_na(t_2m, wind_spd) #leave only complete cases. Measurements missing that time's temperature and/or wind speed data are excluded.

new_combined_data <- mutate(combined_cc, t_difference = t_10m-t_2m) #calculate temperature difference between 10m above sea level and 2m above sea level
new_combined_data_df <- as_tibble(new_combined_data)

new_combined_cc <- new_combined_data %>% 
  mutate(date = make_datetime(year, month, day, hour)) %>% 
  filter(date %in% full_dates)

# Calculate bulk Richardson number from existing parameters and add to dataset
#First find potential temperature using the hydrostatic equation

updated_Richardson_combined_cc <- mutate(new_combined_cc, RichardsonBulk = ((9.81/t_2m)*(t_difference)*(8))/(wind_spd^2))

#calculate variables needed to find the bulk richardson number
t_average_cc <- mutate(updated_Richardson_combined_cc, avg_temp_K = ((t_2m+t_10m)/2)+273.15) #add a variable for average temperature across the vertical column
pressurechange_cc <- mutate(t_average_cc, DeltaP = (-9.81*(p/(8.314*avg_temp_K))*8)) #calculate pressure change over the air column
rm(t_average_cc)
p10_cc <- mutate(pressurechange_cc, p_10m = (DeltaP+p)) #calculate pressure 10 m above sea level
rm(pressurechange_cc)
theta2m_cc <- mutate(p10_cc, PotTemp_2m = (t_2m*(1000/p)^0.286)) #calculate potential temperature at 2 m above sea level
rm(p10_cc)
theta10m_cc <- mutate(theta2m_cc, PotTemp_10m = (t_10m*(1000/p_10m)^0.286)) #calculate potential temperature at 10 m above sea level
rm(theta2m_cc)
deltatheta_cc <- mutate(theta10m_cc, DeltaTheta = PotTemp_10m-PotTemp_2m) #calculate change in potential temperature
rm(theta10m_cc)
updated_Richardson_combined_cc_new <- mutate(deltatheta_cc, RichardsonBulkPotential = ((9.81/t_2m)*(DeltaTheta)*(8))/(wind_spd^2)) #calculate the richardson bulk potential
