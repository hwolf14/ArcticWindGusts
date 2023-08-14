#instead of using gust_time, use the hourly standard deviation of the wind speed

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

minute_summary <- minute_data %>% 
  mutate(wind_spd = ifelse(wind_spd < 0, NA, wind_spd)) %>%
  group_by(site, year, month, day, hour) %>% 
  summarize(wind_mean = mean(wind_spd), 
            gust_time = sum(wind_spd >= 10), 
            wind_sd = sd(wind_spd), 
            wind_max = max(wind_spd), wind_min = min(wind_spd), 
            wind_med = median(wind_spd), nrow = n(), 
            na_ct = sum(is.na(wind_spd)),
            .groups = "drop")

clean_minute_summary <- minute_summary %>% filter(nrow == 60, na_ct == 0) %>% 
  mutate(date = make_datetime(year, month, day, hour))

full_dates <- clean_minute_summary$date

#same as before

combined <- Hourly %>%
  left_join(clean_minute_summary, by = c("site", "year",  "month", "day", "hour"))

combined_cc <- combined %>% drop_na(gust_time, t_2m)

new_combined_cc <- combined_cc %>% 
  mutate(date = make_datetime(year, month, day, hour)) %>% 
  filter(date %in% full_dates)

new_combined_cc <- new_combined_cc %>% 
  left_join(select(clean_minute_summary, site, date, new_gust_time = gust_time, 
                   wind_mean, wind_sd, wind_max, wind_min, wind_med), 
            by = c("site", "date"))

rm(minute_data)
gc()

#Build mgusty1

mgusty1_updated <- quap(
  alist(
    wind_sd.x ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , data=new_combined_cc)

t_2m.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_2m = t_2m.seq )
mu <- link( mgusty1_updated , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.wind_sd <- sim( mgusty1_updated , data=new_combined_cc )
wind_sd.PI <- apply( sim.wind_sd , 2 , PI , prob=0.90 )
rm(sim.wind_sd)
gc()

#plot with mgusty1

plot( wind_sd.x ~ t_2m , new_combined_cc , col=col.alpha(rangi2,0.5) )
lines( t_2m.seq , mu.mean )
shade( mu.PI , t_2m.seq )

mgusty1plot_updated <- ggplot(new_combined_cc, aes(t_2m, wind_sd.x)) +
  geom_smooth()
mgusty1plot_updated + coord_cartesian(xlim = c(-30, 10))

gam_updated <- gam((wind_sd.x) ~ s(t_2m + t_2m^2), data = new_combined_cc)

plot(gam_updated, add = TRUE, ylim = c(0,10), ylab = "Hourly Standard Deviation", xlab = "temperature")
points(wind_sd.x ~ t_2m, data = slice_sample(new_combined_cc, n = 2000), col = alpha(rangi2, 0.5))
lines(t_2m.seq, mu.mean)
shade(mu.PI, t_2m.seq)

precis(mgusty1_updated)

#Build mgusty2

combined_cc2 <- combined %>% drop_na(gust_time, wind_spd)

new_combined_cc2 <- combined_cc2 %>% 
  mutate(date = make_datetime(year, month, day, hour)) %>% 
  filter(date %in% full_dates)

new_combined_cc2 <- new_combined_cc2 %>% 
  left_join(select(clean_minute_summary, site, date, new_gust_time = gust_time, 
                   wind_mean, wind_sd, wind_max, wind_min, wind_med), 
            by = c("site", "date"))

mgusty2_updated <- quap(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * wind_spd^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , data=new_combined_cc2)

#plotting with mgusty2_updated

precis(mgusty2_updated)

wind_spd.seq <- seq( from=0 , to=100 , length.out=30 )
pred_dat <- list( wind_spd_s=wind_spd.seq , wind_spd_s2=wind_spd.seq^2 )
mu <- link( mgusty2_updated , data=new_combined_cc2 )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory
sim.wind_sd <- sim( mgusty2_updated , data=new_combined_cc2 )
wind_sd.PI <- apply( sim.wind_sd , 2 , PI , prob=0.90 )
rm(sim.wind_sd)
gc()

plot( wind_sd ~ wind_spd , new_combined_cc2 , col=col.alpha(rangi2,0.5) )
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

ggplot(new_combined_cc2, aes(wind_spd, wind_sd)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

mgusty2plot_updated <- ggplot(new_combined_cc2, aes(wind_spd, wind_sd)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

mgusty2plot_updated + coord_cartesian(xlim = c(0, 20))

#Build mgusty3

dat_slim_updated <- list(
  wind_sd = new_combined_cc$wind_sd.x,
  t_2m = new_combined_cc$t_2m
)
str(dat_slim_updated)

mgusty3_updated <- ulam(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(0.5, 2),
    b1 ~ dnorm(0, 0.1),
    b2 ~ dnorm(0, 0.1),
    sigma ~ dexp(1)
  ) ,
  data=dat_slim_updated , chains=1 )

mgusty3.1_updated <- ulam(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(1, 1),
    b1 ~ dnorm(0, 0.3),
    b2 ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ) ,
  data=dat_slim_updated , chains=4 , cores=4 , iter = 1000 )

##plot with mgusty3_updated and mgusty3.1_updated

t_2m.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_2m = t_2m.seq )
mu <- link( mgusty3.1_updated , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.wind_sd <- sim( mgusty3.1_updated , data=dat_slim_updated )
wind_sd.PI <- apply( sim.wind_sd , 2 , PI , prob=0.90 )
rm(sim.wind_sd)
gc()

dat_slim_df_updated <- as_tibble(dat_slim_updated)

plot( wind_sd ~ t_2m , dat_slim_df_updated , col=col.alpha(rangi2,0.5) )
lines( t_2m.seq , mu.mean )
shade( mu.PI , t_2m.seq )

ggplot(dat_slim_df_updated, aes(t_2m, wind_sd)) +  geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))
ggplot(dat_slim_df_updated, aes(t_2m, wind_sd)) + geom_point(color = alpha(rangi2, 0.1)) + geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))

plot3.1_updated <- ggplot(dat_slim_df_updated, aes(t_2m, wind_sd)) + geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))
plot3.1_updated + coord_cartesian(ylim = c(0, 10))

#Build new mgusty4

dat_slim_2_updated <- list(
  wind_sd = new_combined_cc2$wind_sd,
  wind_spd = new_combined_cc2$wind_spd
)
str(dat_slim_2_updated)

mgusty4_updated <- ulam(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * wind_spd^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 0.3),
    b2 ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ) , 
  data=dat_slim_2_updated , chains=1 )

mgusty4.1_updated <- ulam(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * wind_spd^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 0.3),
    b2 ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ) , 
  data=dat_slim_2_updated , chains=4 , cores=4 , iter = 1000)

#plotting with both ulam models

library(tidybayes)
library(tidybayes.rethinking)
library(bayesplot)
mcmc_dens(mgusty4.1_updated@stanfit)

post_mgusty3_updated <- extract.samples(mgusty3_updated)
post_mgusty3.1_updated <- extract.samples(mgusty3.1_updated)

dens(post_mgusty3_updated$sigma, lwd = 1)
dens(post_mgusty3.1_updated$sigma, lwd = 1, col = rangi2, add = TRUE)

mgusty3_draws_updated <- tidy_draws(mgusty3_updated) %>% mutate(model = "mgusty3_updated")
mgusty3.1_draws_updated <- tidy_draws(mgusty3.1_updated) %>% mutate(model = "mgusty3.1_updated")

bind_rows(mgusty3_draws_updated, mgusty3.1_draws_updated) %>% 
  ggplot(aes(x = sigma, color = model, fill = model)) + 
  geom_density(size = 1, alpha = 0.3)

post_mgusty4_updated <- extract.samples(mgusty4_updated)
post_mgusty4.1_updated <- extract.samples(mgusty4.1_updated)

dens(post_mgusty4$sigma, lwd = 1)
dens(post_mgusty4.1$sigma, lwd = 1, col = rangi2, add = TRUE)

mgusty4_draws <- tidy_draws(mgusty4_updated) %>% mutate(model = "mgusty4_updated")
mgusty4.1_draws <- tidy_draws(mgusty4.1_updated) %>% mutate(model = "mgusty4.1_updated")

bind_rows(mgusty4_draws_updated, mgusty4.1_draws_updated) %>% 
  ggplot(aes(x = sigma, color = model, fill = model)) + 
  geom_density(size = 1, alpha = 0.3)

wind_spd.seq <- seq( from=0 , to=100 , length.out=30 )
pred_dat <- list( wind_spd_s=wind_spd.seq , wind_spd_s2=wind_spd.seq^2 )
mu <- link( mgusty4.1_updated , data=dat_slim_2_updated )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory
sim.wind_sd <- sim( mgusty4.1_updated , data=dat_slim_2_updated )
wind_sd.PI <- apply( sim.wind_sd , 2 , PI , prob=0.90 )
rm(sim.wind_sd)
gc()

dat_slim_2_df_updated <- as_tibble(dat_slim_2_updated)

ggplot(dat_slim_2_df_updated, aes(wind_spd, wind_sd)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

plot4.1_updated <- ggplot(dat_slim_2_df_updated, aes(wind_spd, wind_sd)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

plot4.1_updated + coord_cartesian(xlim = c(0, 20))

#Build new mgusty5, temperature difference model

newest_combined_cc <- new_combined_cc
print(newest_combined_cc)

newest_combined_cc <- mutate(newest_combined_cc, t_difference = t_10m-t_2m)

mgusty5_updated <- quap(
  alist(
    wind_sd.x ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , data=newest_combined_cc)

#plot with mgusty5 to graph t_difference

t_difference.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_difference = t_difference.seq )
mu <- link( mgusty5_updated , data=newest_combined_cc )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.wind_sd <- sim( mgusty5_updated , data=newest_combined_cc )
wind_sd.PI <- apply( sim.wind_sd , 2 , PI , prob=0.90 )
rm(sim.wind_sd)
gc()

ggplot(newest_combined_cc, aes(t_difference, wind_sd)) +
  geom_smooth()

ggplot(data = newest_combined_cc, aes (x = t_difference)) + geom_density() + xlim(-2, 2)

ggplot(newest_combined_cc, aes(x=t_difference, y=wind_sd.x)) + 
  geom_bin2d(bins = 100, mapping = aes(fill = log(..ndensity..))) +
  theme_bw()

ggplot(newest_combined_cc, aes(x=t_difference, y=wind_sd.x)) + 
  geom_bin2d(bins = c(30,30), mapping = aes(fill = log(..ndensity..))) +
  scale_fill_viridis_c(option = "B") +
  xlim(-2,2) +
  ylim(0,5)
  theme_bw()

#Build new mgusty6

mgusty6_updated <- ulam(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 0.3),
    b2 ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ) ,
  data=dat_slim_updated , chains=4 , cores=4 , iter = 1000 )

# plot with mgusty6_updated

t_difference.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_difference = t_difference.seq )
mu <- link( mgusty6_updated , data=newest_combined_cc )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.wind_sd <- sim( mgusty6_updated , data=newest_combined_cc )
wind_sd.PI <- apply( sim.wind_sd , 2 , PI , prob=0.90 )
rm(sim.wind_sd)
gc()

newest_combined_cc_df <- as_tibble(newest_combined_cc)

plot( wind_sd.x ~ t_difference , newest_combined_cc_df , col=col.alpha(rangi2,0.1))
lines( t_difference.seq , mu.mean )
shade( mu.PI , t_difference.seq )

ggplot(newest_combined_cc_df, aes(t_difference, wind_sd)) +
  geom_smooth() + xlim(-2, 2) + ylim (0, 1)

# code updated mgusty models 3-5 into brms

mgusty3.1prior_updated <- c(
  set_prior("normal(5, 5)", class = "Intercept", lb = 0),
  set_prior("normal(0, 1)", class = "b", coef = "t_2m"),
  set_prior("exponential(0.05)", class="sigma"))
mgusty3.1_fit_formula_updated <- bf(wind_sd.x ~ t_2m + I(t_2m^2))
b3.1_updated <-
  brm(mgusty3.1_fit_formula_updated,
      data = new_combined_cc,
      family = gaussian,
      prior = mgusty3.1prior_updated,
      chains = 4, cores = 4, backend = "cmdstanr" )

mgusty4.1prior_updated <- c(
  set_prior("normal(5, 5)", class = "Intercept", lb = 0),
  set_prior("normal(0, 1)", class = "b", coef = "wind_spd"),
  set_prior("exponential(0.05)", class="sigma"))
mgusty4.1_fit_formula_updated <- bf(wind_sd.x ~ wind_spd + I(wind_spd^2))
b4.1_updated <-
  brm(mgusty4.1_fit_formula_updated,
      data = new_combined_cc,
      family = gaussian,
      prior = mgusty4.1prior,
      chains = 4, cores = 4, backend = "cmdstanr" )

mgusty5.1prior_updated <- c(
  set_prior("normal(5, 5)", class = "Intercept", lb = 0),
  set_prior("normal(0, 1)", class = "b", coef = "t_difference"),
  set_prior("exponential(0.05)", class="sigma"))
mgusty5.1_fit_formula_updated <- bf(wind_sd ~ t_difference + I(t_difference^2))
b5.1_updated <-
  brm(mgusty5.1_fit_formula_updated,
      data = newest_combined_cc_df,
      family = gaussian,
      prior = mgusty5.1prior_updated,
      chains = 4, cores = 4, backend = "cmdstanr" )

# Calculate bulk Richardson number from existing parameters and add to dataset
#First find potential temperature using the hydrostatic equation

Richardson_combined_cc <- mutate(newest_combined_cc, RichardsonBulk = ((9.81/t_2m)*(t_difference)*(8))/(wind_spd^2))

t_average_cc <- mutate(Richardson_combined_cc, avg_temp_K = ((t_2m+t_10m)/2)+273.15)
pressurechange_cc <- mutate(t_average_cc, DeltaP = (-9.81*(p/(8.314*avg_temp_K))*8))
rm(t_average_cc)
p10_cc <- mutate(pressurechange_cc, p_10m = (DeltaP+p))
rm(pressurechange_cc)
theta2m_cc <- mutate(p10_cc, PotTemp_2m = (t_2m*(1000/p)^0.286))
rm(p10_cc)
theta10m_cc <- mutate(theta2m_cc, PotTemp_10m = (t_10m*(1000/p_10m)^0.286))
rm(theta2m_cc)
deltatheta_cc <- mutate(theta10m_cc, DeltaTheta = PotTemp_10m-PotTemp_2m)
rm(theta10m_cc)
Richardson_combined_cc_new <- mutate(deltatheta_cc, BulkRichardsonNumber = ((9.81/t_2m)*(DeltaTheta)*(8))/(wind_spd^2))

# Plot bulk Richardson number
ggplot(Richardson_combined_cc_new, aes(x=BulkRichardsonNumber, y=wind_sd.x)) +
  geom_bin2d(bins = c(30,30), mapping = aes(fill = log(..ndensity..))) +
  scale_fill_viridis_c(option = "B") +
  xlim(-1,1) +
  ylim(0,10) +
  theme_bw()

ggplot(Richardson_combined_cc_new, aes(x=RichardsonBulkPotential, y=wind_sd)) +
  xlim(-1,1) +
  geom_point(alpha = 0.1, na.rm = TRUE)

# Filter out obvious errors in raw data
new_combined_cc <- new_combined_cc[-c(36543, 36544, 42809, 42810, 42954:42976, 47563:47582, 137149:137220), ]
# filter out anything with sd above 5 and any hour with wind max > 35 m/s
new_combined_cc <- new_combined_cc[ new_combined_cc$wind_sd.x <= 5 & new_combined_cc$wind_max.x <= 35, ]
new_combined_cc2 <- new_combined_cc2[ new_combined_cc2$wind_sd <= 5 & new_combined_cc2$wind_max <= 35, ]

#try a multiple regression model

dat_slim_MR <- list(
  wind_sd = new_combined_cc$wind_sd.x,
  t_2m = new_combined_cc$t_2m,
  wind_spd = new_combined_cc$wind_spd
)
str(dat_slim_MR)
dat_slim_MR_df <- as_tibble(dat_slim_MR)
dat_slim_MR_df <- dat_slim_MR_df %>% drop_na(wind_spd)

mgustyMR_quap <- quap(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * t_2m,
    a ~ dnorm(0.5, 0.5),
    b1 ~ dnorm(0, 0.02),
    b2 ~ dnorm(0, 0.05),
    sigma ~ dexp(1)
  ) , data=dat_slim_MR_df )

mgustyMR <- ulam(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * t_2m,
    a ~ dnorm(0.5, 0.5),
    b1 ~ dnorm(0, 0.02),
    b2 ~ dnorm(0, 0.05),
    sigma ~ dexp(1)
  ) , data=dat_slim_MR_df, chains =1 )

#HOLY SHIT IT WORKS. Let's do 4 chains!!!

mgustyMR_4chains <- ulam(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * t_2m,
    a ~ dnorm(0.5, 0.5),
    b1 ~ dnorm(0, 0.02),
    b2 ~ dnorm(0, 0.05),
    sigma ~ dexp(1)
  ) , data=dat_slim_MR_df, chains =4, cores=4, iter=1000 )

postmgustyMR4chains <- extract.samples( mgustyMR_4chains ,  clean = FALSE)
str(postmgustyMR4chains)

#this isn't for the model I've made in ulam... it's A model but..
linear_model <- lm(wind_sd ~ wind_spd + t_2m, data = dat_slim_MR_df)
plot(x=predict(linear_model), y=dat_slim_MR_df$wind_sd,
     xlab='Predicted Values',
     ylab='Actual Values',
     xlim = c(0,2),
     ylim = c(0,5),
     main='Predicted vs. Actual Values')
abline(a = 0,
       b = 1,
       col = "red",
       lwd = 2)

cor(x=predict(my_mod), y=dat_slim_MR_df$wind_sd)

##trying to get ulam model plotted
mgustyMR_4chains %>%
  spread_draws(r_condition[])

mgustyMR_draws <- linpred_draws(
  mgustyMR_4chains,
  dat_slim_MR_df,
  value = ".linpred",
  post = NULL,
  ndraws = 3,
  dpar = FALSE)

mgustyMR_predicteddraws <- predicted_draws(
  mgustyMR_4chains,
  dat_slim_MR_df,
  value = ".prediction",
  ndraws = 5)

plot (wind_sd ~ .prediction , mgustyMR_predicteddraws , col=col.alpha(rangi2,0.5),
      xlab='Predicted Values',
      ylab='Actual Values',
      main='Predicted vs. Actual wind_sd values for Bayesian Model')

cor(x=mgustyMR_predicteddraws$.prediction, y=mgustyMR_predicteddraws$wind_sd)

plot (wind_sd.x ~ month , new_combined_cc, col=col.alpha(rangi2,0.5),
      xlab='Month',
      ylab='Wind SD',
      ylim=c(3,5),
      main='Wind Gusts by Month')

plot (wind_sd.x ~ wind_mean.x , new_combined_cc, col=col.alpha(rangi2,0.5),
      xlab='Mean Wind Speed (m/s)',
      ylab='Wind SD',
      ylim=c(3,5),
      main='Wind Gusts by Mean Wind Speed')

plot (wind_sd.x ~ year , new_combined_cc, col=col.alpha(rangi2,0.5),
      xlab='Year',
      ylab='Wind SD',
      main='Wind Gusts by Year')

new_combined_cc_df <- as_tibble(new_combined_cc)

ggplot(new_combined_cc_df, aes(year, wind_sd.x)) +
  geom_smooth()

ggplot(new_combined_cc_df, aes(month, wind_sd.x)) +
  geom_smooth()

#Improving linear model... or attempting to

Richardson_dat_slim_MR <- list(
  wind_sd = Richardson_combined_cc$wind_sd.x,
  t_2m = Richardson_combined_cc$t_2m,
  wind_spd = Richardson_combined_cc$wind_spd,
  BulkRichardsonNumber = Richardson_combined_cc$RichardsonBulk
)
str(Richardson_dat_slim_MR)
Richardson_dat_slim_MR_df <- as_tibble(Richardson_dat_slim_MR)
Richardson_dat_slim_MR_df <- Richardson_dat_slim_MR_df %>% drop_na(wind_spd)
Richardson_dat_slim_MR_df <- Richardson_dat_slim_MR_df[!is.infinite(rowSums(Richardson_dat_slim_MR_df)),]
Richardson_dat_slim_MR_df_cc <- Richardson_dat_slim_MR_df[complete.cases(Richardson_dat_slim_MR_df), ]

linear_model_2 <- lm(wind_sd ~ wind_spd * wind_spd^2 * t_2m * t_2m^2 * BulkRichardsonNumber * BulkRichardsonNumber^2, data = Richardson_dat_slim_MR_df_cc)
plot(x=predict(linear_model_2), y=Richardson_dat_slim_MR_df_cc$wind_sd,
     xlab='Predicted Values',
     ylab='Actual Values',
     xlim = c(0,2),
     ylim = c(0,5),
     main='Predicted vs. Actual Values')
abline(a = 0,
       b = 1,
       col = "red",
       lwd = 2)

cor(x=predict(linear_model_2), y=Richardson_dat_slim_MR_df_cc$wind_sd)

#try improving ulam model to see what adding squared terms does

mgustyMR_test <- ulam(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * t_2m + b3 * BulkRichardsonNumber,
    a ~ dnorm(0.5, 0.5),
    b1 ~ dnorm(0, 0.02),
    b2 ~ dnorm(0, 0.05),
    b3 ~ dnorm(0, 0.04),
    sigma ~ dexp(1)
  ) , data=Richardson_dat_slim_MR_df_cc, chains =1 )

mgustyMR_4chains_test <- ulam(
  alist(
    wind_sd ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * t_2m + b3 * BulkRichardsonNumber,
    a ~ dnorm(0.5, 0.5),
    b1 ~ dnorm(0, 0.02),
    b2 ~ dnorm(0, 0.05),
    b3 ~ dnorm(0, 0.04),
    sigma ~ dexp(1)
  ) , data=Richardson_dat_slim_MR_df_cc, chains = 4, cores = 4, iter = 1000)

mgustyMR4chains_test_predicteddraws <- predicted_draws(
  mgustyMR_4chains_test,
  Richardson_dat_slim_MR_df_cc,
  value = ".prediction",
  ndraws = 5)

plot (wind_sd ~ .prediction , mgustyMR4chains_test_predicteddraws , col=col.alpha(rangi2,0.5),
      xlab='Predicted Values',
      ylab='Actual Values',
      main='Predicted vs. Actual wind_sd values for Bayesian Model')

cor(x=mgustyMR4chains_test_predicteddraws$.prediction, y=mgustyMR4chains_test_predicteddraws$wind_sd)
