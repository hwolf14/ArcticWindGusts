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

minute_summary <- MinuteNoArtifacts %>% 
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
mgusty1plot_updated + coord_cartesian(xlim = c(-30, 10)) +
ggtitle("Correlation between Temperature and Wind_SD") +
  ylab("Wind_SD") + 
  xlab("t_2m (°C)")

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
    wind_sd.x ~ dnorm( mu , sigma ) ,
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

ggplot(new_combined_cc2, aes(wind_spd, wind_sd.x)) +
  geom_smooth() +
ggtitle("Correlation between Wind Speed and Wind_SD") +
  ylab("Wind_SD") + 
  xlab("Wind Speed (m/s)") +
  xlim(0,20)
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
  wind_sd = new_combined_cc2$wind_sd.x,
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

newest_combined_cc <- subset(newest_combined_cc, newest_combined_cc$wind_sd.x<5)

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
new_combined_cc2 <- new_combined_cc2[ new_combined_cc2$wind_sd.x <= 5 & new_combined_cc2$wind_max.x <= 35, ]

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

cor(x=predict(linear_model), y=dat_slim_MR_df$wind_sd)

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
      lines=predict(lm(wind_sd ~ .prediction)),
      main='Predicted vs. Actual wind_sd values for Multiple Regression Model')

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


plot (x=mgustyMR_predicteddraws$.prediction, y=mgustyMR_predicteddraws$wind_sd,
      xlab='Predicted Values',
      ylab='Actual Values',
      xlim = c(0,2),
      ylim = c(0,5),
      main='Predicted vs. Actual values MR Model')
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

#Look more into high-gust time series beginning with t_2m
plot ( t_2m ~ date , new_combined_cc_df , col=col.alpha(rangi2,0.5),
      xlab='Date',
      ylab='Temperature',
      main='Temperature Over Time in Utqiagvik')

ggplot(new_combined_cc_df, aes(date, t_2m)) +
  geom_smooth()

ggplot(data=new_combined_cc_df, aes(x=date, y=t_2m)) +
  geom_smooth(se=FALSE, color="black", aes(group=1)) +
  theme_classic() +
  ylab("Average air temperature (°C)") + 
  xlab("Sampling date")

new_combined_cc_df <- as_tibble(new_combined_cc)
new_combined_cc_df_dates <- mutate(new_combined_cc_df, date= as.Date(date, format= "%d.%m.%Y"))

# Define Start and end times for the subset as R objects that are the time class
startTime <- as.Date("2014-01-01")
endTime <- as.Date("2016-12-31")

# create a start and end time R object
start.end <- c(startTime,endTime)
start.end

# View data for 2014-2016 only
AirTempDaily_2014to2016 <- ggplot(new_combined_cc_df_dates, aes(date, t_2m)) +
  geom_point(na.rm=TRUE, size=1) + 
  ggtitle("Air Temperature\n 2014-2016\n Utqiagvik, Alaska") +
  xlab("Date") + ylab("Air Temperature (C)")+ 
  (scale_x_date(limits=start.end))

AirTempDaily_2014to2016

ggplot(new_combined_cc_df_dates, aes(date, t_2m)) +
  geom_point(na.rm=TRUE) +
  ggtitle("Air Temperature 1994-2022\n Utqiagvik, Alaska") +
  scale_x_date(c(as.Date("2014-01-01"),as.Date("2016-12-31")))

AirTempTimeSeries<- ggplot(new_combined_cc_df_dates, aes(date, t_2m)) +
  geom_point(na.rm=TRUE) +
  ggtitle("Air Temperature 1994-2022\n Utqiagvik, Alaska")

#now plot pressure over time

plot ( p ~ date , new_combined_cc_df , col=col.alpha(rangi2,0.5),
       xlab='Pressure (hPa)',
       ylab='Date',
       main='Pressure Over Time in Utqiagvik')

ggplot(data=new_combined_cc_df, aes(x=date, y=p)) +
  geom_smooth(se=FALSE, color="black", aes(group=1)) +
  theme_classic() +
  ylab("Average Pressure (hPa)") + 
  xlab("Sampling date")

# View data for 2014-2016 only
PressureDaily_2014to2016 <- ggplot(new_combined_cc_df_dates, aes(date, p)) +
  geom_smooth(na.rm=TRUE, linewidth=1) + 
  ggtitle("Pressure\n 2014-2016\n Utqiagvik, Alaska") +
  xlab("Date") + ylab("Pressure (hPa)")+ 
  (scale_x_date(limits=start.end))

PressureDaily_2014to2016

#wind speed 2014-2016

ggplot(new_combined_cc_df_dates, aes(date, wind_spd)) +
  geom_smooth(na.rm=TRUE, linewidth=1) + 
  ggtitle("Wind Speed\n Utqiagvik, Alaska") +
  xlab("Date") + ylab("Wind Speed (m/s)")


WindSpd_2014to2016 <- ggplot(new_combined_cc_df_dates, aes(date, wind_spd)) +
  geom_point(na.rm=TRUE, size=1) + 
  ggtitle("Wind Speed\n 2014-2016\n Utqiagvik, Alaska") +
  xlab("Date") + ylab("Wind Speed (m/s)")+ 
  (scale_x_date(limits=start.end))

WindSpd_2014to2016

#compare to just wind gusts
ggplot(data=new_combined_cc_df, aes(x=date, y=wind_sd.x)) +
  geom_point(color="black", aes(group=1)) +
  theme_classic() +
  ylab("wind_sd") + 
  xlab("Sampling date")

WindGusts_2014to2016 <- ggplot(new_combined_cc_df_dates, aes(date, wind_sd.x)) +
  geom_point(na.rm=TRUE, size=1) + 
  ggtitle("Wind Gusts\n 2014-2016\n Utqiagvik, Alaska") +
  xlab("Date") + ylab("wind_sd")+ 
  (scale_x_date(limits=start.end))

WindGusts_2014to2016

#Look around the hours where high SD values are, +/- 6 hours or so

library

Minute_dates <- MinuteNoArtifacts %>% mutate(date = make_datetime(year, month, day, hour, minute))

Minute_dates$date=as.POSIXct(date[,1], format="%y/%m/%d %H:%M:%S")


ggplot(data=Minute_dates,aes(x=date, y=t_2m)) + 
  geom_point(color="black", aes(group=1)) +
  ggtitle("Wind Gusts\n 2004-12-01\n Utqiagvik, Alaska") +
  ylab("Temperature (ºC)") + 
  xlab("Time") + 
  ylim (-25,-15) +
  scale_x_datetime(limits = ymd_h(c("2004-12-01 12", "2004-12-02 0")))

ggplot(data=Minute_dates,aes(x=date, y=wind_spd)) + 
  geom_point(color="black", aes(group=1)) +
  ggtitle("Wind Gusts\n 2004-12-01\n Utqiagvik, Alaska") +
  ylab("Wind Speed (m/s)") + 
  xlab("Time") + 
  ylim(0,20) +
  scale_x_datetime(limits = ymd_h(c("2004-12-01 12", "2004-12-02 0")))

ggplot(data=Minute_dates,aes(x=date, y=p)) + 
  geom_point(color="black", aes(group=1)) +
  ggtitle("Wind Gusts\n 2005-02-08\n Utqiagvik, Alaska") +
  ylab("Pressure (hPa)") + 
  xlab("Time") + 
  ylim(990,1010) +
  scale_x_datetime(limits = ymd_h(c("2005-02-08 16", "2005-02-09 4")))

ggplot(data=Minute_dates,aes(x=date, y=wind_dir)) + 
  geom_point(color="black", aes(group=1)) +
  ggtitle("Wind Gusts\n 2015-12-25\n Utqiagvik, Alaska") +
  ylab("Wind Direction (degrees)") + 
  xlab("Time") + 
  ylim(0,360) +
  scale_x_datetime(limits = ymd_h(c("2015-12-25 14", "2015-12-25 20")))

#look at high gusts only

HighGusts_combined_cc_df <- as_tibble(HighGusts_combined_cc) 

ggplot(data=HighGusts_combined_cc_df,aes(x=date, y=wind_sd.x)) +
  geom_point(color="black", aes(group=1)) +
  ggtitle("High Wind Gusts\n Utqiagvik, Alaska") +
  ylab("Wind_SD") + 
  xlab("Time")

HighGusts_combined_cc <- new_combined_cc %>% filter(3<wind_sd.x & 5>wind_sd.x)

mean(HighGusts_combined_cc_df$t_2m, na.rm=TRUE)

#Make month datasets
Month1_combined_cc <- new_combined_cc %>% filter(month==1)
Month2_combined_cc <- new_combined_cc %>% filter(month==2)
Month3_combined_cc <- new_combined_cc %>% filter(month==3)
Month4_combined_cc <- new_combined_cc %>% filter(month==4)
Month5_combined_cc <- new_combined_cc %>% filter(month==5)
Month6_combined_cc <- new_combined_cc %>% filter(month==6)
Month7_combined_cc <- new_combined_cc %>% filter(month==7)
Month8_combined_cc <- new_combined_cc %>% filter(month==8)
Month9_combined_cc <- new_combined_cc %>% filter(month==9)
Month10_combined_cc <- new_combined_cc %>% filter(month==10)
Month11_combined_cc <- new_combined_cc %>% filter(month==11)
Month12_combined_cc <- new_combined_cc %>% filter(month==12)

#view wind gust data over time for a particular month
ggplot(data=Month12_combined_cc,aes(x=date, y=wind_sd.x)) +
  geom_point(color="black", aes(group=1)) +
  ggtitle("High Wind Gusts\n December \n Utqiagvik, Alaska") +
  ylab("Wind_SD") + 
  xlab("Year")

ggplot(data=Month12_combined_cc,aes(x=date, y=wind_sd.x)) +
  geom_smooth(color="black", aes(group=1)) +
  ggtitle("High Wind Gusts\n December \n Utqiagvik, Alaska") +
  ylab("Wind_SD") + 
  xlab("Year")

#view only cases where sd > 2
Month1_high <- Month1_combined_cc %>% filter(wind_sd.x>2)
Month2_high <- Month2_combined_cc %>% filter(wind_sd.x>2)
Month3_high <- Month3_combined_cc %>% filter(wind_sd.x>2)
Month4_high <- Month4_combined_cc %>% filter(wind_sd.x>2)
Month5_high <- Month5_combined_cc %>% filter(wind_sd.x>2)
Month6_high <- Month6_combined_cc %>% filter(wind_sd.x>2)
Month7_high <- Month7_combined_cc %>% filter(wind_sd.x>2)
Month8_high <- Month8_combined_cc %>% filter(wind_sd.x>2)
Month9_high <- Month9_combined_cc %>% filter(wind_sd.x>2)
Month10_high <- Month10_combined_cc %>% filter(wind_sd.x>2)
Month11_high <- Month11_combined_cc %>% filter(wind_sd.x>2)
Month12_high <- Month12_combined_cc %>% filter(wind_sd.x>2)
Month1_count <- Month1_high %>% count(year)
Month2_count <- Month2_high %>% count(year)
Month3_count <- Month3_high %>% count(year)
Month4_count <- Month4_high %>% count(year)
Month5_count <- Month5_high %>% count(year)
Month6_count <- Month6_high %>% count(year)
Month7_count <- Month7_high %>% count(year)
Month8_count <- Month8_high %>% count(year)
Month9_count <- Month9_high %>% count(year)
Month10_count <- Month10_high %>% count(year)
Month11_count <- Month11_high %>% count(year)
Month12_count <- Month12_high %>% count(year)

#add missing years with 0 values
Month1_count <- Month1_count %>% add_row (year=1994:1995, n=0)
Month1_count <- Month1_count %>% add_row (year=1998, n=0)
Month1_count <- Month1_count %>% add_row (year=2004, n=0)
Month1_count <- Month1_count %>% add_row (year=2006:2007, n=0)
Month1_count <- Month1_count %>% add_row (year=2012:2015, n=0)
Month1_count <- Month1_count %>% add_row (year=2018, n=0)
Month1_count <- Month1_count %>% add_row (year=2020:2022, n=0)
Month2_count <- Month2_count %>% add_row (year=1994, n=0)
Month2_count <- Month2_count %>% add_row (year=1995, n=0)
Month2_count <- Month2_count %>% add_row (year=1997:1998, n=0)
Month2_count <- Month2_count %>% add_row (year=2004, n=0)
Month2_count <- Month2_count %>% add_row (year=2010, n=0)
Month2_count <- Month2_count %>% add_row (year=2012, n=0)
Month2_count <- Month2_count %>% add_row (year=2017, n=0)
Month2_count <- Month2_count %>% add_row (year=2020, n=0)
Month3_count <- Month3_count %>% add_row (year=1994, n=0)
Month3_count <- Month3_count %>% add_row (year=1996, n=0)
Month3_count <- Month3_count %>% add_row (year=1997, n=0)
Month3_count <- Month3_count %>% add_row (year=1999:2003, n=0)
Month3_count <- Month3_count %>% add_row (year=2006:2007, n=0)
Month3_count <- Month3_count %>% add_row (year=2009:2012, n=0)
Month3_count <- Month3_count %>% add_row (year=2015, n=0)
Month3_count <- Month3_count %>% add_row (year=2017:2019, n=0)
Month3_count <- Month3_count %>% add_row (year=2022, n=0)
Month4_count <- Month4_count %>% add_row (year=1994:2001, n=0)
Month4_count <- Month4_count %>% add_row (year=2003:2007, n=0)
Month4_count <- Month4_count %>% add_row (year=2009:2013, n=0)
Month4_count <- Month4_count %>% add_row (year=2017:2020, n=0)
Month5_count <- Month5_count %>% add_row (year=1994:2001, n=0)
Month5_count <- Month5_count %>% add_row (year=2005:2006, n=0)
Month5_count <- Month5_count %>% add_row (year=2009, n=0)
Month5_count <- Month5_count %>% add_row (year=2011, n=0)
Month5_count <- Month5_count %>% add_row (year=2016:2018, n=0)
Month5_count <- Month5_count %>% add_row (year=2020:2022, n=0)
Month6_count <- Month6_count %>% add_row (year=1994:1999, n=0)
Month6_count <- Month6_count %>% add_row (year=2001:2003, n=0)
Month6_count <- Month6_count %>% add_row (year=2005, n=0)
Month6_count <- Month6_count %>% add_row (year=2007, n=0)
Month6_count <- Month6_count %>% add_row (year=2009:2013, n=0)
Month6_count <- Month6_count %>% add_row (year=2017:2020, n=0)
Month7_count <- Month7_count %>% add_row (year=1994:2002, n=0)
Month7_count <- Month7_count %>% add_row (year=2004:2005, n=0)
Month7_count <- Month7_count %>% add_row (year=2007:2008, n=0)
Month7_count <- Month7_count %>% add_row (year=2010:2012, n=0)
Month7_count <- Month7_count %>% add_row (year=2015:2018, n=0)
Month8_count <- Month8_count %>% add_row (year=1996, n=0)
Month8_count <- Month8_count %>% add_row (year=1998, n=0)
Month8_count <- Month8_count %>% add_row (year=2001:2003, n=0)
Month8_count <- Month8_count %>% add_row (year=2007:2008, n=0)
Month8_count <- Month8_count %>% add_row (year=2011:2013, n=0)
Month8_count <- Month8_count %>% add_row (year=2015:2016, n=0)
Month8_count <- Month8_count %>% add_row (year=2020, n=0)
Month9_count <- Month9_count %>% add_row (year=1994:1998, n=0)
Month9_count <- Month9_count %>% add_row (year=2001:2002, n=0)
Month9_count <- Month9_count %>% add_row (year=2005:2008, n=0)
Month9_count <- Month9_count %>% add_row (year=2012, n=0)
Month9_count <- Month9_count %>% add_row (year=2014, n=0)
Month9_count <- Month9_count %>% add_row (year=2016, n=0)
Month9_count <- Month9_count %>% add_row (year=2018:2019, n=0)
Month10_count <- Month10_count %>% add_row (year=1994, n=0)
Month10_count <- Month10_count %>% add_row (year=2003:2005, n=0)
Month10_count <- Month10_count %>% add_row (year=2021, n=0)
Month11_count <- Month11_count %>% add_row (year=1994:1995, n=0)
Month11_count <- Month11_count %>% add_row (year=2006, n=0)
Month11_count <- Month11_count %>% add_row (year=2008, n=0)
Month11_count <- Month11_count %>% add_row (year=2016, n=0)
Month11_count <- Month11_count %>% add_row (year=2018, n=0)
Month12_count <- Month12_count %>% add_row (year=1995, n=0)
Month12_count <- Month12_count %>% add_row (year=1997, n=0)
Month12_count <- Month12_count %>% add_row (year=2000, n=0)
Month12_count <- Month12_count %>% add_row (year=2006, n=0)
Month12_count <- Month12_count %>% add_row (year=2009:2012, n=0)

ggplot(data=Month12_count,aes(x=year, y=n)) +
  geom_point(color="black", aes(group=1)) +
  geom_smooth(method=lm) +
  ggtitle("Occurrences of High Wind Gusts\n December") +
  ylab("# of measurements where wind_sd>2") + 
  xlab("Year") +
  xlim (1994, 2022)

#changes in temp over the years by month
ggplot(data=Month12_combined_cc,aes(x=date, y=t_2m)) +
  geom_point(color="black", aes(group=1)) +
  geom_smooth(method=lm) +
  ggtitle("December Temperatures\n Utqiagvik, Alaska") +
  ylab("Temperature (C)") + 
  xlab("Year")

#changes in wind spd over the years by month
ggplot(data=Month12_combined_cc,aes(x=date, y=wind_spd)) +
  geom_point(color="black", aes(group=1)) +
  geom_smooth(method=lm) +
  ggtitle("December Wind Speeds\n Utqiagvik, Alaska") +
  ylab("Wind Speed (m/s)") + 
  xlab("Year")

#is wind speed in general higher in winter?
ggplot(data=new_combined_cc_df,aes(x=month, y=wind_spd)) +
  geom_smooth(color="black", aes(group=1)) +
  ggtitle("Wind Speed by Month\n Utqiagvik, Alaska") +
  ylab("Wind Speed (m/s)") + 
  xlab("Month")

#pressure and wind_sd
plot (wind_sd.x ~ p , new_combined_cc, col=col.alpha(rangi2,0.5),
      xlab='Pressure (hPa)',
      ylab='Wind SD',
      main='Wind Gusts and Pressure')
cor(x=new_combined_cc$wind_sd.x, y=new_combined_cc$p, use = "complete.obs")

#minute plots
ggplot(data=Minute_dates,aes(x=date, y=t_2m)) + 
  geom_point(color="black", aes(group=1)) +
  ggtitle("Wind Gusts\n 2019-08-02\n Utqiagvik, Alaska") +
  ylab("Pressure (hPa)") + 
  xlab("Time") + 
  scale_x_datetime(limits = ymd_h(c("2019-8-11 0", "2019-8-11 23")))

#filter out clear artefacts

Minute <- Minute[Minute$wind_spd >= 0, ]

dates <- as.POSIXct(c())
MinuteRows <- Minute |> mutate(row = seq(n()))

MinuteNoArtifacts <- MinuteRows[-c(5136745:5137573, 5137598:5137689, 5137728:5137798, 5138067:5138184, 5649385, 5649500, 5649501, 5650506, 5650508:5650535, 10512561, 10512562, 10512567, 10512568, 10512806, 10512811, 10512812, 10512826, 10512851, 10512853, 10512870, 10512871, 10512921, 10512930, 10512933, 10512945, 10512947:10512950, 10548981:10548982, 10548985, 10549001:10549010, 10549014, 10549040, 10549060:10549062, 10549087:10549146, 10549294, 10589989:10589990, 10589995:10590013, 10590026, 10590035, 10590355:10590356, 10590359, 10590361, 10590364, 10590366:10590368, 10590380, 11591548, 11591549, 11591552:11591553, 11591555:11591556, 11591559, 11591562, 11591564, 6025143, 6025145:6026406, 8543808, 8543810:8543812, 8543814:8543818), ]

#get slope of regression line for each monthly plot
summary(lm(formula = Month1_combined_cc$wind_spd ~ Month1_combined_cc$year))
summary(lm(formula = Month1_combined_cc$wind_spd ~ Month1_combined_cc$year))
#Plot m values for each month for gust frequency, temp, wind spd

month <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
gust_slope <- c(0.05468, 0.00936, 0.2429, 0.04286, 0.01429, 0.02167, 0.09332, 0.3035, 0.1924, 0.1253, 0.25972, 0.14067)
temp_slope <- c(0.14972, 0.07747, 0.1603, 0.05521, 0.09262, 0.023258, 0.05198, 0.01883, 0.06779, 0.2429, 0.1890, 0.1285)
wind_spd_slope <- c(0.023519, 0.007336, 0.016447, 0.021963, 0.001949, 0.028198, -0.005238, -0.0002307, -0.004906, 0.011290, 0.028367, 0.041099)
mean_st_dev_slope <- c(0.007056, 0.005667, 0.007224, 0.0042022, 0.005018, 0.006747, 0.006214, 0.008362, 0.005922, 0.007642, 0.009847, 0.009264)
mean_st_dev_slope_error <- c(0.0002719, 0.0002755, 0.000261, 0.0002068, 0.000198, 0.0001863, 0.0002128, 0.000248, 0.0002607, 0.0003094, 0.0003064, 0.0002893)
wind_spd_slope_error <- c(0.003036, 0.007336, 0.00294, 0.002473, 0.002174, 0.001928, 0.002148, 0.0023891, 0.002493, 0.002944, 0.003128, 0.003136)
temp_slope_error <-c(0.00615, 0.007286, 0.006139, 0.00576, 0.003992, 0.002591, 0.00307, 0.00299, 0.002683, 0.002429, 0.00189, 0.005753)
slopes_df <- data.frame(month, gust_slope, temp_slope, wind_spd_slope, mean_st_dev_slope, mean_st_dev_slope_error, wind_spd_slope_error, temp_slope_error)

ggplot(data=slopes_df,aes(x=month, y=wind_spd_slope)) + 
  geom_point(color="black", aes(group=1)) +
  ggtitle("Change in Wind Speed by Month From 1994-2022") +
  ylab("slope (wind speed (m/s) / year)") + 
  xlab("Month") +
  xlim(0,12)

#Figure 5 for my thesis is below
ggplot(data=slopes_df,aes(x=month, y=gust_slope)) + 
  geom_point(color="black", aes(group=1)) +
  ggtitle("Change in Gust Frequency by Month From 1994-2022") +
  ylab("slope (change in # of occurrences of sd>2/year)") + 
  xlab("Month") +
  xlim(0,12)

#take out August 1 and see how trend changes, 08/14/2023
year <- c(1994:2021)
n <- c(1,1,0,2,0,1,1,0,0,0,1,1,2,0,0,2,1,0,0,0,3,0,0,1,5,8,0,11)
Month8_count_no0801 <- data.frame(year,n)

ggplot(data=Month8_count_no0801,aes(x=year, y=n)) +
  geom_point(color="black", aes(group=1)) +
  geom_smooth(method=lm) +
  ggtitle("Occurrences of High Wind Gusts\n August Without 0801") +
  ylab("# of measurements where wind_sd>2") + 
  xlab("Year") +
  xlim (1994, 2022)

#Plot mean wind speed over time by year and months
year_windspd_means <- ddply(new_combined_cc_df, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE))
ggplot(data=year_windspd_means,aes(x=year, y=mean_wind_spd)) +
  geom_point(color="black", aes(group=1)) +
  geom_smooth(method=lm) +
  ggtitle("Annual Mean Wind Speed Over Time") +
  ylab("Mean Wind Speed (m/s)") + 
  xlab("Year") +
  xlim (1994, 2022)

Month1_combined_cc_summary <- ddply(Month1_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month2_combined_cc_summary <- ddply(Month2_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month3_combined_cc_summary <- ddply(Month3_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month4_combined_cc_summary <- ddply(Month4_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month5_combined_cc_summary <- ddply(Month5_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month6_combined_cc_summary <- ddply(Month6_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month7_combined_cc_summary <- ddply(Month7_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month8_combined_cc_summary <- ddply(Month8_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month9_combined_cc_summary <- ddply(Month9_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month10_combined_cc_summary <- ddply(Month10_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month11_combined_cc_summary <- ddply(Month11_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))
Month12_combined_cc_summary <- ddply(Month12_combined_cc, .(year), summarize, mean_wind_spd=mean(wind_spd, na.rm=TRUE), mean_p=mean(p, na.rm=TRUE), mean_wind_sd=mean(wind_sd.x, na.rm=TRUE), mean_t_2m=(mean(t_2m, na.rm=TRUE)))

ggplot(data=Month12_combined_cc_summary,aes(x=year, y=mean_wind_sd)) +
  geom_point(color="black", aes(group=1)) +
  geom_smooth(method=lm) +
  ggtitle("December Mean Standard Deviation Over Time") +
  ylab("Mean Standard Deviation") + 
  xlab("Year") +
  xlim (1994, 2022)

#make error bar plot 08/18/2023
summary(lm(formula = Month12_combined_cc$t_2m ~ Month12_combined_cc$year))
#added these values to slopes_df above, now make error bar plots
ggplot(data=slopes_df, aes(x=month,y=temp_slope)) +
  geom_point() +
  ggtitle("Changes in Mean Temperature, 1994-2022") +
  geom_errorbar(aes(ymax = temp_slope + temp_slope_error, ymin = temp_slope - temp_slope_error),
                position = "dodge")

ggplot(data=slopes_df, aes(x=month,y=wind_spd_slope)) +
  geom_point() +
  ggtitle("Changes in Mean Wind Speed, 1994-2022") +
  geom_errorbar(aes(ymax = wind_spd_slope + wind_spd_slope_error, ymin = wind_spd_slope - wind_spd_slope_error),
                position = "dodge")

#regression test 08/23/2023
month12highfit_sd <- lm(formula = Month12_high$wind_sd.x ~ Month12_high$year)
summary(month12highfit_sd)

#regression test 09/04/2023 for wind sd and 09/15/2023 for wind spd and t_2m
Month12meanfit_t_2m <- lm(formula = Month12_combined_cc_summary$mean_t_2m ~ Month12_combined_cc_summary$year)
summary(Month12meanfit_t_2m)

#make monthly mean temperature plots over time
ggplot(data=Month12_combined_cc_summary,aes(x=year, y=mean_t_2m)) +
  geom_point(color="black", aes(group=1)) +
  geom_smooth(method=lm) +
  ggtitle("December Mean Temperature Over Time") +
  ylab("Mean Temperature (Celsius)") + 
  xlab("Year") +
  xlim (1994, 2022)

#Temperature exceedance frequency: how many hours do we have thaw conditions?
Month12_thawcount <- Month12_combined_cc %>% filter(t_2m>-1.8)
Month12_thawcount <- Month12_thawcount %>% group_by(year) %>% tally()

Month12_thawcount <- Month12_thawcount %>% add_row (year=1994:1997, n=0)

ggplot(data=Month12_thawcount,aes(x=year, y=n)) +
  geom_point(color="black", aes(group=1)) +
  geom_smooth(method=lm) +
  ggtitle("December hours with chaw conditions") +
  ylab("# of hours above thaw conditions") + 
  xlab("Year") +
  xlim (1994, 2022)

Month12meanfit_thawcount <- lm(formula = Month12_thawcount$n ~ Month12_thawcount$year)
summary(Month12meanfit_thawcount)

#11/28/2023: remake Mean Wind SD plot
ggplot(data=slopes_df,aes(x=month, y=mean_st_dev_slope)) + 
  geom_point(color="black", aes(group=1)) +
  ggtitle("Changes in Mean Wind SD, 1994-2022") +
  ylab("Rate of Change of Mean Wind_SD)") + 
  xlab("Month") +
  xlim(0,12) +
  geom_errorbar(aes(ymax = mean_st_dev_slope + mean_st_dev_slope_error, ymin = mean_st_dev_slope - mean_st_dev_slope_error),
                position = "dodge")
