#update gust_time parameter so it includes only points 50% above mean wind spd

library(tidyverse)
library(lubridate)

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

#create new combined dataset to include updated gust_time

updated_combined <- Hourly %>%
  left_join(updated_clean_minute_summary, by = c("site", "year",  "month", "day", "hour"))

updated_combined_cc <- updated_combined %>% drop_na(gust_time, t_2m)

updated_new_combined_cc <- updated_combined_cc %>% 
  mutate(date = make_datetime(year, month, day, hour)) %>% 
  filter(date %in% full_dates)

updated_new_combined_cc <- updated_new_combined_cc %>% 
  left_join(select(clean_minute_summary, site, date, new_gust_time = gust_time, 
                   wind_mean, wind_sd, wind_max, wind_min, wind_med), 
            by = c("site", "date"))

rm(minute_data)
gc()

#Build mgusty1

mgusty1_updated <- quap(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , data=updated_combined_cc)

t_2m.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_2m = t_2m.seq )
mu <- link( mgusty1_updated , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.gust_time <- sim( mgusty1_updated , data=updated_combined_cc )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

#plot with mgusty1

plot( gust_time ~ t_2m , updated_new_combined_cc , col=col.alpha(rangi2,0.5) )
lines( t_2m.seq , mu.mean )
shade( mu.PI , t_2m.seq )

mgusty1plot_updated <- ggplot(updated_new_combined_cc, aes(t_2m, gust_time)) +
  geom_smooth()
mgusty1plot_updated + coord_cartesian(xlim = c(-30, 10))

gam_updated <- gam(gust_time ~ s(t_2m + t_2m^2), data = updated_new_combined_cc)

plot(gam, add = TRUE, ylim = c(0,10), ylab = "Gust time", xlab = "temperature")
points(gust_time ~ t_2m, data = slice_sample(updated_new_combined_cc, n = 2000), col = alpha(rangi2, 0.5))
lines(t_2m.seq, mu.mean)
shade(mu.PI, t_2m.seq)

precis(mgusty1_updated)

#Build mgusty2

updated_combined_cc2 <- updated_combined %>% drop_na(gust_time, wind_spd)

updated_new_combined_cc2 <- updated_combined_cc2 %>% 
  mutate(date = make_datetime(year, month, day, hour)) %>% 
  filter(date %in% full_dates)

updated_new_combined_cc2 <- updated_new_combined_cc2 %>% 
  left_join(select(updated_clean_minute_summary, site, date, new_gust_time = gust_time, 
                   wind_mean, wind_sd, wind_max, wind_min, wind_med), 
            by = c("site", "date"))

mgusty2_updated <- quap(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * wind_spd^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , data=updated_combined_cc2)

#plotting with mgusty2_updated

precis(mgusty2_updated)

wind_spd.seq <- seq( from=0 , to=100 , length.out=30 )
pred_dat <- list( wind_spd_s=wind_spd.seq , wind_spd_s2=wind_spd.seq^2 )
mu <- link( mgusty2_updated , data=updated_combined_cc2 )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory
sim.gust_time <- sim( mgusty2_updated , data=updated_combined_cc2 )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

plot( gust_time ~ wind_spd , updated_new_combined_cc2 , col=col.alpha(rangi2,0.5) )
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

ggplot(updated_new_combined_cc2, aes(wind_spd, gust_time)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

mgusty2plot_updated <- ggplot(updated_new_combined_cc2, aes(wind_spd, gust_time)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

mgusty2plot_updated + coord_cartesian(xlim = c(0, 20))

#Build mgusty3

dat_slim_updated <- list(
  gust_time = updated_new_combined_cc$gust_time,
  t_2m = updated_new_combined_cc$t_2m
)
str(dat_slim_updated)

mgusty3_updated <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) ,
  data=dat_slim_updated , chains=1 )

mgusty3.1_updated <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
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

sim.gust_time <- sim( mgusty3.1_updated , data=dat_slim_updated )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

dat_slim_df_updated <- as_tibble(dat_slim_updated)

plot( gust_time ~ t_2m , dat_slim_df_updated , col=col.alpha(rangi2,0.5) )
lines( t_2m.seq , mu.mean )
shade( mu.PI , t_2m.seq )

ggplot(dat_slim_df_updated, aes(t_2m, gust_time)) +  geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))
ggplot(dat_slim_df_updated, aes(t_2m, gust_time)) + geom_point(color = alpha(rangi2, 0.1)) + geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))

plot3.1_updated <- ggplot(dat_slim_df_updated, aes(t_2m, gust_time)) + geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))
plot3.1_updated + coord_cartesian(ylim = c(0, 10))

#Build new mgusty4

dat_slim_2_updated <- list(
  gust_time = updated_new_combined_cc2$gust_time,
  wind_spd = updated_new_combined_cc2$wind_spd
)
str(dat_slim_2_updated)

mgusty4_updated <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * wind_spd^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , 
  data=dat_slim_2_updated , chains=1 )

mgusty4.1_updated <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * wind_spd^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , 
  data=dat_slim_2_updated , chains=4 , cores=4 , iter = 1000)

#plotting with both ulam models

library(tidybayes)
library(tidybayes.rethinking)
library(bayesplot)
mcmc_dens(mgusty4.1_updated@stanfit)

post_mgusty3_updated <- extract.samples(mgusty3)
post_mgusty3.1_updated <- extract.samples(mgusty3.1)

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
sim.gust_time <- sim( mgusty4.1_updated , data=dat_slim_2_updated )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

dat_slim_2_df_updated <- as_tibble(dat_slim_2_updated)

ggplot(dat_slim_2_df_updated, aes(wind_spd, gust_time)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

plot4.1_updated <- ggplot(dat_slim_2_df_updated, aes(wind_spd, gust_time)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

plot4.1_updated + coord_cartesian(xlim = c(0, 20))

#Build new mgusty5, temperature difference model

updated_newest_combined_cc <- updated_new_combined_cc
print(updated_newest_combined_cc)

updateD_newest_combined_cc <- mutate(updated_newest_combined_cc, t_difference = t_10m-t_2m)

mgusty5_updated <- quap(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , data=updated_newest_combined_cc)

#plot with mgusty5 to graph t_difference

t_difference.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_difference = t_difference.seq )
mu <- link( mgusty5_updated , data=updated_newest_combined_cc )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.gust_time <- sim( mgusty5_updated , data=updated_newest_combined_cc )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

ggplot(updated_newest_combined_cc, aes(t_difference, gust_time)) +
  geom_smooth()

ggplot(data = updated_newest_combined_cc, aes (x = t_difference)) + geom_density() + xlim(-2, 2)

ggplot(updated_newest_combined_cc, aes(x=t_difference, y=gust_time)) + 
  geom_bin2d(bins = 100, mapping = aes(fill = log(..ndensity..))) +
  theme_bw()

ggplot(updated_newest_combined_cc, aes(x=t_difference, y=gust_time)) + 
  geom_bin2d(bins = c(30,30), mapping = aes(fill = log(..ndensity..))) +
  scale_fill_viridis_c(option = "B") +
  xlim(-2,2) +
  theme_bw()

#Build new mgusty6

mgusty6_updated <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) ,
  data=dat_slim_updated , chains=4 , cores=4 , iter = 1000 )

# plot with mgusty6_updated

t_difference.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_difference = t_difference.seq )
mu <- link( mgusty6_updated , data=updated_newest_combined_cc )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.gust_time <- sim( mgusty6_updated , data=updated_newest_combined_cc )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

updated_newest_combined_cc_df <- as_tibble(updated_newest_combined_cc)

plot( gust_time ~ t_difference , updated_newest_combined_cc_df , col=col.alpha(rangi2,0.1))
lines( t_difference.seq , mu.mean )
shade( mu.PI , t_difference.seq )

ggplot(updated_newest_combined_cc_df, aes(t_difference, gust_time)) +
  geom_smooth() + xlim(-2, 2) + ylim (0, 61)

# code updated mgusty models 3-5 into brms

mgusty3.1prior_updated <- c(
  set_prior("normal(5, 5)", class = "Intercept", lb = 0),
  set_prior("normal(0, 1)", class = "b", coef = "t_2m"),
  set_prior("exponential(0.05)", class="sigma"))
mgusty3.1_fit_formula_updated <- bf(gust_time ~ t_2m + I(t_2m^2))
b3.1_updated <-
  brm(mgusty3.1_fit_formula_updated,
      data = updated_new_combined_cc,
      family = gaussian,
      prior = mgusty3.1prior_updated,
      chains = 4, cores = 4, backend = "cmdstanr" )

mgusty4.1prior_updated <- c(
  set_prior("normal(5, 5)", class = "Intercept", lb = 0),
  set_prior("normal(0, 1)", class = "b", coef = "wind_spd"),
  set_prior("exponential(0.05)", class="sigma"))
mgusty4.1_fit_formula_updated <- bf(gust_time ~ wind_spd + I(wind_spd^2))
b4.1_updated <-
  brm(mgusty4.1_fit_formula_updated,
      data = updated_new_combined_cc,
      family = gaussian,
      prior = mgusty4.1prior,
      chains = 4, cores = 4, backend = "cmdstanr" )

mgusty5.1prior_updated <- c(
  set_prior("normal(5, 5)", class = "Intercept", lb = 0),
  set_prior("normal(0, 1)", class = "b", coef = "t_difference"),
  set_prior("exponential(0.05)", class="sigma"))
mgusty5.1_fit_formula_updated <- bf(gust_time ~ t_difference + I(t_difference^2))
b5.1_updated <-
  brm(mgusty5.1_fit_formula_updated,
      data = updated_newest_combined_cc_df,
      family = gaussian,
      prior = mgusty5.1prior_updated,
      chains = 4, cores = 4, backend = "cmdstanr" )

# Calculate bulk Richardson number from existing parameters and add to dataset
#First find potential temperature using the hydrostatic equation

updated_Richardson_combined_cc <- mutate(updated_newest_combined_cc, RichardsonBulk = ((9.81/t_2m)*(t_difference)*(8))/(wind_spd^2))

t_average_cc <- mutate(updated_Richardson_combined_cc, avg_temp_K = ((t_2m+t_10m)/2)+273.15)
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
updated_Richardson_combined_cc_new <- mutate(deltatheta_cc, RichardsonBulkPotential = ((9.81/t_2m)*(DeltaTheta)*(8))/(wind_spd^2))

# Plot bulk Richardson number
ggplot(updated_Richardson_combined_cc_new, aes(x=RichardsonBulkPotential, y=gust_time)) +
  geom_bin2d(bins = c(30,30), mapping = aes(fill = log(..ndensity..))) +
  scale_fill_viridis_c(option = "B") +
  xlim(-1,1) +
  theme_bw()

ggplot(updated_Richardson_combined_cc_new, aes(x=RichardsonBulkPotential, y=gust_time)) +
  xlim(-1,1) +
  geom_point(alpha = 0.1, na.rm = TRUE)