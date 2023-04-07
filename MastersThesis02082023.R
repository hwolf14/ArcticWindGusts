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

library(rethinking)
library(mgcv)

combined <- read_rds("data/combined_data.Rds")

combined <- combined %>%
  drop_na(gust_time)

combined_cc <- combined %>% drop_na(gust_time, t_2m)

mgusty1 <- quap(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , data=combined_cc)

t_2m.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_2m = t_2m.seq )
mu <- link( mgusty1 , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory



sim.gust_time <- sim( mgusty1 , data=combined_cc )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

plot( gust_time ~ t_2m , new_combined_cc , col=col.alpha(rangi2,0.5) )
lines( t_2m.seq , mu.mean )
shade( mu.PI , t_2m.seq )

ggplot(new_combined_cc, aes(t_2m, gust_time)) +
  geom_smooth()

mgusty1plot <- ggplot(new_combined_cc, aes(t_2m, gust_time)) +
  geom_smooth()
mgusty1plot + coord_cartesian(ylim = c(0, 10))

gam <- gam(gust_time ~ s(t_2m + t_2m^2), data = new_combined_cc)

plot(gam, add = TRUE, ylim = c(0,60), ylab = "Gust time", xlab = "temperature")
points(gust_time ~ t_2m, data = slice_sample(new_combined_cc, n = 2000), col = alpha(rangi2, 0.5))
lines(t_2m.seq, mu.mean)
shade(mu.PI, t_2m.seq)

precis(mgusty1)

combined_cc_1 <- combined_cc %>% mutate(has_gusts = gust_time > 0)


gusty_zif_model <- 
  "data{
  int <lower=0> N;
  int gust_time[N];
  vector[N] t_2m;
}
transformed data {
  vector[N] t_2m_2;
  t_2m_2 = t_2m^2;
}
parameters{
  real p_zero;
  real a;
  real b1;
  real b2;
  real<lower=0> sigma;
}
model{
  vector[N] mu;
  real mu_temp;
  
  sigma ~ exponential( 0.05 );
  pz ~ beta(3,3);
  b2 ~ normal( 0 , 1 );
  b1 ~ normal( 0 , 1 );
  a ~ normal( 5 , 5 );
    mu_temp = inv_logit(a + b1 * t_2m + b2 * t_2m_2);
  for ( i in 1:N ) {
    if (gust_time == 0) {
      target += log1m(pz) + binomial_lpmf(gust_time[i]|60,mu_temp[i]);
    } else {
      target += log_mix(pz,0,binomial_lpmf(0|1,mu_temp[i]));
    }
}
"

stan_data <- list(N = nrow(combined_cc),
                  t_2m = combined_cc$t_2m,
                  gust_time = as.integer(combined_cc$gust_time))

stan_file <- write_stan_file(gusty_zif_model)
new_model <- cmdstan_model(stan_file)
new_fit <- new_model$sample(data = stan_data, chains = 4)

new_model_draws <- tidy_draws(new_fit)
new_model_pred <- linpred_draws(new_fit, list(N = nrows(pred_dat)))

combined_cc2 <- combined %>% drop_na(gust_time, wind_spd)

new_combined_cc2 <- combined_cc2 %>% 
  mutate(date = make_datetime(year, month, day, hour)) %>% 
  filter(date %in% full_dates)

new_combined_cc2 <- new_combined_cc2 %>% 
  left_join(select(clean_minute_summary, site, date, new_gust_time = gust_time, 
                   wind_mean, wind_sd, wind_max, wind_min, wind_med), 
            by = c("site", "date"))

mgusty2 <- quap(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * wind_spd^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , data=combined_cc2)

precis(mgusty2)

wind_spd.seq <- seq( from=0 , to=100 , length.out=30 )
pred_dat <- list( wind_spd_s=wind_spd.seq , wind_spd_s2=wind_spd.seq^2 )
mu <- link( mgusty2 , data=combined_cc2 )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory
sim.gust_time <- sim( mgusty2 , data=combined_cc2 )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

plot( gust_time ~ wind_spd , new_combined_cc2 , col=col.alpha(rangi2,0.5) )
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

ggplot(new_combined_cc2, aes(wind_spd, gust_time)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

mgusty2plot <- ggplot(new_combined_cc2, aes(wind_spd, gust_time)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

mgusty2plot + coord_cartesian(xlim = c(0, 20))

library(tidyverse)
library(lubridate)

minute_data <- read_rds("data/minute_data.Rds")

minute_summary <- minute_data %>% 
  mutate(wind_spd = ifelse(wind_spd < 0, NA, wind_spd)) %>%
  group_by(site, year, month, day, hour) %>% 
  summarize(gust_time = sum(wind_spd >= 10), wind_mean = mean(wind_spd), 
            wind_sd = sd(wind_spd), 
            wind_max = max(wind_spd), wind_min = min(wind_spd), 
            wind_med = median(wind_spd), nrow = n(), 
            na_ct = sum(is.na(wind_spd)),
            .groups = "drop")

clean_minute_summary <- minute_summary %>% filter(nrow == 60, na_ct == 0) %>% 
  mutate(date = make_datetime(year, month, day, hour))

full_dates <- clean_minute_summary$date

new_combined_cc <- combined_cc %>% 
  mutate(date = make_datetime(year, month, day, hour)) %>% 
  filter(date %in% full_dates)

new_combined_cc <- new_combined_cc %>% 
  left_join(select(clean_minute_summary, site, date, new_gust_time = gust_time, 
                   wind_mean, wind_sd, wind_max, wind_min, wind_med), 
            by = c("site", "date"))

rm(minute_data)
gc()

ggplot(new_combined_cc, aes(x = wind_spd, y = new_gust_time)) + geom_smooth()
ggplot(new_combined_cc, aes(x = wind_spd, y = wind_sd)) + geom_smooth()
ggplot(new_combined_cc, aes(x = wind_spd, y = wind_max)) + geom_smooth()
ggplot(new_combined_cc, aes(x = wind_spd, y = wind_min)) + geom_smooth()

new_combined_cc %>% 
  select(wind_spd, wind_max, wind_min, wind_sd, gust_time) %>%
  pivot_longer(cols = -wind_spd, names_to = "var", values_to = "value") %>%
  ggplot(aes(x = wind_spd, y = value, color = var)) + 
  geom_smooth() +
  scale_color_brewer(palette = "Dark2", name = "Variable") +
  labs(x = "Wind speed", y = "Minute data")

dat_slim <- list(
  gust_time = new_combined_cc$gust_time,
  t_2m = new_combined_cc$t_2m
)
str(dat_slim)

mgusty3 <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) ,
  data=dat_slim , chains=1 )

mgusty3.1 <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) ,
  data=dat_slim , chains=4 , cores=4 , iter = 1000 )

t_2m.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_2m = t_2m.seq )
mu <- link( mgusty3.1 , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.gust_time <- sim( mgusty3.1 , data=dat_slim )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

plot( gust_time ~ t_2m , dat_slim_df , col=col.alpha(rangi2,0.5) )
lines( t_2m.seq , mu.mean )
shade( mu.PI , t_2m.seq )

dat_slim_df <- as_tibble(dat_slim)

ggplot(dat_slim_df, aes(t_2m, gust_time)) +  geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))
ggplot(dat_slim_df, aes(t_2m, gust_time)) + geom_point(color = alpha(rangi2, 0.1)) + geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))

plot3.1 <- ggplot(dat_slim_df, aes(t_2m, gust_time)) + geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))
plot3.1 + coord_cartesian(ylim = c(0, 10))


dat_slim_2 <- list(
  gust_time = new_combined_cc2$gust_time,
  wind_spd = new_combined_cc2$wind_spd
)
str(dat_slim_2)

mgusty4 <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * wind_spd^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , 
  data=dat_slim_2 , chains=1 )

mgusty4.1 <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * wind_spd + b2 * wind_spd^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , 
  data=dat_slim_2 , chains=4 , cores=4 , iter = 1000)

library(tidybayes)
library(tidybayes.rethinking)
library(bayesplot)
mcmc_dens(mgusty4.1@stanfit)

post_mgusty3 <- extract.samples(mgusty3)
post_mgusty3.1 <- extract.samples(mgusty3.1)

dens(post_mgusty3$sigma, lwd = 1)
dens(post_mgusty3.1$sigma, lwd = 1, col = rangi2, add = TRUE)

mgusty3_draws <- tidy_draws(mgusty3) %>% mutate(model = "mgusty3")
mgusty3.1_draws <- tidy_draws(mgusty3.1) %>% mutate(model = "mgusty3.1")

bind_rows(mgusty3_draws, mgusty3.1_draws) %>% 
  ggplot(aes(x = sigma, color = model, fill = model)) + 
  geom_density(size = 1, alpha = 0.3)

post_mgusty4 <- extract.samples(mgusty4)
post_mgusty4.1 <- extract.samples(mgusty4.1)

dens(post_mgusty4$sigma, lwd = 1)
dens(post_mgusty4.1$sigma, lwd = 1, col = rangi2, add = TRUE)

mgusty4_draws <- tidy_draws(mgusty4) %>% mutate(model = "mgusty4")
mgusty4.1_draws <- tidy_draws(mgusty4.1) %>% mutate(model = "mgusty4.1")

bind_rows(mgusty4_draws, mgusty4.1_draws) %>% 
  ggplot(aes(x = sigma, color = model, fill = model)) + 
  geom_density(size = 1, alpha = 0.3)

wind_spd.seq <- seq( from=0 , to=100 , length.out=30 )
pred_dat <- list( wind_spd_s=wind_spd.seq , wind_spd_s2=wind_spd.seq^2 )
mu <- link( mgusty4.1 , data=dat_slim_2 )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory
sim.gust_time <- sim( mgusty4.1 , data=dat_slim_2 )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

dat_slim_2_df <- as_tibble(dat_slim_2)

ggplot(dat_slim_2_df, aes(wind_spd, gust_time)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

plot4.1 <- ggplot(dat_slim_2_df, aes(wind_spd, gust_time)) +
  geom_smooth()
lines( wind_spd.seq , mu.mean )
shade( mu.PI , wind_spd.seq )
shade( height.PI , wind_spd.seq )

plot4.1 + coord_cartesian(xlim = c(0, 20))

newest_combined_cc <- new_combined_cc
print(newest_combined_cc)

newest_combined_cc <- mutate(newest_combined_cc, t_difference = t_10m-t_2m)

mgusty5 <- quap(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) , data=newest_combined_cc)

t_difference.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_difference = t_difference.seq )
mu <- link( mgusty5 , data=newest_combined_cc )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.gust_time <- sim( mgusty5 , data=newest_combined_cc )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

ggplot(newest_combined_cc, aes(t_difference, gust_time)) +
  geom_smooth()

mgusty6 <- ulam(
  alist(
    gust_time ~ dnorm( mu , sigma ) ,
    mu <- a + b1 * t_2m + b2 * t_2m^2,
    a ~ dnorm(5, 5),
    b1 ~ dnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dexp(0.05)
  ) ,
  data=dat_slim , chains=4 , cores=4 , iter = 1000 )

t_difference.seq <- seq( from=-49.0 , to=22.8 , length.out=30 )
pred_dat <- list( t_difference = t_difference.seq )
mu <- link( mgusty6 , data=newest_combined_cc )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.90 )
rm(mu) # delete the huge matrix mu now that we're done with it.
gc() # garbage-collect to free R's unused memory

sim.gust_time <- sim( mgusty6 , data=newest_combined_cc )
gust_time.PI <- apply( sim.gust_time , 2 , PI , prob=0.90 )
rm(sim.gust_time)
gc()

newest_combined_cc_df <- as_tibble(newest_combined_cc)

plot( gust_time ~ t_difference , newest_combined_cc_df , col=col.alpha(rangi2,0.1))
lines( t_difference.seq , mu.mean )
shade( mu.PI , t_difference.seq )

ggplot(newest_combined_cc_df, aes(t_difference, gust_time)) +
  geom_smooth() + xlim(-2, 2) + ylim (0, 61)

library(brms)

mgusty3.1prior <- c(
  set_prior("normal(5, 5)", class = "Intercept", lb = 0),
  set_prior("normal(0, 1)", class = "b", coef = "t_2m"),
  set_prior("exponential(0.05)", class="sigma"))
mgusty3.1_fit_formula <- bf(gust_time ~ t_2m + I(t_2m^2))
b3.1 <-
  brm(mgusty3.1_fit_formula,
    data = new_combined_cc,
      family = gaussian,
      prior = mgusty3.1prior,
      chains = 4, cores = 4, backend = "cmdstanr" )

mgusty4.1prior <- c(
  set_prior("normal(5, 5)", class = "Intercept", lb = 0),
  set_prior("normal(0, 1)", class = "b", coef = "wind_spd"),
  set_prior("exponential(0.05)", class="sigma"))
mgusty4.1_fit_formula <- bf(gust_time ~ wind_spd + I(wind_spd^2))
b4.1 <-
  brm(mgusty4.1_fit_formula,
      data = new_combined_cc,
      family = gaussian,
      prior = mgusty4.1prior,
      chains = 4, cores = 4, backend = "cmdstanr" )

mgusty5.1prior <- c(
  set_prior("normal(5, 5)", class = "Intercept", lb = 0),
  set_prior("normal(0, 1)", class = "b", coef = "t_difference"),
  set_prior("exponential(0.05)", class="sigma"))
mgusty5.1_fit_formula <- bf(gust_time ~ t_difference + I(t_difference^2))
b5.1 <-
  brm(mgusty5.1_fit_formula,
      data = newest_combined_cc_df,
      family = gaussian,
      prior = mgusty5.1prior,
      chains = 4, cores = 4, backend = "cmdstanr" )

ggplot(new_combined_cc, aes(x=t_2m, y=GustTime)) +
  geom_point(size=3) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.2, fill=colourcodes[3]) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.5, fill=colourcodes[3]) +
  geom_line(data=hill_bayes_ribbons, aes(y=Estimate), colour=colourcodes[3], 
            size=1)

ggplot(data = newest_combined_cc, aes (x = t_difference)) + geom_density() + xlim(-2, 2)

ggplot(newest_combined_cc, aes(x=t_difference, y=gust_time)) + 
  geom_bin2d(bins = 100, mapping = aes(fill = log(..ndensity..))) +
  theme_bw()

ggplot(newest_combined_cc, aes(x=t_difference, y=gust_time)) + 
  geom_bin2d(bins = c(30,30), mapping = aes(fill = log(..ndensity..))) +
  scale_fill_viridis_c(option = "B") +
  xlim(-2,2) +
  theme_bw()

Richardson_combined_cc <- mutate(newest_combined_cc, RichardsonBulk = ((9.81/t_2m)*(t_difference)*(8))/(wind_spd^2))

ggplot(Richardson_combined_cc, aes(x=RichardsonBulk, y=gust_time)) +
  geom_bin2d(bins = c(30,30), mapping = aes(fill = log(..ndensity..))) +
  scale_fill_viridis_c(option = "B") +
  xlim(-1,1) +
  theme_bw()

# use hydrostatic equation to find potential temperature
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
Richardson_combined_cc_new <- mutate(deltatheta_cc, RichardsonBulkPotential = ((9.81/t_2m)*(DeltaTheta)*(8))/(wind_spd^2))

ggplot(Richardson_combined_cc_new, aes(x=RichardsonBulkPotential, y=gust_time)) +
  geom_bin2d(bins = c(30,30), mapping = aes(fill = log(..ndensity..))) +
  scale_fill_viridis_c(option = "B") +
  xlim(-1,1) +
  theme_bw()
ggplot(Richardson_combined_cc_new, aes(x=RichardsonBulkPotential, y=gust_time)) +
  xlim(-1,1) +
  geom_point(alpha = 0.1, na.rm = TRUE)
