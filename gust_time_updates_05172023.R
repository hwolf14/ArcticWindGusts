#update gust_time parameter so it includes only points 50% above mean wind spd

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

updated_combined <- read_rds("data/combined_data.Rds")

updated_combined_cc <- updated_combined %>% drop_na(gust_time, t_2m)

updated_new_combined_cc <- updated_combined_cc %>% 
  mutate(date = make_datetime(year, month, day, hour)) %>% 
  filter(date %in% full_dates)

updated_new_combined_cc <- updated_new_combined_cc %>% 
  left_join(select(clean_minute_summary, site, date, new_gust_time = gust_time, 
                   wind_mean, wind_sd, wind_max, wind_min, wind_med), 
            by = c("site", "date"))

rm(minute_data)
gc

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

##this plot shows some gust_times at 60? Should not be possible. Something is wrong.
#a calculation is wrong? combined_data needs to be updated?

ggplot(updated_new_combined_cc, aes(t_2m, gust_time)) +
  geom_smooth()

mgusty1plot_updated <- ggplot(updated_new_combined_cc, aes(t_2m, gust_time)) +
  geom_smooth()
mgusty1plot_updated + coord_cartesian(ylim = c(0, 10))

gam_updated <- gam(gust_time ~ s(t_2m + t_2m^2), data = updated_new_combined_cc)

plot(gam, add = TRUE, ylim = c(0,60), ylab = "Gust time", xlab = "temperature")
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


plot3.1 <- ggplot(dat_slim_df, aes(t_2m, gust_time)) + geom_line(data = dat_smooth, linewidth = 1) + geom_ribbon(data = dat_smooth, aes(ymin = low, ymax = high), fill = alpha("darkgreen", 0.3))
plot3.1 + coord_cartesian(ylim = c(0, 10))