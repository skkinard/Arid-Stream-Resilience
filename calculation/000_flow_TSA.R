# flow_record_2017_to_2020
# Sean Kinard
# 2023-02-15

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
rm(list=ls())

setwd('/home/kinard/Documents/Research/Dissertation/02_Resilience')

library(tidyverse) # all tidyverse packages
library(lubridate) # date managment
library(ggdark) # cool dark ggplot theme
library(modeest)

# -----------------------------------------------------------------------------
# Data Prep
# -----------------------------------------------------------------------------
# load data
site_order <- c("TR", "SF", "AR", "MR", "PD", "PL", "GC", "WM", "EM")
my_colors <- c('#F5191CFF', '#E78200FF', '#E8A117FF',
               '#EABB22FF', '#CBC988FF', '#9FC095FF',
               '#6BB699FF', '#3DAAA6FF', '#3B99B1FF')

d_site <- read_csv('Data/site_data.csv')%>%
  select(Site, AnnualRain)%>%
  unique()%>%
  rename(site_code=Site, annualrain=AnnualRain) %>%
  mutate(site_code=substr(site_code, 1,2))

d_flow <- read_csv('Data/flow_record_allsites_1990_to_2020.csv') %>%
  mutate(week = week(collection_date)) %>%
  left_join(d_site) %>%
  mutate(site_code = fct_relevel(site_code, site_order)) %>%
  mutate(alltime = 'all') %>% # create dummie variable for site: alltime
  mutate(seas = case_when(
    month %in% c(11,12,1,2,3) ~ 'cool_dry',
    month %in% c(4,5,6) ~ 'cool_wet',
    month %in% c(7,8) ~ 'hot_dry',
    month %in% c(9,10) ~ 'hot_wet'),
    ymonth = paste(year, month, sep='-'),
    yqtr = paste(year, qtr, sep='-'),
    yweek = paste(year, week, sep='-'),
    yseas = paste(year, seas, sep='-'))

# remove NA and site-years with less than 24 measurements
d_flow <- d_flow %>%
  na.omit() %>%
  group_by(site_code, year) %>%
  dplyr::summarize(n = length(q)) %>%
  right_join(d_flow%>%na.omit()) %>%
  filter(n > 24) %>%
  select(!n)

# -----------------------------------------------------------------------------
# Calculation Functions
# -----------------------------------------------------------------------------
make_DC <- function(x){
  if(length(x) <= 5) {return(NA)}
  else {
  x_out <- rep(NA,length(x))
  
  for(i in 2:length(x)) {
    x_out[i] <- abs(x[i]-x[i-1])
  }
  
  x_out[1] <- median(x_out, na.rm=T)
  
  return(x_out) }
}
mu_under_q <- function(x, q) {
  
  qx <- quantile(x, probs=c(q), na.rm=T)
  qx_under <- x[c(which(x<=qx))]
  qx_mu <- mean(qx_under)
  
  return(qx_mu)
}
mu_over_q <- function(x, q) {
  
  qx <- quantile(x, probs=c(q), na.rm=T)
  qx_over <- x[c(which(x>=qx))]
  qx_mu <- mean(qx_over)
  
  return(qx_mu)
}
make_stats <- function(my_data, time_scale) {
  
  # ---------------------------------------------------------------------------
  # Summarize Statistics
  my_stats <- my_data %>%
    na.omit() %>%
    rename(my_period = time_scale) %>%
    group_by(site_code, my_period) %>%
    dplyr::summarize(
      XXX_med = median(q, na.rm=T),
      XXX_min = min(q, na.rm=T),
      XXX_max = max(q, na.rm=T),
      XXX_mu = mean(q, na.rm=T),
      XXX_sd = sd(q, na.rm=T),
      XXX_rsd = XXX_sd / XXX_med,
      XXX_q25 = quantile(q, probs=c(.25), na.rm=T),
      XXX_q75 = quantile(q, probs=c(.75), na.rm=T),
      XXX_n = length(q),
      XXX_WCDC = sum((make_DC(q)/XXX_med)^2, na.rm=T)/XXX_n,
      XXX_CDC =  sum((make_DC(q)/XXX_med), na.rm=T)/XXX_n,
      XXX_FS = abs(mu_over_q(q, .85)- XXX_med) / XXX_med,
      XXX_DS = abs(mu_under_q(q, .15)- XXX_med) / XXX_med) %>%
    pivot_longer(cols=contains('XXX'), names_to='x', values_to='y') %>%
    filter(!is.na(y)) %>%
    filter(!is.infinite(y)) %>%
    mutate(x=str_replace_all(x, 'XXX', time_scale)) %>%
    separate(x, into=c('scale', 'stat')) %>%
    pivot_wider(names_from=stat, values_from=y) %>%
    rename({{time_scale}} := my_period)
  
  # ---------------------------------------------------------------------------
  # Pulse Percentages
  PP <- my_stats %>%
    rename(my_period = time_scale) %>%
    select(site_code, my_period, q25, med) %>%
    right_join(d_flow) %>%
    mutate(LF = ifelse(q < q25, 1, 0),
           HF3 = ifelse(q > 3*med, 1, 0),
           HF7 = ifelse(q > 7*med, 1, 0),
           HF15 = ifelse(q > 15*med, 1, 0)) %>%
    group_by(site_code, my_period) %>%
    dplyr::summarize(LFPP = sum(LF)/length(LF),
                     HFPP3 = sum(HF3)/length(LF),
                     HFPP7 = sum(HF7)/length(LF),
                     HFPP15 = sum(HF15)/length(LF)) %>%
    rename({{time_scale}} := my_period)
  
  # ---------------------------------------------------------------------------
  # merge
  
  my_output <- left_join(my_stats, PP)
  
  return(my_output)
  
}
combine_scale_stats <- function(my_data) {
  my_data %>%
    pivot_longer(cols=med:HFPP15, names_to='nam', values_to='val') %>%
    mutate(ID=paste(substr(scale,1,2), nam, sep='_')) %>%
    select(-scale, -nam) %>%
    pivot_wider(names_from=ID, values_from=val)
}

# -----------------------------------------------------------------------------
# Calculations
# -----------------------------------------------------------------------------

fstats_alltime <- make_stats(d_flow, "alltime")
fstats_year <- make_stats(d_flow, "year")
fstats_qtr <- make_stats(d_flow, "qtr")
fstats_yqtr <- make_stats(d_flow, "yqtr")
fstats_month <- make_stats(d_flow, "month")
fstats_ymonth <- make_stats(d_flow, "ymonth")
fstats_week <- make_stats(d_flow, "week")
fstats_yweek <- make_stats(d_flow, "yweek")
fstats_seas <- make_stats(d_flow, "seas") %>%
  mutate(seas = fct_relevel(seas, c("cool_dry", "cool_wet", 
                                    "hot_dry", "hot_wet")))
fstats_yseas <- make_stats(d_flow, "yseas")%>%
  mutate(seas = fct_relevel(seas, c("cool_dry", "cool_wet", 
                                    "hot_dry", "hot_wet")))

# Megaframe
d_mega <- d_flow %>% 
  left_join(fstats_alltime %>% combine_scale_stats()) %>% 
  left_join(fstats_year %>% combine_scale_stats()) %>% 
  left_join(fstats_qtr %>% combine_scale_stats()) %>% 
  left_join(fstats_seas %>% combine_scale_stats()) %>% 
  left_join(fstats_month %>% combine_scale_stats()) %>% 
  left_join(fstats_week %>% combine_scale_stats()) %>%
  left_join(fstats_yqtr %>% combine_scale_stats()) %>%
  left_join(fstats_ymonth %>% combine_scale_stats()) %>%
  left_join(fstats_yweek %>% combine_scale_stats()) %>%
  left_join(fstats_yseas %>% combine_scale_stats())

# -----------------------------------------------------------------------------
# Visualization: Time Series
# -----------------------------------------------------------------------------
ts_base <- function(my_data) {
  my_data %>%
    ggplot(aes(collection_date, y=q_region)) +
    geom_point(alpha=.2, color='skyblue') +
    dark_theme_gray(base_size=14) +
    scale_x_date(date_breaks = '1 year', date_labels = "%y") +
    scale_color_manual(values=my_colors) +
    ylab('Regional Discharge') +
    scale_y_log10() +
    xlab('Year') +
    ggtitle('Regional Discharge Time Series: Q_region = mean( daily_Q / site_Q_mu )')
}

regional_discharge_time_series_longterm <- d_mega %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_mu, na.rm=T)) %>%
  ts_base() +
  geom_smooth(method='loess', se=F, span=.2, color='red3')

region_flow_time_series_shortterm <- d_mega %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_mu, na.rm=T)) %>%
  mutate(time_period = ifelse(collection_date <= ymd('1997-01-01'), '1990-1997',
                       ifelse(collection_date > ymd('1997-01-01') &
                              collection_date <= ymd('2004-01-01'), '1997-2004',
                      ifelse(collection_date >= ymd('2004-01-01') &
                             collection_date < ymd('2011-01-01'), '2004-2011', '2011-2021')))) %>%
  ts_base() +
  facet_wrap(~time_period, scales='free_x', ncol=1) +
  geom_smooth(method='loess', se=F, span=.025, color='red3')
  
# -----------------------------------------------------------------------------
# Histogram
# -----------------------------------------------------------------------------
hist_base <- function(my_data) {
  my_data %>%
    ggplot(aes(log(q), fill=site_code, color = site_code)) +
    geom_density(alpha=.1) +
    scale_color_manual(values=my_colors) +
    scale_fill_manual(values=my_colors) +
    dark_theme_grey() +
    ggtitle('Smoothed Density of Daily flows 1990-2020') +
    ylab(element_blank()) +
    xlab('Log Discharge')
}

flow_density <- d_flow %>%
  hist_base() 

flow_density_facet <- d_flow %>%
  hist_base() +
  facet_wrap(~site_code, ncol=3)

# Visualization: Season
# -----------------------------------------------------------------------------

# Seasonal flow Oscillation
seas_flow_osc <- d_mega %>%
  select(year, week, we_mu, ye_mu, site_code, annualrain) %>%
  unique() %>%
  ggplot(aes(week, y=we_mu/ye_mu, color=site_code)) +
  geom_smooth(method = "loess", se=F, span=.4) +
  dark_theme_gray() +
  ylab('weekly_avg / yr_avg') +
  scale_x_continuous(breaks = seq(from=1, to=53, by=4)) +
  scale_color_manual(values=my_colors) +
  scale_y_log10() +
  ggtitle('Seasonal Oscillations in Daily Discharge')

# new
d_mega %>%
  select(year, week, mo_mu, ye_mu, site_code, annualrain) %>%
  unique() %>%
  ggplot(aes(week, y=we_mu/ye_mu, color=site_code)) +
  geom_smooth(method = "loess", se=F, span=.4) +
  dark_theme_gray() +
  ylab('weekly_avg / yr_avg') +
  scale_x_continuous(breaks = seq(from=1, to=53, by=4)) +
  scale_color_manual(values=my_colors) +
  scale_y_log10() +
  ggtitle('Seasonal Oscillations in Daily Discharge')



pulse_weekly <- d_flow %>%
  select(site_code, wk_LFPP, wk_HFPP3, wk_HFPP7, wk_HFPP15) %>%
  pivot_longer(cols = wk_LFPP:wk_HFPP15, names_to='variable', values_to = 'values') %>%
  unique() %>%
  ggplot(aes(site_code, values)) +
  facet_wrap(~variable, scales='free') +
  geom_boxplot() +
  dark_theme_grey()

# monthly
d_flow <- d_flow %>%
  mutate(mo_LF = ifelse(q < mo_q25, 1, 0),
         mo_HF3 = ifelse(q > 3*mo_mu, 1, 0),
         mo_HF7 = ifelse(q > 7*mo_mu, 1, 0),
         mo_HF15 = ifelse(q > 15*mo_mu, 1, 0)) %>%
  group_by(site_code, month) %>%
  dplyr::summarize(mo_LFPP = sum(mo_LF)/length(mo_LF),
                   mo_HFPP3 = sum(mo_HF3)/length(mo_LF),
                   mo_HFPP7 = sum(mo_HF7)/length(mo_LF),
                   mo_HFPP15 = sum(mo_HF15)/length(mo_LF)) %>%
  left_join(d_flow)  

pulse_monthly <- d_flow %>%
  select(site_code, mo_LFPP, mo_HFPP3, mo_HFPP7, mo_HFPP15) %>%
  pivot_longer(cols = mo_LFPP:mo_HFPP15, names_to='variable', values_to = 'values') %>%
  unique() %>%
  ggplot(aes(site_code, values)) +
  facet_wrap(~variable, scales='free') +
  geom_boxplot() +
  dark_theme_grey()

# site
d_flow <- d_flow %>%
  mutate(si_LF = ifelse(q < si_q25, 1, 0),
         si_HF3 = ifelse(q > 3*si_mu, 1, 0),
         si_HF7 = ifelse(q > 7*si_mu, 1, 0),
         si_HF15 = ifelse(q > 15*si_mu, 1, 0)) %>%
  group_by(site_code) %>%
  dplyr::summarize(si_LFPP = sum(si_LF)/length(si_LF),
                   si_HFPP3 = sum(si_HF3)/length(si_LF),
                   si_HFPP7 = sum(si_HF7)/length(si_LF),
                   si_HFPP15 = sum(si_HF15)/length(si_LF)) %>%
  left_join(d_flow)  

pulse_site <- d_flow %>%
  select(site_code, si_LFPP, si_HFPP3, si_HFPP7, si_HFPP15) %>%
  pivot_longer(cols = si_LFPP:si_HFPP15, names_to='variable', values_to = 'values') %>%
  unique() %>%
  ggplot(aes(site_code, values)) +
  facet_wrap(~variable, scales='free') +
  geom_point(size = 3, color = 'skyblue') +
  dark_theme_grey()

# -----------------------------------------------------------------------------
# Flow Variation (RSD, CDC, WCDC)
# -----------------------------------------------------------------------------

# Relative standard deviation discharge at different scales
p_rsd <- d_flow %>%
  select(site_code, contains('rsd')) %>%
  unique() %>%
  pivot_longer(cols=contains('rsd'), names_to = 'variable', values_to = 'value') %>%
  ggplot(aes(site_code, value)) +
  facet_wrap(~variable) +
  geom_point() +
  scale_y_log10() +
  dark_theme_grey()

# Daily changes
make_DC <- function(x){
  x_out <- rep(NA,length(x))
  
  for(i in 2:length(x)) {
    x_out[i] <- abs(x[i]-x[i-1])
  }
  
  x_out[1] <- median(x_out, na.rm=T)
  
  return(x_out)
}

# Weighted Cumulative Relative Daily Changes per year
# weighting refers to the squared daily changes which amplifies larger swings compared to smaller swings
# CDC is unweighted; so there's no detection for few-large changes vs numerous-small changes
d_flow <- d_flow %>%
  group_by(site_code) %>%
  dplyr::summarize(WCDC = sum((make_DC(q)/si_mu)^2, na.rm=T)/(length(q)/365),
                   CDC =  sum((make_DC(q)/si_mu), na.rm=T)/(length(q)/365)) %>%
  left_join(d_flow)

flow_variability <- d_flow %>%
  select(site_code, annualrain, si_rsd, WCDC, CDC) %>%
  unique() %>%
  mutate(seasonality_weigthed = scale(WCDC),
         seasonality = scale(CDC),
         relative_sd = scale(si_rsd)) %>%
  pivot_longer(cols=c(seasonality_weigthed, seasonality, relative_sd), 
               names_to='variable', values_to='value') %>%
  ggplot(aes(annualrain, value, fill=variable, shape=variable)) +
  geom_line(aes(x=annualrain, y=value),
            linetype=2,color=NA) +
  geom_point(size = 9) +
  dark_theme_grey(base_size=12) +
  scale_fill_manual(values = c('red', 'skyblue', 'blue')) +
  scale_shape_manual(values = c(21,22,23)) +
  ylab(element_blank()) +
  xlab(element_blank())
# seasonality overpredicts variation at the ends of the gradient, where flows are either consistently low or consistenly high, and underpredicts variation in the middle, where flows experience wide swings. This is due to the fact that the sum of changes for few-large swings is equivalent to the sum of changes for numerous-small swings. So, I modified the seasonality calculation to be the cumalitve day-to-day changes squared : sum of ((q_i-q_i-1) / q_mu)^2. Since not all time series were 20 years in length, I divided the seasonality values by the number of days/365. Lastly, in order to make comparisons between indices, I scaled the variables. Here, I also report the scaled, relative standard deviation which should approximate stream variability throughout the year.

# -----------------------------------------------------------------------------
# Typical flood & Drought Strength
# -----------------------------------------------------------------------------
make_Typical_Strengths <- function(my_data, time_scale, 
                                   time_abbreviation) {
  # select columns with my abbreviation
  my_output <- my_data %>%
    select(site_code, collection_date, q,
           year, month, week, alltime,
           contains(time_abbreviation)) %>%
    rename('my_period'=time_scale)
  
  # remove my abbreviation from columns
  colnames(my_output) <- str_replace_all(
    colnames(my_output), 
    paste(time_abbreviation, '_', sep=''),
    'xxx_')
  
  # calculate TFS
  my_output <- my_output %>%
    filter(q>=xxx_q75) %>%
    group_by(site_code, my_period) %>%
    dplyr::summarize(xxx_q75_mu = mean(q, na.rm=T)) %>%
    left_join(my_output) %>%
    mutate(xxx_TFS = abs(xxx_q75_mu - xxx_med) / xxx_med )
  
  # calculate TDS
  my_output <- my_output %>%
    filter(q<=xxx_q25) %>%
    group_by(site_code, my_period) %>%
    dplyr::summarize(xxx_q25_mu = mean(q, na.rm=T)) %>%
    left_join(my_output) %>%
    mutate(xxx_TDS = abs(xxx_q25_mu - xxx_med) / xxx_med )
  
  # colnames
  colnames(my_output) <- str_replace_all(colnames(my_output),
                                         'xxx', time_abbreviation)
  
  # Merge to original data
  my_output <- left_join(my_data %>% 
                           select(
                             !contains(paste(time_abbreviation, '_', sep=''))), 
                         my_output %>% select(!my_period))
  
  return(my_output) }

# Iterate function across scales of comparison

# weekly
d_flow <- make_Typical_Strengths(my_data = d_flow,
                                 time_scale = 'week',
                                 time_abbreviation = 'wk')
# monthly
d_flow <- make_Typical_Strengths(my_data = d_flow,
                                 time_scale = 'month',
                                 time_abbreviation = 'mo')

# Annually
d_flow <- make_Typical_Strengths(my_data = d_flow,
                                 time_scale = 'year',
                                 time_abbreviation = 'yr')

# site
d_flow <- make_Typical_Strengths(my_data = d_flow,
                                 time_scale = 'alltime',
                                 time_abbreviation = 'si')


# visualize outputs
TS_base <- function(my_data) {
  my_data %>%
    ggplot(aes(site_code, value, fill=annualrain)) +
    facet_wrap(~time_period, ncol=1, scales='free') +
    geom_boxplot(aes()) +
    geom_point(aes(y=mu), size=4, shape = 23, color='white', fill='black') +
    dark_theme_grey(base_size=12) +
    paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
    xlab(element_blank()) +
    ylab(element_blank())
}

TS_means <- d_flow %>%
  group_by(site_code) %>%
  dplyr::summarize(si_TFS_mu = mean(si_TFS, na.rm=T),
                   yr_TFS_mu = mean(yr_TFS, na.rm=T),
                   mo_TFS_mu = mean(mo_TFS, na.rm=T),
                   wk_TFS_mu = mean(wk_TFS, na.rm=T),
                   si_TDS_mu = mean(si_TDS, na.rm=T),
                   yr_TDS_mu = mean(yr_TDS, na.rm=T),
                   mo_TDS_mu = mean(mo_TDS, na.rm=T),
                   wk_TDS_mu = mean(wk_TDS, na.rm=T)) %>%
  pivot_longer(cols=contains('mu'), names_to='variable', values_to='value') %>%
  separate(variable, into = c('time_period', 'metric', 'stat'), sep='_') %>%
  select(!stat) %>%
  rename(mu=value)

p_TDS <- d_flow  %>%
  select(site_code, annualrain, contains('TDS')) %>%
  pivot_longer(cols=c('si_TDS', 'yr_TDS', 'mo_TDS', 'wk_TDS'),
               names_to='variable', values_to='value') %>%
  unique() %>%
  separate(variable, into = c('time_period', 'metric'), sep='_') %>%
  left_join(TS_means) %>%
  mutate(time_period=fct_relevel(time_period, 'si', 'yr', 'mo', 'wk')) %>%
  TS_base() +
  ggtitle('Typical Drought Strength = (mean(q<25%ile) - q_mu) / q_mu')


p_TFS <- d_flow  %>%
  select(site_code, annualrain, contains('TFS')) %>%
  pivot_longer(cols=c('si_TFS', 'yr_TFS', 'mo_TFS', 'wk_TFS'),
               names_to='variable', values_to='value') %>%
  unique() %>%
  separate(variable, into = c('time_period', 'metric'), sep='_') %>%
  left_join(TS_means) %>%
  mutate(time_period=fct_relevel(time_period, 'si', 'yr', 'mo', 'wk')) %>%
  TS_base() +
  ggtitle('Typical Flood Strength = (mean(q>75%ile) - q_mu) / q_mu')

# -----------------------------------------------------------------------------
# Flow Metric Summary
# -----------------------------------------------------------------------------

make_stat_mu <- function(my_data, my_abbreviation) {
  colnames(my_data) <- str_replace_all(colnames(my_data),
                                       paste(my_abbreviation, '_',sep=''),
                                       '')
  
  my_output <- my_data %>%
    group_by(site_code) %>%
    dplyr::summarize(XXX_med_mu = mean(med, na.rm=T),
                     XXX_min_mu = mean(min, na.rm=T),
                     XXX_max_mu = mean(max, na.rm=T),
                     XXX_mu_mu = mean(mu, na.rm=T),
                     XXX_sd_mu = mean(sd, na.rm=T),
                     XXX_rsd_mu = XXX_sd_mu / XXX_mu_mu,
                     XXX_q25_mu = mean(q25, probs=c(.25), na.rm=T),
                     XXX_q75_mu = mean(q75, probs=c(.75), na.rm=T))
  
  colnames(my_output) <- str_replace_all(colnames(my_output),
                                         'XXX',
                                         my_abbreviation)
  
  return(my_output) }

site_stats_long <- stats_site %>% 
  left_join(make_stat_mu(stats_week, 'wk')) %>%
  left_join(make_stat_mu(stats_month, 'mo')) %>%
  left_join(make_stat_mu(stats_year%>%filter(!is.na(yr_sd)), 'yr')) %>%
  rename(site=site_code) %>%
  pivot_longer(cols=contains('_'), names_to='metric', values_to='value') %>%
  separate(metric, into=c('time_period', 'stat'), sep='_') %>%
  rename(site_code=site) %>%
  mutate(time_period = fct_relevel(time_period, 'si', 'yr', 'mo', 'wk')) %>%
  full_join(TS_means %>%
  rename(stat=metric, value=mu) ) %>%
  full_join(d_flow %>%
  select(site_code, CDC, WCDC) %>%
  unique() %>%
  pivot_longer(cols=c(CDC,WCDC), names_to='stat', values_to='value') %>%
  mutate(time_period = 'si') ) %>%
  left_join(d_flow%>%select(site_code,annualrain)%>%unique()) %>%
  arrange(site_code, time_period, stat)

  
scatterplot_stats_all_periods <- site_stats_long %>%
  filter(time_period != 'wk') %>%
  filter(!is.na(value)) %>%
  filter(!is.infinite(value)) %>%
  ggplot(aes(annualrain, value, fill=time_period, shape=time_period)) +
  facet_wrap(~stat, scales='free') +
  geom_jitter(size=5, color='black', alpha=.1,
              position = position_dodge(width=5)) +
  geom_jitter(size=5, aes(color=time_period), fill=NA,
              position = position_dodge(width=5)) +
  dark_theme_grey(base_size = 10) +
  scale_shape_manual(values = c(21,22,24)) +
  scale_fill_manual(values = my_colors[c(1,4,9)]) +
  scale_color_manual(values = my_colors[c(1,4,9)]) +
  scale_y_log10()
  
# Choosing my hydro-metrics:
# WCDC:site
# max:yr
# med:yr
# min:yr
# mu:site
# q25:yr
# q75:yr
# 

  
  
  
  pivot_wider(names_from=stat, values_from=value) %>%





# hurricane harvey August 2017
d_flow %>%
  filter(collection_date < ymd('2017-09-15') & 
           collection_date > ymd('2017-08-15')) %>%
  flow_base() +
  annotate("rect", xmin = ymd('2017-08-21'), xmax = ymd('2017-09-01'), 
           ymin = 0, ymax = 1000000, alpha = .2, fill = 'grey80', color = NA) +
  geom_vline(xintercept=ymd('2017-08-26'), color = 'chartreuse', alpha = .5) +
  facet_wrap(~site_code, ncol=3, scales='free') 

# flood metrics: max flow, duration, base flow

# base flow (moving average vs 30 prior to hurricane)
colnames(d_flow)
d_flow %>%
  arrange(site_code, collection_date) %>%
  group_by(site_code) %>%
  mutate(avg_30day = rollmean(q, k=30))


moving_average <- function(series, klags) {
  return(
    lag(
      zoo::rollmean(series, klags, fill = NA), floor(klags / 2)
    )
  )
}

moving_averages <- d_flow %>%
  select(! year:qtr) %>%
  arrange(site_code, collection_date) %>%
  group_by(site_code) %>%
  mutate(
    MA_3 = moving_average(q, 3),
    MA_4 = moving_average(q, 4),
    MA_5 = moving_average(q, 5),
    MA_6 = moving_average(q, 6),
  ) %>%
  ungroup()

moving_averages %>% filter(site_code =='AR')


# max flow [xmin = ymd('2017-08-21'), xmax = ymd('2017-09-01')]


  



# Region flood 2018
d_flow %>%
  filter(collection_date < ymd('2018-08-01') & 
           collection_date > ymd('2018-06-01')) %>%
  flow_base() +
  facet_wrap(~site_code, ncol=3, scales='free') +
  geom_vline(xintercept=ymd('2018-06-18'))

