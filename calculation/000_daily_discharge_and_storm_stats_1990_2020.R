# daily_discharge_and_storm_stats_1990_2020
# Sean Kinard
# 2023-02-24

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
# rm(list=ls())

# setwd('/home/kinard/Documents/Research/Dissertation/02_Resilience')

# library(tidyverse) # all tidyverse packages
# library(lubridate) # date managment
# library(ggdark) # cool dark ggplot theme

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
  dplyr::rename(site_code=Site, annualrain=AnnualRain) %>%
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
    x_out <- rep(NA,length(x))
  
  for(i in 2:length(x)) {
    x_out[i] <- abs(x[i]-x[i-1])
  }
  
  x_out[1] <- median(x_out, na.rm=T)
  
  return(x_out)
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
    dplyr::rename(my_period = time_scale) %>%
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
      XXX_WCDC = if(length(q)<=14) {NA} else{
        sum((make_DC(q)/XXX_med)^2, na.rm=T)/XXX_n},
      XXX_CDC =  if(length(q)<=14) {NA} else{
        sum((make_DC(q)/XXX_med), na.rm=T)/XXX_n},
      XXX_FS = abs(mu_over_q(q, .85)- XXX_med) / XXX_med,
      XXX_DS = abs(mu_under_q(q, .15)- XXX_med) / XXX_med) %>%
    pivot_longer(cols=contains('XXX'), names_to='x', values_to='y') %>%
    filter(!is.na(y)) %>%
    filter(!is.infinite(y)) %>%
    mutate(x=str_replace_all(x, 'XXX', time_scale)) %>%
    separate(x, into=c('scale', 'stat')) %>%
    pivot_wider(names_from=stat, values_from=y) %>%
    dplyr::rename({{time_scale}} := my_period)
  
  # ---------------------------------------------------------------------------
  # Pulse Percentages
  PP <- my_stats %>%
    dplyr::rename(my_period = time_scale) %>%
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
    dplyr::rename({{time_scale}} := my_period)
  
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

# storm_stats
storm_stats <- function(df, start, stop, ID) {
  
  # ---------------------------------------------------------------------------
  # Stats
  my_stats <- df %>%
    filter(collection_date >= start &
             collection_date <= stop) %>%
    group_by(site_code) %>%
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
      XXX_WCDC = if(length(q)<=14) {NA} else{
        sum((make_DC(q)/XXX_med)^2, na.rm=T)/XXX_n},
      XXX_CDC =  if(length(q)<=14) {NA} else{
        sum((make_DC(q)/XXX_med), na.rm=T)/XXX_n},
      XXX_FS = abs(mu_over_q(q, .85)- XXX_med) / XXX_med,
      XXX_DS = abs(mu_under_q(q, .15)- XXX_med) / XXX_med) %>%
    pivot_longer(cols=contains('XXX'), names_to='x', values_to='y') %>%
    filter(!is.na(y)) %>%
    filter(!is.infinite(y)) %>%
    mutate(x=str_replace_all(x, 'XXX', ID)) %>%
    separate(x, into=c('scale', 'stat')) %>%
    pivot_wider(names_from=stat, values_from=y)
  
  # ---------------------------------------------------------------------------
  # Pulse Percentages
  PP <- my_stats %>%
    select(site_code, q25, med) %>%
    right_join(d_flow %>%
                 filter(collection_date >= start &
                          collection_date <= stop)) %>%
    mutate(LF = ifelse(q < q25, 1, 0),
           HF3 = ifelse(q > 3*med, 1, 0),
           HF7 = ifelse(q > 7*med, 1, 0),
           HF15 = ifelse(q > 15*med, 1, 0)) %>%
    group_by(site_code) %>%
    dplyr::summarize(LFPP = sum(LF)/length(LF),
                     HFPP3 = sum(HF3)/length(LF),
                     HFPP7 = sum(HF7)/length(LF),
                     HFPP15 = sum(HF15)/length(LF))
  # ---------------------------------------------------------------------------
  # merge
  
  my_output <- left_join(my_stats, PP)
  
  return(my_output)
}
return_time = function(my_data, flood_max, base_level, start, end) {
  # filter a data window
  d_window <- my_data %>%
    filter(collection_date >= start &
             collection_date <= end)
  
  # find day of max flood
  d_f_max_date <- d_window %>%
    dplyr::rename(my_max = {{flood_max}}) %>%
    mutate(flood_max_date = ifelse(q==my_max, collection_date, NA)) %>%
    filter(!is.na(flood_max_date)) %>%
    mutate(flood_max_date = as_date(flood_max_date)) %>%
    select(site_code, flood_max_date) 
  
  # Find minimum date where flow is at or below base level condition
  d_return_date <- d_f_max_date %>%
    right_join(my_data) %>%
    filter(collection_date >= flood_max_date) %>%
    dplyr::rename(my_base = {{base_level}}) %>%
    na.omit() %>%
    mutate(is_base_level = ifelse(q<=my_base, collection_date, NA)) %>%
    group_by(site_code) %>%
    dplyr::summarize(return_date = min(is_base_level, na.rm=T)) %>%
    mutate(return_date = as_date(return_date))
  
  # Calculate return interval (Flood Duration)
  my_output <- d_f_max_date %>%
    left_join(d_return_date) %>%
    mutate(flood_duration = as.numeric(return_date - flood_max_date)) %>%
    left_join(d_window %>%
                select(site_code, ye_med, {{base_level}}, {{flood_max}})%>%
                unique()) %>%
    dplyr::rename(my_max = {{flood_max}}) %>%
    mutate(flood_ratio = my_max / ye_med) %>%
    select(-my_max) %>%
    left_join(d_window %>%
                select(site_code,{{flood_max}})%>%
                unique())
  
  return(my_output)
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
fstats_yseas <- make_stats(d_flow, "yseas")
fstats_seas <- make_stats(d_flow, "seas")%>%
  mutate(seas = fct_relevel(seas, c("cool_dry", "cool_wet", 
                                    "hot_dry", "hot_wet")))

# year-weeks requires segmented formulation to avoid overloading
fywk1 <- make_stats(d_flow%>%filter(year<1995), "yweek")
fywk2 <- make_stats(d_flow%>%filter(year>=1995 & year<2000), "yweek")
fywk3 <- make_stats(d_flow%>%filter(year>=2000 & year<2005), "yweek")
fywk4 <- make_stats(d_flow%>%filter(year>=2010 & year<2015), "yweek")
fywk5 <- make_stats(d_flow%>%filter(year>=2015 & year<2020), "yweek")
fywk6 <- make_stats(d_flow%>%filter(year>=2020), "yweek")
fstats_yweek <- full_join(fywk1, fywk2) %>%
  full_join(fywk3) %>%
  full_join(fywk4) %>%
  full_join(fywk5) %>%
  full_join(fywk6)
rm(fywk1, fywk2, fywk3, fywk4, fywk5, fywk6)

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
  left_join(fstats_yseas %>% combine_scale_stats()) %>%
  mutate(time_period = case_when(
    year %in% c(1990:1994) ~ '1990-1994',
    year %in% c(1995:1999) ~ '1995-1999',
    year %in% c(2000:2004) ~ '2000-2004',
    year %in% c(2005:2010) ~ '2005-2010',
    year %in% c(2010:2014) ~ '2010-2014',
    year %in% c(2015:2019) ~ '2015-2019',
    year %in% c(2020:2025) ~ '2020-2025'  ))

# Storm stats
hurricane_start <- ymd('2017-08-15')
hurricane_end <- ymd('2017-09-15')
regionflood_start <- ymd('2018-06-15')
regionflood_end <- ymd('2018-07-25')
fstats_hur <- storm_stats(d_mega, hurricane_start, hurricane_end, 'hu') 
fstats_rfl <- storm_stats(d_mega, regionflood_start, regionflood_end, 'rf')
# merge storms
d_mega <- d_mega %>% 
  left_join(fstats_hur%>% combine_scale_stats()) %>%
  left_join(fstats_rfl%>% combine_scale_stats())
# flood duration
RI_hur <- return_time(my_data=d_mega, 
flood_max='hu_max', 
base_level='mo_med',
start=hurricane_start,
end=hurricane_end)
colnames(RI_hur) <- c(colnames(RI_hur)[1], 
                      paste('hu_',colnames(RI_hur)[-1],sep=''))
RI_rfl <- return_time(my_data=d_mega, 
                      flood_max='rf_max', 
                      base_level='mo_med',
                      start=regionflood_start,
                      end=regionflood_end)
colnames(RI_rfl) <- c(colnames(RI_rfl)[1], 
                      paste('rf_',colnames(RI_rfl)[-1],sep=''))
# merge
d_mega <- d_mega %>% 
  left_join(RI_hur) %>%
  left_join(RI_rfl)
  
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
    xlab('Year')
}

regional_discharge_time_series <- d_mega %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_med, na.rm=T)) %>%
  ts_base() +
  geom_smooth(method='loess', se=F, span=.2, color='red3')

regional_discharge_time_series_facet <- d_mega %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_mu, na.rm=T)) %>%
  mutate(time_period = ifelse(collection_date <= ymd('1997-01-01'), '1990-1997',
                       ifelse(collection_date > ymd('1997-01-01') &
                              collection_date <= ymd('2004-01-01'), '1997-2004',
                      ifelse(collection_date >= ymd('2004-01-01') &
                             collection_date < ymd('2011-01-01'), '2004-2011', '2011-2021')))) %>%
  ts_base() +
  facet_wrap(~time_period, scales='free_x', ncol=1) +
  geom_smooth(method='loess', se=F, span=.001, color='red3')
  
# -----------------------------------------------------------------------------
# Visualization: Histogram
# -----------------------------------------------------------------------------
hist_base <- function(my_data) {
  my_data %>%
    ggplot(aes(log(q), fill=site_code, color = site_code)) +
    geom_density(alpha=.1) +
    scale_color_manual(values=my_colors) +
    scale_fill_manual(values=my_colors) +
    dark_theme_grey() +
    ylab(element_blank()) +
    xlab('Log Discharge')
}

flow_density <- d_flow %>%
  hist_base()

flow_density_sites <- d_flow %>%
  hist_base() +
  facet_wrap(~site_code, ncol=3)

# -----------------------------------------------------------------------------
# Visualization: Base Flows
# -----------------------------------------------------------------------------
box_yr_base <- function(my_data) {
  my_data %>%
    ggplot(aes(site_code, val, color=annualrain, fill=annualrain)) +
    geom_boxplot(alpha=.1) +
    geom_boxplot(fill=NA) +
    #geom_point(size=5, shape = 21, alpha=.1) +
    #geom_point(size=5, shape = 21, fill=NA) +
    paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
    paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction=-1) +
    dark_theme_grey(base_size=12) +
    scale_y_log10() +
    xlab(element_blank()) +
    ylab('Discharge (l/s)') +
    theme(legend.position = 'none')
  
}

base_flows <- d_mega %>%
  select(annualrain, ye_med, site_code) %>%
  unique() %>%
  dplyr::rename(val = ye_med) %>%
  box_yr_base()

# -----------------------------------------------------------------------------
# Visualization: Season
# -----------------------------------------------------------------------------

# weekly
seas_weekly <- d_mega %>%
  select(year, week, we_mu, al_mu, site_code, annualrain) %>%
  unique() %>%
  ggplot(aes(week, y=we_mu/al_mu, color=site_code)) +
  geom_smooth(method = "loess", se=F, span=.4, alpha=.3) +
  dark_theme_gray() +
  ylab('weekly_avg / yr_avg') +
  scale_x_continuous(breaks = seq(from=1, to=53, by=4)) +
  scale_color_manual(values=my_colors) +
  scale_y_log10() 

# monthly
seas_monthly <- d_mega %>%
  filter(year>=2000) %>%
  select(year, ymonth, month, ym_mu, ye_med, site_code, annualrain) %>%
  unique() %>%
  ggplot(aes(month, y=ym_mu/ye_med, color=site_code)) +
  geom_smooth(method = "loess", se=F, span=.1, alpha=.3) +
  geom_point(alpha=.1) +
  dark_theme_gray() +
  scale_x_continuous(breaks = seq(from=1, to=12, by=1)) +
  ylab('monthly_mu / yr_med') +
  scale_color_manual(values=my_colors) +
  scale_y_log10()

# monthly_faceted
annual_variation_in_month <- d_mega %>%
  filter(year>=2000) %>%
  select(year, ymonth, month, ym_mu, ye_med, site_code, annualrain) %>%
  unique() %>%
  ggplot(aes(year, y=ym_mu/ye_med, color=site_code)) +
  geom_smooth(method = "loess", se=F, color='purple', span=.1) +
  geom_point(alpha=.3) +
  dark_theme_gray() +
  ylab('weekly_avg / yr_avg') +
  #scale_x_continuous(breaks = seq(from=1, to=60, by=3)) +
  scale_color_manual(values=my_colors) +
  scale_y_log10() +
  facet_wrap(~month, scales='free_x', ncol=2)

# Consistent troughs in February and August
# consistent peaks in May and September

# -----------------------------------------------------------------------------
# Visualization: Variability
# -----------------------------------------------------------------------------
# 1990-2020: measures of intra-annual variation
q_variability <- d_mega %>%
  select(annualrain, site_code, ye_FS, ye_DS, mo_rsd, mo_WCDC, ye_LFPP,
         ye_HFPP3, ye_HFPP7, ye_HFPP15) %>%
  pivot_longer(cols=ye_FS:ye_HFPP15, names_to='var', values_to = 'val') %>%
  box_yr_base() +
  ylab(element_blank()) +
  facet_wrap(~var, scales='free')
  
# -----------------------------------------------------------------------------
# Visualization Hurricane Harvey
# -----------------------------------------------------------------------------
# 2016-2019 regional hydrograph
hu_context <- d_mega %>%
  filter(year %in% 2016:2019) %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_med, na.rm=T)) %>%
  left_join(d_flow %>% select(collection_date, nday, year) %>% unique()) %>%
  select(-site_code) %>%
  unique() %>%
  ggplot(aes(nday, y=q_region)) +
  geom_vline(aes(xintercept=271), alpha = .5, linetype=2) +
  geom_point(alpha=.4, color='skyblue') +
  dark_theme_gray(base_size=14) +
  scale_color_manual(values=my_colors) +
  ylab('Regional Discharge Ratio (daily mean / annual mean)') +
  scale_y_log10() +
  xlab('Day') +
  facet_wrap(~year) +
  geom_smooth(method='loess', se=F, span=.1, color='red3')
  
# 2017 hydrographs faceted by site
hu_sites <- d_mega %>%
  filter(year == 2017) %>%
  filter(nday %% 7 == 0) %>%
  select(site_code, annualrain, yw_mu, yw_med, collection_date) %>%
  unique() %>%
  ggplot(aes(collection_date, y=yw_mu/yw_med, 
             fill=annualrain, color = annualrain)) +
  geom_vline(aes(xintercept=ymd('2017-08-27')), alpha=.5, linetype=2) +
  #geom_point(size=3, shape=21, alpha=.3) +
  #geom_point(size=3, shape=21, fill=NA) +
  dark_theme_gray(base_size=14) +
  scale_x_date(date_breaks = '1 month', date_labels = "%m") +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction=-1) +
  ylab('Discharge Ratio (Weekly Mean / Annual Median)') +
  scale_y_log10() +
  xlab('Year') +
  geom_smooth(method='loess', se=F, span=.1, linewidth=.4, alpha=.7) +
  facet_wrap(~site_code) +
  theme(legend.position='none')

# Visualize hurricane duration
hu_duration <- d_mega %>%
  filter(year == 2017) %>%
  filter(collection_date > hurricane_start) %>%
  ggplot(aes(collection_date, y=q, 
             fill=annualrain, color = annualrain)) +
  geom_vline(aes(xintercept=hu_flood_max_date), alpha=.5, linetype=2) +
  geom_vline(aes(xintercept=hu_return_date), alpha=.5, linetype=2) +
  geom_hline(data = filter(d_mega, collection_date==hurricane_start),
             aes(yintercept=al_med), alpha=.5, linetype=2) +
  geom_line() +
  dark_theme_gray(base_size=14) +
  scale_x_date(date_breaks = '1 month', date_labels = "%m") +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction=-1) +
  ylab('Discharge (l/s)') +
  scale_y_log10() +
  xlab('Month') +
  facet_wrap(~site_code) +
  theme(legend.position='none')

# -----------------------------------------------------------------------------
# Visualization June 2018 Regional Flood Event
# -----------------------------------------------------------------------------
# 2017-2020 regional hydrograph
rf_context <- d_mega %>%
  filter(year %in% 2017:2020) %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_med, na.rm=T)) %>%
  left_join(d_flow %>% select(collection_date, nday, year) %>% unique()) %>%
  select(-site_code) %>%
  unique() %>%
  ggplot(aes(nday, y=q_region)) +
  geom_vline(aes(xintercept=200), alpha=.5, linetype=2) +
geom_point(alpha=.4, color='skyblue') +
  dark_theme_gray(base_size=14) +
  scale_color_manual(values=my_colors) +
  ylab('Regional Discharge Ratio (daily mean / annual mean)') +
  scale_y_log10() +
  xlab('Day') +
  facet_wrap(~year) +
  geom_smooth(method='loess', se=F, span=.1, color='red3')

# 2018 sites faceted
rf_sites <- d_mega %>%
  filter(year == 2018) %>%
  ggplot(aes(collection_date, y=yw_mu/yw_med, 
             fill=annualrain, color = annualrain)) +
  geom_vline(aes(xintercept=rf_flood_max_date), alpha=.5, linetype=2) +
  dark_theme_gray(base_size=14) +
  scale_x_date(date_breaks = '1 month', date_labels = "%m") +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction=-1) +
  ylab('Discharge Ratio (Weekly Mean / Annual Median)') +
  scale_y_log10() +
  xlab('Month') +
  geom_smooth(method='loess', se=F, span=.1, linewidth=.4, alpha=.7) +
  facet_wrap(~site_code) +
  theme(legend.position='none')

# Visualize storm duration
rf_duration <- d_mega %>%
  filter(year == 2018) %>%
  filter(collection_date > regionflood_start) %>%
  select(site_code, annualrain, q, mo_med, 
         collection_date, rf_return_date,rf_flood_max_date) %>%
  unique() %>%
  ggplot(aes(collection_date, y=q, 
             fill=annualrain, color = annualrain)) +
  geom_vline(aes(xintercept=rf_flood_max_date), alpha=.5, linetype=2) +
  geom_vline(aes(xintercept=rf_return_date), alpha=.5, linetype=2) +
  geom_hline(data = filter(d_mega, collection_date==regionflood_start),
             aes(yintercept=al_med), alpha=.5, linetype=2) +
  geom_line() +
  dark_theme_gray(base_size=14) +
  scale_x_date(date_breaks = '1 month', date_labels = "%m") +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction=-1) +
  ylab('Discharge (l/s)') +
  scale_y_log10() +
  xlab('Month') +
  #geom_smooth(method='loess', se=F, span=.1, linewidth=.4, alpha=.7) +
  facet_wrap(~site_code) +
  theme(legend.position='none')

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
# export megaframe
write_csv(d_mega, 'Data/daily_discharge_and_storm_stats_1990_2020.csv')

# Site Summary Table
site_summary_table <- fstats_alltime %>% 
  left_join(RI_hur%>%select(-hu_mo_med)%>%unique()) %>%
  left_join(RI_rfl%>%select(-rf_mo_med)%>%unique())
# export summarytable
write_csv(site_summary_table, 'Data/storm_stats_summary.csv')

ggsave('Figures/base_flows.pdf',
       plot = base_flows,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/flow_density.pdf',
       plot = flow_density,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/flow_density_sites.pdf',
       plot = flow_density_sites,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/hu_duration.pdf',
       plot = hu_duration,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/hu_sites.pdf',
       plot = hu_sites,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/q_variability.pdf',
       plot = q_variability,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/regional_discharge_time_series.pdf',
       plot = regional_discharge_time_series,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/regional_discharge_time_series_facet.pdf',
       plot = regional_discharge_time_series_facet,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/rf_context.pdf',
       plot = rf_context,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/rf_duration.pdf',
       plot = rf_duration,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/rf_sites.pdf',
       plot = rf_sites,
       width = 9,
       height = 9,
       units = c("in"))

# -----------------------------------------------------------------------------

# End daily_discharge_and_storm_stats_1990_2020











