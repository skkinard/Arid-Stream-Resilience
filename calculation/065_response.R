# response_figures
# Sean Kinard
# 1-16-2023

#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------

# Clean Workspace
rm(list=ls())

# set working directory:
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience")

# load packages
library(tidyverse)
library(lubridate)

#------------------------------------------------------------------------------
# Tidy data
#------------------------------------------------------------------------------

# length
d_len <- read_csv('Data/forklength_all_megaframe.csv') %>%
              group_by(site_code, collection_date) %>%
              dplyr::summarize(length_mu = mean(lengthmm))

# density
d_den <- read_csv('Data/fish_d100.csv') %>%
              group_by(site_code, collection_date) %>%
              dplyr::summarize(d100_mu = mean(d100)) 

# biomass
d_bio <- read_csv('Data/fish_biomass_megaframe.csv') %>%
              group_by(site_code, collection_date) %>%
              dplyr::summarize(biomass_mu = mean(biomass)) 

# diversity
d_div <- read_csv('Data/RTC_fish_diversity.csv') %>%
  select( ! c(richness, simpson))

# my environmental variables
my_evars <- c('annualrain', 'flsh', 'lfpp', 'hfpp3',
              'q_2w_max', 'depth_mx', 'silt, gravel',
              'total_algae', 'canopy_density_mid',
              'conductivity', 'do_mg_l',
              'nh4_n', 'no3n', 'ortho_p')

# short_term environment
d_ste <- read_csv('Data/site_daily_nogap.csv') %>%
  select(site_code, collection_date, any_of(my_evars)) %>%
  mutate(cd=collection_date) %>%
  separate(cd, into = c('year', 'month', 'day')) %>%
  mutate(year = as.numeric(year),
         month = as.numeric(month),
         day = as.numeric(day)) %>%
  mutate(qtr=ifelse(month<4, 'Q1',
             ifelse(month>=4 & month<7, 'Q2',
             ifelse(month>=7 & month <10, 'Q3', 
             ifelse(month>=10 & month <13, 'Q4', NA)))))

# long-term environment
d_lte <- read_csv('Data/site_longterm_nogap.csv') %>%
  select(site_code, any_of(my_evars))

# merge
df <- d_len %>%
  left_join(d_den) %>%
  left_join(d_bio) %>%
  left_join(d_div) %>%
  left_join(d_ste) %>%
  left_join(d_lte)

#------------------------------------------------------------------------------
# Calculate baseline
#------------------------------------------------------------------------------
my_bio_vars <- c('lengthmm', 'd100', 'biomass', 'shannon')

# 2020 Annual Average (A20)
A20 <- df %>%
  filter(collection_date > ymd('2020-01-01')) %>%
  group_by(site_code) %>%
  dplyr::summarize(length_A20_mu=mean(length_mu, na.rm=T),
                   length_A20_sd=sd(length_mu, na.rm=T),
                   d100_A20_mu=mean(d100_mu, na.rm=T),
                   d100_A20_sd=sd(d100_mu, na.rm=T),
                   biomass_A20_mu=mean(biomass_mu, na.rm=T),
                   biomass_A20_sd=sd(biomass_mu, na.rm=T),
                   shannon_A20_mu=mean(shannon, na.rm=T),
                   shannon_A20_sd=sd(shannon, na.rm=T))

# 2020 Seasonal Average (S20)
S20 <- df %>%
  filter(collection_date > ymd('2020-01-01')) %>%
  group_by(site_code, qtr) %>%
  dplyr::summarize(length_S20_mu=mean(length_mu, na.rm=T),
                   d100_S20_mu=mean(d100_mu, na.rm=T),
                   biomass_S20_mu=mean(biomass_mu, na.rm=T),
                   shannon_S20_mu=mean(shannon, na.rm=T) ) %>%
full_join(A20 %>% filter(site_code == 'PD') %>% # PD S20 Q4 missing substitute A20
            select(site_code, contains('_mu')) %>%
            mutate(qtr = 'Q4') %>%
            dplyr::rename(length_S20_mu = length_A20_mu,
                   d100_S20_mu = d100_A20_mu,
                   biomass_S20_mu = biomass_A20_mu,
                   shannon_S20_mu = shannon_A20_mu) )

# Baseline = 2020_seasonal_mean, 2020_annual_standard deviation,
df <- A20 %>%
  select(site_code, contains('_sd')) %>%
  right_join(S20) %>%
  right_join(df) %>% 
  select(site_code, collection_date, qtr, everything()) %>%
  arrange(site_code, collection_date)

#------------------------------------------------------------------------------
# Calculate Response Ratio
#------------------------------------------------------------------------------


# Realized Disturbance (RR response ratio)
#  (value-baseline)/baseline 
make_RR <- function(x, baseline_x) { 
  (x-baseline_x) / baseline_x }

df %>%
  mutate(RR_length_com = make_RR(length_mu_com, length_com_S20_mu),
         RR_length_spp = make_RR(length_mu_spp, length_com_S20_mu),
         RR_d100_com = make_RR(d100_mu_com, d100_com_S20_mu),
         RR_d100_spp = make_RR(d100_mu_spp, d100_com_S20_mu),
         RR_biomass_com = make_RR(biomass_mu_com, biomass_com_S20_mu),
         RR_biomass_spp = make_RR(biomass_mu_spp, biomass_com_S20_mu),
         RR_shannon_com = make_RR(shannon, shannon_S20_mu))

#------------------------------------------------------------------------------
# Calculate Recovery interval (RI in days)
#------------------------------------------------------------------------------
# compare value to baseline with t-test, first match = recovery
# or time to 0 LRR

# function determines status relative to baseline
is_normal <- function(x, y_mu, y_sd, y_n) {
  y_margin <- qt(0.95, df=y_n-1)*y_sd / sqrt(y_n)
  lower_CL <- y_mu - y_margin
  upper_CL <- y_mu + y_margin
  
  if(any(x < lower_CL)) { return('under')  }
  else if(any(x > upper_CL)) { return('over')   }
  else if(any(x > lower_CL & x < upper_CL)) { return('ordinary')  }  
  else {return(NA)}
  }


df %>%
  select(site_code, collection_date, 
         d100_mu_com, d100_com_S20_mu, d100_com_A20_sd) %>%
  arrange(site_code, collection_date) %>% 
  filter(! is.na(d100_mu_com)) %>% 
  filter(! is.na(d100_com_S20_mu)) %>% 
  filter(! is.na(d100_com_A20_sd)) %>% 
  unique() %>%
  mutate(y_margin = qt(0.95, df=4-1)*d100_com_A20_sd / sqrt(4),
         LCL = d100_com_S20_mu - y_margin,
         UCL = d100_com_S20_mu + y_margin,
         RS_d100_com = ifelse(d100_mu_com < LCL, 'under',
                       ifelse(d100_mu_com > UCL, 'over',
                       ifelse(d100_mu_com > LCL &
                                d100_mu_com < UCL, 'ordinary', NA))) )
  
  








# Recover Status (RS)
test <- df %>%
  mutate(RS_length_com = 
           is_normal(length_mu_com, length_com_S20_mu, length_com_A20_sd, 4),
         RS_length_spp = 
           is_normal(length_mu_spp, length_spp_S20_mu, length_spp_A20_sd, 4),
         RS_d100_com = 
           is_normal(d100_mu_com, d100_com_S20_mu, d100_com_A20_sd, 4),
         RS_d100_spp = 
           is_normal(d100_mu_spp, d100_spp_S20_mu, d100_spp_A20_sd, 4),
         RS_biomass_com = 
           is_normal(biomass_mu_com, biomass_com_S20_mu, biomass_com_A20_sd, 4),
         RS_biomass_spp = 
           is_normal(biomass_mu_spp, biomass_spp_S20_mu, biomass_spp_A20_sd, 4),
         RS_shannon_com = 
           is_normal(shannon, shannon_S20_mu, shannon_A20_sd, 4) )

# check
test %>% arrange(site_code, collection_date, species) %>%
  select(site_code, collection_date, species, RS_d100_com, 
         d100_mu_com, d100_com_S20_mu, d100_com_A20_sd) %>%
  print(n=100)

is_normal(4.93, 5.54, 4.21, 4)

is_normal(df$d100_mu_com[1],
          df$d100_com_S20_mu[1],
          df$d100_com_A20_sd[1],
          4)



# Recovery Slope
# LRR / RI