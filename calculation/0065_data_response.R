# data_response
# Sean Kinard
# 2023-02-15

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
library(ggpubr)

#------------------------------------------------------------------------------
# Tidy data
#------------------------------------------------------------------------------

# length (mm)
d_length <- read_csv('Data/fish_Length_stats.csv') 

# density (individuals/100m2)
d_density <- read_csv('Data/fish_d100.csv')

# biomass (g/100m2)
d_biomass <- read_csv('Data/fish_B100.csv') 

# diversity (richness, shannon, simpson)
d_diversity <- read_csv('Data/RTC_fish_diversity.csv')

# Hurricane flood information



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

# Baseline = 2020_seasonal_mean, 2020_annual_standard deviation,
df <- A20 %>%
  right_join(df) %>% 
  select(site_code, collection_date, qtr, everything()) %>%
  arrange(site_code, collection_date) %>%
  filter(collection_date < ymd('2020-01-01'))