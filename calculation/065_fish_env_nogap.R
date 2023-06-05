# fish_env_nogap
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
#------------------------------------------------------------------------------
# merge
#------------------------------------------------------------------------------

df <- d_len %>%
  left_join(d_den) %>%
  left_join(d_bio) %>%
  left_join(d_div) %>%
  left_join(d_ste) %>%
  left_join(d_lte)

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(df, 'Data/fish_env_nogap.csv')


# end fish_env_nogap