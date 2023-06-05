# 002_site_daily
# 1-19-23
# Sean Kinard

# Extract a priori environmental variables
# fill gaps with seasonal means (3x 4 month seasons)
# if no seasonal mean, fill with annual mean

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

# load data
d <- read_csv('Data/fish_biomass_megaframe.csv')
d_meta <- read_csv('Data/fish_rapid_terrg_megaframe_metadata.csv')

#------------------------------------------------------------------------------
# RAPID Daily environmental variables (extraction)
#------------------------------------------------------------------------------

devars <- d %>% 
  select(site_code, collection_date, q, q_2w_mn, q_2w_max, q_2w_min,
         diatoms_mean, green_algae_mean, 
         depth_mx_mean, width_mean, gravel_mean, 
         silt_mean, canopy_density_mid_mean,
         conductivity_mean, do_mg_l_mean, 
         nh4_n_mean, no3n_mean, ortho_p_mean , bluegreen_cyano_mean ) %>%
  distinct() %>%
  mutate(qtr = ifelse(collection_date >= ymd('2017-09-01') &
                      collection_date < ymd('2018-01-01'), '2017-Q4',
                ifelse(collection_date >= ymd('2018-01-01') &
                       collection_date < ymd('2018-04-01'),'2018-Q1',
                ifelse(collection_date >= ymd('2018-04-01') &
                       collection_date < ymd('2018-07-01'), '2018-Q2',
                ifelse(collection_date >= ymd('2018-07-01') &
                       collection_date < ymd('2018-10-01'), '2018-Q3',
                ifelse(collection_date >= ymd('2018-10-01') &
                       collection_date < ymd('2019-01-01'), '2018-Q4',
                ifelse(collection_date >= ymd('2018-01-01') &
                       collection_date < ymd('2020-04-01'),'2020-Q1',
                ifelse(collection_date >= ymd('2020-04-01') &
                       collection_date < ymd('2020-07-01'), '2020-Q2',
                ifelse(collection_date >= ymd('2020-07-01') &
                       collection_date < ymd('2020-10-01'), '2020-Q3',
                ifelse(collection_date >= ymd('2020-10-01') &
                       collection_date < ymd('2021-01-01'), '2020-Q4',
                       'other')))))))))) %>%
  mutate(q = ifelse(is.infinite(q), NA, q),
         q_2w_mn = ifelse(is.infinite(q_2w_mn), NA, q_2w_mn),
         q_2w_max = ifelse(is.infinite(q_2w_max), NA, q_2w_max),
         q_2w_min = ifelse(is.infinite(q_2w_min), NA, q_2w_min))

#------------------------------------------------------------------------------
# Daily environmental variables: Replace -Inf & NA with mean per site*season
#------------------------------------------------------------------------------


# creates mean for site*season (variable_y)
devars <- devars %>%
  group_by(site_code, qtr) %>%
  dplyr::summarize(q_x = mean(q, na.rm=T), 
                   q_2w_mn_x = mean(q_2w_mn, na.rm=T),
                   q_2w_max_x = mean(q_2w_max, na.rm=T),
                   q_2w_min_x = mean(q_2w_min, na.rm=T),
                   diatoms_mean_x = mean(diatoms_mean, na.rm=T),
                   green_algae_mean_x = mean(green_algae_mean, na.rm=T),
                   depth_mx_mean_x = mean(depth_mx_mean, na.rm=T),
                   width_mean_x = mean(width_mean, na.rm=T),
                   gravel_mean_x = mean(gravel_mean, na.rm=T),
                   silt_mean_x = mean(silt_mean, na.rm=T),
                   canopy_density_mid_mean_x = mean(canopy_density_mid_mean, na.rm=T),
                   conductivity_mean_x = mean(conductivity_mean, na.rm=T),
                   do_mg_l_mean_x = mean(do_mg_l_mean, na.rm=T),
                   nh4_n_mean_x = mean(nh4_n_mean, na.rm=T),
                   no3n_mean_x = mean(no3n_mean, na.rm=T),
                   ortho_p_mean_x = mean(ortho_p_mean, na.rm=T),
                   bluegreen_cyano_mean_x = mean(bluegreen_cyano_mean, na.rm=T)
  ) %>%
  right_join(devars)

# creates mean for site (variable_y)
devars <- devars %>%
  group_by(site_code) %>%
  dplyr::summarize(q_y = mean(q, na.rm=T), 
                   q_2w_mn_y = mean(q_2w_mn, na.rm=T),
                   q_2w_max_y = mean(q_2w_max, na.rm=T),
                   q_2w_min_y = mean(q_2w_min, na.rm=T),
                   diatoms_mean_y = mean(diatoms_mean, na.rm=T),
                   green_algae_mean_y = mean(green_algae_mean, na.rm=T),
                   depth_mx_mean_y = mean(depth_mx_mean, na.rm=T),
                   width_mean_y = mean(width_mean, na.rm=T),
                   gravel_mean_y = mean(gravel_mean, na.rm=T),
                   silt_mean_y = mean(silt_mean, na.rm=T),
                   canopy_density_mid_mean_y = mean(canopy_density_mid_mean, na.rm=T),
                   conductivity_mean_y = mean(conductivity_mean, na.rm=T),
                   do_mg_l_mean_y = mean(do_mg_l_mean, na.rm=T),
                   nh4_n_mean_y = mean(nh4_n_mean, na.rm=T),
                   no3n_mean_y = mean(no3n_mean, na.rm=T),
                   ortho_p_mean_y = mean(ortho_p_mean, na.rm=T),
                   bluegreen_cyano_mean_y = mean(bluegreen_cyano_mean, na.rm=T)
  ) %>%
  right_join(devars)


# gapfill NAs with season if possible and where no season, annual mean
devars <- devars %>%
  mutate(q = ifelse(is.na(q_x), q_y,
                    ifelse(is.na(q), q_x, q)),
         q_2w_mn = ifelse(is.na(q_2w_mn_x), q_2w_mn_y,
                          ifelse(is.na(q_2w_mn), q_2w_mn_x, q)),
         q_2w_max = ifelse(is.na(q_2w_max_x), q_2w_max_y,
                           ifelse(is.na(q_2w_max), q_2w_max_x, q_2w_max)),
         q_2w_min = ifelse(is.na(q_2w_min_x), q_2w_min_y,
                           ifelse(is.na(q_2w_min), q_2w_min_x, q_2w_min)),
         diatoms_mean = ifelse(is.na(diatoms_mean_x), diatoms_mean_y,
                               ifelse(is.na(diatoms_mean), diatoms_mean_x, diatoms_mean)),
         green_algae_mean = ifelse(is.na(green_algae_mean_x), 
                                   green_algae_mean_y,
                                   ifelse(is.na(green_algae_mean),
                                          green_algae_mean_x, 
                                          green_algae_mean)),
         depth_mx_mean = ifelse(is.na(depth_mx_mean_x), 
                                depth_mx_mean_y,
                                ifelse(is.na(depth_mx_mean),
                                       depth_mx_mean_x, 
                                       depth_mx_mean)),
         width_mean = ifelse(is.na(width_mean_x), 
                             width_mean_y,
                             ifelse(is.na(width_mean),
                                    width_mean_x, 
                                    width_mean)),
         gravel_mean = ifelse(is.na(gravel_mean_x), 
                              gravel_mean_y,
                              ifelse(is.na(gravel_mean),
                                     gravel_mean_x, 
                                     gravel_mean)),
         silt_mean = ifelse(is.na(silt_mean_x), 
                            silt_mean_y,
                            ifelse(is.na(silt_mean),
                                   silt_mean_x, 
                                   silt_mean)),
         canopy_density_mid_mean = ifelse(is.na(canopy_density_mid_mean_x), 
                                          canopy_density_mid_mean_y,
                                          ifelse(is.na(canopy_density_mid_mean),
                                                 canopy_density_mid_mean_x, 
                                                 canopy_density_mid_mean)),
         conductivity_mean = ifelse(is.na(conductivity_mean_x), 
                                    conductivity_mean_y,
                                    ifelse(is.na(conductivity_mean),
                                           conductivity_mean_x, 
                                           conductivity_mean)),
         do_mg_l_mean = ifelse(is.na(do_mg_l_mean_x), 
                               do_mg_l_mean_y,
                               ifelse(is.na(do_mg_l_mean),
                                      do_mg_l_mean_x, 
                                      do_mg_l_mean)),
         nh4_n_mean = ifelse(is.na(nh4_n_mean_x), 
                             nh4_n_mean_y,
                             ifelse(is.na(nh4_n_mean),
                                    nh4_n_mean_x, 
                                    nh4_n_mean)),
         no3n_mean = ifelse(is.na(no3n_mean_x), 
                            no3n_mean_y,
                            ifelse(is.na(no3n_mean),
                                   no3n_mean_x, 
                                   no3n_mean)),
         ortho_p_mean = ifelse(is.na(ortho_p_mean_x), 
                               ortho_p_mean_y,
                               ifelse(is.na(ortho_p_mean),
                                      ortho_p_mean_x, 
                                      ortho_p_mean)),
         bluegreen_cyano_mean = ifelse(is.na(bluegreen_cyano_mean_x), 
                                       bluegreen_cyano_mean_y,
                                       ifelse(is.na(bluegreen_cyano_mean),
                                              bluegreen_cyano_mean_x, 
                                              bluegreen_cyano_mean)) ) %>%
  select(! contains('_x')) %>%
  select(! contains('_y'))

colnames(devars) <- str_replace_all(colnames(devars),'_mean', '')

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------

# write files
write_csv(devars, '/home/kinard/Documents/Research/Dissertation/02_Resilience/Data/site_daily_nogap.csv')


#------------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------------
rm(d, d_meta)

# End 002_site_daily