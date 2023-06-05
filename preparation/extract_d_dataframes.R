# d_dataframes
# Sean Kinard
# 2023-01-27

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------

# load packages
library(tidyverse)
library(lubridate)

#------------------------------------------------------------------------------
# My Objects
#------------------------------------------------------------------------------

# Helper objects
site_order <- c("TR", "SF", "AR", "MR", "PD", "PL", "GC", "WM", "EM")
my_index <- c('site_code', 'collection_date', 'collection_period')
my_taxonomic <- c('order', 'family', 'genus', 'species', 'lowest_taxon', 'commonname')
my_short_evars <- c('q_2wk_max', 'diatoms', 'green_algae', 'bluegreen_cyano',
                    'depth_mx', 'width', 'gravel', 'silt', 'canopy_density_mid', 
                    'conductivity', 'do_mg_l', 'nh4_n', 'no3n', 'ortho_p')
my_algae_vars <- c('diatoms', 'green_algae', 'bluegreen_cyano')
my_long_evars <- c('basin_size')
my_hydro_vars <- c('med', 'rsd', 'WCDC', 'FS', 'DS', 'LFPP', 'HFPP3', 'HFPP7',
                   'hu_flood_duration', 'hu_flood_ratio', 'hu_ye_med', 
                   'rf_flood_duration', 'rf_flood_ratio', 'rf_ye_med')
my_rain_vars <- c('storm_rain')
my_evars <- c('annualrain', my_short_evars, my_long_evars, my_hydro_vars, my_rain_vars)
my_colors <- c('#F5191CFF', '#E78200FF', '#E8A117FF',
               '#EABB22FF', '#CBC988FF', '#9FC095FF',
               '#6BB699FF', '#3DAAA6FF', '#3B99B1FF')

# Helper functions
add_time_vars <- function(my_data) {
  # Requires columns 'collection_date', and 'site_code'
  
  
  # Add collection_period which allows site-site comparisons by date
  my_data <- read_csv('Data/forklength_all_megaframe.csv') %>%
    select(collection_date, collection_period) %>%
    unique() %>%
    right_join(my_data)
  
  # Add time variables derived from colleciton_date
  my_data <- my_data %>%
    mutate(alltime = 'all') %>%
    mutate(cd2 = collection_date) %>%
    separate(cd2, into = c('year', 'month', 'day')) %>%
    mutate(week=week(collection_date)) %>%
    mutate(year=as.numeric(year),
           month=as.numeric(month),
           day=as.numeric(day)) %>%
    mutate(nday=month*30+day) %>%
    mutate(qtr = case_when(
      month < 4 ~ 'Q1',
      month > 3 & month < 7 ~ 'Q2',
      month > 6 & month < 10 ~ 'Q3', 
      month > 9 ~ 'Q4') ) %>%
    mutate(seas = case_when(
      month %in% c(11,12,1,2,3) ~ 'cool_dry',
      month %in% c(4,5,6) ~ 'cool_wet',
      month %in% c(7,8) ~ 'hot_dry',
      month %in% c(9,10) ~ 'hot_wet')) %>%
    mutate(ymonth = paste(as.character(year), as.character(month), sep='-'), 
           yqtr = paste(as.character(year), as.character(qtr), sep='-'), 
           yweek = paste(as.character(year), as.character(week), sep='-'), 
           yseas = paste(as.character(year), as.character(seas), sep='-'))
  
  
  return(my_data) }

fix_site_order <- function(my_data) {
  my_data <- my_data %>%
    mutate(site_code = fct_relevel(site_code, site_order)) 
  return(my_data) }

add_taxonomic <- function(my_data) {
  
  my_taxonomic <- c('order', 'family', 'genus', 
                    'species', 'lowest_taxon')
  
  # Taxanomic info
  d_spp <- read_csv('Data/forklength_all_megaframe.csv') %>%
    select(any_of(my_taxonomic)) %>%
    unique()
  
  my_data <- my_data %>% left_join(d_spp)
  
  return(my_data)
  
}

#------------------------------------------------------------------------------
# Tidy Bio Data
#------------------------------------------------------------------------------

# helper functions
# Sum or mean by community , taxonomic family, taxonomic species

my_sums <- function(my_data, my_variable) {
  
  df <- my_data %>%
    dplyr::rename(x_var = my_variable)
  
  d_com <- df %>%
    group_by(site_code, collection_date) %>%
    dplyr::summarize(XXX_sum_com = sum(x_var)) %>% 
    ungroup()
  
  d_fam <- df %>%
    group_by(site_code, collection_date, family) %>%
    dplyr::summarize(XXX_sum_fam = sum(x_var)) %>% 
    ungroup()
  
  d_spe <- df %>%
    group_by(site_code, collection_date, lowest_taxon) %>%
    dplyr::summarize(XXX_sum_spe = sum(x_var)) %>%
    left_join(select(my_data, family, lowest_taxon) %>% unique()) %>% 
    ungroup()
  
  my_output <- d_spe %>%
    left_join(d_fam) %>%
    left_join(d_com) %>%
    select(contains('XXX'), site_code, collection_date, family, lowest_taxon)
  
  colnames(my_output) <- str_replace_all(
    colnames(my_output), 'XXX', my_variable)
  
  return(my_output) }

my_means <- function(my_data, my_variable) {
  
  df <- my_data %>%
    dplyr::rename(x_var = my_variable)
  
  d_com <- df %>%
    group_by(site_code, collection_date) %>%
    dplyr::summarize(XXX_mu_com = mean(x_var),
                     XXX_n_com = length(x_var),
                     XXX_sd_com = sd(x_var)) %>% ungroup()
  
  d_fam <- df %>%
    group_by(site_code, collection_date, family) %>%
    dplyr::summarize(XXX_mu_fam = mean(x_var),
                     XXX_n_fam = length(x_var),
                     XXX_sd_fam = sd(x_var)) %>% ungroup()
  
  d_spe <- df %>%
    group_by(site_code, collection_date, lowest_taxon) %>%
    dplyr::summarize(XXX_mu_spe = mean(x_var),
                     XXX_n_spe = length(x_var),
                     XXX_sd_spe = sd(x_var)) %>%
    left_join(select(my_data, family, lowest_taxon) %>% unique())
  
  my_output <- d_spe %>%
    left_join(d_fam) %>%
    left_join(d_com) %>%
    select(contains('XXX'), site_code, collection_date, family, lowest_taxon)
  
  colnames(my_output) <- str_replace_all(
    colnames(my_output), 'XXX', my_variable)
  
  return(my_output) }

# taxonomic
d_spp <- read_csv('Data/forklength_all_megaframe.csv') %>%
  select(any_of(my_taxonomic)) %>%
  unique()

# length
d_length <- read_csv('Data/forklength_all_megaframe.csv') %>%
  select(lengthmm, any_of(my_index), any_of(my_taxonomic)) %>%
  dplyr::rename(TL = lengthmm) %>%
  add_taxonomic() %>%
  my_means(my_variable='TL')
  
# density
d_density <- read_csv('Data/fish_d100.csv') %>%
  select(contains('d100'), any_of(my_index), lowest_taxon)  %>%
  add_taxonomic()

# biomass
d_biomass <- read_csv('Data/fish_biomass_megaframe.csv') %>%
  filter(biomass>0) %>%
  mutate(B100 = biomass*100) %>%
  select(B100, any_of(my_index), any_of(my_taxonomic)) %>%
  add_taxonomic() %>%
  my_sums(my_variable='B100')

# diversity
d_diversity <- read_csv('Data/RTC_fish_diversity.csv') %>%
  add_time_vars()

# merge for bioframe
d_biovars <-
  d_length %>%
  left_join(d_density) %>%
  left_join(d_biomass) %>%
  left_join(d_diversity) %>%
  select(site_code, collection_date, contains('_')) %>%
  mutate(lowest_taxon = str_replace_all(lowest_taxon, ' ','')) %>%
  left_join(d_spp %>%
              mutate(lowest_taxon = str_replace_all(lowest_taxon, ' ','')) )

#------------------------------------------------------------------------------
# Tidy Environmental Data
#------------------------------------------------------------------------------

# load environmental data
d_ste <- read_csv('Data/site_daily_nogap.csv') # short_term environment
d_lte <- read_csv('Data/site_longterm_nogap.csv') # long-term environment
d_hydr <- read_csv('Data/storm_stats_summary.csv') # Hydrological Variables
d_rain <- read_csv('Data/HH_Rain.csv') %>%
  dplyr::rename(storm_rain = `storm.rain`) %>%
  mutate(site_code = substr(SITE, 1,2)) %>%
  select(-SITE) # Storm Rain

# Merge environmental data
d_environment <- d_ste %>%
  left_join(d_lte) %>%
  left_join(d_hydr) %>%
  left_join(d_rain)

# select variables
d_environment <- d_environment %>%
  select(any_of(my_index), any_of(my_evars))

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------

write_csv(d_length, 'Data/d_length.csv')
write_csv(d_density, 'Data/d_density.csv')
write_csv(d_biomass, 'Data/d_biomass.csv')
write_csv(d_diversity, 'Data/d_diversity.csv')
write_csv(d_biovars, 'Data/d_biovars.csv')
write_csv(d_environment, 'Data/d_environment.csv')

#------------------------------------------------------------------------------

# End d_dataframes