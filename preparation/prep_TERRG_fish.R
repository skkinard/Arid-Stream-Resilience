# prep_TERRG_fish
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

source('preparation/merge_toolkit.R') # load packages and helper-functions

#------------------------------------------------------------------------------
# Tidy TERRG Fish
#------------------------------------------------------------------------------
fish_terrg <- read_csv("source_data/Strickland-data/terrg_fish_fixed.csv") %>%
  mutate(Site_Code = substr(Site, 1,2)) %>%
  separate(Date, into = c("Month", "Day", "Year"), sep="/") %>%
  mutate(Year = paste('20',Year,sep='')) %>%
  unite("Collection_Date", c("Year", "Month", "Day"), sep='-') %>%
  mutate(Collection_Date = ymd(Collection_Date)) %>%
  create_terrg_period() %>%
  dplyr::select(Site_Code, Collection_Date, Collection_Period, site_period, CommonName, GenusSpecies, Lengthmm, Method)

#------------------------------------------------------------------------------
# Export Clean data
#------------------------------------------------------------------------------
write_csv(fish_terrg, 'clean_data/terrg_fish.csv')

#------------------------------------------------------------------------------
# End prep_TERRG_fish