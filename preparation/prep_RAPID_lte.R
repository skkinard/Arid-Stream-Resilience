# prep_RAPID_lte
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
# Set Workspace
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

# load packages and helper-functions
source('preparation/merge_toolkit.R')
#------------------------------------------------------------------------------
# RAPID long term environmental data
#------------------------------------------------------------------------------
sites_rapid <- read_csv("source_data/TERRG-Drive/RAPID_FINAL_DATA/Site_Data.csv") %>%
  mutate(Site_Code = substr(Site, 1,2)) %>%
  dplyr::rename(Site_name = FullName) %>%
  dplyr::select(Site_Code, Site_name, STAID, Lat, Lon, AnnualRain, everything()) %>%
  dplyr::select(!Site)

# Fix erroneous basin size
sites_rapid <- sites_rapid %>% select( ! BasinSize) %>%
  left_join(read_csv(
    'source_data/TERRG-Drive/RAPID_FINAL_DATA/GAGES_RAPID_Basin_Comparison_Table.csv') %>%
              dplyr::rename(Site_Code = RAPID_Site, 
                            BasinSize = gage_sqkm) %>%
              select(Site_Code, BasinSize) )

# Substitute nearest site to fill missing data
sites_rapid <- sites_rapid %>% 
  arrange(AnnualRain) %>% 
  fill(everything(), .direction = "downup")

#------------------------------------------------------------------------------
# Export Clean data
#------------------------------------------------------------------------------
write_csv(sites_rapid, 'clean_data/RAPID_lte.csv')
#------------------------------------------------------------------------------
# End prep_RAPID_lte