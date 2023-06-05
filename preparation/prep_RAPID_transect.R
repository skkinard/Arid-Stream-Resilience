# prep_RAPID_transect
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
  left_join(read_csv('source_data/TERRG-Drive/RAPID_FINAL_DATA/GAGES_RAPID_Basin_Comparison_Table.csv') %>%
              dplyr::rename(Site_Code = RAPID_Site, 
                            BasinSize = gage_sqkm) %>%
              select(Site_Code, BasinSize) )

#------------------------------------------------------------------------------
# RAPID Transect
#------------------------------------------------------------------------------

# load data, create site_period
transect_rapid <- read_csv("source_data/TERRG-Drive/RAPID_FINAL_DATA/Site_Day_TransectData.csv") %>%
  mutate(Site_Code = substr(Site, 1,2)) %>%
  unite("Collection_Date", c("Year", "Month", "Day"), sep='-') %>%
  mutate(Collection_Date = ymd(Collection_Date)) %>%
  create_site_period() %>%
  dplyr::select(!c(UID3, UID, Site)) %>%
  mutate(UID = paste(Site_Code, Collection_Date, Transect, sep='_')) %>%
  dplyr::select(UID, Site_Code, Collection_Date, Transect, Collection_Period, site_period, everything())

# code friendly column headers
colnames(transect_rapid) <- str_to_lower(colnames(transect_rapid))
colnames(transect_rapid) <- str_replace_all(colnames(transect_rapid), '\\.', '_')
colnames(transect_rapid) <- str_replace_all(colnames(transect_rapid), '__', '_')

# Calculate mean and sd for transect measurements (grouped by site and sampling date)
# n =< 4 for all means and sds
transect_rapid <- transect_rapid %>%
  mutate(depth_mx = ifelse(depth_mx > 3, 
                           0.01*depth_mx, 
                           depth_mx),
         depth_mid = ifelse(depth_mid > 3, 
                            0.01*depth_mid,
                            depth_mid) ) %>%
  pivot_longer(cols = bluegreen_cyano : ortho_p,
               names_to = "measure",
               values_to = "value") %>%
  group_by(site_period, measure) %>%
  dplyr::summarize(xmean = mean(value, na.rm=T)) %>%
  mutate(measure = paste(measure, 'mean', sep='_')) %>%
  pivot_wider(names_from = measure, values_from = xmean) %>%
  left_join( transect_rapid %>%
               pivot_longer(cols = bluegreen_cyano : ortho_p,
                            names_to = "measure",
                            values_to = "value") %>%
               group_by(site_period, measure) %>%
               dplyr::summarize(xsd = sd(value, na.rm=T) ) %>%
               mutate(measure = paste(measure, 'sd', sep = '_')) %>%
               pivot_wider(names_from = measure, values_from = xsd) )

#------------------------------------------------------------------------------
# Export Clean data
#------------------------------------------------------------------------------
write_csv(transect_rapid, 'clean_data/RAPID_transect.csv')
#------------------------------------------------------------------------------
# End prep_RAPID_transect