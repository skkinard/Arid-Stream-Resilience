# prep_terrg_transect
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

source('preparation/merge_toolkit.R') # load packages and helper-functions

#------------------------------------------------------------------------------
# TERRG Transect
#------------------------------------------------------------------------------

transect_terrg <- read_csv("source_data/TERRG-Drive/TERRG environmental data/mlevel_fielddata_TERRG.csv") %>%
  mutate(Site_Code = substr(Site, 1,2)) %>%
  unite("Collection_Date", c("Year", "Month", "Day"), sep='-') %>%
  mutate(Collection_Date = ymd(Collection_Date)) %>%
  create_terrg_period() %>%
  dplyr::select(!c(UID3, UID, Site)) %>%
  mutate(UID = paste(Site_Code, Collection_Date, Transect, sep='_')) %>%
  dplyr::select(UID, Site_Code, Collection_Date, Transect, Collection_Period, site_period, Cycle, everything()) %>%
  mutate(Depth.Mx = as.numeric(Depth.Mx)) %>% # fix outliers where depth was recorded in cm instead of m
  mutate(Depth.Mx = ifelse(Depth.Mx > 3, Depth.Mx*0.01, Depth.Mx)) 

# code friendly column headers
colnames(transect_terrg) <- str_to_lower(colnames(transect_terrg))
colnames(transect_terrg) <- str_replace_all(colnames(transect_terrg), '\\.', '_')
colnames(transect_terrg) <- str_replace_all(colnames(transect_terrg), '__', '_')

# convert character vectors to numeric vectors
transect_terrg <- transect_terrg %>%
  mutate(do_sat = as.numeric(do_sat),
         depth_right = as.numeric(depth_right),
         depth_mn = as.numeric(depth_mn),
         depth_mx = as.numeric(depth_mx),
         bank_slope_left = as.numeric(bank_slope_left),
         bank_slope_right = as.numeric(bank_slope_right),
         canopy_total = as.numeric(canopy_total),
         canopy_density_right = as.numeric(canopy_density_right),
         canopy_density_mid = as.numeric(canopy_density_mid),
         canopy_density_left = as.numeric(canopy_density_left) )

# Calculate mean and sd for transect measurements (grouped by site and sampling date)
# n =< 4 for all means and sds
transect_terrg <- transect_terrg %>%
  pivot_longer(cols = bluegreen_cyano : bank_slope_right,
               names_to = "measure",
               values_to = "value") %>%
  group_by(site_period, measure) %>%
  dplyr::summarize(xmean = mean(value, na.rm=T)) %>%
  mutate(measure = paste(measure, 'mean', sep='_')) %>%
  pivot_wider(names_from = measure, values_from = xmean) %>%
  left_join( transect_terrg %>%
               pivot_longer(cols = bluegreen_cyano : bank_slope_right,
                            names_to = "measure",
                            values_to = "value") %>%
               group_by(site_period, measure) %>%
               dplyr::summarize(xsd = sd(value, na.rm=T) ) %>%
               mutate(measure = paste(measure, 'sd', sep = '_')) %>%
               pivot_wider(names_from = measure, values_from = xsd) )

#------------------------------------------------------------------------------
# Fill missing data with nearest
#------------------------------------------------------------------------------
transect_terrg <- transect_terrg %>%
  separate(site_period, into=c('site', 'period'), sep='_')

all_periods <- transect_terrg %>% pull(period) %>% unique()
all_sites <- transect_terrg %>% pull(site) %>% unique()
d_cross <- crossing(all_sites, all_periods) %>% print(n=20) %>%
  rename(site=all_sites, period=all_periods) # all site-period combinations

transect_terrg <- left_join(d_cross, transect_terrg) %>% # expand df with gaps
  group_by(site) %>%
  fill(everything(), .direction = "downup") %>% # fill gaps with row below or above
  unite('site_period', site:period, sep='_')

#------------------------------------------------------------------------------
# Export clean data
#------------------------------------------------------------------------------
write_csv(transect_terrg, 'clean_data/terrg_transect.csv')

#------------------------------------------------------------------------------
# End prep_terrg_transect