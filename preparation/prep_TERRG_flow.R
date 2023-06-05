# prep_TERRG_flow
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

source('preparation/merge_toolkit.R') # load packages and helper-functions

#------------------------------------------------------------------------------
# TERRG Flow
#------------------------------------------------------------------------------
terrg_flow <- read_csv("source_data/TERRG-Drive/TERRG_DailyFlowData/flow_data_TERRG.csv") %>%
  mutate(Site_Code = substr(Site, 1,2)) %>%
  mutate(Collection_Date = ymd(date)) %>% 
  dplyr::select(Site_Code, Collection_Date, Q, Q.2w.med, Q.2w.mn, Q.2w.max, Q.2w.min ) %>%
  create_terrg_period() %>%
  select(-Site_Code, -Collection_Date)

#------------------------------------------------------------------------------
# Fill missing data with nearest
#------------------------------------------------------------------------------
terrg_flow <- terrg_flow %>%
  separate(site_period, into=c('site', 'period'), sep='_')

all_periods <- terrg_flow %>% pull(period) %>% unique()
all_sites <- terrg_flow %>% pull(site) %>% unique()
d_cross <- crossing(all_sites, all_periods) %>% print(n=20) %>%
  rename(site=all_sites, period=all_periods) # all site-period combinations

terrg_flow <- left_join(d_cross, terrg_flow) %>% # expand df with gaps
  group_by(site) %>%
  fill(everything(), .direction = "downup") %>% # fill gaps with row below or above
  unite('site_period', site:period, sep='_')

#------------------------------------------------------------------------------
# Export Merged dataframe
#------------------------------------------------------------------------------
write_csv(terrg_flow, 'clean_data/terrg_flow.csv')

#------------------------------------------------------------------------------
# End prep_TERRG_flow