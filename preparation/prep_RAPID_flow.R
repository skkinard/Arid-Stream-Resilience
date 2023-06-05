# prep_RAPID_flow
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
# Set Workspace
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

# load packages and helper-functions
source('preparation/merge_toolkit.R')

fish_rapid <- read_csv('clean_data/RAPID_fish.csv')

#------------------------------------------------------------------------------
# RAPID Flow data
#------------------------------------------------------------------------------
flow_rapid <-  read_csv("source_data/TERRG-Drive/RAPID_FINAL_DATA/Site_Day_Flow.csv") %>%
  mutate(Site_Code = substr(Site, 1,2)) %>%
  unite(col = Collection_Date, "Year":"Day", sep = "-") %>%
  mutate(Collection_Date = ymd(Collection_Date)) %>%
  dplyr::select(!c(UID, Cycle, Site)) %>%
  dplyr::select(Site_Code, Collection_Date, everything()) 

#  Collection_Dates in fish data that don't match flow data
fish_rapid %>%
  anti_join(flow_rapid) %>%
  dplyr::select(Site_Code, Collection_Date) %>%
  unique()

#  <chr>      <date>         
# 1 AR        2017-09-16     
# 2 GC        2018-11-15     
# 3 TR        2018-02-02

#  Collection_Dates in flow data that don't match fish data
flow_rapid %>%
  dplyr::select(Site_Code, Collection_Date) %>%
  unique() %>%
  filter(Site_Code %in% c('AR', 'GC', 'TR')) %>%
  mutate(id = seq(1:51)) %>%
  pivot_wider(names_from=Site_Code, values_from = Collection_Date) %>% 
  print(n=51)
# AR 2017-09-12
# GC no flow data from november
# TR 2018-02-03

# match flow dates within several days
flow_rapid <- flow_rapid %>%
  mutate(Collection_Date = as.character(Collection_Date)) %>%
  mutate(Collection_Date = ifelse(
    Site_Code =='AR' & Collection_Date == '2017-09-12', '2017-09-16',
    ifelse(Site_Code == 'TR' & Collection_Date == '2018-02-03',
           '2018-02-02', Collection_Date)) ) %>%
  mutate(Collection_Date = ymd(Collection_Date))

#  Collection_Dates in fish data that don't match flow data
fish_rapid %>%
  anti_join(flow_rapid) %>%
  dplyr::select(Site_Code, Collection_Date) %>%
  unique()
# 2 GC        2018-11-15 

# more than 2 site-dates for flow data?
flow_rapid %>%
  group_by(Site_Code, Collection_Date) %>%
  dplyr::summarize(n=length(Q)) %>%
  filter(n>1)
# AR        2018-04-07          2

flow_rapid %>%
  distinct() %>% filter(Site_Code == 'AR' & Collection_Date == ymd('2018-04-07'))

flow_rapid <- flow_rapid %>% distinct()
#------------------------------------------------------------------------------
# Export Clean Data
#------------------------------------------------------------------------------
write_csv(flow_rapid, 'clean_data/RAPID_flow.csv')
#------------------------------------------------------------------------------
# End prep_RAPID_flow