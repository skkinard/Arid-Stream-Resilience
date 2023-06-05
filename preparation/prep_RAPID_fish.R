# prep_RAPID_fish
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
# RAPID Fish
#------------------------------------------------------------------------------
fish_rapid <- read_csv("source_data/Strickland-data/rapid fish.csv") %>%
  mutate(Site_Code = substr(Site, 1,2)) %>%
  separate(Date, into = c("Month", "Day", "Year"), sep="/") %>%
  mutate(Year = paste('20',Year,sep='')) %>%
  unite("Collection_Date", c("Year", "Month", "Day"), sep='-') %>%
  mutate(Collection_Date = ymd(Collection_Date)) %>%
  create_site_period() %>%
  dplyr::select(Site_Code, Collection_Date, Collection_Period, site_period, CommonName, GenusSpecies, Lengthmm, Method) %>%
  mutate(GenusSpecies = str_replace_all(GenusSpecies, '  ', ' '))

#------------------------------------------------------------------------------
# Convert RAPID Forklength to Total Length
#------------------------------------------------------------------------------
# load conversion data
FL_to_TL_data <- read_csv('source_data/Strickland-data/FishSpeciesRED.csv') %>%
  select(CommonName, LLa, LLb)

# Convert Length data
fish_rapid <- fish_rapid %>%
  left_join(FL_to_TL_data) %>%
  mutate(Lengthmm = LLa + LLb*Lengthmm) %>% 
  select( ! c(LLa, LLb))

# Error check
filter(fish_rapid, is.na(Lengthmm)) %>% print(n=200)

#------------------------------------------------------------------------------
# Redbreast 'Lepomis auritius' spelling error
#------------------------------------------------------------------------------
fish_rapid %>% arrange(GenusSpecies) %>% pull(GenusSpecies) %>% unique()
# 2 levels "Lepomis auritius"         "Lepomis auritus"

# Fix spelling error
fish_rapid <- fish_rapid %>%
  mutate(GenusSpecies = str_replace_all(GenusSpecies, 
                                        "Lepomis auritius", "Lepomis auritus"))
#------------------------------------------------------------------------------
# Export Clean data
#------------------------------------------------------------------------------
write_csv(fish_rapid, 'clean_data/RAPID_fish.csv')
#------------------------------------------------------------------------------
# End prep_RAPID_fish