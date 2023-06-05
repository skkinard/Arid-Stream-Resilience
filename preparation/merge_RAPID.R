# merge_RAPID
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
# Set Workspace
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

# load packages and helper-functions
source('preparation/merge_toolkit.R')

# load clean data
rapid_fish <- read_csv('clean_data/RAPID_fish.csv')
fish_species <- read_csv('clean_data/fish_species.csv') 
rapid_lte <- read_csv('clean_data/RAPID_lte.csv')
rapid_flow <- read_csv('clean_data/RAPID_flow.csv')
rapid_transect <- read_csv('clean_data/RAPID_transect.csv')

#------------------------------------------------------------------------------
# RAPID Merge datasets
#------------------------------------------------------------------------------
m_rapid <- rapid_fish %>%
  dplyr::rename(genus_species = GenusSpecies) %>%
  left_join(fish_species) %>%
  left_join(rapid_lte) %>%
  left_join(rapid_flow) %>%
  left_join(rapid_transect)

# formatting
m_rapid <- m_rapid %>% r_friendly_colnames() %>% fix_fish_spelling()

#------------------------------------------------------------------------------
# Export Merged dataframe
#------------------------------------------------------------------------------
write_csv(m_rapid, 'clean_data/m_RAPID.csv')

#------------------------------------------------------------------------------
# End merge_RAPID