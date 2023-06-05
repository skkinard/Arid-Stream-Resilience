# merge_TERRG_length
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

source('preparation/merge_toolkit.R') # load packages and helper-functions

fish_terrg <- read_csv('clean_data/terrg_fish.csv')
flow_terrg <- read_csv('clean_data/terrg_flow.csv')
transect_terrg <- read_csv('clean_data/terrg_transect.csv')
nutrient_terrg <- read_csv('clean_data/terrg_nutrient.csv')
sites_rapid <- read_csv('clean_data/RAPID_lte.csv')
fish_species <- read_csv('clean_data/fish_species.csv')

#------------------------------------------------------------------------------
# TERRG length date matching
#------------------------------------------------------------------------------
# check date matches
fish_terrg_slim <- fish_terrg %>% 
  dplyr::select(site_period) %>% 
  distinct() # 35 long

# check flow merge to terrg fish
flow_terrg_slim <- flow_terrg %>% 
  dplyr::select(site_period, Q) %>% distinct() # 127 long

left_join(fish_terrg_slim, flow_terrg_slim) #  35
left_join(fish_terrg_slim, flow_terrg_slim) %>% filter(is.na(Q)) # check is good

# check transect merge to terrg fish
transect_terrg_slim <- transect_terrg %>% 
  dplyr::select(site_period, green_algae_mean) %>% distinct() # 88 long
left_join(fish_terrg_slim, transect_terrg_slim) # length is 35
left_join(fish_terrg_slim, transect_terrg_slim) %>% 
  filter(is.na(green_algae_mean)) # check is good

# check nutrients merge to terrg fish
nutrient_terrg_slim <- nutrient_terrg %>% dplyr::select(site_period,
                                                        no3n_mean) %>% 
  distinct() # 131 long
left_join(fish_terrg_slim, nutrient_terrg_slim)
left_join(fish_terrg_slim, nutrient_terrg_slim) %>% 
  filter(is.na(no3n_mean)) # check is good

#------------------------------------------------------------------------------
# merge terrg length data
#------------------------------------------------------------------------------
m_terrg <- fish_terrg %>%
  dplyr::rename(genus_species = GenusSpecies) %>%
  left_join(fish_species) %>%
  left_join(sites_rapid) %>%
  left_join(flow_terrg) %>%
  left_join(nutrient_terrg) %>%
  left_join(transect_terrg)

# formatting
m_terrg <- m_terrg %>% r_friendly_colnames() %>% fix_fish_spelling()

#------------------------------------------------------------------------------
# Export Clean Data
#------------------------------------------------------------------------------
write_csv(m_terrg, 'clean_data/m_terrg_length.csv')

#------------------------------------------------------------------------------
# End merge_TERRG_length