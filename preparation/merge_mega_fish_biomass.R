# merge_mega_fish_biomass
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

source('preparation/merge_toolkit.R') # load packages and helper-functions

fish_terrg <- read_csv('clean_data/terrg_fish.csv')
mega_length <- read_csv('clean_data/mega_fish_length.csv')
production <- read_csv('clean_data/terrg_fish_production.csv')

#------------------------------------------------------------------------------
# TERRG production merge
#------------------------------------------------------------------------------
# Extract minimum terrg taxonomic info
extract_tax <- fish_terrg %>% 
  dplyr::rename(genus_species = GenusSpecies,
                commonname = CommonName) %>%
  dplyr::select(commonname, genus_species) %>%
  mutate(gs=genus_species) %>%
  separate(gs, into = c('genus', 'species')) %>%
  distinct()

# Extract merged environmental data
fish_vars <- c( "commonname", "genus_species", "lengthmm", "method", "order", 
                "family", "genus", "species", "lowest_taxon", "comments" )

extract_env <- mega_length %>% 
  dplyr::select(!any_of(fish_vars)) %>% 
  distinct()

# merge production to extracts
d_production <- production %>%
  left_join(extract_tax) %>%
  select(! project) %>%
  arrange(site_code, collection_date) %>%
  left_join(extract_env)

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(d_production, 'clean_data/mega_fish_biomass.csv')

#------------------------------------------------------------------------------
# End merge_mega_fish_biomass


