# diversity_biomass_calculate
# Sean Kinard
# 1-30-2023

#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
# Clean Workspace
rm(list=ls())

# set working directory:
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience")

# load packages
library(iNEXT)
library(tidyverse)
library(lubridate)

#------------------------------------------------------------------------------
# Tidy Data
#------------------------------------------------------------------------------
# load data
bm <- read_csv('Data/fish_biomass_megaframe.csv')

# format: rows = species,  column = 'site_code collection_date'
df <- bm %>%
  mutate(biomass = ceiling(biomass*100)) %>% # biomass (g per 100m2)
  mutate(UID=paste(site_code, collection_date, sep=' ')) %>%
  select(genus_species, UID, biomass) %>%
  distinct() %>%
  pivot_wider(names_from = UID,
              values_from = biomass,
              values_fill = 0) %>%
  column_to_rownames(var = "genus_species") %>%
  as.data.frame()

#------------------------------------------------------------------------------
# Calculate Hill Numbers
#------------------------------------------------------------------------------
# supported data formats (abundance, incidence_raw, incidence_frequency)
# 'abundance' -> Individual‐based abundance data (datatype="abundance"): Input data for each assemblage/site include species abundances in an empirical sample of n individuals (“reference sample”). When there are N assemblages, input data consist of an S by N abundance matrix, or N lists of species abundances.

# Create iNEXt objects:
f_next <- iNEXT(df)

# Extract Hill numbers from iNEXT objects
f_hill <- f_next$AsyEst %>%
  separate(Assemblage, c('site_code', 'collection_date'), sep = ' ') %>%
  as_tibble()

f_div <- f_hill %>%
  select(site_code:Observed) %>%
  pivot_wider(names_from = Diversity, values_from = Observed) %>%
  dplyr::rename(richness = `Species richness`,
         shannon = `Shannon diversity`,
         simpson = `Simpson diversity`)

#------------------------------------------------------------------------------
# Export site-date-Diversity-BiomassBased
#------------------------------------------------------------------------------
write_csv(f_hill, 'Data/RTC_BiomassBased_fhill.csv')
write_csv(f_div, 'Data/RTC_fish_diversity_BiomassBased.csv')

#------------------------------------------------------------------------------
# Clean up
rm(list=ls())
#------------------------------------------------------------------------------

# End diversity_Biomass_calculate