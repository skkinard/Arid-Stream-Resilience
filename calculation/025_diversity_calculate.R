# diversity_calculate
# Sean Kinard
# 1-25-2023

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
d <- read_csv('Data/fish_d100.csv')

# format: rows = species,  column = 'site_code collection_date'
df <- d %>%
  mutate(UID=paste(site_code, collection_date, sep=' ')) %>%
  select(lowest_taxon, UID, d100_spe) %>%
  pivot_wider(names_from=UID, values_from = d100_spe, values_fill = 0) %>%
  column_to_rownames(var = "lowest_taxon") %>%
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
# Export site-date-Diversity
#------------------------------------------------------------------------------
write_csv(f_hill, 'Data/RTC_fhill.csv')
write_csv(f_div, 'Data/RTC_fish_diversity.csv')

#------------------------------------------------------------------------------
# Clean up
rm(list=ls())
#------------------------------------------------------------------------------

# End diversity_calculate