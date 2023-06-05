# prep_fish_species
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
# Load Packages
library(tidyverse)
library(lubridate)

# Set Workspace
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

# load data
fish_species <- read_csv("source_data/Kinard_data/my_fish_species.csv") 
#------------------------------------------------------------------------------
# Fix spelling errors
fish_species <- fish_species %>%
  mutate(order = str_replace_all(order,
                                 'Cyprinidontiformes', 
                                 'Cyprinodontiformes'),
         order = str_replace_all(order, 
                                 'Centrachiformes', 
                                 'Centrarchiformes'))
#------------------------------------------------------------------------------
# export clean data
write_csv(fish_species, 'clean_data/fish_species.csv')
#------------------------------------------------------------------------------
# end prep_fish_species