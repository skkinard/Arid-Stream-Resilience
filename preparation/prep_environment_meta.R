# prep_environment_meta
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
# Load Packages
library(tidyverse)

# Set Workspace
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

# load data
environment_meta <- read_csv("source_data/TERRG-Drive/RAPID_FINAL_DATA/meta_data.csv") 
#------------------------------------------------------------------------------
# R-friendly column headers
environment_meta <- environment_meta %>%
  dplyr::rename(data_sheet = `Data Sheet`, 
                column_name = `Column Name`, 
                units = Units, 
                description = Description)
#------------------------------------------------------------------------------
# export clean data
write_csv(environment_meta, 'clean_data/environment_meta.csv')
#------------------------------------------------------------------------------
# end prep_environment_meta