# merge_toolkit
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
# Load Packages
library(tidyverse)
library(lubridate)

# site_period = lumping sampling events + or - 2 weeks to beginning of month
create_site_period <- function(d) {
  my_d <- d %>%
    mutate(my_date = Collection_Date) %>%
    separate(my_date, c("Year", "Month", "Day"), sep = '-') %>%
    mutate(Day = as.numeric(Day),
           Month = as.numeric(Month)) %>%
    mutate(Collection_Month = ifelse(Day <= 20,
                                     Month, as.character(Month + 1))) %>%
    mutate(Collection_Period = ym(paste(Year, as.character(Collection_Month), sep = '-'))) %>%
    mutate(site_period = paste(Site_Code, Collection_Period, sep = '_')) %>%
    dplyr::select( - Year, -Month, -Day, -Collection_Month) 
  
  return(my_d)  }

create_terrg_period <- function(d) {
  my_d <- d %>%
    mutate(my_date = Collection_Date) %>%
    separate(my_date, c("Year", "Month", "Day"), sep = '-') %>%
    mutate(Day = as.numeric(Day),
           Month = as.numeric(Month)) %>%
    mutate(Collection_Month = Month) %>%
    mutate(Collection_Period = ym(paste(Year, 
                                        as.character(Collection_Month), 
                                        sep = '-'))) %>%
    mutate(site_period = paste(Site_Code, Collection_Period, sep = '_')) %>%
    dplyr::select( - Year, -Month, -Day, -Collection_Month) 
  
  return(my_d) }

fix_fish_spelling <- function(d) {
  d %>%
    mutate(genus = str_replace_all(genus,'Menidia', 'Minidia' ),
           species = str_replace_all(species, 'cyanoguttatum', 'cyanoguttatus')) %>%
    mutate(family = str_replace_all(family, 'Pecidae', 'Percidae')) %>%
    mutate(commonname = str_replace_all(commonname, 'formosa', 'latipinna'), 
           genus_species = str_replace_all(genus_species, 'formosa', 'latipinna'), 
           species = str_replace_all(species, 'formosa', 'latipinna'), 
           lowest_taxon = str_replace_all(lowest_taxon, 'formosa', 'latipinna')) }

r_friendly_colnames <- function(d) {
  output <- d
  colnames(output) <- str_to_lower(colnames(output))
  colnames(output) <- str_replace_all(colnames(output), '\\.', '_') 
  return(output) }

#------------------------------------------------------------------------------
# End merge_toolkit