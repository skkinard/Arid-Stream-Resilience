# 001_site_longterm
# 1-19-23
# Sean Kinard

# Extract a priori environmental variables
# fill gaps with seasonal means (3x 4 month seasons)
# if no seasonal mean, fill with annual mean

#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------

# Clean Workspace
rm(list=ls())

# set working directory:
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience")

# load packages
library(tidyverse)
library(lubridate)

# load data
d <- read_csv('Data/fish_biomass_megaframe.csv')
d_meta <- read_csv('Data/fish_rapid_terrg_megaframe_metadata.csv')

#------------------------------------------------------------------------------
# Long Term environmental variables: data extraction
#------------------------------------------------------------------------------

ltevars <- d %>%
  select(site_code,annualrain, developedland:otherland, hfpp7:vardf,
         bfi:medianflow, basinsize, flsh:seasonality ) %>%
  distinct()

#------------------------------------------------------------------------------
# TR missing flow metrics
#------------------------------------------------------------------------------
# Fernando Thesis has HFPP7 for TR = 0.049
ltevars <- ltevars %>%
  mutate(hfpp7 = ifelse(site_code == 'TR', 0.049, hfpp7))

# Estimate hfpp3 for TR

est_hfpp3 <- function(my_data, my_site_code) {
  
  hfpp3_lm <- summary(lm(hfpp3 ~ hfpp7 + annualrain, data = my_data))
  
  yint <- hfpp3_lm$coefficients[1]
  m1 <- hfpp3_lm$coefficients[2]
  m2 <- hfpp3_lm$coefficients[3]
  x1 <- my_data %>% filter(site_code == my_site_code) %>% pull(hfpp7)
  x2 <- my_data %>% filter(site_code == my_site_code) %>% pull(annualrain)
  
  est_hfpp3 <- yint + m1*x1 +m2*x2
  
  return(est_hfpp3)
}

# replace NA with modeled hfpp3 for TR

ltevars <- ltevars %>%
  mutate(hfpp3 = ifelse(site_code == 'TR', est_hfpp3(ltevars, 'TR'), hfpp3))

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------

# write files
write_csv(ltevars, '/home/kinard/Documents/Research/Dissertation/02_Resilience/Data/site_longterm_nogap.csv')

#------------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------------
rm(d, d_meta)

# End 001_site_longterm