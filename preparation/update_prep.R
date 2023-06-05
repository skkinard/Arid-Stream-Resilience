# update_prep
# runs directory containing all data preparation files with starting with 'prep_'
# Sean Kinard
# Last Edit: 2023-06-02
#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

my_scripts <- list.files("preparation")
my_scripts <- my_scripts[substr(my_scripts, 1,4) == 'prep']

#------------------------------------------------------------------------------
# Source scripts
#------------------------------------------------------------------------------
for(i in 1:length(my_scripts)) {
  source(paste("preparation/", my_scripts[i], sep=''))
  rm(list=ls())
  my_scripts <- list.files("preparation")
  my_scripts <- my_scripts[substr(my_scripts, 1,4) == 'prep']   }

#------------------------------------------------------------------------------
# End update_R_cal