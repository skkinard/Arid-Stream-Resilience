# merge_mega_fish_length
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

source('preparation/merge_toolkit.R') # load packages and helper-functions

m_rapid <- read_csv('clean_data/m_RAPID.csv')
m_terrg <- read_csv('clean_data/m_terrg_length.csv')

#------------------------------------------------------------------------------
# merge length data
#------------------------------------------------------------------------------
d_forklength <- full_join(m_rapid %>% mutate(project='RAPID'),
                          m_terrg%>% mutate(project='TERRG'))

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(d_forklength, 'clean_data/mega_fish_length.csv')

#------------------------------------------------------------------------------
# End merge_mega_fish_length

