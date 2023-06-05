# prep_TERRG_nutrient
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

source('preparation/merge_toolkit.R') # load packages and helper-functions

#------------------------------------------------------------------------------
# TERRG Nutrients
#------------------------------------------------------------------------------

nutrient_terrg <- read_csv('source_data/TERRG-Drive/Nutrients_CLB.csv') %>%
  mutate(Site_Code = substr(Site, 1,2)) %>%
  dplyr::rename(Collection_Date = Date,
                Ortho_P = OrthoP,
                NO3N = NO3_N) %>%
  create_terrg_period %>%
  dplyr::select(Site_Code, Collection_Date, site_period, 
                Transect, NO3N, NH4_N, Ortho_P, DOC)

colnames(nutrient_terrg) <- str_to_lower(colnames(nutrient_terrg))

# Calculate mean and sd for transect measurements (grouped by site and sampling date)
# n =< 4 for all means and sds
nutrient_terrg <- nutrient_terrg %>%
  pivot_longer(cols = no3n : doc,
               names_to = "measure",
               values_to = "value") %>%
  group_by(site_period, measure) %>%
  dplyr::summarize(xmean = mean(value, na.rm=T)) %>%
  mutate(measure = paste(measure, 'mean', sep='_')) %>%
  pivot_wider(names_from = measure, values_from = xmean) %>%
  left_join( nutrient_terrg %>%
               pivot_longer(cols = no3n : doc,
                            names_to = "measure",
                            values_to = "value") %>%
               group_by(site_period, measure) %>%
               dplyr::summarize(xsd = sd(value, na.rm=T) ) %>%
               mutate(measure = paste(measure, 'sd', sep = '_')) %>%
               pivot_wider(names_from = measure, values_from = xsd) )

#------------------------------------------------------------------------------
# Fill missing data with nearest
#------------------------------------------------------------------------------
nutrient_terrg <- nutrient_terrg %>%
  separate(site_period, into=c('site', 'period'), sep='_')

all_periods <- nutrient_terrg %>% pull(period) %>% unique()
all_sites <- nutrient_terrg %>% pull(site) %>% unique()
d_cross <- crossing(all_sites, all_periods) %>% print(n=20) %>%
  rename(site=all_sites, period=all_periods) # all site-period combinations

nutrient_terrg <- left_join(d_cross, nutrient_terrg) %>% # expand df with gaps
  group_by(site) %>%
  fill(everything(), .direction = "downup") %>% # fill gaps with row below or above
  unite('site_period', site:period, sep='_')

#------------------------------------------------------------------------------
# Export Clean Data
#------------------------------------------------------------------------------
write_csv(nutrient_terrg, 'clean_data/terrg_nutrient.csv')

#------------------------------------------------------------------------------
# End prep_TERRG_nutrient