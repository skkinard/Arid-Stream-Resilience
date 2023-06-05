# Log Response Ratio Calculation
# Sean Kinard
# 3-16-2023
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

setwd('/home/kinard/Documents/Research/Dissertation/02_Resilience')

library(tidyverse) # all tidyverse packages
library(lubridate) # date managment

# load data
source('R_Scripts/047_ordination.R')
source('R_Scripts/000_biological_scales_comparison.R')

# -----------------------------------------------------------------------------
# Baseline (2020 qtrly x site): Biomass
# -----------------------------------------------------------------------------
bas_biomass <- d_biomass %>%
  filter(collection_date > ymd('2020-01-01')) %>%
  dplyr::rename(bas_B100_com = B100_sum_com,
                bas_B100_fam = B100_sum_fam,
                bas_B100_spe = B100_sum_spe) %>%
  add_time_vars() %>%
  select(contains('bas'), site_code, lowest_taxon, qtr)
# Biomass: annual mean species
mu20_biomass <- bas_biomass %>%
  group_by(site_code, lowest_taxon) %>%
  dplyr::summarize(mu2020_B100_spe = mean(bas_B100_spe, na.rm=T)) 
# Biomass: annual mean family
mu20_biomass <- bas_biomass %>%
  add_taxonomic() %>%
  group_by(site_code, family) %>%
  dplyr::summarize(mu2020_B100_fam = mean(bas_B100_fam, na.rm=T)) %>%
  right_join(mu20_biomass)
# Biomass: annual mean community
mu20_biomass <- bas_biomass %>%
  group_by(site_code) %>%
  dplyr::summarize(mu2020_B100_com = mean(bas_B100_com, na.rm=T)) %>%
  right_join(mu20_biomass)

# -----------------------------------------------------------------------------
# Baseline (2020 qtrly x site): Abundance
# -----------------------------------------------------------------------------
bas_density <- d_density %>%
  filter(collection_date > ymd('2020-01-01')) %>%
  dplyr::rename(bas_d100_com = d100_com,
         bas_d100_fam = d100_fam,
         bas_d100_spe = d100_spe) %>%
  add_time_vars() %>%
  select(contains('bas'), site_code, lowest_taxon, qtr)
# Abundance: annual mean species
mu20_density <- bas_density %>%
  group_by(site_code, lowest_taxon) %>%
  dplyr::summarize(mu2020_d100_spe = mean(bas_d100_spe, na.rm=T))
# Abundance: annual mean family
mu20_density <- bas_density %>%
  add_taxonomic() %>%
  group_by(site_code, family) %>%
  dplyr::summarize(mu2020_d100_fam = mean(bas_d100_fam, na.rm=T)) %>%
  right_join(mu20_density)
# Abundance: annual mean community
mu20_density <- bas_density %>%
  group_by(site_code) %>%
  dplyr::summarize(mu2020_d100_com = mean(bas_d100_com, na.rm=T)) %>%
  right_join(mu20_density)

# -----------------------------------------------------------------------------
# Baseline (2020 qtrly x site): Diversity
# -----------------------------------------------------------------------------
bas_diversity <- d_diversity %>%
  filter(collection_date > ymd('2020-01-01')) %>%
  dplyr::rename(bas_shannon = shannon,
         bas_simpson = simpson,
         bas_richness = richness) %>%
  add_time_vars() %>%
  select(contains('bas'), site_code, qtr)
# Diveristy: annual average
mu20_diversity <- bas_diversity %>%
  group_by(site_code) %>%
  dplyr::summarize(mu2020_shannon = mean(bas_shannon, na.rm=T),
                   mu2020_simpson = mean(bas_simpson, na.rm=T),
                   mu2020_richness = mean(bas_richness, na.rm=T))

# -----------------------------------------------------------------------------
# Baseline (2020 qtrly x site): Composition
# -----------------------------------------------------------------------------
bas_axis <- rda_all$sites_long %>%
  as_tibble() %>%
  select(site_code, collection_date, axis1, axis2) %>%
  filter(collection_date > ymd('2020-01-01')) %>%
  dplyr::rename(bas_axis1 = axis1,
         bas_axis2 = axis2) %>%
  mutate(collection_date = as_date(collection_date)) %>%
  add_time_vars() %>%
  select(contains('bas'), site_code, qtr)

mu2020_centroid <- bas_axis %>%
  group_by(site_code) %>%
  dplyr::summarize(mu2020_axis1_mu = mean(bas_axis1, na.rm=T),
                   mu2020_axis2_mu = mean(bas_axis2, na.rm=T))

bas_ctrdist <- bas_axis %>%
  left_join(mu2020_centroid) %>%
  mutate(bas_ctrdist = sqrt((bas_axis1-mu2020_axis1_mu)^2 + 
                                         (bas_axis2-mu2020_axis2_mu)^2)) %>%
  select(contains('bas'), site_code, qtr)

mu2020_ctrdist <- bas_ctrdist %>%
  group_by(site_code) %>%
  dplyr::summarize(mu2020_ctrdist = mean(bas_ctrdist, na.rm=T))
# -----------------------------------------------------------------------------
# Log Response Ratio Calculations: Biomass
# -----------------------------------------------------------------------------

RR_biomass <- d_biomass %>%
  add_time_vars() %>%
  select(contains('B100'), site_code, family, lowest_taxon, collection_date, qtr) %>%
  left_join(mu20_biomass) %>%
  mutate(mu2020_B100_com = ifelse(is.na(mu2020_B100_com), 0, mu2020_B100_com),
         mu2020_B100_fam = ifelse(is.na(mu2020_B100_fam), 0, mu2020_B100_fam),
         mu2020_B100_spe = ifelse(is.na(mu2020_B100_spe), 0, mu2020_B100_spe) ) %>%
  left_join(bas_biomass) %>%
  mutate( # gapfill baselines with 2020 mean
    bas_biomass_com = ifelse(is.na(bas_B100_com), 
                             mu2020_B100_com, 
                             bas_B100_com),
    bas_biomass_fam = ifelse(is.na(bas_B100_fam), 
                             mu2020_B100_fam,
                             bas_B100_fam),
    bas_biomass_spe = ifelse(is.na(bas_B100_spe), 
                             mu2020_B100_spe, 
                             bas_B100_spe)) %>%
  mutate( # calculate response ratio (RR)
    RR_biomass_com = (B100_sum_com-bas_biomass_com ) / 
                             (bas_biomass_com + .001),
    RR_biomass_fam = (B100_sum_fam-bas_biomass_fam) / 
                             (bas_biomass_fam + .001),
    RR_biomass_spe = (B100_sum_spe-bas_biomass_spe) / 
                             (bas_biomass_spe + .001)) %>%
  mutate( # calculate log response ratio (LRR)
    LRR_biomass_com = log(abs(RR_biomass_com)),
    LRR_biomass_fam = log(abs(RR_biomass_com)),
    LRR_biomass_spe = log(abs(RR_biomass_spe))) %>%
  mutate( # replace infinity with -100000 LRR
    LRR_biomass_com = ifelse(
      is.infinite(LRR_biomass_com), -100000, LRR_biomass_com),
    LRR_biomass_fam = ifelse(
      is.infinite(LRR_biomass_fam), -100000, LRR_biomass_fam),
    LRR_biomass_spe = ifelse(
      is.infinite(LRR_biomass_spe), -100000, LRR_biomass_spe)) %>%
  filter(collection_date < ('2020-01-01')) %>%
  ungroup() %>%
  fix_site_order()

# Visualize Response Ratios
RR_biomass %>%
  select(site_code, family, collection_date, RR_biomass_fam) %>%
  ggplot(aes(x=collection_date, y = RR_biomass_fam)) +
  facet_wrap(~family, scales='free') +
  geom_smooth(aes(color=site_code),
              method = "loess", se=F, span=.4, alpha=.1, linewidth=.3) +
  geom_point(aes(fill=site_code), 
             shape=21, size=3, alpha=.6) +
  geom_point(aes(color=site_code), 
             shape=21, size=3) +
  dark_theme_gray() +
  scale_color_manual(values = my_colors, breaks = site_order) +
  scale_fill_manual(values = my_colors, breaks = site_order) +
  dark_theme_grey(base_size = 12) +
  scale_y_log10()

# -----------------------------------------------------------------------------
# Log Response Ratio Calculations: Abundance
# -----------------------------------------------------------------------------
RR_density <- d_density %>%
  add_time_vars() %>%
  select(contains('d100'), site_code, family, lowest_taxon, collection_date, qtr) %>%
  left_join(mu20_density) %>%
  mutate(mu2020_d100_com = ifelse(is.na(mu2020_d100_com), 0, mu2020_d100_com),
         mu2020_d100_fam = ifelse(is.na(mu2020_d100_fam), 0, mu2020_d100_fam),
         mu2020_d100_spe = ifelse(is.na(mu2020_d100_spe), 0, mu2020_d100_spe) ) %>%
  left_join(bas_density) %>%
  mutate( # gapfill baselines with 2020 mean
    bas_density_com = ifelse(is.na(bas_d100_com), 
                             mu2020_d100_com, 
                             bas_d100_com),
    bas_density_fam = ifelse(is.na(bas_d100_fam), 
                             mu2020_d100_fam,
                             bas_d100_fam),
    bas_density_spe = ifelse(is.na(bas_d100_spe), 
                             mu2020_d100_spe, 
                             bas_d100_spe)) %>%
  mutate( # calculate response ratio (RR)
    RR_density_com = (d100_com-bas_density_com ) / 
      (bas_density_com + .001),
    RR_density_fam = (d100_fam-bas_density_fam) / 
      (bas_density_fam + .001),
    RR_density_spe = (d100_spe-bas_density_spe) / 
      (bas_density_spe + .001)) %>%
  mutate( # calculate log response ratio (LRR)
    LRR_density_com = log(abs(RR_density_com)),
    LRR_density_fam = log(abs(RR_density_com)),
    LRR_density_spe = log(abs(RR_density_spe))) %>%
  mutate( # replace infinity with -100000 LRR
    LRR_density_com = ifelse(
      is.infinite(LRR_density_com), -100000, LRR_density_com),
    LRR_density_fam = ifelse(
      is.infinite(LRR_density_fam), -100000, LRR_density_fam),
    LRR_density_spe = ifelse(
      is.infinite(LRR_density_spe), -100000, LRR_density_spe)) %>%
  filter(collection_date < ('2020-01-01')) %>%
  ungroup() %>%
  fix_site_order()

# -----------------------------------------------------------------------------
# Log Response Ratio Calculations: Diversity
# -----------------------------------------------------------------------------
RR_diversity <- d_diversity %>%
  add_time_vars() %>%
  select(shannon, simpson, richness, site_code, collection_date, qtr) %>%
  left_join(mu20_diversity) %>% # merge 2020 annual means
  mutate(mu2020_shannon = ifelse(is.na(mu2020_shannon), 0, mu2020_shannon),
         mu2020_simpson = ifelse(is.na(mu2020_simpson), 0, mu2020_simpson),
         mu2020_richness = ifelse(is.na(mu2020_richness), 0, mu2020_richness) ) %>%
  left_join(bas_diversity) %>% # merge 2020 quarterly baslines
  mutate( # gapfill quarterly baselines with 2020 means
    bas_shannon = ifelse(is.na(bas_shannon), 
                             mu2020_shannon, 
                             bas_shannon),
    bas_simpson = ifelse(is.na(bas_simpson), 
                             mu2020_simpson,
                             bas_simpson),
    bas_richness = ifelse(is.na(bas_richness), 
                             mu2020_richness, 
                             bas_richness)) %>%
  mutate( # calculate response ratio (RR)
    RR_shannon = (shannon-bas_shannon ) / 
      (bas_shannon + .001),
    RR_simpson = (simpson-bas_simpson) / 
      (bas_simpson + .001),
    RR_richness = (richness-bas_richness) / 
      (bas_richness + .001)) %>%
  mutate( # calculate log response ratio (LRR)
    LRR_shannon = log(abs(RR_shannon)),
    LRR_simpson = log(abs(RR_shannon)),
    LRR_richness = log(abs(RR_richness))) %>%
  mutate( # replace infinity with -100000 LRR
    LRR_shannon = ifelse(
      is.infinite(LRR_shannon), -100000, LRR_shannon),
    LRR_simpson = ifelse(
      is.infinite(LRR_simpson), -100000, LRR_simpson),
    LRR_richness = ifelse(
      is.infinite(LRR_richness), -100000, LRR_richness)) %>%
  filter(collection_date < ('2020-01-01')) %>%
  ungroup() %>%
  fix_site_order()
# -----------------------------------------------------------------------------
# Log Response Ratio Calculations: Composition
# -----------------------------------------------------------------------------
RR_composition <- d_axis %>%
  left_join(mu2020_centroid) %>%
  add_time_vars() %>%
  select(contains('axis'), site_code, collection_date, qtr) %>%
  left_join(bas_ctrdist) %>%
  left_join(mu2020_ctrdist) %>%
  mutate(ctrdist = sqrt((axis1-mu2020_axis1_mu)^2 + 
                          (axis2-mu2020_axis2_mu)^2)) %>%
  mutate( # gapfill quarterly baselines with 2020 means
    bas_ctrdist = ifelse(is.na(bas_ctrdist), 
                         mu2020_ctrdist, 
                         bas_ctrdist)) %>%
  mutate( # calculate response ratio (RR)
    RR_ctrdist = (ctrdist-bas_ctrdist ) / 
      (bas_ctrdist + .001)) %>%
  mutate( # calculate log response ratio (LRR)
    LRR_ctrdist = log(abs(RR_ctrdist))) %>%
  mutate( # replace infinity with -100000 LRR
    LRR_ctrdist = ifelse(is.infinite(LRR_ctrdist), -100000, LRR_ctrdist)) %>%
  filter(collection_date < ('2020-01-01')) %>%
  ungroup() %>%
  fix_site_order()

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
slim_data <- function(my_data) {
  my_columns <- c('site_code', 'collection_date', 'lowest_taxon')
  my_data %>%
    select(any_of(my_columns), contains('RR')) }

write_csv(RR_biomass %>% slim_data(), 'Data/RR_biomass.csv')
write_csv(RR_density %>% slim_data(), 'Data/RR_density.csv')
write_csv(RR_diversity %>% slim_data(), 'Data/RR_diversity.csv')
write_csv(RR_composition %>% slim_data(), 'Data/RR_composition.csv')

RR_combined <- RR_biomass %>% slim_data() %>%
  left_join(RR_density%>%slim_data()) %>%
  left_join(RR_diversity%>%slim_data()) %>%
  left_join(RR_composition%>%slim_data())

write_csv(RR_combined, 'Data/RR_combined.csv')

# -----------------------------------------------------------------------------
# End Log Response Ratio Calculation