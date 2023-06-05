# Resistance_Resilience_Calculations
# Sean Kinard
# 2023-03-16

# -----------------------------------------------------------------------------
# Strickland Methods:
# -----------------------------------------------------------------------------
# Resistance was calculated as a log-response ratio (LRR), the natural log of the maximum change created by the storm divided by the baseline value. LRR revealed the magnitude of impact or general resistance of the system for each community metric. 

# Resilience was calculated as return time, the number of days until values resembled the baseline values.

# -----------------------------------------------------------------------------
# setup
# -----------------------------------------------------------------------------

setwd('/home/kinard/Documents/Research/Dissertation/02_Resilience')

library(tidyverse) # all tidyverse packages
library(lubridate) # date managment

source('R_Scripts/067_LRR_calc.R')

d_storm <- read_csv('Data/storm_stats_summary.csv')

# -----------------------------------------------------------------------------
# Maximum change after hurricane
# -----------------------------------------------------------------------------
LRR_max_spe <- RR_combined %>%
  filter(collection_date < ymd('2018-01-01')) %>%
  select(any_of(my_index), lowest_taxon, contains('spe')) %>%
  left_join(d_storm%>%select(any_of(my_index), contains('hu_'))) %>%
  pivot_longer(
    cols=contains('RR'), names_to='LRR_metric', values_to='LRR_value') %>%
  group_by(site_code, lowest_taxon, LRR_metric) %>%
  dplyr::summarize(LRR_max_spe = max(LRR_value))

LRR_max_fam <- RR_combined %>%
  filter(collection_date < ymd('2018-01-01')) %>%
  add_taxonomic() %>%
  select(any_of(my_index), family, contains('fam')) %>%
  left_join(d_storm%>%select(any_of(my_index), contains('hu_'))) %>%
  pivot_longer(
    cols=contains('RR'), names_to='LRR_metric', values_to='LRR_value') %>%
  group_by(site_code, family, LRR_metric) %>%
  dplyr::summarize(LRR_max_fam = max(LRR_value))

LRR_max_com <- RR_combined %>%
  filter(collection_date < ymd('2018-01-01')) %>%
  add_taxonomic() %>%
  select( ! contains('spe')) %>%
  select( ! contains('fam')) %>%
  select( ! contains('LRR')) %>%
  left_join(d_storm%>%select(any_of(my_index), contains('hu_'))) %>%
  pivot_longer(
    cols=contains('RR'), names_to='RR_metric', values_to='RR_value') %>%
  group_by(site_code, RR_metric) %>%
  dplyr::summarize(RR_max_com = max(RR_value))

# Visualize
LRR_max_com %>%
  left_join(d_environment) %>%
  fix_site_order() %>%
  ggplot(aes(x=annualrain, y=RR_max_com)) +
  facet_wrap(~RR_metric, scales='free') +
  geom_smooth(method='loess', se=F, color = 'white', linewidth=.5, linetype=2) +
  geom_point(aes(fill=site_code),
             shape=21, size=3, alpha=.3) +
  geom_point(aes(color=site_code),
             shape=21, size=3)+
  dark_theme_gray() +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors)

# -----------------------------------------------------------------------------
# Return Interval after hurricane
# -----------------------------------------------------------------------------
LRR_RI_spe <- RR_combined %>%
  filter(collection_date < ymd('2018-01-01')) %>%
  select(any_of(my_index), lowest_taxon, contains('spe')) %>%
  left_join(d_storm%>%select(any_of(my_index), contains('hu_'))) %>%
  pivot_longer(
    cols=contains('RR'), names_to='LRR_metric', values_to='LRR_value') %>%
  group_by(site_code, lowest_taxon, LRR_metric) %>%
  dplyr::summarize(LRR_max_spe = max(LRR_value))

LRR_RI_fam <- RR_combined %>%
  filter(collection_date < ymd('2018-01-01')) %>%
  add_taxonomic() %>%
  select(any_of(my_index), family, contains('fam')) %>%
  left_join(d_storm%>%select(any_of(my_index), contains('hu_'))) %>%
  pivot_longer(
    cols=contains('RR'), names_to='LRR_metric', values_to='LRR_value') %>%
  group_by(site_code, family, LRR_metric) %>%
  dplyr::summarize(LRR_max_fam = max(LRR_value))

LRR_RI_com <- RR_combined %>%
  filter(collection_date < ymd('2018-01-01')) %>%
  add_taxonomic() %>%
  select( ! contains('spe')) %>%
  select( ! contains('fam')) %>%
  select( ! contains('LRR')) %>%
  left_join(d_storm%>%select(any_of(my_index), contains('hu_'))) %>%
  pivot_longer(
    cols=contains('RR'), names_to='RR_metric', values_to='RR_value') %>%
  group_by(site_code, RR_metric) %>%
  dplyr::summarize(RR_max_com = max(RR_value))

# Visualize
LRR_max_com %>%
  left_join(d_environment) %>%
  fix_site_order() %>%
  ggplot(aes(x=annualrain, y=RR_max_com)) +
  facet_wrap(~RR_metric, scales='free') +
  geom_smooth(method='loess', se=F, color = 'white', linewidth=.5, linetype=2) +
  geom_point(aes(fill=site_code),
             shape=21, size=3, alpha=.3) +
  geom_point(aes(color=site_code),
             shape=21, size=3)+
  dark_theme_gray() +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors)

# visualise RR vs baseline
RR_biomass %>%
  select(site_code, collection_date, bas_biomass_com, B100_sum_com) %>%
  unique() %>%
  filter(collection_date < ymd('2020-01-01')) %>%
  left_join(d_environment) %>%
  fix_site_order() %>%
  ggplot(aes(x=collection_date, y=B100_sum_com)) +
  facet_wrap(~site_code, scales='free', ncol=3) +
  geom_point(aes(fill=site_code, color=site_code),
             shape=21, size=3, alpha=.3) +
  geom_point(aes(y= bas_biomass_com, color=site_code),
             shape=22, size=3) +
  dark_theme_gray() +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors)
