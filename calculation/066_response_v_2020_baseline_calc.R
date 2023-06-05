# response_v_2020_baseline_calc
# Sean Kinard
# 1-31-2023

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
library(ggpubr)

#------------------------------------------------------------------------------
# Tidy data
#------------------------------------------------------------------------------

# length
d_len <- read_csv('Data/forklength_all_megaframe.csv') %>%
              group_by(site_code, collection_date) %>%
              dplyr::summarize(length_mu = mean(lengthmm))

# density
d_den <- read_csv('Data/fish_d100.csv') %>%
              group_by(site_code, collection_date) %>%
              dplyr::summarize(d100_mu = mean(d100)) 

# biomass
d_bio <- read_csv('Data/fish_biomass_megaframe.csv') %>%
              group_by(site_code, collection_date) %>%
              dplyr::summarize(biomass_mu = mean(biomass)) 

# diversity
d_div <- read_csv('Data/RTC_fish_diversity.csv') %>%
  select( ! c(richness, simpson))

# my environmental variables
my_evars <- c('annualrain', 'flsh', 'lfpp', 'hfpp3',
              'q_2w_max', 'depth_mx', 'silt, gravel',
              'total_algae', 'canopy_density_mid',
              'conductivity', 'do_mg_l',
              'nh4_n', 'no3n', 'ortho_p')

# short_term environment
d_ste <- read_csv('Data/site_daily_nogap.csv') %>%
  select(site_code, collection_date, any_of(my_evars)) %>%
  mutate(cd=collection_date) %>%
  separate(cd, into = c('year', 'month', 'day')) %>%
  mutate(year = as.numeric(year),
         month = as.numeric(month),
         day = as.numeric(day)) %>%
  mutate(qtr=ifelse(month<4, 'Q1',
             ifelse(month>=4 & month<7, 'Q2',
             ifelse(month>=7 & month <10, 'Q3', 
             ifelse(month>=10 & month <13, 'Q4', NA)))))

# long-term environment
d_lte <- read_csv('Data/site_longterm_nogap.csv') %>%
  select(site_code, any_of(my_evars))

# merge
df <- d_len %>%
  left_join(d_den) %>%
  left_join(d_bio) %>%
  left_join(d_div) %>%
  left_join(d_ste) %>%
  left_join(d_lte)

#------------------------------------------------------------------------------
# Calculate baseline
#------------------------------------------------------------------------------

# 2020 Annual Average (A20)
A20 <- df %>%
  filter(collection_date > ymd('2020-01-01')) %>%
  group_by(site_code) %>%
  dplyr::summarize(length_A20_mu=mean(length_mu, na.rm=T),
                   length_A20_sd=sd(length_mu, na.rm=T),
                   d100_A20_mu=mean(d100_mu, na.rm=T),
                   d100_A20_sd=sd(d100_mu, na.rm=T),
                   biomass_A20_mu=mean(biomass_mu, na.rm=T),
                   biomass_A20_sd=sd(biomass_mu, na.rm=T),
                   shannon_A20_mu=mean(shannon, na.rm=T),
                   shannon_A20_sd=sd(shannon, na.rm=T))

# Baseline = 2020_seasonal_mean, 2020_annual_standard deviation,
df <- A20 %>%
  right_join(df) %>% 
  select(site_code, collection_date, qtr, everything()) %>%
  arrange(site_code, collection_date) %>%
  filter(collection_date < ymd('2020-01-01'))

#------------------------------------------------------------------------------
# Calculate Response Ratio
#------------------------------------------------------------------------------

# Realized Disturbance (RR response ratio)
df <- df %>%
  mutate(RR_length = (length_mu-length_A20_mu) / length_A20_mu ,
         RR_d100 = (d100_mu-d100_A20_mu) / d100_A20_mu ,
         RR_biomass = (biomass_mu-biomass_A20_mu) / biomass_A20_mu ,
         RR_shannon = (shannon-shannon_A20_mu) / shannon_A20_mu)

# Log absolute value Realized Disturbance (LRR response ratio)
df <- df %>%
  mutate(LRR_length = log(abs((length_mu-length_A20_mu) / length_A20_mu)) ,
         LRR_d100 = log(abs((d100_mu-d100_A20_mu) / d100_A20_mu)) ,
         LRR_biomass = log(abs((biomass_mu-biomass_A20_mu) / biomass_A20_mu)) ,
         LRR_shannon = log(abs((shannon-shannon_A20_mu) / shannon_A20_mu)))

#------------------------------------------------------------------------------
# Calculate Recover Status (RS)
#------------------------------------------------------------------------------
# compare value to baseline with t-test, first match = recovery
# or time to 0 LRR

# Recover Status (RS)
df <- df %>% arrange(site_code, collection_date) %>% # density
  mutate(y_margin = qt(0.95, df=4-1)*d100_A20_sd / sqrt(4),
         LCL_d100 = d100_A20_mu - y_margin,
         UCL_d100 = d100_A20_mu + y_margin,
         RS_d100 = ifelse(d100_mu < LCL_d100, 'under',
                          ifelse(d100_mu > UCL_d100, 'over',
                                 ifelse(d100_mu > LCL_d100 &
                                          d100_mu < UCL_d100, 'ordinary', NA))) )
df <- df %>% arrange(site_code, collection_date) %>% # length
  mutate(y_margin = qt(0.95, df=4-1)*length_A20_sd / sqrt(4),
         LCL_length = length_A20_mu - y_margin,
         UCL_length = length_A20_mu + y_margin,
         RS_length = ifelse(length_mu < LCL_length, 'under',
                          ifelse(length_mu > UCL_length, 'over',
                                 ifelse(length_mu > LCL_length &
                                          length_mu < UCL_length, 'ordinary', NA))) )
df <- df %>% arrange(site_code, collection_date) %>% # biomass
  mutate(y_margin = qt(0.95, df=4-1)*biomass_A20_sd / sqrt(4),
         LCL_biomass = biomass_A20_mu - y_margin,
         UCL_biomass = biomass_A20_mu + y_margin,
         RS_biomass = ifelse(biomass_mu < LCL_biomass, 'under',
                          ifelse(biomass_mu > UCL_biomass, 'over',
                                 ifelse(biomass_mu > LCL_biomass &
                                          biomass_mu < UCL_biomass, 'ordinary', NA))) )
df <- df %>% arrange(site_code, collection_date) %>% # shannon
  mutate(y_margin = qt(0.95, df=4-1)*shannon_A20_sd / sqrt(4),
         LCL_shannon = shannon_A20_mu - y_margin,
         UCL_shannon = shannon_A20_mu + y_margin,
         RS_shannon = ifelse(shannon < LCL_shannon, 'under',
                       ifelse(shannon > UCL_shannon, 'over',
                       ifelse(shannon > LCL_shannon &
                                shannon < UCL_shannon, 'ordinary', NA))) )

#------------------------------------------------------------------------------
# Calculate Recovery Interval (RI days)
#------------------------------------------------------------------------------

df <- df %>% # length
  filter(RS_length == 'ordinary') %>%
  group_by(site_code) %>%
  dplyr::summarize(
    RI_length = as.numeric(min(collection_date) - ymd('2017-08-27'))) %>%
  right_join(df)

df <- df %>% # density
  filter(RS_d100 == 'ordinary') %>%
  group_by(site_code) %>%
  dplyr::summarize(
    RI_d100 = as.numeric(min(collection_date) - ymd('2017-08-27'))) %>%
  right_join(df)

df <- df %>% # biomass
  filter(RS_biomass == 'ordinary') %>%
  group_by(site_code) %>%
  dplyr::summarize(
    RI_biomass = as.numeric(min(collection_date) - ymd('2017-08-27'))) %>%
  right_join(df)

df <- df %>% # shannon
  filter(RS_shannon == 'ordinary') %>%
  group_by(site_code) %>%
  dplyr::summarize(RI_shannon = 
                     as.numeric(min(collection_date) - ymd('2017-08-27'))) %>%
  right_join(df)

#------------------------------------------------------------------------------
# visualize baseline relative to site values
#------------------------------------------------------------------------------

# length
length_v_base20 <- df %>%
  mutate(annualrain = paste(round(annualrain,0), 'cm/yr', sep=' ')) %>%
  ggplot(aes(x=collection_date, y=length_mu, fill=RS_length)) +
  facet_wrap(~annualrain, scales='free') +
  geom_ribbon(aes(x=collection_date, ymin=LCL_length, ymax=UCL_length),
              fill = 'grey50', alpha=.2) +
  geom_smooth(aes(y=length_A20_mu), method = 'lm', 
              color = 'black', linetype = 2, alpha = .5, show.legend = F) +
  geom_point(shape = 21, size = 4) +
  ylab('Total Length (mm)') +
  xlab(element_blank()) +
  labs(fill = element_blank()) +
  ggtitle('Average Fish Length: baseline = 2020 annual average with 90% CL')

# d100
d100_v_base20 <- df %>%
  mutate(annualrain = paste(round(annualrain,0), 'cm/yr', sep=' ')) %>%
  ggplot(aes(x=collection_date, y=d100_mu, fill=RS_d100)) +
  facet_wrap(~annualrain, scales='free') +
  geom_ribbon(aes(x=collection_date, ymin=LCL_d100, ymax=UCL_d100),
              fill = 'grey50', alpha=.2) +
  geom_smooth(aes(y=d100_A20_mu), method = 'lm', 
              color = 'black', linetype = 2, alpha = .5, show.legend = F) +
  geom_point(shape = 21, size = 4) +
  ylab('Density (fish / 100m2)') +
  xlab(element_blank()) +
  labs(fill = element_blank()) +
  ggtitle('Total density (fish/100m2): baseline = 2020 annual average with 90% CL')

# biomass
biomass_v_base20 <- df %>%
  mutate(annualrain = paste(round(annualrain,0), 'cm/yr', sep=' ')) %>%
  ggplot(aes(x=collection_date, y=biomass_mu, fill=RS_biomass)) +
  facet_wrap(~annualrain, scales='free') +
  geom_ribbon(aes(x=collection_date, ymin=LCL_biomass, ymax=UCL_biomass),
              fill = 'grey50', alpha=.2) +
  geom_smooth(aes(y=biomass_A20_mu), method = 'lm', 
              color = 'black', linetype = 2, alpha = .5, show.legend = F) +
  geom_point(shape = 21, size = 4) +
  ylab('Biomass (g/m2)') +
  xlab(element_blank()) +
  labs(fill = element_blank()) +
  ggtitle('Total biomass (g/m2): baseline = 2020 annual average with 90% CL')

# shannon
shannon_v_base20 <- df %>%
  mutate(annualrain = paste(round(annualrain,0), 'cm/yr', sep=' ')) %>%
  ggplot(aes(x=collection_date, y=shannon, fill=RS_shannon)) +
  facet_wrap(~annualrain, scales='free') +
  geom_ribbon(aes(x=collection_date, ymin=LCL_shannon, ymax=UCL_shannon),
              fill = 'grey50', alpha=.2) +
  geom_smooth(aes(y=shannon_A20_mu), method = 'lm', 
              color = 'black', linetype = 2, alpha = .5, show.legend = F) +
  geom_point(shape = 21, size = 4) +
  ylab('Shannon Diversity') +
  xlab(element_blank()) +
  labs(fill = element_blank()) +
  ggtitle('Diversity: baseline = 2020 annual average with 90% CL')

#------------------------------------------------------------------------------
# visualize RR and RI
#------------------------------------------------------------------------------

RR_plot_20base <- df %>% 
  filter(collection_date < ymd('2017-12-31')) %>%
  group_by(annualrain) %>%
  dplyr::summarize(length_mi = min(RR_length),
                   length_mx = max(RR_length),
                   d100_mi = min(RR_d100),
                   d100_mx = max(RR_d100),
                   biomass_mi = min(RR_biomass),
                   biomass_mx = max(RR_biomass),
                   shannon_mi = min(RR_shannon),
                   shannon_mx = max(RR_shannon) ) %>%
  mutate(absmax_length = ifelse(abs(length_mi) > abs(length_mx),
                                length_mi, length_mx),
         absmax_d100 = ifelse(abs(d100_mi) > abs(d100_mx),
                              d100_mi, d100_mx),
         absmax_biomass = ifelse(abs(biomass_mi) > abs(biomass_mx),
                                 biomass_mi, biomass_mx),
         absmax_shannon = ifelse(abs(shannon_mi) > abs(shannon_mx),
                                 shannon_mi, shannon_mx) ) %>%
  select(annualrain, contains('absmax')) %>%
  dplyr::rename(Biomass = absmax_biomass,
                Density = absmax_d100,
                Length = absmax_length,
                Diversity = absmax_shannon) %>%
  pivot_longer(cols=c(Biomass, Density, Length, Diversity),
               names_to='metric', values_to = 'y_val') %>%
  ggplot(aes(x=annualrain, y=y_val)) +
  facet_wrap(~metric,scales='free') +
  geom_smooth(method = "lm", se=FALSE, color="pink3", linetype=2,
              formula =  y ~ poly(x, 2)) +
  geom_smooth(method = "lm", se=FALSE, color="red", linetype=1,
              formula =  y ~ x) +
  geom_point(shape = 21, size = 4, alpha = .5, fill='grey50') +
  stat_cor(label.y = 6, size = 5)+ 
  stat_regline_equation(label.y = 7, size = 5) +
  ylab('Response Ratio') +
  xlab('Annual Rainfall (cm)') +
  ggtitle('Response Ratio: baseline = 2020 annual average')

LRR_plot_20base <- df %>% 
  filter(collection_date < ymd('2017-12-01')) %>%
  group_by(annualrain) %>%
  dplyr::summarize(Biomass = max(abs(LRR_biomass)),
                Density = max(abs(LRR_d100)),
                Length = max(abs(LRR_length)),
                Diversity = max(abs(LRR_shannon))) %>%
  pivot_longer(cols=c(Biomass, Density, Length, Diversity),
               names_to='metric', values_to = 'y_val') %>%
  ggplot(aes(x=annualrain, y=y_val)) +
  facet_wrap(~metric) +
  geom_smooth(method = "lm", se=FALSE, color="pink3", linetype=2,
              formula =  y ~ poly(x, 2)) +
  geom_smooth(method = "lm", se=FALSE, color="red", linetype=1,
              formula =  y ~ x) +
  geom_point(shape = 21, size = 4, alpha = .5) +
  stat_cor(label.y = 6, size = 5)+ 
  stat_regline_equation(label.y = 7, size = 5) +
  ylab('Max Response Ratio') +
  xlab('Annual Rainfall (cm)') +
  ggtitle('Max AbsoluteV Log Response Ratio: baseline = 2020 annual average')

RI_plot_20base <- df %>% 
  select(annualrain, contains('RI')) %>%
  dplyr::rename(Biomass = RI_biomass,
                Density = RI_d100,
                Length = RI_length,
                Diversity = RI_shannon) %>%
  pivot_longer(cols=c(Biomass, Density, Length, Diversity),
               names_to='metric', values_to = 'y_val') %>%
  ggplot(aes(x=annualrain, y=y_val)) +
  facet_wrap(~metric) +
  geom_smooth(method = "lm", se=FALSE, color="pink3", linetype=2,
              formula =  y ~ poly(x, 2)) +
  geom_smooth(method = "lm", se=FALSE, color="red", linetype=1,
              formula =  y ~ x) +
  geom_point(shape = 21, size = 4, alpha = .5, fill='grey50') +
  stat_cor(label.y = 140, size = 5)+ 
  stat_regline_equation(label.y = 150, size = 5) +
  ylim(0,175) +
  ylab('Days to Recovery') +
  xlab('Annual Rainfall (cm)') +
  ggtitle('Recovery Interval: days to 90% CL of baseline (2020 annual average)')

#------------------------------------------------------------------------------
# export
#------------------------------------------------------------------------------

write_csv(df, 'Data/response_v_2020_baseline.csv')

ggsave('Figures/Response_Ratios_2020base.pdf',
       plot = RR_plot_20base,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/LRR_2020base.pdf',
       plot = LRR_plot_20base,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/Recovery_Intervals_2020base.pdf',
       plot = RI_plot_20base,
       width = 9,
       height = 9,
       units = c("in"))

# End response_v_2020_baseline_calc


