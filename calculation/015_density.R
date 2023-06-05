# density_figures
# Sean Kinard
# 1-16-2023

#------------------------------------------------------------------------------
# Calculations: d100 (fish/100m2)
#------------------------------------------------------------------------------

# density = total number of individuals per m2 on this date
d <- read_csv('Data/fish_biomass_megaframe.csv') %>% 
  select(site_code, collection_date,
         genus, species, density) %>%
  unique()

# species
my_sp <- read_csv('Data/my_fish_species.csv')

# make d100 = number of individuals (rounded up) per 100 m2 on this date
d <- d %>%
  mutate(d100 = ceiling(density*100)) %>%
  select( ! density)

# Add taxonomic information
d <- d %>%
  unite('genus_species', genus:species, sep=' ') %>%
  left_join(my_sp %>% select(! comments))

# Sum d100 by community , taxonomic family, taxonomic species

my_sums <- function(my_data) {
  
d_com <- my_data %>%
  group_by(site_code, collection_date) %>%
  dplyr::summarize(d100_com = sum(d100)) %>% ungroup()

d_fam <- my_data %>%
  group_by(site_code, collection_date, family) %>%
  dplyr::summarize(d100_fam = sum(d100)) %>% ungroup()

d_spe <- my_data %>%
  group_by(site_code, collection_date, lowest_taxon) %>%
  dplyr::summarize(d100_spe = sum(d100)) %>%
  filter(d100_spe>0) %>%
  left_join(select(my_data, family, lowest_taxon) %>% unique())

my_output <- d_spe %>%
  left_join(d_fam) %>%
  left_join(d_com) %>%
  select(site_code, collection_date, family, 
         lowest_taxon, d100_com, d100_fam, d100_spe)

return(my_output) }

d <- my_sums(d)

#------------------------------------------------------------------------------
# Export fish_d100
#------------------------------------------------------------------------------

write_csv(d, 'Data/fish_d100.csv')
  
#------------------------------------------------------------------------------
# Prep data for figures
#------------------------------------------------------------------------------

# Format time variables
my_date_vars <- function(my_data) {
  my_data %>%
  mutate(cd2 = collection_date) %>%
  separate(cd2, into = c('year', 'month', 'day')) %>%
  mutate(year=as.numeric(year),
         month=as.numeric(month),
         day=as.numeric(day)) %>%
  mutate(nday=month*30+day) %>%
  mutate(qtr = ifelse(month < 4, 'Q1',
               ifelse(month > 3 & month < 7, 'Q2',
               ifelse(month > 6 & month < 10, 'Q3', 'Q4')))) }

d <- d %>% my_date_vars()

# Add some environmental predictors
d <- d %>%
  left_join(read_csv('Data/site_longterm_nogap.csv') ) %>%
  left_join(read_csv('Data/site_daily_nogap.csv') %>% select( ! qtr) )

#------------------------------------------------------------------------------
# total density vs site
#------------------------------------------------------------------------------

density_com_vs_site <-  d %>%
  mutate(site_code = as.factor(site_code)) %>%
  mutate(site_code = fct_relevel(site_code, site_order)) %>%
  ggplot(aes(x = site_code, y = d100_com, fill = annualrain)) +
  facet_wrap(~year, ncol=1) +
  geom_boxplot() +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  theme_dark(base_size=14) +
  ylab('Density per 100m2') +
  scale_y_log10() +
  labs(fill = 'Rain (cm/yr)') +
  xlab(element_blank())

#------------------------------------------------------------------------------
# total density vs time
#------------------------------------------------------------------------------
base_density <- function(my_data, my_x, my_y) {
  
  my_data$x_var <- my_data %>% pull(my_x)
  my_data$y_var <- my_data %>% pull(my_y)
  
  my_data %>%
    mutate(site_code = as.factor(site_code)) %>%
    mutate(site_code = fct_relevel(site_code, site_order)) %>%
    ggplot(aes(x = x_var, y = y_var, fill = annualrain)) +
    geom_smooth(aes(color = site_code), method = "loess", se=FALSE,
                alpha=.6, span=1.5) +
    geom_point(aes(fill=annualrain), shape = 21, size = 4, alpha = .8) +
    geom_smooth(method = "loess", se=FALSE, color = 'grey30',  alpha=.5, 
                span=1.5) +
    scale_color_manual(values = my_colors) +
    paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
    theme_dark(base_size=14) +
    ylab('Density per 100m2') +
    scale_y_log10() +
    labs(fill = 'Rain\n(cm/yr)',
         color = element_blank())
}

density_com_vs_time <- d %>% 
  base_density('nday', 'd100_com') +
  xlab('Day of the Year') +
  facet_wrap(~year, ncol=1)

density_com_vs_qmax <- d %>% 
  base_density('q_2w_max', 'd100_com') +
  xlab('Max Flow 2 Weeks Prior') +
  scale_x_log10() +
  facet_wrap(~year, ncol=1)

density_com_vs_cond <- d %>% 
  base_density('conductivity', 'd100_com') +
  xlab('Conductivity') +
  scale_x_log10() +
  facet_wrap(~year, ncol=1)

density_com_vs_nitrate <- d %>% 
  base_density('no3n', 'd100_com') +
  xlab('Nitrate') +
  scale_x_log10() +
  facet_wrap(~year, ncol=1)

density_com_vs_green_algae <- d %>% 
  base_density('green_algae', 'd100_com') +
  xlab('green_algae') +
  scale_x_log10() +
  facet_wrap(~year, ncol=1)

density_com_vs_diatoms <- d %>% 
  base_density('diatoms', 'd100_com') +
  xlab('diatoms') +
  scale_x_log10() +
  facet_wrap(~year, ncol=1)

#------------------------------------------------------------------------------
# Taxa-specific densities vs time
#------------------------------------------------------------------------------
my_families <- c('Centrarchidae', 'Chichlidae', 'Lepisosteidae', 'Leuciscidae', 'Poeciliidae')

my_centrarchids <- c('L. auritus',  "L. macrochirus", "L. megalotis", 
                     "L. cyanellus", "L. gulosus")

small_gape <- c('H. cyanoguttatum', 'L. auritus',  "L. macrochirus", "L. megalotis")

large_gape <- c("L. cyanellus", "L. gulosus", "M. salmoides")

density_fam_vs_time <- d %>% 
  filter(family %in% my_families) %>%
  base_density('nday', 'd100_fam') +
  xlab('Time') +
  facet_wrap(~family+year, ncol=3, scales= "free_y") +
  theme(legend.position = 'none') +
  scale_y_continuous()

density_poec_vs_time <- d %>% 
  filter(family == 'Poeciliidae') %>%
  base_density('nday', 'd100_spe') +
  xlab('Time') +
  facet_wrap(~lowest_taxon+year, ncol=3) +
  theme(legend.position = 'none') +
  scale_y_continuous() +
  ylim(c(0,415))

density_sgap_vs_time <- d %>% 
  filter(lowest_taxon %in% small_gape) %>%
  base_density('nday', 'd100_spe') +
  xlab('Time') +
  facet_wrap(~lowest_taxon+year, ncol=3) +
  theme(legend.position = 'none')

density_lgap_vs_time <- d %>% 
  filter(lowest_taxon %in% large_gape) %>%
  base_density('nday', 'd100_spe') +
  xlab('Time') +
  facet_wrap(~lowest_taxon+year, ncol=3) +
  theme(legend.position = 'none')

#------------------------------------------------------------------------------
# Annual Density (AD)
#------------------------------------------------------------------------------

# base plot function
D_base <- function(mydata) {
  
  mydata %>%
    ggplot(aes(x=annualrain, y = y_var)) +
    geom_smooth(method = "lm", se=FALSE, color="grey20", linetype=2,
                formula =  y ~ poly(x, 2)) +
    geom_smooth(method = "lm", se=FALSE, color="grey20", linetype=1,
                formula =  y ~ x) +
    geom_jitter(shape = 21, size = 4, fill = 'white', alpha = .5)+
    stat_cor(label.y = 2.75, size = 2.5)+ #this means at label.yth unit in the y axis, the r squared and p value will be shown
    stat_regline_equation(label.y = 3, size = 2.5) + #this means at label.yth unit regresion line equation will be shown
    scale_y_log10() +
    ylab('Fish / 100m2') +
    xlab('Annual Rainfall (cm)') +
    theme_dark(base_size = 10) }

density_com_vs_rainfall_annual <- d %>%
  select(! c(d100_fam, d100_spe, family, lowest_taxon)) %>%
  unique() %>%
  dplyr::rename(y_var = d100_com) %>%
  D_base() + facet_wrap(~year)

density_fam_vs_rainfall_annual <- d %>%
  select(! c(d100_spe,lowest_taxon)) %>%
  unique() %>%
  filter(family %in% my_families) %>%
  dplyr::rename(y_var = d100_fam) %>%
  D_base(mydata = ) + facet_wrap(~year+family, ncol=4)

density_cent_vs_rainfall_annual <- d %>%
  filter(lowest_taxon %in% my_centrarchids) %>%
  select(! c(d100_fam,family)) %>%
  unique() %>%
  dplyr::rename(y_var = d100_spe) %>%
  D_base(mydata = ) + facet_wrap(~year+lowest_taxon, ncol=5)

density_poec_vs_rainfall_annual <- d %>%
  filter(family == 'Poeciliidae') %>%
  select(! c(d100_fam,family)) %>%
  unique() %>%
  dplyr::rename(y_var = d100_spe) %>%
  D_base(mydata = ) + facet_wrap(~year+lowest_taxon, ncol=2)

caption_density <- 'Fish abundance (indviduals/100m2) plotted against Rainfall (cm/yr). The vertical axis is log-transformed. Data comes from monthly  sampling in 2017 and 2018, and quarterly sampling in 2020.'

#------------------------------------------------------------------------------
# Quarterly Density (QD)
#------------------------------------------------------------------------------
density_com_vs_rainfall_quarterly <- d %>%
  select(! c(d100_fam, d100_spe, family, lowest_taxon)) %>%
  unique() %>%
  D_base() +
  facet_wrap(~year+qtr)



#------------------------------------------------------------------------------
# Clean-up
#------------------------------------------------------------------------------

rm(d, my_sp)

#------------------------------------------------------------------------------
# End density_figures