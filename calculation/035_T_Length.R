# T_Length_figures
# Sean Kinard
# 2023-02-09

#------------------------------------------------------------------------------
# Calculations: L_mu (average length)
#------------------------------------------------------------------------------

d <- read_csv('Data/forklength_all_megaframe.csv') %>% 
  select(site_code, collection_date,
         genus, species, lengthmm) %>%
  filter(lengthmm>0) %>%
  unique()

# species
my_sp <- read_csv('Data/my_fish_species.csv')

# Add taxonomic information
d <- d %>%
  unite('genus_species', genus:species, sep=' ') %>%
  left_join(my_sp %>% select(! comments))

# make L_mu for community , taxonomic family, taxonomic species

my_avg <- function(my_data) {
  
  L_com <- my_data %>%
    group_by(site_code, collection_date) %>%
    dplyr::summarize(L_mu_com = mean(lengthmm),
                     L_me_com = median(lengthmm),
                     L_sd_com = sd(lengthmm),
                     L_n_com = length(lengthmm)) %>% ungroup()
  
  L_fam <- my_data %>%
    group_by(site_code, collection_date, family) %>%
    dplyr::summarize(L_mu_fam = mean(lengthmm),
                     L_me_fam = median(lengthmm),
                     L_sd_fam = sd(lengthmm),
                     L_n_fam = length(lengthmm)) %>% ungroup()
  
  L_spe <- my_data %>%
    group_by(site_code, collection_date, lowest_taxon) %>%
    dplyr::summarize(L_mu_spe = mean(lengthmm),
                     L_me_spe = median(lengthmm),
                     L_sd_spe = sd(lengthmm),
                     L_n_spe = length(lengthmm)) %>%
    filter(L_mu_spe>0) %>%
    left_join(select(my_data, family, lowest_taxon) %>% unique())
  
  my_output <- L_spe %>%
    left_join(L_fam) %>%
    left_join(L_com)
  
  return(my_output) }

d <- my_avg(d)

#------------------------------------------------------------------------------
# Export fish_L_mu
#------------------------------------------------------------------------------

write_csv(d, 'Data/fish_Length_stats.csv')

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
# T_Length vs site
#------------------------------------------------------------------------------

T_Length_com_vs_site <-  d %>%
  mutate(site_code = as.factor(site_code)) %>%
  mutate(site_code = fct_relevel(site_code, site_order)) %>%
  ggplot(aes(x = site_code, y = L_mu_com, fill = annualrain)) +
  facet_wrap(~year, ncol=1) +
  geom_boxplot() +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  theme_dark(base_size=14) +
  ylab('T_Length (mm)') +
  labs(fill = 'Rain (cm/yr)') +
  xlab(element_blank())

#------------------------------------------------------------------------------
# T_Length vs time
#------------------------------------------------------------------------------
base_T_Length <- function(my_data, my_x, my_y) {
  
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
    ylab('T_Length (mm)') +
    labs(fill = 'Rain\n(cm/yr)',
         color = element_blank())
}

T_Length_com_vs_time <- d %>% 
  base_T_Length('nday', 'L_mu_com') +
  xlab('Day of the Year') +
  facet_wrap(~year, ncol=1)

T_Length_com_vs_qmax <- d %>% 
  base_T_Length('q_2w_max', 'L_mu_com') +
  xlab('Max Flow 2 Weeks Prior') +
  scale_x_log10() +
  facet_wrap(~year, ncol=1)

T_Length_com_vs_cond <- d %>% 
  base_T_Length('conductivity', 'L_mu_com') +
  xlab('Conductivity') +
  scale_x_log10() +
  facet_wrap(~year, ncol=1)

T_Length_com_vs_nitrate <- d %>% 
  base_T_Length('no3n', 'L_mu_com') +
  xlab('Nitrate') +
  scale_x_log10() +
  facet_wrap(~year, ncol=1)

T_Length_com_vs_green_algae <- d %>% 
  base_T_Length('green_algae', 'L_mu_com') +
  xlab('green_algae') +
  scale_x_log10() +
  facet_wrap(~year, ncol=1)

T_Length_com_vs_diatoms <- d %>% 
  base_T_Length('diatoms', 'L_mu_com') +
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

T_Length_fam_vs_time <- d %>% 
  filter(family %in% my_families) %>%
  base_T_Length('nday', 'L_mu_fam') +
  xlab('Time') +
  facet_wrap(~family+year, ncol=3, scales= "free_y") +
  theme(legend.position = 'none') +
  scale_y_continuous()

T_Length_poec_vs_time <- d %>% 
  filter(family == 'Poeciliidae') %>%
  mutate(lowest_taxon = fct_relevel(
    lowest_taxon, 'G. affinis', 'P. latipinna', 'P. formosa')) %>%
  base_T_Length('nday', 'L_mu_spe') +
  xlab('Time') +
  facet_wrap(~lowest_taxon+year, ncol=3) +
  theme(legend.position = 'none') +
  scale_y_continuous() 

T_Length_sgap_vs_time <- d %>% 
  filter(lowest_taxon %in% small_gape) %>%
  base_T_Length('nday', 'L_mu_spe') +
  xlab('Time') +
  facet_wrap(~lowest_taxon+year, ncol=3) +
  theme(legend.position = 'none')

T_Length_lgap_vs_time <- d %>% 
  filter(lowest_taxon %in% large_gape) %>%
  base_T_Length('nday', 'L_mu_spe') +
  xlab('Time') +
  facet_wrap(~lowest_taxon+year, ncol=3) +
  theme(legend.position = 'none')

#------------------------------------------------------------------------------
# Annual T_Length (AD)
#------------------------------------------------------------------------------

# base plot function
L_base <- function(mydata) {
  
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
    ylab('T_length (mm)') +
    xlab('Annual Rainfall (cm)') +
    theme_dark(base_size = 10) }

T_Length_com_vs_rainfall_annual <- d %>%
  select(! c(L_mu_fam, L_mu_spe, family, lowest_taxon)) %>%
  unique() %>%
  dplyr::rename(y_var = L_mu_com) %>%
  L_base() + facet_wrap(~year)

T_Length_fam_vs_rainfall_annual <- d %>%
  select(! c(L_mu_spe,lowest_taxon)) %>%
  unique() %>%
  filter(family %in% my_families) %>%
  dplyr::rename(y_var = L_mu_fam) %>%
  L_base() + facet_wrap(~year+family, ncol=4)

T_Length_cent_vs_rainfall_annual <- d %>%
  filter(lowest_taxon %in% my_centrarchids) %>%
  select(! c(L_mu_fam,family)) %>%
  unique() %>%
  dplyr::rename(y_var = L_mu_spe) %>%
  L_base() + facet_wrap(~year+lowest_taxon, ncol=5)

T_Length_poec_vs_rainfall_annual <- d %>%
  filter(family == 'Poeciliidae') %>%
  select(! c(L_mu_fam,family)) %>%
  unique() %>%
  dplyr::rename(y_var = L_mu_spe) %>%
  L_base() + facet_wrap(~year+lowest_taxon, ncol=2)

caption_T_Length <- 'Fish total length (mm) plotted against Rainfall (cm/yr). The vertical axis is log-transformed. Data comes from monthly sampling in 2017 and 2018, and quarterly sampling in 2020.'

#------------------------------------------------------------------------------
# Quarterly T_Length (QD)
#------------------------------------------------------------------------------
T_Length_com_vs_rainfall_quarterly <- d %>%
  select(! c(L_mu_fam, L_mu_spe, family, lowest_taxon)) %>%
  unique() %>%
  dplyr::rename(y_var = L_mu_com) %>%
  L_base() +
  facet_wrap(~year+qtr)

#------------------------------------------------------------------------------
# Clean-up
#------------------------------------------------------------------------------

rm(d, my_sp)

#------------------------------------------------------------------------------
# End T_Length_figures