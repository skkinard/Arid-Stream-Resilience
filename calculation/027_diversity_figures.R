# diversity_figures
# Sean Kinard
# 1-26-2023

#------------------------------------------------------------------------------
# Prep data for figures
#------------------------------------------------------------------------------

# load data
d <- read_csv('Data/RTC_fish_diversity.csv')

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

my_evars <- c('site_code', 'collection_date', 'annualrain', 'q_2w_max',
              'diatoms', 'green_algae', 'no3n', 'ortho_p', 'conductivity', 'depth_mx')

d <- d %>%
  left_join(read_csv('Data/site_longterm_nogap.csv')%>%select(any_of(my_evars)) ) %>%
  left_join(read_csv('Data/site_daily_nogap.csv')%>%select(any_of(my_evars)) )

# pivot longer for diversity metrics
d <- d %>%
  pivot_longer(cols=richness:simpson, names_to = 'd_type', 
               values_to = 'd_value') %>%
  mutate(d_type = str_to_title(d_type),
         d_type = fct_relevel(d_type, 'Richness', 'Shannon', 'Simpson'))

# Change factor levels of site_code to annual rainfall gradient
d <- d %>%
  mutate(site_code = as.factor(site_code)) %>%
  mutate(site_code = fct_relevel(site_code, site_order))

#------------------------------------------------------------------------------
# diversity vs site
#------------------------------------------------------------------------------

div_vs_site <-  d %>%
  ggplot(aes(x = site_code, y = d_value, fill = annualrain)) +
  facet_wrap(~year+d_type, ncol=3) +
  geom_boxplot() +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  theme_dark(base_size=10) +
  ylab('Diversity') +
  labs(fill = 'Rain (cm/yr)') +
  xlab(element_blank())

#------------------------------------------------------------------------------
# diversity vs time
#------------------------------------------------------------------------------
base_div_value <- function(my_data, my_x, my_y) {
  
  my_data$x_var <- my_data %>% pull(my_x)
  my_data$y_var <- my_data %>% pull(my_y)
  
  my_data %>%
    mutate(site_code = as.factor(site_code)) %>%
    mutate(site_code = fct_relevel(site_code, site_order)) %>%
    ggplot(aes(x = x_var, y = y_var, fill = annualrain)) +
    geom_smooth(aes(color = site_code), method = "loess", se=FALSE,
                alpha=.6, span=1.2) +
    geom_point(aes(fill=annualrain), shape = 21, size = 2, alpha = .8) +
    geom_smooth(method = "loess", se=FALSE, color = 'grey30',  alpha=.5,
                span=1.5) +
    scale_color_manual(values = my_colors) +
    paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
    theme_dark(base_size=10) +
    labs(fill = 'Rain\n(cm/yr)',
         color = element_blank()) +
    facet_wrap(~d_type+year, ncol=3) +
    ylab('Diversity')
    
}

div_vs_time <- d %>%
  base_div_value('nday', 'd_value') +
  xlab('Day of the Year')

div_vs_qmax <- d %>%
  base_div_value('q_2w_max', 'd_value') +
  xlab('Max Flow 2 Weeks Prior') +
  scale_x_log10()

div_vs_cond <- d %>%
  base_div_value('conductivity', 'd_value') +
  xlab('Conductivity') +
  scale_x_log10()

div_vs_nitrate <- d %>%
  base_div_value('no3n', 'd_value') +
  xlab('Nitrate') +
  scale_x_log10()

div_vs_green_algae <- d %>%
  base_div_value('green_algae', 'd_value') +
  xlab('Green Algae') +
  scale_x_log10() 

div_vs_diatoms <- d %>%
  base_div_value('diatoms', 'd_value') +
  xlab('Diatoms') +
  scale_x_log10() 

#------------------------------------------------------------------------------
# Annual d_value (AD)
#------------------------------------------------------------------------------

# base plot function

div_vs_rainfall_annual <- d %>%
  ggplot(aes(x=annualrain, y = d_value)) +
  geom_smooth(method = "lm", se=FALSE, color="grey20", linetype=2,
              formula =  y ~ poly(x, 2)) +
  geom_smooth(method = "lm", se=FALSE, color="grey20", linetype=1,
              formula =  y ~ x) +
  geom_jitter(shape = 21, size = 2, fill = 'white', alpha = .5) +
  stat_cor(size = 3.5) + 
  #stat_regline_equation(size = 2.5) +
  ylab('Diversity') +
  xlab('Annual Rainfall (cm)') +
  theme_dark(base_size = 10) +
  facet_wrap(~d_type+year, ncol=3, scales='free')

caption_div <- 'Shannon-Wiener Diversity plotted against time. RAPID data comes from monthly sampling following hurricane Harvey from September 2017 through September 2018. TERRG data includes quarterly sampling (January, May, September, and November) in 2020.'

#------------------------------------------------------------------------------
# Quarterly d_value (QD)
#------------------------------------------------------------------------------
d_qtrly <- function(mydata) {
  mydata %>%
  ggplot(aes(x=annualrain, y = d_value)) +
    geom_smooth(method = "lm", se=FALSE, color="grey20", linetype=2,
                formula =  y ~ poly(x, 2)) +
    geom_smooth(method = "lm", se=FALSE, color="grey20", linetype=1,
                formula =  y ~ x) +
    geom_jitter(shape = 21, size = 2, fill = 'white', alpha = .5) +
    stat_cor(size = 3.5) + 
    #stat_regline_equation(size = 2.5) +
    ylab('Shannon Diversity') +
    xlab('Annual Rainfall (cm)') +
    theme_dark(base_size = 10) +
    facet_wrap(~year+qtr, ncol=4, scales='free')
}


richness_vs_rainfall_quarterly <- d %>%
  filter(d_type == 'Richness') %>%
  d_qtrly()

shannon_vs_rainfall_quarterly <- d %>%
  filter(d_type == 'Shannon') %>%
  d_qtrly()
  

#------------------------------------------------------------------------------
# Clean-up
#------------------------------------------------------------------------------

rm(d)
