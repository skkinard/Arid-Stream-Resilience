# flow_record_2017_to_2020
# Sean Kinard
# 2023-02-15

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
rm(list=ls())

setwd('/home/kinard/Documents/Research/Dissertation/02_Resilience')

library(waterData)
library(tidyverse)
library(lubridate)

# load site data
d_site <- read_csv('Data/site_data.csv') %>%
  mutate(STAID = paste('0', STAID, sep=''),
         Site = substr(Site,1,2))
colnames(d_site) <- str_to_lower(colnames(d_site))

site_order <- c("TR", "SF", "AR", "MR", "PD", "PL", "GC", "WM", "EM")

my_colors <- c('#F5191CFF', '#E78200FF', '#E8A117FF',
               '#EABB22FF', '#CBC988FF', '#9FC095FF',
               '#6BB699FF', '#3DAAA6FF', '#3B99B1FF')
# -----------------------------------------------------------------------------
# Extract flow data from 2017, 2018, and 2020 for sites
# -----------------------------------------------------------------------------

my_staids <- pull(d_site, staid)

# create empty tibble
d_flow <- tibble(
  staid = character(),
  val = numeric(),
  dates = date(),
  site = character() ) %>%
  mutate(dates = as_date(dates))

for (i in 1:length(my_staids)) {
  # https://help.waterdata.usgs.gov/codes-and-parameters/parameters
  
  d_flow <- importDVs(my_staids[i], code = "00060", stat = "00003", 
                   sdate = "2017-01-01", edate = '2020-12-31') %>%
    cleanUp( task = "fix", replace = 0.01) %>% # replace negatives & zeros
    fillMiss(block = 30, pmiss = 40, model = "trend", smooth = TRUE) %>% #fillgaps
    as_tibble() %>% # formatting
    select(! qualcode) %>%
    left_join(d_site%>%select(staid, site)) %>% # adding site names
    mutate(dates = as_date(dates)) %>%
    full_join(d_flow) # merging
  } 

# Formatting and units
d_flow <- d_flow %>% 
  select(! staid) %>%
  mutate(val = val*28.31685) %>% # convert ft3/s to liter/s
  rename(collection_date = dates,
         q = val,
         site_code = site) %>%
  mutate(site_code = fct_relevel(site_code, site_order))

# Date vectors
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

d_flow <- d_flow %>% my_date_vars()

# -----------------------------------------------------------------------------
# Flow Record Figures
# -----------------------------------------------------------------------------
flow_base <- function(my_data) {
  my_data %>%
    ggplot(aes(collection_date, q)) +
    geom_line(color = 'skyblue2') +
    theme_dark(base_size = 12) +
    scale_x_date(date_labels = "%b") +
    scale_y_log10() +
    ylab('Average Daily Discharge (liters/s') +
    xlab(element_blank())
}

# All flow records
flow_record_all_2017_to_2020 <- d_flow %>%
  flow_base() +
  facet_wrap(~site_code * year, ncol=4, scales='free')

# year specific plots
my_years = c(2017, 2018, 2019, 2020)

p <- list()

for(i in 1:length(my_years)) {
  p[[i]] <- d_flow %>%
    filter(near(year, my_years[i])) %>%
    flow_base() +
    facet_wrap(~site_code, ncol=1, scales='free') }

flow_record_2017 <- p[[1]]
flow_record_2018 <- p[[2]]
flow_record_2019 <- p[[3]]
flow_record_2020 <- p[[4]]
# -----------------------------------------------------------------------------
# Export flow records
# -----------------------------------------------------------------------------
write_csv(d_flow, 'Data/flow_record_allsites_2017_to_2020.csv')

ggsave('Figures/flow_record_all_2017_to_2020.pdf',
       plot = flow_record_all_2017_to_2020,
       width = 9,
       height = 18,
       units = c("in"))

ggsave('Figures/flow_record_2017.pdf',
       plot = flow_record_2017,
       width = 9,
       height = 18,
       units = c("in"))

ggsave('Figures/flow_record_2018.pdf',
       plot = flow_record_2018,
       width = 9,
       height = 18,
       units = c("in"))

ggsave('Figures/flow_record_2019.pdf',
       plot = flow_record_2019,
       width = 9,
       height = 18,
       units = c("in"))

ggsave('Figures/flow_record_2020.pdf',
       plot = flow_record_2020,
       width = 9,
       height = 18,
       units = c("in"))

# -----------------------------------------------------------------------------

# End flow_record_2017_to_2020