# site_explore
# Sean Kinard
# 1-19-2023

# -----------------------------------------------------------------------------
# setup
# -----------------------------------------------------------------------------

# load data
d <- read_csv('Data/fish_biomass_megaframe.csv')
d_meta <- read_csv('Data/fish_rapid_terrg_megaframe_metadata.csv')

my_colors <- c('#F5191CFF', '#E78200FF', '#E8A117FF',
               '#EABB22FF', '#CBC988FF', '#9FC095FF',
               '#6BB699FF', '#3DAAA6FF', '#3B99B1FF')

site_order <- d %>% select(annualrain, site_code) %>% arrange(annualrain) %>% unique() %>% pull(site_code)
colnames(d) <- str_replace_all(colnames(d), '_mean' , '')

# -----------------------------------------------------------------------------
# Data Extraction
# -----------------------------------------------------------------------------

# extract columns
site_info <- d %>% 
  select(site_code, site_name, staid, lat, lon) %>% 
  distinct()

lte_watershed <- d %>% 
  select(site_code,
         basinsize, developedland, forestland, cropland, otherland) %>% 
  distinct()

lte_flow <- d %>% 
  select(site_code,
         flsh, hfpp3, lfpp, meanflow, medianflow, seasonality) %>% 
  distinct()

daily_flow <- d %>% 
  select(site_code, annualrain, collection_date, 
         q, q_2w_mn, q_2w_max, q_2w_min) %>% 
  distinct()

daily_algae <- d %>% 
  select(site_code, annualrain, collection_date, 
         bluegreen_cyano, diatoms, green_algae) %>% 
  distinct()

daily_channel <- d %>% 
  select(site_code, annualrain, collection_date,
         depth_mx, width, gravel, 
         silt, canopy_density_mid) %>% 
  distinct()

daily_chemistry <- d %>% 
  select(site_code, annualrain, collection_date, 
         conductivity, do_mg_l, 
         nh4_n, no3n, ortho_p) %>% 
  distinct()

# -----------------------------------------------------------------------------
# lte scatterplots
# -----------------------------------------------------------------------------

scatter_plot <- function(my_data) {

  my_data %>%
    mutate(site_code = fct_relevel(site_code, site_order)) %>%
    pivot_longer(cols = 2:length(colnames(my_data)),
                 names_to = "variable",
                 values_to = "value") %>%
    left_join(d %>% select(site_code, annualrain) %>% distinct()) %>%
    mutate(site_code = fct_relevel(site_code, site_order)) %>%
    ggplot(aes(x=annualrain, y=value)) +
    facet_wrap(~variable, scale = 'free', ncol=2) +
    geom_smooth(method = "loess", se=FALSE, color="grey30", linetype=2,
                size=1.2, span=.8) +
    geom_point(aes(fill=site_code), size = 6, shape=21) +
    scale_fill_manual(values = my_colors) +
    labs(fill='Site') +
    xlab('Rain (cm/yr)') +
    ylab(element_blank()) +
    theme_dark(base_size=14) }

scatter_watershed <- scatter_plot(lte_watershed)
scatter_flow <- scatter_plot(lte_flow)

scatter_lte <- scatter_watershed / 
  scatter_flow + 
  plot_layout(guides='collect', heights = c(2,3))

# -----------------------------------------------------------------------------
# daily boxplots
# -----------------------------------------------------------------------------

explot2 <- function(my_data) {
  
  my_data %>%
    mutate(site_code = fct_relevel(site_code, site_order)) %>%
    pivot_longer(cols = 4:length(colnames(my_data)),
                 names_to = "variable",
                 values_to = "value") %>%
    ggplot(aes(x=site_code, y=value, fill=site_code)) +
    facet_wrap(~variable, scale = 'free', ncol=2) +
    geom_boxplot(show.legend = F) +
    ylab(element_blank())+
    scale_fill_manual(values = my_colors) +
    labs(fill='Site') +
    xlab(element_blank()) +
    ylab(element_blank()) +
    theme_dark(base_size=14) }

boxplot_flow <- explot2(daily_flow) + scale_y_log10()
boxplot_channel <- explot2(daily_channel) 
boxplot_algae <- explot2(daily_algae) 
boxplot_chemistry <- explot2(daily_chemistry)

# -----------------------------------------------------------------------------
# time series
# -----------------------------------------------------------------------------

time_series <- function(my_data) {
  
  d_new <- my_data %>%
    mutate(site_code = fct_relevel(site_code, site_order)) %>%
    pivot_longer(cols = 4:length(colnames(my_data)),
                 names_to = "variable",
                 values_to = "value")
  
  ggplot(data = d_new %>% filter(collection_date < ymd('2020-01-01')),
         aes(x=collection_date, y=value, fill= annualrain)) +
    facet_wrap(~variable, scale = 'free', ncol=2) +
    geom_point(size=3, color = 'black', shape =21, alpha=.5,
               show.legend = F) +
    geom_smooth(aes(color = site_code),
                method = "loess", size = 1.3, se=FALSE, linetype=1, span=1.1,
                show.legend = F) +
    scale_color_manual(values = my_colors) +
    paletteer::scale_fill_paletteer_c("grDevices::Zissou 1",
                                      direction=-1) +
    ylab(element_blank()) +
    xlab(element_blank()) +
    labs(fill='Rain\n(cm/yr)') +
    theme_dark(base_size=14) +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))
}

timeline_flow <- time_series(daily_flow) + 
  scale_y_log10()

timeline_channel <- time_series(daily_channel) 
timeline_algae <- time_series(daily_algae)  + scale_y_log10()
timeline_chemistry <- time_series(daily_chemistry) 

# -----------------------------------------------------------------------------
# RAPID time series (SCALED)
# -----------------------------------------------------------------------------

scaling <- function(my_data) {
  
  x <- my_data %>%
    filter(collection_date < ymd('2020-01-01')) %>%
    select(! annualrain)
  
  old_cols <- colnames(x)
  
  for (i in 3:length(colnames(x))) {
    
    new_col <- paste(colnames(x[i]), 
                     'scaled',
                     sep='_')
    
    x <- x %>%
      select(1,2,i) %>%
      pivot_wider(names_from = site_code, 
                  values_from = 3) %>%
      mutate(AR = scale(AR),
             EM = scale(EM),
             GC = scale(GC),
             MR = scale(MR),
             PD = scale(PD),
             PL = scale(PL),
             SF = scale(SF),
             TR = scale(TR),
             WM = scale(WM) ) %>%
      pivot_longer(cols= AR:WM,
                   names_to = "site_code",
                   values_to = new_col,
                   values_drop_na = T) %>%
      right_join(x) %>%
      select(all_of(old_cols), new_col)
    
    old_cols <- colnames(x)
  }
  
  x <- select(x, site_code, collection_date, contains('scaled'))
  
  x <- my_data %>% select(site_code, annualrain) %>% distinct() %>%
    right_join(x)
  
  return(x) } 

timeline_flow_scaled <- time_series(daily_flow %>% scaling() )
timeline_channel_scaled <- time_series(daily_channel %>% scaling()) 
timeline_algae_scaled <- time_series(daily_algae %>% scaling())
timeline_chemistry_scaled <- time_series(daily_chemistry %>% scaling()) 

# -----------------------------------------------------------------------------
# Algae vs flow or Nutrients, 
# -----------------------------------------------------------------------------
base_algae <- function(my_data, my_variable) {
  
  my_data$my_column <- my_data %>% pull(my_variable)
  
  my_data %>%
    right_join(daily_algae) %>%
    pivot_longer(cols = c(bluegreen_cyano, diatoms, green_algae),
                 names_to = 'algae_type',
                 values_to = 'algae_value') %>%
    filter(collection_date < ymd('2020-01-01')) %>%
    mutate(site_code = fct_relevel(site_code, site_order)) %>%
    ggplot(aes(my_column, algae_value)) +
    facet_wrap(~algae_type, scale = 'free', ncol=1) +
    geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
    geom_point(aes(fill=annualrain), shape = 21, size = 4, alpha = .8)+
    stat_cor(label.y = 1.5, size = 3)+ 
    stat_regline_equation(label.y = 2, size = 3) + 
    scale_color_manual(values = my_colors) +
    paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
    scale_x_log10() +
    ylab('Chlorophyll') +
    theme_dark(base_size=14) 
    
}

algae_v_flow <- base_algae(daily_flow, 'q_2w_max') +
  xlab('Max Discharge in cfs (2 weeks prior)')
  
algae_v_nitrate <- base_algae(daily_chemistry, 'no3n')  +
  xlab('Nitrate')

algae_v_phosphate <- base_algae(daily_chemistry, 'ortho_p')  +
  xlab('Phosphate')


# -----------------------------------------------------------------------------
# Clean-up
# -----------------------------------------------------------------------------
rm(d, d_meta, daily_algae, daily_channel, daily_chemistry, daily_flow,
    lte_flow, lte_watershed, site_info)
