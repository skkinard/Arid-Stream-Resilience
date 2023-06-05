# Log Response Ratio Calculation
# Sean Kinard
# 3-16-2023
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
# load data
source('R_Scripts/047_ordination.R')
source('R_Scripts/000_biological_scales_comparison.R')

# -----------------------------------------------------------------------------
# Function
# -----------------------------------------------------------------------------
bas20 <- function(my_data, my_var, my_level) { 
  # ------------------------------------
  # data prep 
  # ------------------------------------
  df <- my_data %>%
    add_time_vars() %>%
    mutate(com='community') %>%
    dplyr::rename(XX={{my_var}}) %>% # XX = variable
    dplyr::rename(YY={{my_level}}) %>% # YY = level of comparison
    select(any_of(c('site_code', 'collection_date', 'qtr',
                    'XX', 'YY'))) %>%
    group_by(site_code, qtr, collection_date, YY) %>%
    dplyr::summarize(XX=sum(XX,na.rm=T)) %>%
    unique() %>%
    ungroup()
  
  # ------------------------------------
  # 2020 quarterly mean
  # ------------------------------------
  bas20 <- df %>%
    filter(collection_date > ymd('2020-01-01')) %>%
    group_by(site_code, qtr, YY) %>%
    dplyr::summarize(XX_bas20=mean(XX,na.rm=T)) %>%
    unique() %>%
    ungroup()

  # ------------------------------------
  # 2020 annual mean
  # ------------------------------------
  mu20 <- df %>%
    filter(collection_date > ymd('2020-01-01')) %>%
    group_by(site_code, YY) %>%
    dplyr::summarize(XX_mu20 = mean(XX, na.rm=T),
                     XX_sd20 = sd(XX, na.rm=T)) %>%
    ungroup() 
 
  # ------------------------------------
  # Response Ratios = (x-x_base) / (x_base)
  # ------------------------------------
  RR <- df %>%
    filter(collection_date < ('2020-01-01')) %>%
    left_join(mu20) %>%
    mutate(XX_mu20 = ifelse(is.na(XX_mu20), 0, XX_mu20)) %>%
    left_join(bas20) %>%
    mutate( # gapfill baselines with 2020 mean
      XX_bas20_fill = ifelse(is.na(XX_bas20), XX_mu20, XX_bas20)) %>%
    mutate( # calculate response ratio (RR)
      XX_RR_qt20 = (XX-XX_bas20_fill) / (XX_bas20_fill + .01),
      XX_RR_mu20 = (XX-XX_mu20) / (XX_mu20 + .01)) %>%
    mutate( # calculate log response ratio (LRR)
      XX_LRR_qt20 = log(abs(XX_RR_qt20)),
      XX_LRR_mu20 = log(abs(XX_RR_mu20))) %>%
    mutate( # replace infinity with -100000 LRR
      XX_LRR_qt20 = ifelse(is.infinite(XX_LRR_qt20), -100000, XX_LRR_qt20),
      XX_LRR_mu20 = ifelse(is.infinite(XX_LRR_mu20), -100000, XX_LRR_mu20)) %>%
    ungroup() %>%
    fix_site_order()
  
  # ------------------------------------
  # Revert Column names
  # ------------------------------------
  colnames(RR) <- str_replace_all(colnames(RR), 'XX', substr(my_var,1,4))
  colnames(RR) <- str_replace_all(colnames(RR), 'YY', my_level)
  
  return(RR) }

# -----------------------------------------------------------------------------
# Calculate Response Ratios
# -----------------------------------------------------------------------------

# Biomass
RR_B100_spe <- d_biomass %>%
  bas20(my_var = 'B100_sum_spe', 
      my_level = 'lowest_taxon')

RR_B100_fam <- d_biomass %>% 
  bas20(my_var = 'B100_sum_spe', 
        my_level = 'family')

RR_B100_com <- d_biomass %>% 
  bas20(my_var = 'B100_sum_spe', 
        my_level = 'com')

# Abundance
RR_d100_spe <- d_density %>%
  bas20(my_var = 'd100_spe', 
        my_level = 'lowest_taxon')

RR_d100_fam <- d_density %>% 
  bas20(my_var = 'd100_spe', 
        my_level = 'family')

RR_d100_com <- d_density %>% 
  bas20(my_var = 'd100_spe', 
        my_level = 'com')

# Diversity
RR_shan_com <- d_diversity %>% 
  bas20(my_var = 'shannon', 
        my_level = 'com')

RR_simp_com <- d_diversity %>% 
  bas20(my_var = 'simpson', 
        my_level = 'com')

RR_rich_com <- d_diversity %>% 
  bas20(my_var = 'richness', 
        my_level = 'com')

# Composition

RR_cdis_com <- d_axis %>%
  filter(collection_date > ymd('2020-01-01')) %>%
  group_by(site_code) %>%
  dplyr::summarize(
    axis1_mu20 = mean(axis1, na.rm=T),
    axis2_mu20 = mean(axis2, na.rm=T) ) %>%
  ungroup() %>%
  right_join(d_axis) %>%
  mutate(
    cdis = sqrt((axis1-axis1_mu20)^2 + 
                  (axis2-axis2_mu20)^2) ) %>%
  bas20(my_var = 'cdis', 
        my_level = 'com')

# -----------------------------------------------------------------------------
# Visualization: Response Vs Baseline (Community)
# -----------------------------------------------------------------------------
p_RR_RvB <- function(my_data, my_variable) {
  
  colnames(my_data) <- str_replace_all(
    colnames(my_data), my_variable, 'X')
  
  categorical_colors <- c( 'cyan', 'red', 'grey90')
  
  my_data %>%
    fix_site_order() %>%
    mutate(upper_lim = X_bas20_fill + X_sd20,
           lower_lim = X_bas20_fill - X_sd20,
           Comparison=case_when(
             X > upper_lim ~ 'Above',
             X < lower_lim ~ 'Below',
             X <= upper_lim & X >= lower_lim ~ 'Inside') ) %>%
    ggplot() +
    facet_wrap(~site_code, scales='free_y') +
    geom_point(aes(x=collection_date, y = X, fill= Comparison),
               shape=21, size = 3, color = 'black', alpha=.3) +
    geom_point(aes(x=collection_date, y = X, color=Comparison),
               shape=21, size = 3) +
    geom_line(aes(x=collection_date, y = upper_lim),
               linetype=2, linewidth=.4, color = 'grey40') +
    geom_line(aes(x=collection_date, y = lower_lim),
                linetype=2, linewidth=.4, color = 'grey40') +
    dark_theme_gray(base_size=12) +
    scale_fill_manual(values=categorical_colors) +
    scale_color_manual(values=categorical_colors) +
    scale_x_date(date_labels = "%m") +
    ylab(my_variable) +
    xlab('Month') +
    theme(legend.position = 'none')
  }

RvB_com_B100 <- p_RR_RvB(RR_B100_com, 'B100') +
  ggtitle('Community Biomass Deviation From 2020 Quarter')

RvB_com_d100 <- p_RR_RvB(RR_d100_com,'d100') +
  ggtitle('Community Abundance Deviation From 2020 Quarter')

RvB_com_shan <- p_RR_RvB(RR_shan_com, 'shan') +
  ggtitle('Shannon Diversity Deviation From 2020 Quarterly')

RvB_com_rich <-p_RR_RvB(RR_rich_com, 'rich') +
  ggtitle('Species Richness Deviation From 2020 Quarterly')

RvB_com_cdis <-p_RR_RvB(RR_cdis_com, 'cdis') +
  ggtitle('RDA Centroid Distance Deviation From 2020 Quarterly')

# -----------------------------------------------------------------------------
# Visualization: Response Vs Baseline (Family)
# -----------------------------------------------------------------------------
common_families <- c("Leuciscidae", "Poeciliidae", "Ictaluridae", "Centrarchidae")

RvB_leuc_d100 <- p_RR_RvB(
  my_data = RR_d100_fam %>% filter(family=='Leuciscidae'),
  my_variable = 'd100') +
  ggtitle('Leuciscidae Density Deviation From 2020 Quarter')

RvB_poec_d100 <- p_RR_RvB(
  my_data = RR_d100_fam %>% filter(family=='Poeciliidae'),
  my_variable = 'd100') +
  ggtitle('Poeciliidae Density Deviation From 2020 Quarter')

RvB_icta_d100 <- p_RR_RvB(
  my_data = RR_d100_fam %>% filter(family=='Ictaluridae'),
  my_variable = 'd100') +
  ggtitle('Ictaluridae Density Deviation From 2020 Quarter')

RvB_cent_d100 <- p_RR_RvB(
  my_data = RR_d100_fam %>% filter(family=='Centrarchidae'),
  my_variable = 'd100') +
  ggtitle('Centrarchidae Density Deviation From 2020 Quarter')

# -----------------------------------------------------------------------------
# Visualization: variable Vs Baseline (Species)
# -----------------------------------------------------------------------------
species_unlimited

RvB_Gaff_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="G. affinis"),
  my_variable = 'd100') +
  ggtitle('G. affinis Density Deviation From 2020 Quarter')

RvB_Lcya_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. cyanellus"),
  my_variable = 'd100') +
  ggtitle('L. cyanellus Density Deviation From 2020 Quarter')

RvB_Lmac_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. macrochirus"),
  my_variable = 'd100') +
  ggtitle('L. macrochirus Density Deviation From 2020 Quarter')

RvB_Lmeg_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. megalotis"),
  my_variable = 'd100') +
  ggtitle('L. megalotis Density Deviation From 2020 Quarter')

RvB_Laur_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. auritus"),
  my_variable = 'd100') +
  ggtitle('L. auritus Density Deviation From 2020 Quarter')

RvB_Lgul_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. gulosus"),
  my_variable = 'd100') +
  ggtitle('L. gulosus Density Deviation From 2020 Quarter')

RvB_Msal_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="M. salmoides"),
  my_variable = 'd100') +
  ggtitle('M. salmoides Density Deviation From 2020 Quarter')

RvB_Locu_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. oculatus"),
  my_variable = 'd100') +
  ggtitle('L. oculatus Density Deviation From 2020 Quarter')

RvB_Anat_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="A. natalis"),
  my_variable = 'd100') +
  ggtitle('A. natalis Density Deviation From 2020 Quarter')

RvB_Ngyr_d100 <- p_RR_RvB(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="N. gyrinus"),
  my_variable = 'd100') +
  ggtitle('N. gyrinus Density Deviation From 2020 Quarter')

# -----------------------------------------------------------------------------
# Visualization: Response Ratio
# -----------------------------------------------------------------------------

p_RR <- function(my_data, my_variable) {
  
  colnames(my_data) <- str_replace_all(
    colnames(my_data), my_variable, 'X')
  
  categorical_colors <- c( 'cyan', 'red', 'grey90')
  
  my_data %>%
    fix_site_order() %>%
    mutate(upper_lim = X_bas20_fill + X_sd20,
           lower_lim = X_bas20_fill - X_sd20,
           Comparison=case_when(
             X > upper_lim ~ 'Above',
             X < lower_lim ~ 'Below',
             X <= upper_lim & X >= lower_lim ~ 'Inside') ) %>%
    ggplot() +
    facet_wrap(~site_code, scale='fixed') +
    geom_smooth(aes(x=collection_date, y = X_RR_qt20),
                method = 'loess', span=.8,
                se=F, linetype=2, color='grey45', linewidth=.5) +
    geom_point(aes(x=collection_date, y = X_RR_qt20, fill= Comparison),
               shape=21, size = 3, color = 'black', alpha=.3) +
    geom_point(aes(x=collection_date, y = X_RR_qt20, color=Comparison),
               shape=21, size = 3) +
    dark_theme_gray(base_size=12) +
    scale_fill_manual(values=categorical_colors) +
    scale_color_manual(values=categorical_colors) +
    scale_x_date(date_labels = "%m") +
    ylab(my_variable) +
    xlab('Time') +
    theme(legend.position = 'none')
}

p_RR_com_d100 <- p_RR(my_data = RR_d100_com,
     my_variable = 'd100') +
  ggtitle('Community Abundance Response Ratio')

p_RR_com_rich <- p_RR(my_data = RR_rich_com,
     my_variable = 'rich') +
  ggtitle('Species Richness Response Ratio')

p_RR_com_cdist <- p_RR(my_data = RR_cdis_com,
     my_variable = 'cdis') +
  ggtitle('RDA Centroid Distance Response Ratio')

p_RR_leuc <- p_RR(
  my_data = RR_d100_fam %>% filter(family=='Leuciscidae'),
  my_variable = 'd100') +
  ggtitle('Leuciscidae Density Response Ratio') +
  facet_wrap(~site_code, scale='free_y')

p_RR_poec <- p_RR(
  my_data = RR_d100_fam %>% filter(family=='Poeciliidae'),
  my_variable = 'd100') +
  ggtitle('Poeciliidae Density Response Ratio') +
  facet_wrap(~site_code, scale='free_y')

p_RR_icta <- p_RR(
  my_data = RR_d100_fam %>% filter(family=='Ictaluridae'),
  my_variable = 'd100') +
  ggtitle('Ictaluridae Density Response Ratio') +
  facet_wrap(~site_code, scale='free_y')

p_RR_cent <- p_RR(
  my_data = RR_d100_fam %>% filter(family=='Centrarchidae'),
  my_variable = 'd100') +
  ggtitle('Centrarchidae Density Response Ratio') +
  facet_wrap(~site_code, scale='free_y')

p_RR_Gaff <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="G. affinis"),
  my_variable = 'd100') +
  ggtitle('G. affinis Density Response Ratio')

p_RR_Lcya <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. cyanellus"),
  my_variable = 'd100') +
  ggtitle('L. cyanellus Density Response Ratio')

p_RR_Lmac <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. macrochirus"),
  my_variable = 'd100') +
  ggtitle('L. macrochirus Density Response Ratio')

p_RR_Lmeg <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. megalotis"),
  my_variable = 'd100') +
  ggtitle('L. megalotis Density Response Ratio') +
facet_wrap(~site_code, scale='free_y')

p_RR_Laur <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. auritus"),
  my_variable = 'd100') +
  ggtitle('L. auritus Density Response Ratio')+
  facet_wrap(~site_code, scale='free_y')

p_RR_Lgul <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. gulosus"),
  my_variable = 'd100') +
  ggtitle('L. gulosus Density Response Ratio')

p_RR_Msal <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="M. salmoides"),
  my_variable = 'd100') +
  ggtitle('M. salmoides Density Response Ratio')

p_RR_Locu <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="L. oculatus"),
  my_variable = 'd100') +
  ggtitle('L. oculatus Density Response Ratio')

p_RR_Anat <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="A. natalis"),
  my_variable = 'd100') +
  ggtitle('A. natalis Density Response Ratio')

p_RR_Ngyr <- p_RR(
  my_data = RR_d100_spe %>% filter(lowest_taxon=="N. gyrinus"),
  my_variable = 'd100') +
  ggtitle('N. gyrinus Density Response Ratio')
