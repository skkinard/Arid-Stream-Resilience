# site_table
# Sean Kinard
# 1-19-2023

###############################################################################
# setup
###############################################################################

# load data
d <- read_csv('Data/fish_biomass_megaframe.csv')
d_meta <- read_csv('Data/fish_rapid_terrg_megaframe_metadata.csv')
ltevars <- read_csv("Data/site_longterm_nogap.csv")
devars <- read_csv("Data/site_daily_nogap.csv")

site_order <- d %>% select(site_code, annualrain) %>% arrange(annualrain) %>% (unique) %>% pull(site_code)

###############################################################################
# Table LTE
###############################################################################

# Add elevation from Fernando thesis:
site_code <- c('TR', 'SF', 'AR', 'MR', 'PD', 'PL', 'GC', 'WM', 'EM')
Elevation <- c(18,62,47,14,50,17,20,20,20)
ltevars <- tibble(site_code, Elevation) %>%
  right_join(ltevars)
Units = 'm'
Column = 'Elevation'
d_meta <- tibble(Units, Column) %>%
  full_join(d_meta)

# renaming
site_info <- d %>% 
  select(site_code, site_name, staid, lat, lon) %>%
  mutate(site_name = str_replace_all(site_name, ' Creek', '')) %>%
  mutate(site_name = str_replace_all(site_name, 'iver', '')) %>%
  mutate(site_name = str_replace_all(site_name, 'West ', 'W.')) %>%
  mutate(site_name = str_replace_all(site_name, 'East ', 'E.')) %>%
  mutate(site_name = str_replace_all(site_name, 'San ', 'S.')) %>%
  distinct()

# set variable order
v_order = c('Sta.ID', 'Latitude', 'Longitude', 'Elevation', 'Rainfall',
            'Air Temp.', 'Basin', 'Developed', 'Crops', 'Flow',
            'Flashiness', 'HFPP', 'LFPP', 'Season')

# make table renaming variables for publishing
table_LTE <- ltevars %>%
  left_join(site_info) %>%
  mutate(site_code = site_name) %>%
  select(! site_name) %>%
  pivot_longer(cols = Elevation:lon,
               names_to='Variable', values_to = 'y') %>%
  mutate(y = round(y, 2)) %>%
  left_join(d_meta %>% select( ! Description) %>%
              dplyr::rename('Variable' = Column)) %>%
  pivot_wider(names_from=site_code, values_from=y) %>%
  mutate(Variable = str_replace_all(Variable, 'annualrain', 'Rainfall')) %>%
  mutate(Variable = str_replace_all(Variable, 'annualtemp', 'Air Temp.')) %>%
  mutate(Variable = str_replace_all(Variable, 'basinsize', 'Basin')) %>%
  mutate(Variable = str_replace_all(Variable, 'developedland', 'Developed')) %>%
  mutate(Variable = str_replace_all(Variable, 'cropland', 'Crops')) %>%
  mutate(Variable = str_replace_all(Variable, 'flsh', 'Flashiness')) %>%
  mutate(Variable = str_replace_all(Variable, 'hfpp3', 'HFPP')) %>%
  mutate(Variable = str_replace_all(Variable, 'lfpp', 'LFPP')) %>%
  mutate(Variable = str_replace_all(Variable, 'meanflow', 'Flow')) %>%
  mutate(Variable = str_replace_all(Variable, 'seasonality', 'Season')) %>%
  mutate(Variable = str_replace_all(Variable, 'staid', 'Sta.ID')) %>%
  mutate(Variable = str_replace_all(Variable, 'lat', 'Latitude')) %>%
  mutate(Variable = str_replace_all(Variable, 'lon', 'Longitude')) %>%
  mutate(Units = str_replace_all(Units, 'degrees', 'o')) %>%
  mutate(Units = str_replace_all(Units, 'character', '-')) %>%
  mutate(Variable = fct_relevel(Variable, v_order)) %>%
  arrange(Variable)

caption_table_LTE <- 'Table of long-term environmental variables for sampling locations. Values represent 20 year averages. Flow is average annual discharge, HFPP is the proportion of the annual discharge that is 3x higher than the Flow, LFPP is the proportion of discharge below the 25th percentile, Flashiness is the cumulative changes in day to day discharge divided by cumulative annual discharge, Season approximates the degree to which the flow varies during the course of a single year.'

knitr::kable(table_LTE, align=c('l','c','c','c','c','c','c','c','c','c','c'))

###############################################################################
# Table Daily Variables quarterly
###############################################################################
colnames(devars)

dvar_x <- devars %>%
  group_by(site_code, qtr) %>%
  dplyr::summarize(q = mean(q),
                   q_2w_mn = mean(q_2w_mn),
                   q_2w_mn = mean(q_2w_mn),
                   q_2w_min = mean(q_2w_min),
                   diatoms = mean(diatoms),
                   green_algae = mean(green_algae),
                   depth_mx = mean(depth_mx),
                   width = mean(width),
                   silt = mean(silt),
                   canopy_density_mid = mean(canopy_density_mid),
                   conductivity = mean(conductivity),
                   do_mg_l = mean(do_mg_l),
                   nh4_n = mean(nh4_n),
                   no3n = mean(no3n),
                   ortho_p = mean(ortho_p),
                   bluegreen_cyano = mean(bluegreen_cyano),
  ) %>%
  pivot_longer(cols = q:bluegreen_cyano, names_to = 'x', values_to = 'mu')
  
dvar_sd <- devars %>%
  group_by(site_code, qtr) %>%
  dplyr::summarize(q = sd(q),
                   q_2w_mn = sd(q_2w_mn),
                   q_2w_mn = sd(q_2w_mn),
                   q_2w_min = sd(q_2w_min),
                   diatoms = sd(diatoms),
                   green_algae = sd(green_algae),
                   depth_mx = sd(depth_mx),
                   width = sd(width),
                   silt = sd(silt),
                   canopy_density_mid = sd(canopy_density_mid),
                   conductivity = sd(conductivity),
                   do_mg_l = sd(do_mg_l),
                   nh4_n = sd(nh4_n),
                   no3n = sd(no3n),
                   ortho_p = sd(ortho_p),
                   bluegreen_cyano = sd(bluegreen_cyano),
  ) %>%
  pivot_longer(cols = q:bluegreen_cyano, names_to = 'x', values_to = 'sd')

table_devar_qtrly <- left_join(dvar_x, dvar_sd) %>%
  mutate(z=paste(round(mu,1), "\u00b1", round(sd,1), sep=' ')) %>%
  left_join(d_meta %>% select( ! Description) %>%
              dplyr::mutate(Column = str_replace_all(Column, '_mean', '')) %>%
              dplyr::rename('x' = Column)) %>%
  dplyr::mutate(x = str_replace_all(x, 'q', 'D.x')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'q_2w_mn', 'D.mu')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'q_2w_max', 'D.max')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'q_2w_min', 'D.min')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'diatoms', 'Algae.D')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'green_algae', 'Algae.G')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'depth_mx', 'Depth')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'width', 'Width')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'gravel', 'Gravel')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'silt', 'Silt')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'canopy_density_mid', 'Canopy')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'conductivity', 'Cond.')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'do_mg_l', 'D.O')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'nh4_n', 'Amm.')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'no3n', 'Nitr.')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'ortho_p', 'Phos.')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'bluegreen_cyano', 'Algae_BGC')) %>%
  mutate(Variable = paste(x, ' (', Units, ')', sep='')) %>%
  select(site_code, qtr, Variable, z) %>%
  ungroup() %>%
  mutate(site_code = fct_relevel(site_code, site_order)) %>%
  arrange(site_code, Variable, qtr) %>%
  pivot_wider(names_from = site_code, values_from = z)
 

###############################################################################
# Table Daily Variables site
###############################################################################
dvarq_x <- devars %>%
  group_by(site_code) %>%
  dplyr::summarize(q = mean(q),
                   q_2w_mn = mean(q_2w_mn),
                   q_2w_mn = mean(q_2w_mn),
                   q_2w_min = mean(q_2w_min),
                   diatoms = mean(diatoms),
                   green_algae = mean(green_algae),
                   depth_mx = mean(depth_mx),
                   width = mean(width),
                   silt = mean(silt),
                   canopy_density_mid = mean(canopy_density_mid),
                   conductivity = mean(conductivity),
                   do_mg_l = mean(do_mg_l),
                   nh4_n = mean(nh4_n),
                   no3n = mean(no3n),
                   ortho_p = mean(ortho_p),
                   bluegreen_cyano = mean(bluegreen_cyano),
  ) %>%
  pivot_longer(cols = q:bluegreen_cyano, names_to = 'x', values_to = 'mu')

dvarq_sd <- devars %>%
  group_by(site_code) %>%
  dplyr::summarize(q = sd(q),
                   q_2w_mn = sd(q_2w_mn),
                   q_2w_mn = sd(q_2w_mn),
                   q_2w_min = sd(q_2w_min),
                   diatoms = sd(diatoms),
                   green_algae = sd(green_algae),
                   depth_mx = sd(depth_mx),
                   width = sd(width),
                   silt = sd(silt),
                   canopy_density_mid = sd(canopy_density_mid),
                   conductivity = sd(conductivity),
                   do_mg_l = sd(do_mg_l),
                   nh4_n = sd(nh4_n),
                   no3n = sd(no3n),
                   ortho_p = sd(ortho_p),
                   bluegreen_cyano = sd(bluegreen_cyano),
  ) %>%
  pivot_longer(cols = q:bluegreen_cyano, names_to = 'x', values_to = 'sd')

table_devar_site <- left_join(dvarq_x, dvarq_sd) %>%
  mutate(z=paste(round(mu,1), "\u00b1", round(sd,1), sep=' ')) %>%
  left_join(d_meta %>% select( ! Description) %>%
              dplyr::mutate(Column = str_replace_all(Column, '_mean', '')) %>%
              dplyr::rename('x' = Column)) %>%
  dplyr::mutate(x = str_replace_all(x, 'q', 'D.x')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'q_2w_mn', 'D.mu')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'q_2w_max', 'D.max')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'q_2w_min', 'D.min')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'diatoms', 'Algae.D')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'green_algae', 'Algae.G')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'depth_mx', 'Depth')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'width', 'Width')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'gravel', 'Gravel')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'silt', 'Silt')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'canopy_density_mid', 'Canopy')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'conductivity', 'Cond.')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'do_mg_l', 'D.O')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'nh4_n', 'Amm.')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'no3n', 'Nitr.')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'ortho_p', 'Phos.')) %>% 
  dplyr::mutate(x = str_replace_all(x, 'bluegreen_cyano', 'Algae_BGC')) %>%
  mutate(Variable = paste(x, ' (', Units, ')', sep='')) %>%
  select(site_code, Variable, z) %>%
  mutate(site_code = fct_relevel(site_code, site_order)) %>%
  arrange(site_code, Variable) %>%
  pivot_wider(names_from = site_code, values_from = z)

caption_table_devar_site <- "Principal Table of long-term environmental variables for sampling locations. Values represent 20 year averages. Flow is average annual discharge, HFPP is the proportion of the annual discharge that is 3x higher than the Flow, LFPP is the proportion of discharge below the 25th percentile, Flashiness is the cumulative changes in day to day discharge divided by cumulative annual discharge.  Horizontally, sites roughly order according to precipitation regime. Sites with drier climate like Tranquitas and Perdido have higher silt, conductivity, and green algae. Wetter sites like West Mustang, Garcitas, and Mission River have elevated channel depth and channel widths as well as greater discharge and ammonia conncentrations. San Fernando and Aransas are distinguished along the vertical axis with elevated nitrates, phophates, canopy density, gravel, and blue-green cyano bacteria. East Mustang and Placedo appear to share characteristics with many or the other streams."


###############################################################################
# Clean-up
###############################################################################
rm(d, d_meta, devars, dvar_sd, dvar_x, dvarq_sd, dvarq_x, ltevars, site_info,
    Column, Elevation, site_code, Units, v_order)
