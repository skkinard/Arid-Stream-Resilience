# site_PCA
# site_table
# Sean Kinard
# 1-19-2023

# use PCA identify patterns in environmental variables among sample sites

#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------

# load data
ltevars <- read_csv("Data/site_longterm_nogap.csv")
devars <- read_csv("Data/site_daily_nogap.csv")

#------------------------------------------------------------------------------
# Long Term environmental variables: PCA
#------------------------------------------------------------------------------

# Data prep
ltevars <- ltevars %>%
  dplyr::rename(Rainfall = annualrain,
                'Basin' = basinsize,
                'L.Dev.' = developedland,
                'L.Other' = otherland,
                'L.Forest' = forestland,
                'L.Crops' = cropland,
                'F.Flash' = flsh,
                'F.Base' = bfi,
                'F.Var' = vardf,
                'F.Med' = medianflow,
                'F.LFF' = frlfsp,
                'F.HPP3' = hfpp7,
                'F.HPP7' = hfpp3,
                'F.LPP' = lfpp,
                'F.Avg' = meanflow,
                'Season' = seasonality  )

# PCA and plots
PCA_out <- prcomp(ltevars[,-1], scale = TRUE)
summary(PCA_out)
PCA_lte_importance <- as.data.frame(summary(PCA_out)$importance)
PCA_lte_importance$metric <- rownames(PCA_lte_importance)
(PCA_lte_importance <- PCA_lte_importance[,c(8,1:7)])

# diagnostic scree plot
plot(PCA_out, type = "l")

# PCA plot
PCA_lte <- ggbiplot(PCA_out, obs.scale=1, var.scale=1,
         labels = ltevars$site_code,
         labels.size = 4,
         varname.size = 4,
         varname.adjust = 1.2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

caption_plot_PCA_lte <- "Principal Component Analysis of long-term environmental variables for sampling locations. Variables were scaled to have unit variance before the analysis. landuse categorizations are denoted by 'L.' and were estimated using satellite imagery in 2006. Other variables represent 20 year averages. Average annual rainfall is labeled 'Rainfall'. Flow metrics (denoted with 'F.')  include annual average, annual median, base flow index (avg of lowest 7 day stretch / average daily flow), low pulse percentage (% days below 25 percentile of discharge), high flow pulse percentages (% above 3x and 7x average flow), low flow freqency (# times where daily discharge drops below the 5th percentile) and flashiness (the cumulative changes in day to day discharge divided by cumulative annual discharge)."

# generating tables
m_env_lte <- cbind(ltevars[,-1], PCA_out$x[,1:2])

# correlations with PCA1
table_pca_lte <- tibble(cor(m_env_lte[,1:16],m_env_lte[,17:18]) %>% 
         as.data.frame %>%
         rownames_to_column('variable')) %>%
  arrange(desc(abs(PC1)))

caption_PCA_table_cor_de <- "Table containing correlations between the first two axes from principal component analysis and long-term environmental variables."

#------------------------------------------------------------------------------
# Daily environmental variables: PCA
#------------------------------------------------------------------------------

# data prep
site_order <- ltevars %>% arrange(Rainfall) %>% pull(site_code)

devars <- devars %>%
  mutate(site_code = fct_relevel(site_code, site_order)) %>%
  arrange(site_code, collection_date) %>%
  select( ! c(q, q_2w_mn)) %>%
  dplyr::rename(D.max = q_2w_max,
                D.min = q_2w_min,
                A.Dia = diatoms,
                A.FGr = green_algae,
                Depth = depth_mx,
                Width = width,
                Gravel = gravel,
                Silt = silt,
                Canopy = canopy_density_mid,
                Cond. = conductivity,
                Oxy. = do_mg_l,
                Amm. = nh4_n,
                Nitr. = no3n,
                Phos.=ortho_p,
                A.BGC =bluegreen_cyano)

devars20 <- devars %>%
  filter(collection_date > ymd('2020-01-01')) 

devars18<- devars %>%
  filter(collection_date > ymd('2018-01-01') &
           collection_date < ymd('2019-01-01')) 

devars17 <- devars %>% # outlier flow data at WMC 2017-09-30
  filter(collection_date < ymd('2018-01-01') &
           collection_date > ymd('2017-01-01')) 

# PCA and plots

my_pca <- function(my_data) {
  
  # PCA and plots
  PCA_de <- prcomp(my_data[,-c(1:3)], scale = TRUE)
  summary(PCA_de)
  PCA_de_importance <- as.data.frame(summary(PCA_de)$importance)
  PCA_de_importance$metric <- rownames(PCA_de_importance)
  (PCA_de_importance <- PCA_de_importance[,c(8,1:7)])
  
  # diagnostic scree plot
  plot(PCA_de, type = "l")
  
  # PCA plot
  plot_PCA_de <- ggbiplot(PCA_de, obs.scale=1, var.scale=1,
                          labels = my_data$site_code,
                          labels.size = 4,
                          groups = my_data$site_code,
                          linewidth=2,
                          ellipse=T,
                          varname.size = 4,
                          varname.adjust = 3) +
    scale_color_manual(values = my_colors) +
    theme_dark(base_size = 12) +
    theme(legend.position = "none")
  
  # generating tables
  m_env_de <- cbind(my_data[,-c(1:3)], PCA_de$x[,1:2])
  
  # correlations with PCA1
  PCA_table_cor_de <- tibble(cor(m_env_de[,1:15],m_env_de[,16:17]) %>% 
                               as.data.frame %>%
                               rownames_to_column('variable'))
  
  my_output <- list(plot_PCA_de,PCA_table_cor_de)
  
  return(my_output) }

# plots

PCA_ste_2017 <- my_pca(devars17)[[1]]
PCA_ste_2018 <- my_pca(devars18)[[1]]
PCA_ste_2020 <- my_pca(devars20)[[1]]

PCA_ste_all <- my_pca(devars)[[1]] + ggtitle('2017-2020') 

caption_plot_PCA_de <- 'Principal Component Analysis of environmental variables measured during sample events. Variables were scaled to have unit variance before the analysis. Flow is average annual discharge, HFPP is the proportion of the annual discharge that is 3x higher than the Flow, LFPP is the proportion of discharge below the 25th percentile, Flashiness is the cumulative changes in day to day discharge divided by cumulative annual discharge.  Horizontally, sites roughly order according to precipitation regime. Sites with drier climate like Tranquitas and Perdido have higher silt, conductivity, and green algae. Wetter sites like West Mustang, Garcitas, and Mission River have elevated channel depth and channel widths as well as greater discharge and ammonia conncentrations. San Fernando and Aransas are distinguished along the vertical axis with elevated nitrates, phophates, canopy density, gravel, and blue-green cyano bacteria. East Mustang and Placedo appear to share characteristics with many or the other streams.'

# tables
table_pca_ste_17 <- my_pca(devars17)[[2]]
table_pca_ste_18 <- my_pca(devars18)[[2]]
table_pca_ste_20 <- my_pca(devars20)[[2]]
table_pca_ste_all <- my_pca(devars)[[2]]

caption_PCA_table_cor_de <- "Table containing correlations between the first two axes from principal component analysis and daily environmental variables."

#------------------------------------------------------------------------------
# Export Figures
#------------------------------------------------------------------------------
ggsave('Figures/PCA_lte.pdf',
       plot = PCA_lte,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/pca_ste_all.pdf',
       plot = PCA_ste_all,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('Figures/pca_ste_2017.pdf',
       plot = PCA_ste_2017,
       width = 9,
       height = 18,
       units = c("in"))

ggsave('Figures/pca_ste_2018.pdf',
       plot = PCA_ste_2017,
       width = 9,
       height = 18,
       units = c("in"))

ggsave('Figures/pca_ste_2020.pdf',
       plot = PCA_ste_2017,
       width = 9,
       height = 18,
       units = c("in"))

#------------------------------------------------------------------------------
# Clean-up
#------------------------------------------------------------------------------
rm(devars, ltevars, m_env_lte, PCA_out,
   PCA_lte_importance, devars17, devars18, devars20)

