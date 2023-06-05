# terrg_fish_mislabel_fix
# Sean Kinard
# 2023_01-16

# noticed some mislabeled Lepomis microlophus

# Clean Workspace
rm(list=ls())

# Load Packages
library(tidyverse)

# load data
fish_terrg <- read_csv("/home/kinard/Documents/Research/Dissertation/02_Resilience/Data/TERRG-RAPID_MegaFrame/Strickland-data/terrg fish.csv")

# TERRG Fish fix
error1 <- fish_terrg %>% filter(CommonName == "Redspotted Sunfish" & GenusSpecies == "Lepomis microlophus")
error2 <- fish_terrg %>% filter(Site == "GC" & GenusSpecies == "Lepomis microlophus")
error3 <- fish_terrg %>% filter(Site == "GC" & CommonName == "Redspotted Sunfish")
# CommonName Redspotted has two different GenusSpecies labels while Redear only has one

# Solution: Relabel Redspotted to Redear
fish_terrg <- fish_terrg %>%
  mutate(CommonName = ifelse(GenusSpecies == "Lepomis microlophus", "Redear Sunfish", CommonName))

# export
write_csv(fish_terrg, '/home/kinard/Documents/Research/Dissertation/02_Resilience/Data/TERRG-RAPID_MegaFrame/Strickland-data/terrg_fish_fixed.csv')

# end terrg_fish_mislabel_fix