# prep_TERRG_fish_production
# Sean Kinard
# last update: 2023-06-02

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
setwd("/home/kinard/Documents/Research/Dissertation/02_Resilience/data")

source('preparation/merge_toolkit.R') # load packages and helper-functions

fish_terrg <- read_csv('clean_data/terrg_fish.csv')

#------------------------------------------------------------------------------
# TERRG Production
#------------------------------------------------------------------------------
production_meta <- read_csv("source_data/Strickland-data/MetadataSecProdFinalEstimates.csv")

# Create key for production file
key <- read_csv("source_data/Strickland-data/uidnamesNEW.csv") %>%
  mutate(Cycle = as.character(Cycle)) %>%
  separate(UID, into=c('site2', 'year', 'month', 'day'), sep="\\.") %>%
  mutate(Site_Code = substr(Site, 1,2),
         project = ifelse(year > 2018, 'TERRG', 'RAPID')) %>%
  dplyr::select(!c(site2, Site)) %>%
  unite("Collection_Date", c("year", "month", "day"), sep='-') %>%
  mutate(Collection_Date = ymd(Collection_Date)) %>%
  create_site_period() %>%
  full_join(fish_terrg %>% 
              dplyr::select(Site_Code, Collection_Date, site_period) %>%
              distinct() %>%
              mutate(cdate = Collection_Date) %>%
              separate(cdate, into = c('year','Cycle', 'day'), sep = '-') %>%
              mutate(project = 'TERRG',
                     Cycle = as.numeric(Cycle)) %>%
              mutate(Cycle = as.character(Cycle)) %>%
              create_terrg_period() %>%
              dplyr::select(! c(year, day)) )

# pull fish
commonfish <- pull(fish_terrg, CommonName) %>% unique()

production <- read_csv("source_data/Strickland-data/SecProdFinalEstimates.csv") %>%
  filter(Taxon %in% commonfish) %>%
  mutate(Site_Code = substr(Site, 1,2)) %>%
  dplyr::select(Site_Code, Taxon, X3B:NovN) %>%
  pivot_longer(cols = X3B:NovN, names_to="cycle_code", values_drop_na=T) %>%
  mutate(cycle_code = ifelse(substr(cycle_code,1,1) == 'X', 
                             str_replace_all(cycle_code, "X", "RAPID_"),
                             paste("TERRG", cycle_code, sep='_')) ) %>%
  mutate(cycle_code = str_replace(cycle_code, "Jan", "1")) %>%
  mutate(cycle_code = str_replace(cycle_code, "Mar", "3")) %>%
  mutate(cycle_code = str_replace(cycle_code, "May", "5")) %>%
  mutate(cycle_code = str_replace(cycle_code, "Jul", "7")) %>%
  mutate(cycle_code = str_replace(cycle_code, "Sep", "9")) %>%
  mutate(cycle_code = str_replace(cycle_code, "Nov", "11")) %>%
  mutate(cycle_code = str_replace(cycle_code, "B", "_mass")) %>%
  mutate(cycle_code = str_replace(cycle_code, "N", "_count")) %>%
  separate(cycle_code, into=c("project", "Cycle", "estimate_type"), sep="_") %>%
  mutate(id=seq(1:3894)) %>%
  pivot_wider(names_from=estimate_type, values_from=value) %>%
  dplyr::rename(CommonName = Taxon)

# merge rows
production <- full_join(filter(production, ! is.na(mass)) %>% 
                          dplyr::select(-c(id, count)),
                        filter(production, ! is.na(count)) %>% 
                          dplyr::select(-c(id, mass)) ) %>%
  dplyr::select(project, Site_Code, Cycle, everything()) %>%
  left_join(key) %>%
  dplyr::rename(biomass=mass, density = count)

# formatting
production <- production %>% r_friendly_colnames()

#------------------------------------------------------------------------------
# Export Clean Data
#------------------------------------------------------------------------------
write_csv(production, 'clean_data/terrg_fish_production.csv')

#------------------------------------------------------------------------------
# End prep_TERRG_fish_production