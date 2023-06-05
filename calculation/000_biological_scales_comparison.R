# Biological Scales of Comparison
# Sean Kinard
# 2023-03-01

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
# load data
source('R_Scripts/000_d_dataframes.R')

# -----------------------------------------------------------------------------
# Taxonomic groups
# -----------------------------------------------------------------------------
my_taxonomic <- c('order', 'family', 'genus', 
                  'species', 'lowest_taxon')

# Taxanomic info
d_spp <- read_csv('Data/forklength_all_megaframe.csv') %>%
  select(any_of(my_taxonomic)) %>%
  unique()

# -----------------------------------------------------------------------------
# Trophic Groups
# -----------------------------------------------------------------------------
d_trophic <- read_csv('Data/fbase_trophic_level.csv') %>%
  dplyr::rename(trophic=TL_mu_fb) %>%
  select(genus, species, trophic) %>%
  mutate(trophic_category = case_when(
    trophic < 2.5 ~ 'Herbivore',
    trophic  >= 2.5 & trophic < 3.5 ~ 'Omnivore',
    trophic >= 3.5 ~ 'Piscivore'))

# -----------------------------------------------------------------------------
# Functional Groups
# -----------------------------------------------------------------------------
d_function <- read_csv('Data/species_traits.csv') %>%
  separate(genus.species, into=c('genus', 'species'), sep='\\.')

colnames(d_function) <- str_replace_all(colnames(d_function), ' ', '_')

d_function <- d_function %>%
  mutate(rep_strat = case_when(
    Broadcaster > 0 ~ 'broadcast',
    Simple_Nest > 0 ~ 'simple_nest',
    Complex_nest > 0 ~ 'complext_nest',
    Bearer > 0 ~ 'livebearer' )) %>%
  mutate(rep_strat = ifelse(is.na(rep_strat), 'other', rep_strat)) %>%
  mutate(oxy_sens = case_when(
    Sensitive > 0 ~ 'sensitive',
    Mesotolerant > 0 ~ 'mesotolerant',
    Tolerant > 0 ~ 'tolerant')) %>%
  mutate(oxy_sens = ifelse(is.na(oxy_sens), 'other', oxy_sens)) %>%
  select(genus, species, rep_strat, oxy_sens) %>%
  mutate(genus = str_to_title(genus))

# -----------------------------------------------------------------------------
# Endemism (Geographic Distribution)
# -----------------------------------------------------------------------------
# presence or absence
d_present <- read_csv('Data/forklength_all_megaframe.csv') %>%
  select(lengthmm, any_of(my_taxonomic), any_of(my_index)) %>%
  group_by(site_code, collection_date, lowest_taxon) %>%
  dplyr::summarize(present = ifelse(sum(lengthmm)>0, 1, NA)) %>%
  pivot_wider(names_from=lowest_taxon, 
              values_from=present,
              values_fill=0) %>%
  pivot_longer(cols=contains(' '), names_to='lowest_taxon', values_to='presence') 

# Probability of Occurance 
prob_occur <- d_present %>%
  group_by(lowest_taxon, site_code) %>%
  dplyr::summarize(prob_present = sum(presence / length(presence)))

# Regional Species
d_endemic <- prob_occur %>%
  mutate(p_or_a = ifelse(prob_present>0, 1, 0)) %>%
  group_by(lowest_taxon) %>%
  dplyr::summarize(n_sites_occur = sum(p_or_a)) %>%
  arrange(desc(n_sites_occur)) %>%
  mutate(geographic_range = case_when(
    n_sites_occur > 6 ~ 'Unlimited',
    n_sites_occur %in% 4:6 ~ 'Restricted',
    n_sites_occur %in% 1:3 ~ 'Local' ))

# Rare species
species_rare <- d_present %>% 
  group_by(lowest_taxon) %>%
  dplyr::summarize(total_occurances = sum(presence)) %>%
  arrange(total_occurances) %>% 
  filter(total_occurances < 10) %>%
  pull(lowest_taxon)

species_unlimited <- d_endemic %>%
  filter(geographic_range == 'Unlimited') %>%
  pull(lowest_taxon)

species_restricted <- d_endemic %>%
  filter(geographic_range == 'Restricted') %>%
  pull(lowest_taxon)

species_Local <- d_endemic %>%
  filter(geographic_range == 'Local') %>%
  pull(lowest_taxon)

# -----------------------------------------------------------------------------
# Combine Biological Categories
# -----------------------------------------------------------------------------

d_bio_categories <- d_spp %>%
  left_join(d_trophic) %>%
  left_join(d_function) %>%
  left_join(d_endemic)
# -----------------------------------------------------------------------------
# End Biological Scales of Comparison
