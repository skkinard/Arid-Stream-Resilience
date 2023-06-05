# Response v2
# Sean Kinard
# 03-01-23

# -----------------------------------------------------------------------------
# Count and Length vs Time by groups
# -----------------------------------------------------------------------------
# Attempt multiple scaled metrics in the same plot
group_scale <- function(df, grp, grp_x, my_var)
{ # filters by a group and group_x and then scales the aggregated value by site, including instance of zero catch
  df <- df %>% 
    dplyr::rename(g = {{grp}})
  df <- df %>%
    dplyr::rename(v = {{my_var}})
  
  my_output <- df %>%
    filter(g == grp_x) %>%
    group_by(site_code, collection_date) %>%
    dplyr::summarize(value = sum(v)) %>%
    ungroup() %>%
    right_join(d_density%>%select(site_code, collection_date)%>% unique()) %>% 
    mutate(value = ifelse(is.na(value), 0 ,value)) %>%
    arrange(site_code, collection_date) %>%
    pivot_wider(values_from=value, names_from = site_code) %>%
    mutate(AR = scale(AR),
           EM = scale(EM),
           GC = scale(GC),
           MR = scale(MR),
           PD = scale(PD),
           PL = scale(PL),
           SF = scale(SF),
           TR = scale(TR),
           WM = scale(WM)) %>%
    pivot_longer(cols=AR:WM, names_to='site_code', values_to='value') %>%
    na.omit() 
  
  return(my_output)
}

group_scale_plot <- function(grp, grp_x) {
  
d_count_s <- d_length %>%
  add_taxonomic() %>%
  group_by(site_code, collection_date, genus, species) %>%
  dplyr::summarize(cnt=length(TL_mu_spe)) %>%
  ungroup() %>%
  left_join(d_bio_categories) %>%
  group_scale(grp=grp, grp_x = grp_x, my_var='cnt') %>%
  dplyr::rename(count_s=value)

d_size_s <- d_length %>%
  add_taxonomic() %>%
  group_by(site_code, collection_date, genus, species) %>%
  dplyr::summarize(mu=mean(TL_mu_spe)) %>%
  ungroup() %>%
  left_join(d_bio_categories) %>%
  group_scale(grp=grp, grp_x = grp_x, my_var='mu') %>%
  dplyr::rename(size_s=value)
  
d_count_s %>%
  left_join(d_size_s) %>%
  add_time_vars() %>%
  mutate(site_code = fct_relevel(site_code, site_order)) %>%
  pivot_longer(cols=contains('_s'), 
               names_to='var', values_to='value') %>%
  mutate(var = str_replace_all(var, '_s','')) %>%
  ggplot(aes(nday, value)) +
  geom_smooth(aes(color=var),
              method='loess', se=F, linewidth=.25, linetype=2, 
              span=1, alpha=.2) +
  #geom_point(aes(fill=var), shape=21, size = 6, alpha = .1) +
  geom_point(aes(color=var), size = 1.2) +
  #geom_line(aes(color = var)) +
  dark_theme_grey(base_size=12) +
  xlab('Day of the Year') +
  ylab(element_blank()) +
  theme(legend.position='top') +
  scale_fill_manual(grp_x, values=c('darkgoldenrod1', 'skyblue')) +
  scale_color_manual(grp_x, values=c('darkgoldenrod1', 'skyblue')) +
  facet_grid(site_code~year) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

# Family
#plot_nlength_Poeciliidae <- group_scale_plot('family', 'Poeciliidae')
#plot_nlength_Centrarchidae <- group_scale_plot('family', 'Centrarchidae')
#plot_nlength_Ictaluridae <- group_scale_plot('family', 'Ictaluridae')
#plot_nlength_Lepisosteidae <- group_scale_plot('family', 'Lepisosteidae')
#plot_nlength_Leuciscidae <- group_scale_plot('family', 'Leuciscidae')
#plot_nlength_Cichlidae <- group_scale_plot('family', 'Cichlidae')
#plot_nlength_Cyprinodontidae <- group_scale_plot('family', 'Cyprinodontidae')

# Species
#plot_nlength_macrochirus <- group_scale_plot('lowest_taxon', 'L. macrochirus')
#plot_nlength_megalotis <- group_scale_plot('lowest_taxon', 'L. megalotis')
#plot_nlength_gulosus <- group_scale_plot('lowest_taxon', 'L. gulosus')
#plot_nlength_cyanellus <- group_scale_plot('lowest_taxon', 'L. cyanellus')
#plot_nlength_cyanoguttatum <- 
#  group_scale_plot('lowest_taxon', 'H. cyanoguttatum')
#plot_nlength_latipinna <- group_scale_plot('lowest_taxon', 'P. latipinna')
#plot_nlength_variegatus <- group_scale_plot('lowest_taxon', 'C. variegatus')
#plot_nlength_vigilax <- group_scale_plot('lowest_taxon', 'P. vigilax')

# Trophic Category
#plot_nlength_Herbivore <- group_scale_plot('trophic_category', 'Herbivore')
#plot_nlength_Omnivore <- group_scale_plot('trophic_category', 'Omnivore')
#plot_nlength_Piscivore <- group_scale_plot('trophic_category', 'Piscivore')

# Oxygen Sensitivity
# group_scale_plot('oxy_sens', 'sensitive')
# group_scale_plot('oxy_sens', 'mesotolerant')
# group_scale_plot('oxy_sens', 'tolerant')

# Reproductive Strategy
# group_scale_plot('rep_strat', 'simple_nest')
# group_scale_plot('rep_strat', 'complext_nest')
# group_scale_plot('rep_strat', 'livebearer')
# group_scale_plot('rep_strat', 'broadcast')
# group_scale_plot('rep_strat', 'other')

# Endemism
# group_scale_plot('geographic_range', 'Local')
# group_scale_plot('geographic_range', 'Restricted')
# group_scale_plot('geographic_range', 'Unlimited')

# -----------------------------------------------------------------------------