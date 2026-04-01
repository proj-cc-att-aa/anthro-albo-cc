# ================== Future mobility covariate construction ==================

setwd(".../VectorNet_ 2_original")  

# Load supporting functions, future population data, and the reference map.
source("functions_general.R")
source("functions_plotting.R")
source("functions_risk_from_proximity.R")
source("functions_spatial_networks.R")
source("functions_matrix_operations.R")
source("functions_plot_article.R")



folder_data <- "Data"

out <- load(".../Data/climate/Final_covariates_data/Future_data/Final_future_ssp126_population_data_yearly.RData")    
df_population <- population

out <- attach("Data/data_INLA.RData")
map_inla <- map_inla
detach("file:Data/data_INLA.RData")

year_vec_pop <- sort(unique(df_population$year_idx))    # list of all the unique years in population data set

# Compute yearly mobility fluxes from population and derive region-level summaries.
df_mobility <- NULL             
for (yr in year_vec_pop) {
  
  df_pop_year <- df_population %>%
    filter(year_idx == yr) %>%
    select(geo, population)
  
  map_inla_year <- map_inla %>%
    mutate(year_idx = yr)
  
  map_joined <- map_inla_year %>%
    left_join(df_pop_year, by = "geo")
  
  
  
  mobilty_flux_yearly <- mobility_from_population(map_joined)  # 
  
  
  mobility_plot_yearly <- plot_mobility_graph(mobilty_flux_yearly, nodes = "net_flux_abs")     
  

  
  # Derive net outward flux for each region.
  mobilty_flux_yearly$tb_regions <- mobilty_flux_yearly$tb_regions %>%
    mutate(net_flux = flux_out - flux_in)    ## absolute value of net flux ( outgoing flux - incoming flux) as nodes as a column in mobility$tb_regions, the values of net_flux can be negative too
  
  
  df_mobility <-  df_mobility %>% 
    bind_rows(tibble(geo = mobilty_flux_yearly$tb_regions$geo,                      
                     mobility_idx = mobilty_flux_yearly$tb_regions$mobility_idx,    
                     year_idx = yr,                                                 
                     net_flux_out =  mobilty_flux_yearly$tb_regions$net_flux,       
                     total_flux_in = mobilty_flux_yearly$tb_regions$flux_in))      
  
}
  

# Remove inherited vector names from numeric columns.
df_mobility <- df_mobility |>
  mutate(total_flux_in = unname(total_flux_in))

df_mobility <- df_mobility |>
  mutate(net_flux_out = unname(net_flux_out))

# Save the final yearly mobility covariate table.
save(df_mobility, file=".../Data/climate/Final_covariates_data/Future_data/Final_future_ssp126_mobility_data_yearly.RData")
