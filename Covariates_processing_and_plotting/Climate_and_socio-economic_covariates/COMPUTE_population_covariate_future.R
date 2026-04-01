# ================== Future population covariate harmonisation ==================

setwd(".../VectorNet_ 2_original")       ## change the directory as needed

source("functions_general.R")
source("functions_plotting.R")
source("functions_plot_article.R")

# Load the reference map and supporting functions.

out <- attach("Data/data_INLA.RData")
map_inla <- map_inla                           ## attaching the map
detach("file:Data/data_INLA.RData")

# Initialize storage for yearly population tables.
tb_list <- NULL
count <- 0

# Read annual population data and assemble a yearly covariate table.
for (year in 2021:2100) {                   
  
  
  print(sprintf("Processing year: %d", year))     
  
  
  population <- 
    read.csv(paste0(".../Data/climate/Final_run_data/SSPs&Pi_Future_climate_data/Socioeconomic_data/population/SSP5/population/population_", year, ".csv"))%>% 
    as_tibble()
  
  # Remove an index column if it is present in the exported CSV.
if ("X" %in% colnames(population)) {
    population <- population %>% select(-X)
  }
  
  print(population)
  
  
  population <- population %>%
    select(Location, `total.population`)      ## names are with a . in between (like total.population, rural.population etc.)
  
  
  
  #Harmonise location names with the model map.

  
  vals <- tibble(Location=population$Location, 
                 LocationNa=map_inla$LocationNa) %>% 
    mutate(id=1:n())                                   
  
  vals$is_equal <- NA           ## add a NA column to the above column
  
  for (j in 1:nrow(vals)) {
    vals$is_equal[j] <- identical(vals$Location[j],vals$LocationNa[j])
  }
  
  vals %>% filter(!is_equal)       
  
  population <- population %>% mutate(
    Location = case_when(
      Location == "Cote-d?Or" ~ "Cote-d’Or",
      Location == "Cotes-d?Armor" ~ "Cotes-d’Armor",
      TRUE ~ Location
    )
  )
  
  
  
  
  # Store the annual pop
  tb_population <- tibble(
    geo=map_inla$geo,                        ## location code
    LocationNa=population$Location,       ## location name column
    population = population$total.population
  )
  
  count <- count + 1
  tb_list[[count]] <- tb_population    ## each list contain each year total population data
  
}
  

# Combine yearly tables and assign a sequential year index.
tb_list <- lapply(seq_along(tb_list), function(i) {
  tb_list[[i]] %>% mutate(year_idx = i)  # Add year index
})

combined_df <- bind_rows(tb_list)

combined_df <- combined_df %>%
  relocate(year_idx, .after = geo)

print(combined_df)

unique(combined_df$year_idx)

combined_df[combined_df$year_idx == 18,]
tb_list[[18]]

combined_df <- combined_df %>%
  mutate(year = 2021 + year_idx - 1) %>%  # Calculate year based on year_idx (starting from 2021)
  relocate(year, .after = year_idx)      

unique(combined_df$year)     ## check the years 

# Compute population density using NUTS3 area.
population <- combined_df

map_inla <- map_inla %>%
  mutate(area_km_2 = area / 1e6)       ## area to km^2 from m^2

nuts_area <- map_inla %>%
  st_drop_geometry() %>%
  select(geo, area_km_2)

population <- population %>%
  left_join(nuts_area, by = "geo") %>%
  mutate(pop_density = population / area_km_2)

population

# Save the final yearly population covariate table.
save("population",
     file=paste0(".../Data/climate/Final_covariates_data/Future_data/Final_future_ssp585_population_data_yearly.RData"))
