# ================== Future land-use covariate harmonisation ==================

setwd(".../VectorNet_ 2_original")     ## change the directory as needed

source("functions_general.R")
source("functions_plotting.R")

# Load the reference map and supporting functions.

out <- attach("Data/data_INLA.RData")
map_inla <- map_inla                           ## attaching the map
detach("file:Data/data_INLA.RData")

# Initialise storage for yearly land-use tables.
tb_list <- NULL
count <- 0

# Read annual land-use data and assemble a yearly covariate table.
for (year in 2021:2100) {                   
  
  print(sprintf("Processing year: %d", year))      ## to know which year is being processed
  
  
  land_use <- 
    read.csv(paste0(".../Data/climate/Final_run_data/SSPs&Pi_Future_climate_data/Socioeconomic_data/Land_use/SSP1/land_use/land_use_", year, ".csv"))%>% 
    as_tibble()
  
  
  # Remove an index column if it is present in the exported CSV.
  if ("X" %in% colnames(land_use)) {
    land_use <- land_use %>% select(-X)
  }
  
  print(land_use)

  # Harmonise location names with the model map.
  vals <- tibble(Location=land_use$Location, 
                 LocationNa=map_inla$LocationNa) %>% 
    mutate(id=1:n())                                    
  
  vals$is_equal <- NA           ## add a NA column to the above column
  
  for (j in 1:nrow(vals)) {
    vals$is_equal[j] <- identical(vals$Location[j],vals$LocationNa[j])
  }
  
  vals %>% filter(!is_equal)      
  land_use <- land_use %>% mutate(
    Location = case_when(
      Location == "Cote-d?Or" ~ "Cote-d’Or",
      Location == "Cotes-d?Armor" ~ "Cotes-d’Armor",
      TRUE ~ Location
    )
  )
  
  
  
  
  # Store the annual land-use data 
  tb_land_use <- tibble(
    geo=map_inla$geo,                        ## location code
    LocationNa=land_use$Location,       ## location name column
    primary_forested_area = land_use$primary.foreseted.area,
    primary_non_forested_area = land_use$primary.non.foreseted.area,      ## all covariates
    secondary_forested_area = land_use$secondary.foreseted.area,
    croplands = land_use$croplands,
    pasture_land = land_use$pasture.land,
    rangelands = land_use$rangelands,
    urban_land = land_use$urban.land
  )
  
  count <- count + 1
  tb_list[[count]] <- tb_land_use    ## each list contain each year total land_use data
  
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


## to check
combined_df[combined_df$year_idx == 18,]
tb_list[[18]]

combined_df <- combined_df %>%
  mutate(year = 2021 + year_idx - 1) %>%  
  relocate(year, .after = year_idx)     

unique(combined_df$year)     ## check the years and 

# Save the final yearly land-use covariate table.
land_use <- combined_df

save("land_use",
     file=paste0(".../Data/climate/Final_covariates_data/Future_data/Final_future_ssp126_land_use_data_yearly.RData"))
