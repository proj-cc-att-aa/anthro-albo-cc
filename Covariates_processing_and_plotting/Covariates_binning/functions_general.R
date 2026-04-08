# ================== General spatial data utilities ==================

packages <- c("tidyverse", "stringr", "eurostat", "sf", "sp", "spdep")

lapply(packages, library, character.only = TRUE)


# Load map NUTS.
load_map_NUTS <- function(restrict_to_Europe=TRUE) {
  
  map_NUTS <- get_eurostat_geospatial(nuts_level="3")
  
  map_NUTS_2010 <- get_eurostat_geospatial(nuts_level="3", year = "2010")
  
  ind_common <- which(map_NUTS_2010$geo %in% map_NUTS$geo)
  
  map_NUTS <- rbind(map_NUTS[,names(map_NUTS_2010)],
                    map_NUTS_2010[-ind_common,])
  
  shapes <- as_Spatial(map_NUTS$geometry)
  coords <- coordinates(shapes)
  map_NUTS <- map_NUTS %>%
    add_column(longitude = coords[,1], latitude = coords[,2])
  
  if (restrict_to_Europe)
    map_NUTS <- restrict_map_to_Europe(map_NUTS)
  
  map_NUTS <- as_tibble(map_NUTS)
  
  return(map_NUTS)
}




# Restrict map to Europe.
restrict_map_to_Europe <- function(map) {
  cat("Restricting map to region around Europe (longitude: [-10, 31], latitude: [35 75])\n")
  coord_lims <- list(longitude = c(-10,31), latitude = c(35,75))
  coords_restricted_idx <- which(map$longitude >= coord_lims$longitude[1] &
                                   map$longitude <= coord_lims$longitude[2] &
                                   map$latitude >= coord_lims$latitude[1] &
                                   map$latitude <= coord_lims$latitude[2])
  
  map <- map[coords_restricted_idx,]
  
  return(map)
}



# Add missing regions from NUTS map.
add_missing_regions_from_NUTS_map <- function(map_ECDC, mosq_dist, restrict_to_Europe=TRUE) {
  
  geo_not_in_map_ECDC <- unique(setdiff(
    mosq_dist %>% filter(observation_category_str != "no data") %>% pull(geo),
    map_ECDC$geo
  ))
  
  map_NUTS <- load_map_NUTS(restrict_to_Europe=FALSE)
  
  idx_found <- which(geo_not_in_map_ECDC %in% map_NUTS$geo)
  geo_to_add <- geo_not_in_map_ECDC[idx_found]
  cat("Regions that are missing in the ECDC map but which are in the NUTS map:")
  print(geo_to_add)
  
  idx_not_found <- which(!geo_not_in_map_ECDC %in% map_NUTS$geo)
  cat("Regions that are not found in NUTS map (all outside Europe):")
  print(geo_not_in_map_ECDC[idx_not_found])
  
  map_NUTS_to_add <- map_NUTS %>% filter(geo %in% geo_to_add)
  
  map_NUTS_to_add <- map_NUTS_to_add %>%
    select(-id, -NUTS_ID, -LEVL_CODE, -NAME_LATN, -FID) %>%
    rename(CountryISO=CNTR_CODE,
           LocationNa=NUTS_NAME)
  
  map_full <- bind_rows(map_ECDC %>%
                          select(geo, CountryISO, CountryNam, LocationNa, geometry, longitude, latitude) %>%
                          add_column(map="ECDC"),
                        map_NUTS_to_add %>% add_column(map="NUTS"))
  
  if (restrict_to_Europe)
    map_full <- restrict_map_to_Europe(map_full)
  
  return(map_full)
}




# Compute intersections.
compute_intersections <- function(map) {
  
  intersections <- st_intersects(map$geometry)
  
  intersection_size <- intersections %>% lapply(length) %>% unlist()
  intersections <- intersections[which(intersection_size != 1)]
  
  return(intersections)
}




# Find connected components.
find_connected_components <- function(map, geo_neighs=NULL, input_regions=NULL) {
  
  if (is.null(input_regions)) {
    input_regions <- 1:nrow(map)
  }
  if (is.null(geo_neighs)) {
    geo_neighs <- geo_to_neighs(map, use_snap=FALSE)
  }
  
  neighs_tb <- neighbor_list_from_graph(geo_neighs$nb_map)
  
  connected_components <- list()
  for (region_idx in input_regions) {
    
    if (region_idx %in% unlist(connected_components))
      next
    
    idx_connected_component <-
      get_indices_of_connected_component(geo_neighs=geo_neighs, neighs_tb=neighs_tb,
                                         i_start=region_idx)
    
    connected_components <- c(connected_components, list(idx_connected_component))
  }
  
  size <- unlist(lapply(connected_components, length))
  size_order <- order(size, decreasing = TRUE)
  connected_components <- connected_components[size_order]
  size <- size[size_order]
  
  map <- add_column(map, conn_comp=0, size_comp=0)
  for (j in seq_along(connected_components)) {
    
    map$conn_comp[connected_components[[j]]] <- j
    map$size_comp[connected_components[[j]]] <- size[j]
  }
  
  return(map)
}





# Get indices of connected component.
get_indices_of_connected_component <- function(geo_neighs, neighs_tb=NULL, i_start = NULL) {
  
  if (is.null(neighs_tb)) {
    
    neighs_tb <- neighbor_list_from_graph(geo_neighs$nb_map)
  }
  
  neighs <- neighs_tb$neighs
  n_neighs <- neighs_tb$n_neighs
  
  if (is_null(i_start))
    i_start <- which(n_neighs == max(n_neighs, na.rm=TRUE))
  if (!(is_integer(i_start) || is_double(i_start)))
    stop("i_start error")
  if (length(i_start) > 1)
    i_start <- i_start[1]
  
  connected_component <- c(i_start)
  i_curr <- unlist(neighs[[i_start]])
  i_active <- setdiff(i_curr, connected_component)
  while (length(i_active) > 0) {
    
    i_curr <- vector()
    for (i in i_active) {
      i_curr <- union(i_curr, unlist(neighs[[i]]))
    }
    
    connected_component <- union(connected_component, i_active) %>% sort()
    
    i_active <- setdiff(i_curr, connected_component)
  }
  
  return(connected_component)
}




# Neighbor list from graph.
neighbor_list_from_graph <- function(nb_map) {
  
  neighs_tb <- tibble(has_neigh = rep(FALSE, length(nb_map)),
                      neighs = NA,
                      n_neighs = NA)
  
  for (j in seq_along(nb_map)) {
    
    ns <- nb_map[[j]]
    
    if (!all(ns == 0)) {
      
      neighs_tb$has_neigh[j] <- TRUE
      neighs_tb$neighs[[j]] <- list(ns)
      neighs_tb$n_neighs[j] <- length(ns)
    }
  }
  
  return(neighs_tb)
}


# Write to excel.
write_to_excel <- function(tb, file_excel, n_round=NULL, append=FALSE,
                           sheet_name="0") {
  if (!is.null(n_round)) {
    tb <- tb %>% mutate(across(where(is.numeric), ~ round(.x, n_round)))
  }
  write.xlsx(tb, file=file_excel, append=append, sheetName=sheet_name)
}
