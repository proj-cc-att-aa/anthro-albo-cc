# ================== Historical population extraction to CSV ===================

library(ncdf4)
library(ggplot2)
library(terra)
library(leaflet)
library(lubridate)
library(raster)
library(reshape2)
library(dplyr)
library(sf)
library(exactextractr)

setwd(".../hist_pop_data")

for(year in 2010:2022){

  nc_file_path <- paste0("Netcdf_file/pop_count_global_", year, "_30_arc_sec_or_approx_1_km.nc")

  nc_file <- nc_open(nc_file_path)

  rast_file <- rast(nc_file_path)
  print(rast_file)
  plot(rast_file)

  shapf <- load(".../map_selected.Rdata")
  shapf
  geometry = map_selected$geometry
  sf_map <- st_sf(map_selected, crs = 4326)
  shapefile <- st_transform(sf_map, crs = "+proj=longlat +datum=WGS84")

  crs(rast_file) <- "EPSG:4326"

  rast_EU <- crop(rast_file, extent(shapefile))
  plot(rast_EU)

  crs(rast_EU) <- crs(shapefile)

  # Population totals are obtained by summing the overlap-weighted cell counts.
  output_data = exact_extract(
    rast_EU,
    shapefile,
    fun = function(values, coverage_fraction) {
      sum(values * coverage_fraction, na.rm = TRUE)    ## we will keep the function as sum as we want the total population not the aggregated population
    },
    stack_apply = TRUE
  )

  extracted_df <- as.data.frame(output_data)

  extracted_df$region <- shapefile$LocationNa

  final_df <- data.frame(Location = shapefile$LocationNa, extracted_df)

  # The exported table retains total population at the regional level.
  final_df <- final_df %>%
    select(
      Location = Location,
      `total population` = output_data
    )

  final_df

  output_file <- paste0("population/population_", year, ".csv")

  write.csv(final_df, output_file)

  cat("year",year, "done")

}


