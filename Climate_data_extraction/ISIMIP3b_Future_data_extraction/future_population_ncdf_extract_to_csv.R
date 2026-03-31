# ================== Future population NetCDF extraction to regional CSV files ==================

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

setwd(".../Final_run_data")

for(year in 2015:2100){             ### take care of time values for 2010, origin is 2001-01-01, for 2011 to 2020 (origin is 2011-01-01) and from 2021 to 2024 (origin is 2021-01-01)

  nc_file_path <- paste0("Netcdf_file/Global_ssp585_population_", year, ".nc")

  nc_file <- nc_open(nc_file_path)

  rast_file <- rast(nc_file_path)
  print(rast_file)
  plot(rast_file)

  ## loading the shapefile
  shapf <- load(".../map_selected.Rdata")
  shapf
  geometry = map_selected$geometry
  sf_map <- st_sf(map_selected, crs = 4326)
  shapefile <- st_transform(sf_map, crs = "+proj=longlat +datum=WGS84")

  crs(rast_file) <- "EPSG:4326"

  #Cropping the data so that it matches with area covered by shapefile
  rast_EU <- crop(rast_file, extent(shapefile))
  plot(rast_EU)

  crs(rast_EU) <- crs(shapefile)

  # Aggregate by summing the overlapped share of each cell’s count.
  output_data = exact_extract(
    rast_EU,        # people per cell
    shapefile,                   # polygons
    fun = function(values, coverage_fraction) {
      sum(values * coverage_fraction, na.rm = TRUE) # people, ## we will keep the function as sum as we want the total population not the aggregated population
    },
    stack_apply = TRUE
  )

  extracted_df <- as.data.frame(output_data)

  # Add region names or IDs from the shapefile to the dataframe
  extracted_df$region <- shapefile$LocationNa         # Adjust for the column name in shapefile

  final_df <- data.frame(Location = shapefile$LocationNa, extracted_df)


  # Rename columns (to remove the weighted_mean. term)
  final_df <- final_df %>%
    dplyr::select(
      Location = Location,
      `rural population` = fun.rural.population,
      `total population` = fun.total.population,
      `urban population` = fun.urban.population
    )

  final_df

  output_file <- paste0("population/population_", year, ".csv")

  write.csv(final_df, output_file)

  cat("year",year, "done")

}

