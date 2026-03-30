# ================== Counterfactual relative humidity extraction from NetCDF to CSV ==================

library(ncdf4)
library(ggplot2)
library(terra)
library(leaflet)
library(lubridate)
library(terra)
library(raster)
library(terra)
library(reshape2)
library(dplyr)
library(sf)
library(exactextractr)

setwd(".../20CRV-ERA5")

for(year in 2001:2021){
  # Read the daily NetCDF file for each year
  nc_file_path <- paste0("Netcdf_file/20crv3-era5_counterclim_hurs_global_daily_", year, ".nc")

  nc_file <- nc_open(nc_file_path)
  time_var <- ncvar_get(nc_file, "time")
  time_units <- ncatt_get(nc_file, "time", "units")$value
  origin_date <- gsub("days since ", "", time_units)
  origin_date <- strsplit(origin_date, " ")[[1]][1]
  time_values = as_datetime(c(time_var*60*60*24),origin= origin_date)

  # Load the raster and study-area shapefile
  rast_file <- rast(nc_file_path)
  print(rast_file)

  shapf <- load(".../map_selected.Rdata")
  shapf
  geometry = map_selected$geometry
  sf_map <- st_sf(map_selected, crs = 4326)
  shapefile <- st_transform(sf_map, crs = "+proj=longlat +datum=WGS84")

  crs(rast_file) <- "EPSG:4326"
  rast_EU <- crop(rast_file, extent(shapefile))
  crs(rast_EU) <- crs(shapefile)

  # Extract area-weighted means for each region
  cell_area <- cellSize(rast_EU, unit = "km")
  output_data <- exact_extract(rast_EU, shapefile,
                               function(values, coverage_frac, weights) {
                                 weighted.mean(values, coverage_frac * weights, na.rm = TRUE)
                               },
                               weights = cell_area,
                               stack_apply = TRUE)

  extracted_df <- as.data.frame(output_data)
  extracted_df$region <- shapefile$LocationNa

  # Store regions in rows and dates in columns
  final_df <- data.frame(Location = shapefile$LocationNa)
  for (i in seq_along(time_values)) {
    final_df[, as.character(time_values[i])] <- output_data[, i]
  }

  output_file <- paste0("relative_humidity/relative_humidity_mean_", year, ".csv")
  write.csv(final_df, output_file)

  cat("year",year, "done")
}
