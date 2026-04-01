# ================== Historical and scenario-based land-use extraction to CSV ===================

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

setwd(".../Land_use")

for(year in 2006:2023){

  # Historical data are used through 2014, followed by SSP370 files.
  if (year <= 2014) {
    scenario_prefix <- "historical"
  } else if (year >= 2015) {
    scenario_prefix <- "SSP370"
  }

  nc_file_path <- paste0("Netcdf_file/Global_", scenario_prefix, "_land_use_", year, ".nc")

  nc_file <- nc_open(nc_file_path)

  rast_file <- rast(nc_file_path)
  print(rast_file)
  plot(rast_file)

  # Land-use rasters are aggregated 
  rast_file1 <- aggregate(rast_file, fact = 2, fun = mean, na.rm = TRUE)

  print(rast_file1)
  plot(rast_file1)

  shapf <- load(".../map_selected.Rdata")
  shapf
  geometry = map_selected$geometry
  sf_map <- st_sf(map_selected, crs = 4326)
  shapefile <- st_transform(sf_map, crs = "+proj=longlat +datum=WGS84")

  crs(rast_file1) <- "EPSG:4326"

  rast_EU <- crop(rast_file1, extent(shapefile))
  print(rast_EU)
  plot(rast_EU)

  crs(rast_EU) <- crs(shapefile)

  # Area-weighted means are computed while retaining partially covered cells.
  cell_area <- cellSize(rast_EU, unit = "km")

  output_data <- exact_extract(
    rast_EU, shapefile,
    function(values, coverage_frac, weights) {                     ## coverage fraction is the fraction of cell or grid inside polygon
      weighted.mean(values, coverage_frac * weights, na.rm = TRUE)
    },
    weights = cell_area,
    stack_apply = TRUE
  )

  extracted_df <- as.data.frame(output_data)

  extracted_df$region <- shapefile$LocationNa

  final_df <- data.frame(Location = shapefile$LocationNa, extracted_df)

  # Cropland is defined as the mean across crop-related land-use classes.
  final_df <- final_df %>%
    mutate(
      `croplands` = rowMeans(
        select(., 7:11)          # Select columns 7 to 11
      )
    )

  # Only the land-use variables retained for modelling are exported.
  final_df <- final_df %>%
    select(
      Location = Location,
      `primary foreseted area` = fun.primf,
      `primary non foreseted area` = fun.primn,
      `secondary foreseted area` = fun.secdf,
      `croplands` = croplands,
      `pasture land` = fun.pastr,
      `rangelands` = fun.range,
      `urban land` = fun.urban
    )

  final_df

  output_file <- paste0("land_use/land_use_", year, ".csv")

  write.csv(final_df, output_file)

  cat("year",year, "done")

}

