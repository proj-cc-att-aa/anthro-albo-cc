# ================== Future land-use NetCDF extraction to regional CSV files ==================

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

for(year in 2021:2100){             ### take care of time values for 2010, origin is 2001-01-01, for 2011 to 2020 (origin is 2011-01-01) and from 2021 to 2024 (origin is 2021-01-01)

  nc_file_path <- paste0("Netcdf_file/Global_ssp370_land_use_", year, ".nc")

  nc_file <- nc_open(nc_file_path)

  rast_file <- rast(nc_file_path)
  print(rast_file)
  plot(rast_file)

  # Aggregate the raster to 0.5° resolution
  rast_file1 <- aggregate(rast_file, fact = 2, fun = mean, na.rm = TRUE)    ## mean of all grids of 0.25 while aggregating it to 0.5

  print(rast_file1)
  plot(rast_file1)
  ## loading the shapefile
  shapf <- load(".../map_selected.Rdata")
  shapf
  geometry = map_selected$geometry
  sf_map <- st_sf(map_selected, crs = 4326)
  shapefile <- st_transform(sf_map, crs = "+proj=longlat +datum=WGS84")

  crs(rast_file1) <- "EPSG:4326"

  #Cropping the data so that it matches with area covered by shapefile
  rast_EU <- crop(rast_file1, extent(shapefile))
  print(rast_EU)                 ## print to check the resolution

  plot(rast_EU)

  crs(rast_EU) <- crs(shapefile)

  # Cell-area raster (units as per the choice; they cancel in the weighted mean)
  cell_area <- cellSize(rast_EU, unit = "km")  # same grid/CRS as rast_EU

  output_data <- exact_extract(rast_EU, shapefile,
                               function(values, coverage_frac, weights) {
                                 weighted.mean(values, coverage_frac * weights, na.rm = TRUE)     ## coverage fraction is the fraction of cell or grid inside polygon
                               },                                 ## weights are now, coverage fraction of cell inside polygon (default weight works without any secondary weightage)  * area of cell
                               weights = cell_area,                 ## Calculate area of each cell explicitly depending on latitude and longitude according to coordinate system reference (self calculated by exact_extract function by itself)
                               stack_apply = TRUE)             

  extracted_df <- as.data.frame(output_data)

  # Add region names or IDs from the shapefile to the dataframe
  extracted_df$region <- shapefile$LocationNa         # Adjust for the column name in the shapefile

  final_df <- data.frame(Location = shapefile$LocationNa, extracted_df)


  final_df <- final_df %>%
    mutate(
      `croplands` = rowMeans(
        dplyr::select(., 7:11) # Select columns 7 to 11
      )
    )

  # Rename columns (to remove the weighted_mean. term) except croplands
  final_df <- final_df %>%
    dplyr::select(
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
