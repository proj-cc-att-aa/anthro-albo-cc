# ================== Future relative humidity NetCDF extraction to regional CSV files ==================

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

for(year in 2021:2080){             ### take care of time values for 2010, origin is 2001-01-01, for 2011 to 2020 (origin is 2011-01-01) and from 2021 to 2024 (origin is 2021-01-01)

  nc_file_path <- paste0("Netcdf_file/gfdl-esm4_ssp126_hurs_global_daily_", year, ".nc")     ## for ssp370

  nc_file <- nc_open(nc_file_path)
  time_var <- ncvar_get(nc_file, "time")  # Extract time variable
  time_units <- ncatt_get(nc_file, "time", "units")$value

  origin_date <- gsub("days since ", "", time_units)  # Remove the "days since" part
  origin_date <- strsplit(origin_date, " ")[[1]][1]   # Extract the date portion

  time_values = as_datetime(c(time_var*60*60*24),origin= origin_date)

  rast_file <- rast(nc_file_path)
  print(rast_file)

  ## loading the shapefile
  shapf <- load(".../map_selected.Rdata")
  shapf
  geometry = map_selected$geometry
  sf_map <- st_sf(map_selected, crs = 4326)
  shapefile <- st_transform(sf_map, crs = "+proj=longlat +datum=WGS84")

  crs(rast_file) <- "EPSG:4326"

  #Cropping the data so that it matches with area covered by shapefile 
  rast_EU <- crop(rast_file, extent(shapefile))

  crs(rast_EU) <- crs(shapefile)

  # Cell-area raster (units as per choice; they cancel in the weighted mean)
  cell_area <- cellSize(rast_EU, unit = "km")  # same grid/CRS as rast_EU

  output_data <- exact_extract(rast_EU, shapefile,
                               function(values, coverage_frac, weights) {
                                 weighted.mean(values, coverage_frac * weights, na.rm = TRUE)     ## coverage fraction is the fraction of cell or grid inside polygon
                               },                                  ## wights are now , coverage fraction of cell inside polygon (default weight works without any secondary weightage)  * area if cell
                               weights = cell_area,                 ## Calculate area of each cell explicitly depending on latitude and longitude according to coordinate system reference (self calculated by exact_extract function by itself)
                               stack_apply = TRUE)              

  extracted_df <- as.data.frame(output_data)

  # Add region names or IDs from the shapefile to the dataframe
  extracted_df$region <- shapefile$LocationNa         # Adjust for the column name in shapefile

  final_df <- data.frame(Location = shapefile$LocationNa)

  # Add temperature data as columns with dates
  for (i in seq_along(time_values)) {
    final_df[, as.character(time_values[i])] <- output_data[, i]
  }

  output_file <- paste0("relative_humidity/relative_humidity_mean_", year, ".csv")

  write.csv(final_df, output_file)

  cat("year",year, "done")

}





############ COMBINED WORKFLOW FOR ALLL SCENARIO AND ALL MODEL. #############################
root_dir <- ".../SSPs&Pi_Future_climate_data"
map_rdata_path <- ".../map_selected.Rdata"

years_ssp       <- 2015:2100
years_picontrol <- 1995:2100

# scenarios + models (as per our requirement)
scenarios <- c("Pi_control", "SSP1", "SSP3", "SSP5")
models    <- c("GFDL_model", "IPSL_model", "MPI_model", "MRI_model", "UKESM_model")

# NOTE: folder inside each model_dir where we save CSVs = "relative_humidity"
jobs <- list(
  list(var = "hurs", out_subdir = "relative_humidity", file_prefix = "relative_humidity_mean")
)

#  LOAD & PREP SHAPEFILE ONCE
load(map_rdata_path)  # expects object: map_selected
if (!exists("map_selected")) stop("map_selected not found in the RData file.")

poly_sf <- map_selected
if (!inherits(poly_sf, "sf")) poly_sf <- st_as_sf(poly_sf)

# ensure valid geometry + WGS84
poly_sf <- st_make_valid(poly_sf)
poly_sf <- st_transform(poly_sf, 4326)

if (!("LocationNa" %in% names(poly_sf))) {
  stop("Column 'LocationNa' not found in map_selected. Update the ID column name in the code.")
}

# HELPERS

# Parse NetCDF time -> Date vector using the file's own units (no hardcoded origins!)
read_nc_dates <- function(nc_path, time_var = "time") {

  nc <- nc_open(nc_path)
  on.exit(nc_close(nc), add = TRUE)

  tvals  <- ncvar_get(nc, time_var)               ### time variable number 
  tunits <- ncatt_get(nc, time_var, "units")$value  # e.g. "days since 2015-1-1 00:00:00"

  origin_str  <- sub("^days since\\s+", "", tunits)
  origin_str  <- trimws(origin_str)
  origin_date <- strsplit(origin_str, "\\s+")[[1]][1]
  origin_date <- ymd(origin_date)                      ## setting the origin date

  if (is.na(origin_date)) stop("Could not parse origin date from time units: ", tunits)

  is_integer_days <- all(abs(tvals - round(tvals)) < 1e-10)

  if (is_integer_days) {
    return(as.Date(origin_date + days(as.integer(round(tvals)))))      ##  ## now origin date is added to t_vals(file data days since origin date of the file)
  } else {
    origin_dt <- as.POSIXct(origin_date, tz = "UTC")
    return(as.Date(origin_dt + tvals * 86400))
  }
}

# Read variable unit and fill value (so we can convert + set NAflag robustly)
read_nc_var_meta <- function(nc_path, var) {

  nc <- nc_open(nc_path)
  on.exit(nc_close(nc), add = TRUE)

  units <- ncatt_get(nc, var, "units")$value
  fill1 <- ncatt_get(nc, var, "_FillValue")$value
  fill2 <- ncatt_get(nc, var, "missing_value")$value

  fill <- NA_real_
  if (!is.null(fill1) && !is.na(fill1)) fill <- as.numeric(fill1)
  if (is.na(fill) && !is.null(fill2) && !is.na(fill2)) fill <- as.numeric(fill2)

  list(units = units, fill = fill)
}

## Detect which variable name to read (sometimes file could have only 1 var)
detect_precip_varname <- function(nc_path, preferred = "hurs") {

  nc <- nc_open(nc_path)
  on.exit(nc_close(nc), add = TRUE)

  vnames <- names(nc$var)
  if (preferred %in% vnames) return(preferred)

  if (length(vnames) == 0) stop("No variables found in NetCDF: ", nc_path)
  vnames[1]
}

## Find correct hurs NetCDF inside Netcdf_file/
find_nc_for_pr <- function(nc_dir) {

  c2 <- list.files(nc_dir, pattern = "_hurs_.*\\.nc$", full.names = TRUE, ignore.case = TRUE)
  if (length(c2) > 0) return(c2[1])

  c3 <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE, ignore.case = TRUE)
  if (length(c3) > 0) return(c3[1])

  NA_character_
}

# CORE EXTRACTOR: one NetCDF -> yearly CSVs
extract_pr_one_nc_yearly <- function(nc_path, poly_sf, out_dir, years_keep, file_prefix) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  varname <- detect_precip_varname(nc_path, preferred = "hurs")

  dates <- read_nc_dates(nc_path, time_var = "time")
  meta  <- read_nc_var_meta(nc_path, var = varname)

  r <- terra::rast(nc_path)

  # set NAflag using fill value
  if (!is.na(meta$fill)) terra::NAflag(r) <- meta$fill

  if (is.na(terra::crs(r))) terra::crs(r) <- "EPSG:4326"

  # cell area weights for weighted mean
  cell_area <- terra::cellSize(r, unit = "km")

  # map each layer to its year
  years_all <- as.integer(format(dates, "%Y"))

  for (yr in years_keep) {

    idx <- which(years_all == yr)
    if (length(idx) == 0) next

    r_y <- terra::subset(r, idx)
    d_y <- dates[idx]

    out_mat <- exact_extract(
      r_y, poly_sf,
      function(values, coverage_frac, weights) {
        weighted.mean(values, coverage_frac * weights, na.rm = TRUE)
      },
      weights     = cell_area,
      stack_apply = TRUE
    )

    colnames(out_mat) <- as.character(d_y)

    # final df: Location + daily rh columns
    out_df <- cbind(Location = poly_sf$LocationNa, as.data.frame(out_mat))

    out_file <- file.path(out_dir, sprintf("%s_%d.csv", file_prefix, yr))
    write.csv(out_df, out_file, row.names = FALSE)

    message("✅ ", basename(nc_path), " | hurs| year ", yr, " saved -> ", out_file)
  }
}

# MAIN LOOP: ALL scenarios × models (Relative humidity ONLY)
for (sc in scenarios) {

  years_keep <- if (sc == "Pi_control") years_picontrol else years_ssp

  for (md in models) {

    model_dir <- file.path(root_dir, sc, md)
    nc_dir    <- file.path(model_dir, "Netcdf_file")

    if (!dir.exists(model_dir) || !dir.exists(nc_dir)) {
      message("⚠️ Skipping (missing dirs): ", model_dir)
      next
    }

    # find rh nc for this scenario × model
    nc_path <- find_nc_for_pr(nc_dir)

    if (is.na(nc_path) || !file.exists(nc_path)) {
      message("⚠️ No relative hmidity NetCDF found in: ", nc_dir)
      next
    }

    out_dir <- file.path(model_dir, "relative_humidity")

    message("\n--- Running relative_humidity: scenario=", sc, " | model=", md, " ---")
    extract_pr_one_nc_yearly(
      nc_path     = nc_path,
      poly_sf     = poly_sf,
      out_dir     = out_dir,
      years_keep  = years_keep,
      file_prefix = "relative_humidity_mean"
    )
  }
}

message("\n🎉 All extraction done.")
