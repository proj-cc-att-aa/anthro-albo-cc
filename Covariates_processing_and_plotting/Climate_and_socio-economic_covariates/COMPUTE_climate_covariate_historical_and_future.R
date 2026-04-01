# ================== Climate covariate aggregation and preparation ==================

setwd(".../VectorNet_ 2_original")

# Load functions and supporting objects
source("functions_general.R")
source("functions_plotting.R")

# Define helper functions
compute_moving_average <- function(covariate, n_moving_average=30) {

  ma <- function(x, n){stats::filter(x, rep(1 / n, n), sides = 2)}

  d <- dim(covariate)
  covariate_ma <- matrix(NA, nrow=d[1], ncol=d[2]-1)

  for (j in 1:d[1]) {

    series <- unlist(covariate[j,2:d[2]])

    covariate_ma[j,] <- ma(series, n=n_moving_average)
  }

  return(list(moving_average_matrix = covariate_ma,
              moving_average_summ =
                tibble(
                       moving_average_mean=apply(covariate_ma, 1, function(x) {mean(x, na.rm=TRUE)}))))
}

do_plot <- FALSE

out <- attach(".../Data/data_INLA.RData")
map_inla <- map_inla
detach(".../file:Data/data_INLA.RData")

winter_months_pattern  <- "X.....(12|0[12])"

spring_months_pattern  <- "X.....0[3-5]"

summer_months_pattern  <- "X.....0[6-8]"

autumn_months_pattern  <- "X.....(09|1[01])"

# Aggregate yearly climate summaries
tb_list <- NULL
count <- 0

for (year in 2010:2021) {

  print(sprintf("Processing year: %d", year))

  precipitation <-

    read.csv(paste0(".../precipitation_mean_", year, ".csv")) %>%
    as_tibble()

  if ("X" %in% colnames(precipitation)) {
    precipitation <- precipitation %>% dplyr::select(-X)
  }

  print(precipitation)

  vals <- tibble(Location=precipitation$Location,
                 LocationNa=map_inla$LocationNa) %>%
    mutate(id=1:n())

  vals$is_equal <- NA

  for (j in 1:nrow(vals)) {
    vals$is_equal[j] <- identical(vals$Location[j],vals$LocationNa[j])
  }

  vals %>% filter(!is_equal)

  precipitation <- precipitation %>% mutate(
    Location = case_when(
      Location == "Cote-d?Or" ~ "Cote-d’Or",
      Location == "Cotes-d?Armor" ~ "Cotes-d’Armor",
      TRUE ~ Location
    )
  )
  
  #--- Subset based on season (Winter / Spring / Summer / Autumn)
  winter_months <- names(precipitation) %>% str_subset(winter_months_pattern)
  spring_months <- names(precipitation) %>% str_subset(spring_months_pattern)
  summer_months <- names(precipitation) %>% str_subset(summer_months_pattern)
  autumn_months <- names(precipitation) %>% str_subset(autumn_months_pattern)

  precipitation_winter <- precipitation %>% dplyr::select(Location, matches(winter_months))
  precipitation_spring <- precipitation %>% dplyr::select(Location, matches(spring_months))
  precipitation_summer <- precipitation %>% dplyr::select(Location, matches(summer_months))
  precipitation_autumn <- precipitation %>% dplyr::select(Location, matches(autumn_months))

  # Annual computation
  out <- compute_moving_average(precipitation)

  tb_precipitation_summ <- out$moving_average_summ
  rm(out)

  d <- dim(precipitation)
  precipitation_matrix <- as.matrix(precipitation[,2:d[2]])

  tb_precipitation <- tibble(
    geo=map_inla$geo,
    LocationNa=precipitation$Location,
    precipitation_mean = apply(precipitation_matrix, 1, function(x) {mean(x, na.rm=TRUE)})
  )

  # WINTER season
  out <- compute_moving_average(precipitation_winter)
  tb_precipitation_winter_summ <- out$moving_average_summ
  rm(out)

  d <- dim(precipitation_winter)
  precipitation_winter_matrix <- as.matrix(precipitation_winter[, 2:d[2]])

  tb_precipitation <- bind_cols(
    tb_precipitation,
    tibble(
      precipitation_winter_mean    = apply(precipitation_winter_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )
 
  # SPRING season
  out <- compute_moving_average(precipitation_spring)
  tb_precipitation_spring_summ <- out$moving_average_summ
  rm(out)

  d <- dim(precipitation_spring)
  precipitation_spring_matrix <- as.matrix(precipitation_spring[, 2:d[2]])

  tb_precipitation <- bind_cols(
    tb_precipitation,
    tibble(
      precipitation_spring_mean    = apply(precipitation_spring_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  #SUMMER season
  out <- compute_moving_average(precipitation_summer)
  tb_precipitation_summer_summ <- out$moving_average_summ
  rm(out)

  d <- dim(precipitation_summer)
  precipitation_summer_matrix <- as.matrix(precipitation_summer[, 2:d[2]])

  tb_precipitation <- bind_cols(
    tb_precipitation,
    tibble(
      precipitation_summer_mean    = apply(precipitation_summer_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  #AUTUMN season
  out <- compute_moving_average(precipitation_autumn)
  tb_precipitation_autumn_summ <- out$moving_average_summ
  rm(out)

  d <- dim(precipitation_autumn)
  precipitation_autumn_matrix <- as.matrix(precipitation_autumn[, 2:d[2]])

  tb_precipitation <- bind_cols(
    tb_precipitation,
    tibble(
      precipitation_autumn_mean    = apply(precipitation_autumn_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  #---Plot if TRUE
  # if (do_plot) {
  #   plot_sf(left_join(map_inla, tb_precipitation), "precipitation_ma_median",
  #           color_transformation=power_transformation(0.25))
  #   plot_sf(left_join(map_inla, tb_precipitation), "precipitation_summer_ma_median",
  #           color_transformation=power_transformation(0.25))
  #   plot_sf(left_join(map_inla, tb_precipitation), "precipitation_warm_ma_median",
  #           color_transformation=power_transformation(0.25))
  #   plot_sf(left_join(map_inla, tb_precipitation), "precipitation_cold_ma_median",
  #           color_transformation=power_transformation(0.25))
  # 
  #   plot_sf(left_join(map_inla, tb_precipitation), "precipitation_ma_min",
  #           color_transformation=power_transformation(0.25))
  #   plot_sf(left_join(map_inla, tb_precipitation), "precipitation_summer_ma_min",
  #           color_transformation=power_transformation(0.25))
  #   plot_sf(left_join(map_inla, tb_precipitation), "precipitation_warm_ma_min",
  #           color_transformation=power_transformation(0.25))
  #   plot_sf(left_join(map_inla, tb_precipitation), "precipitation_cold_ma_min",
  #           color_transformation=power_transformation(0.25))
  # }

  ## Similarly for temperature
  temperature <-

    read.csv(paste0(".../temperature_mean_", year, ".csv")) %>%
    as_tibble()

  if ("X" %in% colnames(temperature)) {
    temperature <- temperature %>% dplyr::select(-X)
  }

  temperature

  vals <- tibble(Location=temperature$Location,
                 LocationNa=map_inla$LocationNa) %>%
    mutate(id=1:n())

  vals$is_equal <- NA
  for (j in 1:nrow(vals)) {
    vals$is_equal[j] <- identical(vals$Location[j],vals$LocationNa[j])
  }

  vals %>% filter(!is_equal)

  temperature <- temperature %>% mutate(
    Location = case_when(
      Location == "Cote-d?Or" ~ "Cote-d’Or",
      Location == "Cotes-d?Armor" ~ "Cotes-d’Armor",
      TRUE ~ Location
    )
  )

  winter_months <- names(temperature) %>% str_subset(winter_months_pattern)
  spring_months <- names(temperature) %>% str_subset(spring_months_pattern)
  summer_months <- names(temperature) %>% str_subset(summer_months_pattern)
  autumn_months <- names(temperature) %>% str_subset(autumn_months_pattern)

  temperature_winter <- temperature %>% dplyr::select(Location, matches(winter_months))
  temperature_spring <- temperature %>% dplyr::select(Location, matches(spring_months))
  temperature_summer <- temperature %>% dplyr::select(Location, matches(summer_months))
  temperature_autumn <- temperature %>% dplyr::select(Location, matches(autumn_months))

  out <- compute_moving_average(temperature)
  tb_temperature_summ <- out$moving_average_summ
  rm(out)

  d <- dim(temperature)
  temperature_matrix <- as.matrix(temperature[,2:d[2]])

  tb_temperature <- tibble(
    geo=map_inla$geo,
    LocationNa=temperature$Location,
    temperature_mean=apply(temperature_matrix, 1, function(x) {mean(x, na.rm=TRUE)})
  )

  out <- compute_moving_average(temperature_winter)
  tb_temperature_winter_summ <- out$moving_average_summ
  rm(out)

  d <- dim(temperature_winter)
  temperature_winter_matrix <- as.matrix(temperature_winter[, 2:d[2]])

  tb_temperature <- bind_cols(
    tb_temperature,
    tibble(
      temperature_winter_mean    = apply(temperature_winter_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  out <- compute_moving_average(temperature_spring)
  tb_temperature_spring_summ <- out$moving_average_summ
  rm(out)

  d <- dim(temperature_spring)
  temperature_spring_matrix <- as.matrix(temperature_spring[, 2:d[2]])

  tb_temperature <- bind_cols(
    tb_temperature,
    tibble(
      temperature_spring_mean    = apply(temperature_spring_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  out <- compute_moving_average(temperature_summer)
  tb_temperature_summer_summ <- out$moving_average_summ
  rm(out)

  d <- dim(temperature_summer)
  temperature_summer_matrix <- as.matrix(temperature_summer[, 2:d[2]])

  tb_temperature <- bind_cols(
    tb_temperature,
    tibble(
      temperature_summer_mean    = apply(temperature_summer_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  out <- compute_moving_average(temperature_autumn)
  tb_temperature_autumn_summ <- out$moving_average_summ
  rm(out)

  d <- dim(temperature_autumn)
  temperature_autumn_matrix <- as.matrix(temperature_autumn[, 2:d[2]])

  tb_temperature <- bind_cols(
    tb_temperature,
    tibble(
      temperature_autumn_mean    = apply(temperature_autumn_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  # if (do_plot) {
  #   cap_temp <- 18
  #   cap_temp_cold <- -3
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_ma_max_capped=
  #                                ifelse(temperature_ma_max<cap_temp,
  #                                       cap_temp,
  #                                       temperature_ma_max))),
  #           "temperature_ma_max_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_summer_ma_max_capped=
  #                                ifelse(temperature_summer_ma_max<cap_temp,
  #                                       cap_temp,
  #                                       temperature_summer_ma_max))),
  #           "temperature_summer_ma_max_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_warm_ma_max_capped=
  #                                ifelse(temperature_warm_ma_max<cap_temp,
  #                                       cap_temp,
  #                                       temperature_warm_ma_max))),
  #           "temperature_warm_ma_max_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_cold_ma_max_capped=
  #                                ifelse(temperature_cold_ma_max<cap_temp_cold,
  #                                       cap_temp_cold,
  #                                       temperature_cold_ma_max))),
  #           "temperature_cold_ma_max_capped",
  #           color_transformation=power_transformation(1))
  # 
  #   cap_temp <- 10
  #   cap_temp_cold <- -3
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_ma_median_capped=
  #                                ifelse(temperature_ma_median<cap_temp,
  #                                       cap_temp,
  #                                       temperature_ma_median))),
  #           "temperature_ma_median_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_summer_ma_median_capped=
  #                                ifelse(temperature_summer_ma_median<cap_temp,
  #                                       cap_temp,
  #                                       temperature_summer_ma_median))),
  #           "temperature_summer_ma_median_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_warm_ma_median_capped=
  #                                ifelse(temperature_warm_ma_median<cap_temp,
  #                                       cap_temp,
  #                                       temperature_warm_ma_median))),
  #           "temperature_warm_ma_median_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_cold_ma_median_capped=
  #                                ifelse(temperature_cold_ma_median<cap_temp_cold,
  #                                       cap_temp_cold,
  #                                       temperature_cold_ma_median))),
  #           "temperature_cold_ma_median_capped",
  #           color_transformation=power_transformation(1))
  # 
  #   cap_temp <- 3
  #   cap_temp_cold <- -3
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_ma_min_capped=
  #                                ifelse(temperature_ma_min<cap_temp_cold,
  #                                       cap_temp_cold,
  #                                       temperature_ma_min))),
  #           "temperature_ma_min_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_summer_ma_min_capped=
  #                                ifelse(temperature_summer_ma_min<cap_temp,
  #                                       cap_temp,
  #                                       temperature_summer_ma_min))),
  #           "temperature_summer_ma_min_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_warm_ma_min_capped=
  #                                ifelse(temperature_warm_ma_min<cap_temp,
  #                                       cap_temp,
  #                                       temperature_warm_ma_min))),
  #           "temperature_warm_ma_min_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_temperature %>%
  #                       mutate(temperature_cold_ma_min_capped=
  #                                ifelse(temperature_cold_ma_min<cap_temp_cold,
  #                                       cap_temp_cold,
  #                                       temperature_cold_ma_min))),
  #           "temperature_cold_ma_min_capped",
  #           color_transformation=power_transformation(1))
  # }

  
  
  ## Similarly for relative humidity
  relative_humidity <-

    read.csv(paste0(".../relative_humidity_mean_", year, ".csv")) %>%
    as_tibble()

  if ("X" %in% colnames(relative_humidity)) {
    relative_humidity <- relative_humidity %>% dplyr::select(-X)
  }

  relative_humidity

  vals <- tibble(Location=relative_humidity$Location,
                 LocationNa=map_inla$LocationNa) %>%
    mutate(id=1:n())

  vals$is_equal <- NA
  for (j in 1:nrow(vals)) {
    vals$is_equal[j] <- identical(vals$Location[j],vals$LocationNa[j])
  }

  vals %>% filter(!is_equal)

  relative_humidity <- relative_humidity %>% mutate(
    Location = case_when(
      Location == "Cote-d?Or" ~ "Cote-d’Or",
      Location == "Cotes-d?Armor" ~ "Cotes-d’Armor",
      TRUE ~ Location
    )
  )

  winter_months <- names(relative_humidity) %>% str_subset(winter_months_pattern)
  spring_months <- names(relative_humidity) %>% str_subset(spring_months_pattern)
  summer_months <- names(relative_humidity) %>% str_subset(summer_months_pattern)
  autumn_months <- names(relative_humidity) %>% str_subset(autumn_months_pattern)

  relative_humidity_winter <- relative_humidity %>% dplyr::select(Location, matches(winter_months))
  relative_humidity_spring <- relative_humidity %>% dplyr::select(Location, matches(spring_months))
  relative_humidity_summer <- relative_humidity %>% dplyr::select(Location, matches(summer_months))
  relative_humidity_autumn <- relative_humidity %>% dplyr::select(Location, matches(autumn_months))

  out <- compute_moving_average(relative_humidity)
  tb_relative_humidity_summ <- out$moving_average_summ
  rm(out)

  d <- dim(relative_humidity)
  relative_humidity_matrix <- as.matrix(relative_humidity[,2:d[2]])

  tb_relative_humidity <- tibble(
    geo=map_inla$geo,
    LocationNa=relative_humidity$Location,
    relative_humidity_mean = apply(relative_humidity_matrix, 1, function(x) {mean(x, na.rm=TRUE)})
  )

  out <- compute_moving_average(relative_humidity_winter)
  tb_relative_humidity_winter_summ <- out$moving_average_summ
  rm(out)

  d <- dim(relative_humidity_winter)
  relative_humidity_winter_matrix <- as.matrix(relative_humidity_winter[, 2:d[2]])

  tb_relative_humidity <- bind_cols(
    tb_relative_humidity,
    tibble(
      relative_humidity_winter_mean    = apply(relative_humidity_winter_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  out <- compute_moving_average(relative_humidity_spring)
  tb_relative_humidity_spring_summ <- out$moving_average_summ
  rm(out)

  d <- dim(relative_humidity_spring)
  relative_humidity_spring_matrix <- as.matrix(relative_humidity_spring[, 2:d[2]])

  tb_relative_humidity <- bind_cols(
    tb_relative_humidity,
    tibble(
      relative_humidity_spring_mean    = apply(relative_humidity_spring_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  out <- compute_moving_average(relative_humidity_summer)
  tb_relative_humidity_summer_summ <- out$moving_average_summ
  rm(out)

  d <- dim(relative_humidity_summer)
  relative_humidity_summer_matrix <- as.matrix(relative_humidity_summer[, 2:d[2]])

  tb_relative_humidity <- bind_cols(
    tb_relative_humidity,
    tibble(
      relative_humidity_summer_mean    = apply(relative_humidity_summer_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  out <- compute_moving_average(relative_humidity_autumn)
  tb_relative_humidity_autumn_summ <- out$moving_average_summ
  rm(out)

  d <- dim(relative_humidity_autumn)
  relative_humidity_autumn_matrix <- as.matrix(relative_humidity_autumn[, 2:d[2]])

  tb_relative_humidity <- bind_cols(
    tb_relative_humidity,
    tibble(
      relative_humidity_autumn_mean    = apply(relative_humidity_autumn_matrix, 1, function(x) { mean(x, na.rm = TRUE) })
    )
  )

  # if (do_plot) {
  #   cap_humidity <- 0
  #   cap_humidity_cold <- 0
  #   plot_sf(left_join(map_inla, tb_relative_humidity %>%
  #                       mutate(relative_humidity_ma_median_capped=
  #                                ifelse(relative_humidity_ma_median<cap_humidity,
  #                                       cap_humidity,
  #                                       relative_humidity_ma_median))),
  #           "relative_humidity_ma_median_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_relative_humidity %>%
  #                       mutate(relative_humidity_summer_ma_median_capped=
  #                                ifelse(relative_humidity_summer_ma_median<cap_humidity,
  #                                       cap_humidity,
  #                                       relative_humidity_summer_ma_median))),
  #           "relative_humidity_summer_ma_median_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_relative_humidity %>%
  #                       mutate(relative_humidity_warm_ma_median_capped=
  #                                ifelse(relative_humidity_warm_ma_median<cap_humidity,
  #                                       cap_humidity,
  #                                       relative_humidity_warm_ma_median))),
  #           "relative_humidity_warm_ma_median_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_relative_humidity %>%
  #                       mutate(relative_humidity_cold_ma_median_capped=
  #                                ifelse(relative_humidity_cold_ma_median<cap_humidity_cold,
  #                                       cap_humidity_cold,
  #                                       relative_humidity_cold_ma_median))),
  #           "relative_humidity_cold_ma_median_capped",
  #           color_transformation=power_transformation(1))
  # 
  #   plot_sf(left_join(map_inla, tb_relative_humidity %>%
  #                       mutate(relative_humidity_ma_min_capped=
  #                                ifelse(relative_humidity_ma_min<cap_humidity,
  #                                       cap_humidity,
  #                                       relative_humidity_ma_min))),
  #           "relative_humidity_ma_min_capped",
  #           color_transformation=power_transformation(1))
  #   plot_sf(left_join(map_inla, tb_relative_humidity %>%
  #                       mutate(relative_humidity_warm_ma_min_capped=
  #                                ifelse(relative_humidity_warm_ma_min<cap_humidity,
  #                                       cap_humidity,
  #                                       relative_humidity_warm_ma_min))),
  #           "relative_humidity_warm_ma_min_capped",
  #           color_transformation=power_transformation(1))
  #}

  # Collect everything
  climate_summaries_region <- tb_precipitation %>% left_join(tb_temperature) %>%
    left_join(tb_relative_humidity)

  count <- count + 1
  tb_list[[count]] <- climate_summaries_region
}



# Collect all years
n_tables <- count
n_rows <- nrow(climate_summaries_region)

tb <- climate_summaries_region[,1:2]

# Check for missing values
na_columns <- colnames(climate_summaries_region)[colSums(is.na(climate_summaries_region)) > 0]
print(na_columns)

# Prepare yearly output
climate_summaries_region <- tb

tb_list <- lapply(seq_along(tb_list), function(i) {
  tb_list[[i]] %>% mutate(year_idx = i)
})

combined_df <- bind_rows(tb_list)

combined_df <- combined_df %>%
  relocate(year_idx, .after = geo)

print(combined_df)

unique(combined_df$year_idx)

#checking the dataset created 
combined_df[combined_df$year_idx == 14,]
tb_list[[14]]

climate_summaries_region <- combined_df

save("climate_summaries_region",

     file=paste0(".../Future_data/Final_GFDL-ESM4_ssp126_climate_data_yearly.RData"))







# Aggregate covariates for all scenario–model (ALL scenarios × ALL models) combinations for future climate dataset
options(stringsAsFactors = FALSE)

input_root  <- ".../SSPs&Pi_Future_climate_data"

output_root <- ".../Future_data"

map_inla_path <- ".../Data/data_INLA.RData"

scenarios <- c("Pi_control", "SSP1", "SSP3", "SSP5")
models    <- c("GFDL_model", "IPSL_model", "MPI_model", "MRI_model", "UKESM_model")

years_keep <- 2015:2100

load(map_inla_path)
if (!exists("map_inla")) stop("map_inla not found in: ", map_inla_path)

# Season patterns 
winter_months_pattern  <- "X.....(12|0[12])"
spring_months_pattern  <- "X.....0[3-5]"
summer_months_pattern  <- "X.....0[6-8]"
autumn_months_pattern  <- "X.....(09|1[01])"


# Some helper functions

# Define fast I/O and aggregation helpers
read_csv_fast <- function(path) {

  if (requireNamespace("data.table", quietly = TRUE)) {

    df <- data.table::fread(
      path,
      header      = FALSE,          # important: match read.csv behavior
      data.table  = FALSE,
      check.names = FALSE
    )

  } else {

    df <- read.csv(
      path,
      header           = FALSE,
      check.names      = FALSE,
      stringsAsFactors = FALSE
    )

  }

  tibble::as_tibble(df)
}



# Clean CSV: remove index column if present, and fix encoding of known problematic names
clean_daily_csv <- function(df) {

  if ("Location" %in% names(df)) {

  } else {

    first_cell <- as.character(df[[1]][1])

    if (!identical(first_cell, "Location")) {
      stop(
        "CSV does not contain 'Location' column and first cell is not 'Location'.\n",
        "First cell = ", first_cell, "\n",
        "So this file is neither normal-header nor header-in-first-row."
      )
    }

    new_names <- as.character(df[1, , drop = TRUE])
    df <- df[-1, , drop = FALSE]

    names(df) <- make.names(new_names, unique = TRUE)
  }

  # Remove index column if present
  if ("X" %in% colnames(df)) {
    df <- df %>% dplyr::select(-X)
  }

  if (!("Location" %in% names(df))) {
    stop("Header fix failed: still no 'Location' column after cleaning.")
  }

  num_cols <- setdiff(names(df), "Location")
  df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

  df %>%
    dplyr::mutate(
      Location = dplyr::case_when(
        Location == "Cote-d?Or"     ~ "Cote-d’Or",
        Location == "Cotes-d?Armor" ~ "Cotes-d’Armor",
        TRUE ~ Location
      )
    )
}


# select the days (for rowmeans later) according to the patterns of seasons
select_season <- function(df, pattern) {
  cols <- names(df) %>% str_subset(pattern)
  if (length(cols) == 0) {
    stop("No columns matched season pattern: ", pattern, " (check column name format)")
  }
  df %>% dplyr::select(Location, all_of(cols))
}


# Summarise one variable (precip / temp / rh): annual mean + seasonal means
summarise_variable <- function(df_daily, var_prefix, map_inla) {

  idx <- match(map_inla$LocationNa, df_daily$Location)
  if (any(is.na(idx))) {
    warning("Some Location names in map_inla not found in CSV for ", var_prefix,
            ". Proceeding without reordering; please verify Location alignment.")
  } else {
    df_daily <- df_daily[idx, , drop = FALSE]
  }

  d  <- dim(df_daily)
  mat_annual <- as.matrix(df_daily[, 2:d[2], drop = FALSE])
  annual_mean <- rowMeans(mat_annual, na.rm = TRUE)

  df_winter <- select_season(df_daily, winter_months_pattern)
  df_spring <- select_season(df_daily, spring_months_pattern)
  df_summer <- select_season(df_daily, summer_months_pattern)
  df_autumn <- select_season(df_daily, autumn_months_pattern)

  winter_mean <- rowMeans(as.matrix(df_winter[, -1, drop = FALSE]), na.rm = TRUE)
  spring_mean <- rowMeans(as.matrix(df_spring[, -1, drop = FALSE]), na.rm = TRUE)
  summer_mean <- rowMeans(as.matrix(df_summer[, -1, drop = FALSE]), na.rm = TRUE)
  autumn_mean <- rowMeans(as.matrix(df_autumn[, -1, drop = FALSE]), na.rm = TRUE)

  tibble(
    geo        = map_inla$geo,
    LocationNa = df_daily$Location,
    !!paste0(var_prefix, "_mean")         := annual_mean,
    !!paste0(var_prefix, "_winter_mean")  := winter_mean,
    !!paste0(var_prefix, "_spring_mean")  := spring_mean,
    !!paste0(var_prefix, "_summer_mean")  := summer_mean,
    !!paste0(var_prefix, "_autumn_mean")  := autumn_mean
  )
}


# Build full climate table for ONE year (precip + temp + rh)
process_one_year <- function(year, base_dir, map_inla) {

  prec_dir <- file.path(base_dir, "precipitation")
  temp_dir <- file.path(base_dir, "temperature_mean")
  rh_dir   <- file.path(base_dir, "relative_humidity")

  f_prec <- file.path(prec_dir, sprintf("precipitation_mean_%d.csv", year))
  f_temp <- file.path(temp_dir, sprintf("temperature_mean_%d.csv", year))
  f_rh   <- file.path(rh_dir,   sprintf("relative_humidity_mean_%d.csv", year))

  if (!file.exists(f_prec)) stop("Missing precipitation CSV: ", f_prec)
  if (!file.exists(f_temp)) stop("Missing temperature CSV: ", f_temp)
  if (!file.exists(f_rh))   stop("Missing relative humidity CSV: ", f_rh)

  precipitation      <- clean_daily_csv(read_csv_fast(f_prec))
  temperature        <- clean_daily_csv(read_csv_fast(f_temp))
  relative_humidity  <- clean_daily_csv(read_csv_fast(f_rh))

  tb_precipitation <- summarise_variable(precipitation,     "precipitation",     map_inla)
  tb_temperature   <- summarise_variable(temperature,       "temperature",       map_inla)
  tb_rh            <- summarise_variable(relative_humidity, "relative_humidity", map_inla)

  climate_summaries_region_year <- tb_precipitation %>%
    left_join(tb_temperature,   by = c("geo", "LocationNa")) %>%
    left_join(tb_rh,            by = c("geo", "LocationNa"))

  climate_summaries_region_year
}


# Run ONE scenario × model ( for all future years 2015..2100) and save
run_scenario_model <- function(sc, md) {

  message("\n============================================================")
  message("Running aggregation for scenario = ", sc, " | model = ", md)
  message("Years: ", min(years_keep), " - ", max(years_keep))
  message("============================================================")

  base_dir <- file.path(input_root, sc, md)

  if (!dir.exists(base_dir)) {
    warning("Skipping (missing base_dir): ", base_dir)
    return(invisible(NULL))
  }
  
  # Output folder (output directory creation, recursive directory)
  out_dir <- file.path(output_root, sc, md)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  out_file <- file.path(
    out_dir,
    sprintf("Final_%s_%s_climate_data_yearly_%d_%d.RData",
            sc, md, min(years_keep), max(years_keep))
  )

  tb_list <- vector("list", length(years_keep))

  for (i in seq_along(years_keep)) {
    year <- years_keep[i]
    message("  -> Processing year: ", year)

    tb_list[[i]] <- tryCatch(
      process_one_year(year, base_dir, map_inla),
      error = function(e) {
        message("  ⚠️ Year failed: ", year, " | ", conditionMessage(e))
        return(NULL)
      }
    )

    # free the memory
    if (i %% 5 == 0) gc()
  }

  tb_list <- tb_list[!vapply(tb_list, is.null, logical(1))]
  if (length(tb_list) == 0) {
    warning("No years succeeded for scenario=", sc, " model=", md)
    return(invisible(NULL))
  }

  tb_list <- lapply(seq_along(tb_list), function(i) {
    tb_list[[i]] %>% mutate(year_idx = i)
  })

  combined_df <- bind_rows(tb_list) %>%
    relocate(year_idx, .after = geo)

  climate_summaries_region <- combined_df

  save("climate_summaries_region", file = out_file)
  message("✅ Saved: ", out_file)

  invisible(NULL)
}



# Run the full scenario–model workflow
for (sc in scenarios) {
  for (md in models) {
    run_scenario_model(sc, md)
  }
}



message("\n🎉 All scenario×model aggregations finished.")








#########  Plot climate maps for factual and counterfactual scenario
source("functions_plot_article.R")



out = load(".../Current_data/Final_20CRV3-ERA5_factual_climate_data_yearly.RData")
combined_factual_data <- climate_summaries_region


#Since climate data counterfactual is from year 2001 to 2021 (year_idx is from 1 to 21), so taking only the data from year 2011 to 2021 
combined_factual_data <- combined_factual_data %>%
  filter(year_idx > 10) %>%
  mutate(year_idx = year_idx - 10)

unique(combined_factual_data$year_idx)

combined_factual_data <- combined_factual_data %>%
  filter(year_idx <= 11)

unique(combined_factual_data$year_idx)

combined_factual_data <- combined_factual_data %>%
  mutate(year = 2011 + year_idx - 1) %>%
  relocate(year, .after = geo)

unique(combined_factual_data$year)

out = load(".../Current_data/Final_20CRV3-ERA5_counterfactual_climate_data_yearly.RData")
combined_counterfactual_data <- climate_summaries_region

unique(combined_counterfactual_data$year_idx)

combined_counterfactual_data <- combined_counterfactual_data %>%
  filter(year_idx > 10) %>%
  mutate(year_idx = year_idx - 10)

unique(combined_counterfactual_data$year_idx)

combined_counterfactual_data <- combined_counterfactual_data %>%
  mutate(year = 2011 + year_idx - 1) %>%     # map 1→2011, 2→2012, ..., 11→2021
  relocate(year, .after = geo)

unique(combined_counterfactual_data$year)


# ──  choices ───────────────────────────────────────────────────────────────
# <- numeric year to plot
#year      <- 2017   ## for any particular year
# year <- NULL  for aggregation across all years
year <- NULL                           


## variable
#variable  <- "temperature"             # <- prefix or full name ("temperature", precipitation, etc.)
#variable  <- "precipitation"  
variable <- "relative_humidity"

filename_stub_hist <- paste0(".../map_plot_covariates", variable,"/",variable)

filename_stub_diff <- paste0(".../diff_(factual-counterfactual)_", variable,"/",variable)


# RELATIVE HUMIDITY (change as per the variable)
palette_name         <- "BuGn"
palette_name_diff     <- "PuOr"
palette_reverse_rh_diff  <- FALSE

## FOR TEMPERATURE
#palette_name         = "YlOrRd"  # historical (sequential)
#palette_name_diff    = "PuBuGn"  # diff (sequential; hist- pi_control)
#palette_reverse_precip_diff <- FALSE


## FOR PRECIPITAION 
# palette_name      <- "Blues"
# palette_name_diff <- "BrBG"
# palette_reverse_precip_diff <- FALSE   # if we want +ΔP wetter = blue


plot_adjusted_and_delta_maps(
  adjusted_df = combined_factual_data,
  pi_df       = combined_counterfactual_data,
  map_inla    = map_inla,
  year        = year,
  variable    = variable,

  filename_stub = filename_stub_hist,
  ext           = "tiff",
  width_in      = 3.5, height_in = 4,
  dpi           = 600,
  compression   = "lzw",
  bg            = "transparent",

  hist_legend_ext         = "tiff",
  hist_legend_dpi         = 600,
  hist_legend_filename    = NULL,
  hist_legend_width_in    = 1.4,
  hist_legend_height_in   = 5.0,
  hist_legend_compression = "lzw",
  hist_legend_bg          = "transparent",

  filename_stub_diff      = filename_stub_diff,

  diff_legend_ext         = "tiff",
  diff_legend_dpi         = 600,
  diff_legend_filename    = NULL,
  diff_legend_width_in    = 1.4,
  diff_legend_height_in   = 5.0,
  diff_legend_compression = "lzw",
  diff_legend_bg          = "transparent",

  na_fill            = "grey95",
  border_col         = NA,
  border_size        = NA,
  outer_boundary_col = "black",
  outer_boundary_size= 0.1,

  palette_name         = palette_name,
  palette_n            = 100,
  palette_reverse      = FALSE,
  palette_name_diff    = palette_name_diff,
  palette_n_diff            = 100,
  palette_reverse_diff = palette_reverse_precip_diff
)
# If a warning appears that means we are getting different legends for each seasonal aggregation of climate variable.
