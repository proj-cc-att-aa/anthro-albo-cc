# ================== Load and process future data ==================


# Loading and arranging the training data i.e. ERA5 land
setwd(".../Data")

# Loading some package and functions 
source(".../functions_general.R")
source(".../functions_plotting.R")

# Load the function bin+covariates 
source(".../bin_covariates.R")




# Loading the map first
out <- attach("data_INLA.RData")
map_inla <- map_inla
detach(paste0("file:", "data_INLA.RData"), character.only = TRUE) # # dettaching to save memory


# Loading the mosquito data
out <- load("mosq_dist_until_2023.RData")
mosq_dist <- mosq_dist %>% filter(date_start>=2010) %>%
  mutate(year_idx = match(date_start, sort(unique(date_start))))
tb_inla <- mosq_dist %>% filter(geo %in% map_inla$geo)
rm(mosq_dist)

data_spde_mosq <- tb_inla
map_spde <- map_inla %>% filter(geo %in% data_spde_mosq$geo)


## Load data (first historical as its bin ranges will be used for future) ----------------

# Specifying the location of climate covariates
covariates_folder_data <- "climate/Final_covariates_data/Current_data"

# Loading population 
out <- load(paste0(covariates_folder_data, "/Final_hist_population_data_yearly.RData"))
df_population <- population

# Loading mobility covariates 
out <- load(paste0(covariates_folder_data, "/Final_hist_mobility_data_yearly.RData"))
df_mobility <- df_mobility


## taking out the scale factor of mobility so that we can use the same scaling factor for future prediction data too
max_flux_in_train <- df_mobility %>%
  dplyr::summarise(m = max(total_flux_in, na.rm = TRUE)) %>%
  dplyr::pull(m)

df_mobility <- df_mobility %>%
  mutate(total_flux_in = {
    m <- suppressWarnings(max(total_flux_in, na.rm = TRUE))
    if (is.finite(m) && m > 0) total_flux_in / m else total_flux_in
  })


# Loading climate covariates 
out <- load(paste0(covariates_folder_data, "/Final_20CRV3-ERA5_factual_climate_data_yearly.RData"))
df_climate <- climate_summaries_region

unique(df_climate$year_idx)

df_climate <- df_climate %>%
  filter(year_idx > 9) %>% # drop 1..9
  mutate(year_idx = year_idx - 9) # reindex the year_idx column of left data to start at 1

unique(df_climate$year_idx)
# 1:14


# Loading land use (historical is used for training)
out <- load(paste0(covariates_folder_data, "/Final_hist_land_use_scaled_data_yearly.RData"))
df_landuse <- land_use


# Rationale: Drop `croplands` because LUH cropland categories are noisy/uncertain and consists of multiple categories
stopifnot("croplands" %in% names(df_landuse))

df_landuse <- df_landuse %>%
  dplyr::select(-croplands)


# Loading proximity 
out <- load(paste0(covariates_folder_data, "/Final_hist_proximity_data_yearly.RData"))
df_proximity <- df_proximity


# we will first save the proximity for year 2024, to use in prediction data starting in 2024 and then after that we construct recurring proximity in prediction data using posterior samples one by one
# Get the last year in df_proximity
last_year <- max(df_proximity$year_idx, na.rm = TRUE)


proximity_for_year_2024_for_prediction <- df_proximity %>%
  filter(year_idx == last_year)

# Define output path
out_path <- ".../Output/A-Final_output/Proximity_for_some years_for_prediction/Proximity_for_year_2024_for_prediction.RData"

# Save as .RData
save(proximity_for_year_2024_for_prediction, file = out_path)


# IMPORTANT - save the proximity of year 2022, in case we want to start prediction from there (since our counterfactual data ends in 2021)
# Get the year 2022 in df_proximity
year_2022_idx <- 13

proximity_for_year_2022_for_prediction <- df_proximity %>%
  filter(year_idx == year_2022_idx)

# Define output path
out_path2 <- ".../Output/A-Final_output/Proximity_for_some years_for_prediction/Proximity_for_year_2022_for_prediction.RData"

# Save as .RData
save(proximity_for_year_2022_for_prediction, file = out_path2)
# IMPORTANT


# Dropping the last year’s proximity rows
df_proximity <- df_proximity %>%
  group_by(geo) %>%
  filter(year_idx < max(year_idx)) %>% # keeps all but the last index per geo
  ungroup()

unique(df_proximity$year_idx)
# 1:14



# joining all the covariates

coerce_keys <- function(x) {
  x %>% # # checking that geo and year_idx columns exist
    mutate(
      geo = as.character(geo), # # convert all geo columns to charcter
      year_idx = as.integer(year_idx) # # year_idx columns to interger
    )
}

df_population  <- coerce_keys(df_population)  %>% distinct(geo, year_idx, .keep_all = TRUE)
df_mobility    <- coerce_keys(df_mobility)    %>% distinct(geo, year_idx, .keep_all = TRUE)
df_climate     <- coerce_keys(df_climate)     %>% distinct(geo, year_idx, .keep_all = TRUE)
df_landuse     <- coerce_keys(df_landuse)     %>% distinct(geo, year_idx, .keep_all = TRUE)
df_proximity     <- coerce_keys(df_proximity)     %>% distinct(geo, year_idx, .keep_all = TRUE)

df_population_sel <- df_population %>%
  dplyr::select(geo, year_idx, year, LocationNa, population, area, population_density)

df_mobility_sel <- df_mobility %>%
  dplyr::select(geo, year_idx, mobility_idx, net_flux_out, total_flux_in)

df_climate_sel <- df_climate %>%
  dplyr::select(-LocationNa) %>% # already in population
  dplyr::select(geo, year_idx, everything())

df_landuse_sel <- df_landuse %>%
  dplyr::select(-LocationNa, -year) %>% # already in population
  dplyr::select(geo, year_idx, everything())

df_proximity_sel <- df_proximity %>%
  dplyr::select(geo, year_idx, everything()) # already aligned keys; no LocationNa/year here

# Sanity checks on keys 
stopifnot(nrow(df_population_sel) == nrow(distinct(df_population_sel, geo, year_idx)))
stopifnot(nrow(df_mobility_sel)   == nrow(distinct(df_mobility_sel,   geo, year_idx)))
stopifnot(nrow(df_climate_sel)    == nrow(distinct(df_climate_sel,    geo, year_idx)))
stopifnot(nrow(df_landuse_sel)    == nrow(distinct(df_landuse_sel,    geo, year_idx)))
stopifnot(nrow(df_proximity_sel)  == nrow(distinct(df_proximity_sel,  geo, year_idx)))

# Join (population is the base)
covariates_all_2010_2023 <- df_population_sel %>%
  left_join(df_mobility_sel,  by = c("geo", "year_idx")) %>%
  left_join(df_climate_sel,   by = c("geo", "year_idx")) %>%
  left_join(df_landuse_sel,   by = c("geo", "year_idx")) %>%
  left_join(df_proximity_sel, by = c("geo", "year_idx")) %>%
  clean_names()

n0 <- nrow(df_population_sel)
n_all <- nrow(covariates_all_2010_2023)

message(sprintf("Rows after join: %s (base had %s)", n_all, n0)) # # total of 18690 rows for Sandra map

# Extra double check step
miss_mob   <- df_population_sel %>% anti_join(df_mobility_sel,  by = c("geo","year_idx")) %>% nrow()
miss_clim  <- df_population_sel %>% anti_join(df_climate_sel,   by = c("geo","year_idx")) %>% nrow()
miss_land  <- df_population_sel %>% anti_join(df_landuse_sel,   by = c("geo","year_idx")) %>% nrow()
miss_prox  <- df_population_sel %>% anti_join(df_proximity_sel, by = c("geo","year_idx")) %>% nrow()

message(sprintf("Unmatched rows vs mobility: %d | climate: %d | land-use: %d | proximity: %d",
                miss_mob, miss_clim, miss_land, miss_prox))

# Check the number of rows
rows_before <- nrow(covariates_all_2010_2023)

# Drop first year covariate data as it is not used 
covariates_all <- covariates_all_2010_2023 %>%
  filter(!(year_idx == 1 | year == 2010))

rows_after <- nrow(covariates_all)
message(sprintf("Dropped %d rows (first year). Remaining: %d",
                rows_before - rows_after, rows_after))

covariates_all <- covariates_all %>%
  mutate(year_idx = year_idx - 1L)

# Drop first year from mosquito data too as it is not used
data_spde1 <- data_spde_mosq %>%
  filter(!(date_start == 2010))

data_spde1 <- data_spde1 %>%
  mutate(year_idx = year_idx - 1L)

# Some adjustment of map first
map_spde$region_id <- 1:nrow(map_spde)
map_spde$spatial_idx <- map_spde$region_id # adding a unique region id

data_spde1 <- left_join(data_spde1 %>% dplyr::select(-matches("region_id")),
                        as_tibble(map_spde) %>% dplyr::select(geo, region_id)) %>%
  dplyr::arrange(year_idx, region_id) %>% relocate(geo, region_id)


# Pull the two vectors
g_spde <- data_spde1$geo
g_cov  <- covariates_all$geo

# Same length?
cat("Lengths — data_spde:", length(g_spde), "| covariates_all:", length(g_cov), "\n")

# Same unique sets?
miss_in_cov  <- setdiff(g_spde, g_cov)
miss_in_spde <- setdiff(g_cov,  g_spde)
cat("Missing in covariates_all:", length(miss_in_cov), "\n")
cat("Missing in data_spde:",      length(miss_in_spde), "\n")
if (length(miss_in_cov))  cat("Examples (in spde not in cov):",  head(miss_in_cov),  "\n")
if (length(miss_in_spde)) cat("Examples (in cov not in spde):", head(miss_in_spde), "\n")

data_spde1 <- data_spde1 %>% mutate(data_idx=1:n()) %>% # add unique data_idx for each row
  relocate(data_idx, region_id, year_idx)



## Removing some extra columns like population (not density), area, mobility_idx, net_flux_out (not total flux in) and climate annual mean (used for plottings, only seasonal aggregates used as covariates )
col_combined_cov <- colnames(covariates_all)

cols_to_drop <- c(
  "population",
  "area",
  "mobility_idx",
  "net_flux_out",
  "precipitation_mean",
  "temperature_mean",
  "relative_humidity_mean",
  "y_proximity_exp100",
  "y_proximity_exp200",
  "y_proximity_exp050"
)

# Remove them (safely) and keep everything else
covariates_all <- covariates_all %>% dplyr::select(-dplyr::any_of(cols_to_drop))

info_cols <- c("geo","year_idx","year","location_na") # first 4 columns (NOT covariates)
stopifnot(all(info_cols %in% names(covariates_all)))

# Split info vs covariates
info_df <- covariates_all %>% dplyr::select(all_of(info_cols))
cov_df  <- covariates_all %>% dplyr::select(-all_of(info_cols)) %>% dplyr::select(where(is.numeric))

# Drop rows with ANY NA in covariates
rows_before <- nrow(covariates_all)
good_idx    <- complete.cases(cov_df)
if (!any(good_idx)) stop("All rows have NA in covariates; nothing to analyze.")
dropped_rows <- rows_before - sum(good_idx) # # total number of dropped rows from entire data

# Dropping from both info_df as well as cov_df (combined)
covariates_all_clean <- dplyr::bind_cols(
  info_df[good_idx, , drop = FALSE],
  cov_df [good_idx, , drop = FALSE]
)

message(sprintf("Rows before NA filtering: %d | After: %d | Dropped: %d",
                rows_before, nrow(covariates_all_clean), dropped_rows))



# Check - same number of rows and same ordering as before
stopifnot(nrow(data_spde1) == nrow(covariates_all))
stopifnot(identical(data_spde1$geo, covariates_all$geo))
stopifnot(identical(data_spde1$year_idx,  covariates_all$year_idx))

keep_keys <- covariates_all_clean %>% dplyr::select(geo, year_idx) %>% distinct()

data_spde <- data_spde1 %>% semi_join(keep_keys, by = c("geo","year_idx"))

# checks - new data_spde with clean covairates
stopifnot(nrow(data_spde) == nrow(covariates_all_clean))
stopifnot(setequal(paste(data_spde$geo, data_spde$year_idx),
                   paste(covariates_all_clean$geo, covariates_all_clean$year_idx)))





#########CORRELATION##############
## we repeat the correlation computation in future data preperation too so that we get the same training covariates and then use this covariate list to extract the data for future covariates too

geo_order <- covariates_all_clean %>%
  distinct(geo) %>%
  mutate(.geo_order = row_number())

#### aggregating for all years (one value for each geos)
covariates_by_region <- covariates_all_clean %>%
  group_by(geo) %>%
  summarise(
    location_na = first(location_na),
    n_years     = n_distinct(year_idx),
    across(
      where(is.numeric) & !any_of(c("year_idx", "year")),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  
  left_join(geo_order, by = "geo") %>%
  arrange(.geo_order) %>%
  dplyr::select(-.geo_order)

print(covariates_by_region) 

## Now new info columns for aggregated data of covariates
info_cols2 <- c("geo","location_na","n_years", "presence_previous_year")

info_df_cleaned_agg <- covariates_by_region %>% dplyr::select(all_of(info_cols2)) # info cols defined earlier
cov_df_cleaned_agg  <- covariates_by_region %>% dplyr::select(-all_of(info_cols2)) %>% dplyr::select(where(is.numeric))



# Assess covariate correlation
C_full <- cor(cov_df_cleaned_agg, use = "pairwise.complete.obs", method = "pearson")

# Heatmap BEFORE pruning
ord_full <- hclust(as.dist(1 - abs(C_full)))$order
hm_full <- pheatmap::pheatmap(
  abs(C_full)[ord_full, ord_full],
  clustering_distance_rows = as.dist(1 - abs(C_full)),
  clustering_distance_cols = as.dist(1 - abs(C_full)),
  clustering_method = "complete",
  legend = TRUE,
  main = "",
  display_numbers = FALSE,
  silent = TRUE
)

hm_full
# Grid::grid.newpage(); grid::grid.draw(hm_full$gtable)
# Dev.off()


## Now applying the cutoff
#  Greedy pruning: drop variable with larger mean |r| (any |r| > cutoff) in the most-correlated pair and, 
## compute each correlation pair and then remove the one with low overall mean correlation |r| mean

prune_cor_matrix <- function(C, cutoff, keep = character(0)) {
  keep <- intersect(keep, colnames(C)) 
  drop_vars <- character(0)
  current_C <- C
  
  over_cutoff <- function(M, cutoff) {
    if (ncol(M) <= 1) return(FALSE)
    any(abs(M[upper.tri(M, diag = FALSE)]) > cutoff, na.rm = TRUE)
  }
  
  while (over_cutoff(current_C, cutoff)) {
    A <- abs(current_C)
    A[lower.tri(A, diag = TRUE)] <- NA
    max_idx <- which(A == max(A, na.rm = TRUE), arr.ind = TRUE)[1, ]
    v1 <- colnames(current_C)[max_idx[1]]
    v2 <- colnames(current_C)[max_idx[2]]
    
    mean_r <- function(v) mean(abs(current_C[v, setdiff(colnames(current_C), v)]), na.rm = TRUE)
    
    # Protect "keep" variables
    if (v1 %in% keep && v2 %in% keep) {
      warning(sprintf("Pair (%s, %s) exceeds cutoff but both are in `keep`; stopping pruning.", v1, v2))
      break
    } else if (v1 %in% keep) {
      drop <- v2
    } else if (v2 %in% keep) {
      drop <- v1
    } else {
      drop <- if (mean_r(v1) >= mean_r(v2)) v1 else v2
    }
    
    drop_vars <- c(drop_vars, drop)
    keep_cols <- setdiff(colnames(current_C), drop)
    current_C <- current_C[keep_cols, keep_cols, drop = FALSE]
  }
  
  list(kept = colnames(current_C), dropped = drop_vars, C_kept = current_C)
}

## keep the variables protected
protected_keep <- c("temperature_summer_mean", "population_density")
pruned <- prune_cor_matrix(C_full, cutoff = 0.7, keep = protected_keep)  ## standard cutoff

cat("\n---- Pruning summary (cutoff = 0.70, with protection) ----\n")
cat(sprintf("Initial covariates: %d\n", ncol(C_full)))
cat(sprintf("Dropped: %d\n", length(pruned$dropped)))
if (length(pruned$dropped)) cat("Dropped vars:\n  ", paste(pruned$dropped, collapse = ", "), "\n")
cat(sprintf("Kept: %d\n", length(pruned$kept)))
cat("Kept vars:\n  ", paste(pruned$kept, collapse = ", "), "\n\n")


# Again re-split our original datset...i.e.  Split info vs covariates ----
info_df_cleaned <- covariates_all_clean %>% dplyr::select(all_of(info_cols))
cov_df_cleaned  <- covariates_all_clean %>% dplyr::select(-all_of(info_cols))

#prev_col <- ("presence_previous_year")      ## no need to keep this column or covariate
cols_to_keep <- pruned$kept

cov_df_pruned <- cov_df_cleaned %>%
  dplyr::select(dplyr::any_of(cols_to_keep))

covariates_all_pruned <- dplyr::bind_cols(info_df_cleaned, cov_df_pruned)

# Changing the name to final
final_covariates <- covariates_all_pruned

# Verify the final columns in training data
colnames(final_covariates)







####################### NOW LOADING THE FUTURE PREDICTION DATA also and preparing it ###############


## Above training data will be used to create training bins for climate covariates and then we assign the prediction bins based on the same range created by training bins

## Important######
## Now since SSPs and picontrol will be compared (for any climate model), the socioeconomic covariates will be same for both ssp (any one) and picontrol data
## so we will create first data with picontrol and with corresponding ssp's socioeconomic condition and compare with actual ssp climate with that ssp's socioeconomic condition
## In this way we allow socioeconomic condition to vary according to SSP's but change the climate to picontrol (with same ssp) and ssp climate (with same ssp covariates) and compare both of them to take out only the effects of climate change 





################ FIRST PREPARING PREDICTION DATA FOR SSP126 SCENARIO #################################


# ## loading climate covariates pi control climate  for ssp126 socioeconomic merging
# out <- load(paste0(future_covariates_location, "/Final_ensemble_Pi_control_delta_corrected_climate_data_yearly_2022_2100.RData"))
# 
# df_climate_pred_2022_2080 <- climate_summaries_region_pi_control   ## contain all  prediction climate variable (precipitation, temperature and humidity) mean and seasonal annual aggregates in of NUTS3 regions in map
# 
# unique(df_climate_pred_2022_2080$year_idx)   


# loading ssp126 climate covariates for  ssp126 socioeconomic merging
future_covariates_location <- ".../Data/climate/Final_covariates_data/Future_data"

# Unique(df_climate_pred_2022_2080$year_idx)


out <- load(paste0(future_covariates_location, "/Final_ensemble_SSP1_delta_corrected_climate_data_yearly_2022_2100.RData"))

df_climate_pred_2022_2080 <- climate_summaries_region_stitched
unique(df_climate_pred_2022_2080$year_idx)

# Unique(df_climate_pred_2022_2080$year_idx)





# Loading population ssp1 future
out_pop <- load(paste0(future_covariates_location, "/Final_future_ssp126_population_data_yearly.RData"))

df_pop_pred_2021_2080 <- population

## start the prediction from year 2023
df_pop_pred_2022_2080 <- df_pop_pred_2021_2080 %>%
  filter(year_idx > 1) %>% # drop 1
  mutate(year_idx = year_idx - 1) # reindex the year_idx column of left data to start at 1

unique(df_pop_pred_2022_2080$year_idx)

df_pop_pred_2022_2080 <- df_pop_pred_2022_2080 %>%
  dplyr::rename(population_density = pop_density)

# Loading land_use ssp1 future
out_land_use <- load(paste0(future_covariates_location, "/Final_future_ssp126_land_use_data_yearly.RData"))

df_landuse_pred_2021_2080 <- land_use


df_landuse_pred_2022_2080 <- df_landuse_pred_2021_2080 %>%
  filter(year_idx > 1) %>% # drop 1
  mutate(year_idx = year_idx - 1) 

unique(df_landuse_pred_2022_2080$year_idx)

# Loading mobility ssp1 future
out_mob <- load(paste0(future_covariates_location, "/Final_future_ssp126_mobility_data_yearly.RData"))

df_mobility_pred_2021_2080 <- df_mobility


df_mobility_pred_2022_2080 <- df_mobility_pred_2021_2080 %>%
  filter(year_idx > 1) %>% # drop 1
  mutate(year_idx = year_idx - 1) # reindex the year_idx column of left data to start at 1

unique(df_mobility_pred_2022_2080$year_idx)


# Scale future total_flux_in to match the training scale
df_mobility_pred_2022_2080 <- df_mobility_pred_2022_2080 %>%
  dplyr::mutate(total_flux_in = total_flux_in / max_flux_in_train)



coerce_keys <- function(x) {
  x %>% #
    mutate(
      geo = as.character(geo), 
      year_idx = as.integer(year_idx) 
    )
}

df_population_pred  <- coerce_keys(df_pop_pred_2022_2080)  %>% distinct(geo, year_idx, .keep_all = TRUE)
df_mobility_pred    <- coerce_keys(df_mobility_pred_2022_2080)    %>% distinct(geo, year_idx, .keep_all = TRUE)
df_climate_pred     <- coerce_keys(df_climate_pred_2022_2080)     %>% distinct(geo, year_idx, .keep_all = TRUE) 
df_landuse_pred     <- coerce_keys(df_landuse_pred_2022_2080)     %>% distinct(geo, year_idx, .keep_all = TRUE)

df_population_sel <- df_population_pred %>%
  dplyr::select(geo, year_idx, year, LocationNa, population, area_km_2, population_density)

df_mobility_sel <- df_mobility_pred %>%
  dplyr::select(geo, year_idx, mobility_idx, net_flux_out, total_flux_in)

df_climate_pred_2022_2080_sel <- df_climate_pred   %>% 
  dplyr::select(-LocationNa) %>% 
  dplyr::select(geo, year_idx, everything())

df_landuse_sel <- df_landuse_pred %>%
  dplyr::select(-LocationNa, -year) %>% 
  dplyr::select(geo, year_idx, everything())

# checks on keys 
stopifnot(nrow(df_population_sel) == nrow(distinct(df_population_sel, geo, year_idx)))
stopifnot(nrow(df_mobility_sel)   == nrow(distinct(df_mobility_sel,   geo, year_idx)))
stopifnot(nrow(df_climate_pred_2022_2080_sel)    == nrow(distinct(df_climate_pred_2022_2080_sel,    geo, year_idx)))
stopifnot(nrow(df_landuse_sel)    == nrow(distinct(df_landuse_sel,    geo, year_idx)))

# Join (population is the base)
prediction_covariates_all_2022_2080 <- df_population_sel %>%
  left_join(df_mobility_sel,  by = c("geo", "year_idx")) %>%
  left_join(df_climate_pred_2022_2080_sel,   by = c("geo", "year_idx")) %>%
  left_join(df_landuse_sel,   by = c("geo", "year_idx")) %>%
  clean_names()

n0 <- nrow(df_population_sel)
n_all <- nrow(prediction_covariates_all_2022_2080)

message(sprintf("Rows after join: %s (base had %s)", n_all, n0)) # 

# Extra double check step
miss_mob   <- df_population_sel %>% anti_join(df_mobility_sel,  by = c("geo","year_idx")) %>% nrow()
miss_clim_counter  <- df_population_sel %>% anti_join(df_climate_pred_2022_2080_sel,   by = c("geo","year_idx")) %>% nrow()
miss_land  <- df_population_sel %>% anti_join(df_landuse_sel,   by = c("geo","year_idx")) %>% nrow()

message(sprintf("Unmatched rows (of population) vs mobility: %d | climate: %d | land-use: %d ",
                miss_mob, miss_clim_counter, miss_land))



# Check the number of rows
rows_before <- nrow(prediction_covariates_all_2022_2080)
rows_before
unique(prediction_covariates_all_2022_2080$year_idx)

covariates_all_pred <- prediction_covariates_all_2022_2080 %>% dplyr::select(-dplyr::any_of(cols_to_drop))
colnames(covariates_all_pred)

# Desired training columns (and order)
desired_cols <- names(final_covariates)
desired_cols

# Clean desired_cols first
desired_cols_pred <- desired_cols[desired_cols %in% names(covariates_all_pred)] # # only those desired cols which are in pred data
# The proximity is gone now as it is not in pred data

covariates_all_pred_final <- covariates_all_pred %>%
  dplyr::select(all_of(desired_cols_pred))

colnames(covariates_all_pred_final)

id_cols <- c("geo","year_idx","year","location_na")
cols_to_check <- setdiff(names(covariates_all_pred_final), id_cols)

# Row-wise completeness + quick summary
good_idx <- complete.cases(covariates_all_pred_final[cols_to_check])

n_total <- nrow(covariates_all_pred_final)
n_good  <- sum(good_idx)
n_bad   <- n_total - n_good

if (n_bad == 0) {
  message("✅ No missing or NA values in selected columns (", length(cols_to_check), " columns; ", n_total, " rows).")
} else {
  message("⚠️ Found ", n_bad, " / ", n_total, " rows with at least one NA in selected columns.")
  message("First few problematic row indices: ", paste(head(which(!good_idx), 10), collapse = ", "))
}



# Now final counter climates
final_counter_covariates <- covariates_all_pred_final


# we will use the same range and bin number for binning as in training or counterfactual data binning 
# But instead of going through each binning_cols like temerature or precipitation etc. we just use a loop to go through all the columns at once
# CAUTION - No proximity binning for future covariates
binning_cols <- names(final_covariates)[
  grepl("temperature|precipitation|humidity|population", names(final_covariates))
  # No proximity binning for future covareates
]
binning_cols

final_covariates_binned         <- final_covariates
final_counter_covariates_binned <- final_counter_covariates


##### Merging bins after inspecting log odds of fitted random effects
get_merge_spec <- function(var) {
  if (var == "population_density") {
    return(list("2:4", "6:13"))
  } else if (var == "temperature_winter_mean") {
    return(list("1:6","13:15"))
  } else if (var == "temperature_summer_mean") {
    return(list("1:3","14:15"))
  } else if (var == "relative_humidity_winter_mean") {
    return(list("1:4", "13:15"))
  } else if (var == "relative_humidity_spring_mean") {
    return(list("1:4", "10:11"))
  } else if (var == "y_proximity_exp300") {
    return(list("8:9", "13:15"))
  } else {
    return(NULL)
  }
}



get_n_bins <- function(var) {
  if (var == "population_density") 15 else 15     ## same initial bins for population too
}

# Decide if merging actually happened
was_merged <- function(out) {
  if (!("bin_idx_merged" %in% names(out$train_bin))) return(FALSE)
  !isTRUE(all(out$train_bin$bin_idx_merged == out$train_bin$bin_idx))
}




#------------------------------------------------------------
# CORE: attach a SINGLE canonical bin index + keys
#    - If merging happened: use merged indices and merged keys
#    - Else: use unmerged indices and unmerged keys
# ------------------------------------------------------------
attach_bins <- function(train_df, pred_df, out, col) {
  merged <- was_merged(out)
  
  # Column names to write/join
  idx_col  <- paste0("bin_idx_",   col)
  rng_col  <- paste0("bin_range_", col) # readable range label
  mean_col <- paste0("bin_mean_",  col) # bin mean of the raw variable
  
  if (merged) {
    # Use MERGED indices
    train_df <- dplyr::mutate(train_df, !!idx_col := out$train_bin$bin_idx_merged) # joining train bins if merging is done
    pred_df  <- dplyr::mutate(pred_df,  !!idx_col := out$pred_bin$bin_idx_merged) # joining pred bins if merging is done
    
    keys_tr <- out$keys_train_merged %>%
      dplyr::rename(
        !!rng_col  := dplyr::starts_with("bin_range_"),
        !!mean_col := dplyr::starts_with("bin_mean_")
      )
    
    keys_pr <- out$keys_pred_merged %>%
      dplyr::rename(
        !!rng_col  := dplyr::starts_with("bin_range_"),
        !!mean_col := dplyr::starts_with("bin_mean_")
      )
    
    msg <- "MERGED (since merging has been done)"
  } else {
    # Use UNMERGED indices
    train_df <- dplyr::mutate(train_df, !!idx_col := out$train_bin$bin_idx) # joining train bins if no merging is done
    pred_df  <- dplyr::mutate(pred_df,  !!idx_col := out$pred_bin$bin_idx) #  joining pred bins if no merging is done
    
    # Use UNMERGED key tables
    keys_tr <- out$keys_train %>%
      dplyr::rename(
        !!rng_col  := dplyr::starts_with("bin_range_"),
        !!mean_col := dplyr::starts_with("bin_mean_")
      )
    
    keys_pr <- out$keys_pred %>%
      dplyr::rename(
        !!rng_col  := dplyr::starts_with("bin_range_"),
        !!mean_col := dplyr::starts_with("bin_mean_")
      )
    
    msg <- "UNMERGED (since no merging has been done)"
  }
  
  train_df <- train_df %>%
    dplyr::left_join(
      keys_tr %>% dplyr::select(bin_idx, !!rng_col, !!mean_col),
      by = dplyr::join_by(!!rlang::sym(idx_col) == bin_idx)
    )
  
  pred_df  <- pred_df %>%
    dplyr::left_join(
      keys_pr %>% dplyr::select(bin_idx, !!rng_col, !!mean_col),
      by = dplyr::join_by(!!rlang::sym(idx_col) == bin_idx)
    )
  
  message(sprintf("[%s] Using %s bins; wrote %s and joined %s/%s.",
                  col, msg, idx_col, rng_col, mean_col))
  
  list(train_df = train_df, pred_df = pred_df)
}


# Extracting the bin details columns from res and writing back the three new columns from res
write_back_by_keys <- function(orig_df, res_df, col, keys = c("geo","year_idx")) {
  idx_col  <- paste0("bin_idx_",   col)
  rng_col  <- paste0("bin_range_", col)
  mean_col <- paste0("bin_mean_",  col)
  
  out_cols <- c(keys, idx_col, rng_col, mean_col)
  
  add_df <- res_df %>% dplyr::select(dplyr::all_of(out_cols))
  
  orig_df %>%
    dplyr::select(-dplyr::any_of(c(idx_col, rng_col, mean_col))) %>%
    dplyr::left_join(add_df, by = keys)
}



# LOOP OVER ALL COLUMNS IN `binning_cols`
for (col in binning_cols) {
  message("------------------------------------------------------------")
  message(sprintf("Processing binning for variable: %s", col))
  message("------------------------------------------------------------")
  
  merge_spec <- get_merge_spec(col)
  n_bins     <- get_n_bins(col)
  
  # FUNCTION FOR NO MERGING AND MERGING BOTH
  out <- bin_and_merge_one(
    train_df   = final_covariates,
    pred_df    = final_counter_covariates,
    col        = col,
    n_bins     = n_bins, 
    right = TRUE, #  right-closed intervals: (a, b] for bins
    include_lowest = TRUE,
    merge_spec = merge_spec, # if NULL → unmerged, otherwise merging as per merge_spec
    key_fmt    = "[%.1f,%.1f]"
  )
  
  
  # A USING TRAINING BINS TO COMPUTE PROXIMITY
  merged_bin_scenario <- out # also contain information about original bins too
  
  # TRUE if merging is done (FALSE if not)
  was_merged(merged_bin_scenario)
  
  # Attach ONE canonical bin index + readable keys
  res <- attach_bins(
    train_df = final_covariates_binned,
    pred_df  = final_counter_covariates_binned,
    out      = merged_bin_scenario,
    col      = col
  )
  
  
  final_covariates_binned         <- write_back_by_keys(final_covariates_binned,         res$train_df, col)
  final_counter_covariates_binned <- write_back_by_keys(final_counter_covariates_binned, res$pred_df,  col)
  
  message(sprintf("Finished binning for %s", col))
}



# After loop:
# Training data is already saved in load_and_process_data

# Paths (already set)
out_dir <- ".../Data/climate/Final_covariates_data/Future_data/Binned_final_data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# Renaming the data to make things clearer

# # first for picontrol climate with ssp126 socioeconomic

# final_GFDL_ESM4_picontrol_climate_with_ssp126_socioeconomic_binned_data <- final_counter_covariates_binned
# 
# ## chage the variable and file name accordingly
# save(final_GFDL_ESM4_picontrol_climate_with_ssp126_socioeconomic_binned_data,
#      file = file.path(out_dir, "final_GFDL-ESM4_picontrol_climate_with_ssp126_socioeconomic_binned_data.RData"))


# Next for ssp126 climate with ssp126 socioeconomic
final_GFDL_ESM4_ssp126_climate_with_ssp126_socioeconomic_binned_data <- final_counter_covariates_binned

# Chage the variable and file name accordingly
save(final_GFDL_ESM4_ssp126_climate_with_ssp126_socioeconomic_binned_data,
     file = file.path(out_dir, "final_GFDL-ESM4_ssp126_climate_with_ssp126_socioeconomic_binned_data.RData"))

############### END FOR SSP126 future data preparation ##############################################



######################################################################################################
## Repeat the same process as above for SSP370/picontrol or SSP585/picontrol climate data with SSP3 or SSP5 socioeconomic variables
## save them seperately and use it for future predictions
######################################################################################################



