# ================== Historical and counterfactual data loading and preprocessing ==================

setwd(".../Data")

# Source project functions
source(".../functions_general.R")
source(".../functions_plotting.R")
source(".../bin_covariates.R")


# Load input data
out <- attach("data_INLA.RData")
map_inla <- map_inla
detach(paste0("file:", "data_INLA.RData"), character.only = TRUE)


# Load data -------
out <- load("mosq_dist_until_2023.RData")
mosq_dist <- mosq_dist %>% filter(date_start>=2010) %>%
  mutate(year_idx = match(date_start, sort(unique(date_start))))
tb_inla <- mosq_dist %>% filter(geo %in% map_inla$geo)
rm(mosq_dist)

data_spde_mosq <- tb_inla
map_spde <- map_inla %>% filter(geo %in% data_spde_mosq$geo)


# Define data locations
covariates_folder_data <- "climate/Final_covariates_data/Current_data"

out <- load(paste0(covariates_folder_data, "/Final_hist_population_data_yearly.RData"))

# Load and prepare covariates
df_population <- population



## mobility 
out <- load(paste0(covariates_folder_data, "/Final_hist_mobility_data_yearly.RData"))
df_mobility <- df_mobility

# scale the total_flux in column as we are going to use this as as fixed effect and it should be comparable to land use covariates which are between 0 and 1
df_mobility <- df_mobility %>%
  mutate(total_flux_in = {
    m <- suppressWarnings(max(total_flux_in, na.rm = TRUE))
    if (is.finite(m) && m > 0) total_flux_in / m else total_flux_in
  })



out <- load(paste0(covariates_folder_data, "/Final_20CRV3-ERA5_factual_climate_data_yearly.RData"))
df_climate <- climate_summaries_region

unique(df_climate$year_idx)


### Since climate data is from year 2001 to 2023 (year_idx is from 1 to 23), so taking only the data from year 2010 to 2023 
df_climate <- df_climate %>%
  filter(year_idx > 9) %>%
  mutate(year_idx = year_idx - 9)

unique(df_climate$year_idx)

out <- load(paste0(covariates_folder_data, "/Final_hist_land_use_scaled_data_yearly.RData"))
df_landuse <- land_use


# Rationale: Drop `croplands` because LUH cropland categories are noisy/uncertain and distributed across many crops type categories.
stopifnot("croplands" %in% names(df_landuse))

df_landuse <- df_landuse %>%
  dplyr::select(-croplands)


## proximity covariatw
out <- load(paste0(covariates_folder_data, "/Final_hist_proximity_data_yearly.RData"))
df_proximity <- df_proximity

# dropping the last year’s proximity rows
df_proximity <- df_proximity %>%
  group_by(geo) %>%
  filter(year_idx < max(year_idx)) %>%
  ungroup()

unique(df_proximity$year_idx)



## joining all the covariates 
coerce_keys <- function(x) {
  x %>%
    mutate(
      geo = as.character(geo),
      year_idx = as.integer(year_idx)
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
  dplyr::select(-LocationNa) %>%
  dplyr::select(geo, year_idx, everything())

df_landuse_sel <- df_landuse %>%
  dplyr::select(-LocationNa, -year) %>%
  dplyr::select(geo, year_idx, everything())

df_proximity_sel <- df_proximity %>%
  dplyr::select(geo, year_idx, everything())

stopifnot(nrow(df_population_sel) == nrow(distinct(df_population_sel, geo, year_idx)))
stopifnot(nrow(df_mobility_sel)   == nrow(distinct(df_mobility_sel,   geo, year_idx)))
stopifnot(nrow(df_climate_sel)    == nrow(distinct(df_climate_sel,    geo, year_idx)))
stopifnot(nrow(df_landuse_sel)    == nrow(distinct(df_landuse_sel,    geo, year_idx)))
stopifnot(nrow(df_proximity_sel)  == nrow(distinct(df_proximity_sel,  geo, year_idx)))

# Join covariate tables
covariates_all_2010_2023 <- df_population_sel %>%
  left_join(df_mobility_sel,  by = c("geo", "year_idx")) %>%
  left_join(df_climate_sel,   by = c("geo", "year_idx")) %>%
  left_join(df_landuse_sel,   by = c("geo", "year_idx")) %>%
  left_join(df_proximity_sel, by = c("geo", "year_idx")) %>%
  clean_names()

n0 <- nrow(df_population_sel)
n_all <- nrow(covariates_all_2010_2023)

message(sprintf("Rows after join: %s (base had %s)", n_all, n0))

miss_mob   <- df_population_sel %>% anti_join(df_mobility_sel,  by = c("geo","year_idx")) %>% nrow()
miss_clim  <- df_population_sel %>% anti_join(df_climate_sel,   by = c("geo","year_idx")) %>% nrow()
miss_land  <- df_population_sel %>% anti_join(df_landuse_sel,   by = c("geo","year_idx")) %>% nrow()
miss_prox  <- df_population_sel %>% anti_join(df_proximity_sel, by = c("geo","year_idx")) %>% nrow()

message(sprintf("Unmatched rows vs mobility: %d | climate: %d | land-use: %d | proximity: %d",
                miss_mob, miss_clim, miss_land, miss_prox))

rows_before <- nrow(covariates_all_2010_2023)


### Drop first year covariate data as it is not used in training (use either condition; both included for safety)
covariates_all <- covariates_all_2010_2023 %>%
  filter(!(year_idx == 1 | year == 2010))

rows_after <- nrow(covariates_all)
message(sprintf("Dropped %d rows (first year). Remaining: %d",
                rows_before - rows_after, rows_after))

covariates_all <- covariates_all %>%
  mutate(year_idx = year_idx - 1L)

### Drop first year from mosuito data too as it is not used in training (use either condition; both included for safety)
data_spde1 <- data_spde_mosq %>%
  filter(!(date_start == 2010))

data_spde1 <- data_spde1 %>%
  mutate(year_idx = year_idx - 1L)

map_spde$region_id <- 1:nrow(map_spde)
map_spde$spatial_idx <- map_spde$region_id

data_spde1 <- left_join(data_spde1 %>% dplyr::select(-matches("region_id")),
                        as_tibble(map_spde) %>% dplyr::select(geo, region_id)) %>%
  dplyr::arrange(year_idx, region_id) %>% relocate(geo, region_id)

g_spde <- data_spde1$geo
g_cov  <- covariates_all$geo


#  Same length?
cat("Lengths — data_spde:", length(g_spde), "| covariates_all:", length(g_cov), "\n")

#  Same unique sets?
miss_in_cov  <- setdiff(g_spde, g_cov)
miss_in_spde <- setdiff(g_cov,  g_spde)
cat("Missing in covariates_all:", length(miss_in_cov), "\n")
cat("Missing in data_spde:",      length(miss_in_spde), "\n")
if (length(miss_in_cov))  cat("Examples (in spde not in cov):",  head(miss_in_cov),  "\n")
if (length(miss_in_spde)) cat("Examples (in cov not in spde):", head(miss_in_spde), "\n")




data_spde1 <- data_spde1 %>% mutate(data_idx=1:n()) %>%
  relocate(data_idx, region_id, year_idx)

col_combined_cov <- colnames(covariates_all)

## Removing some extra columns like population (not density), area, mobility_idx, net_flux_out(not total flux in) and climate annual mean (used for plotting, only seasonal aggregates used as covariates )
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

covariates_all <- covariates_all %>% dplyr::select(-dplyr::any_of(cols_to_drop))

info_cols <- c("geo","year_idx","year","location_na")
stopifnot(all(info_cols %in% names(covariates_all)))

info_df <- covariates_all %>% dplyr::select(all_of(info_cols))
cov_df  <- covariates_all %>% dplyr::select(-all_of(info_cols)) %>% dplyr::select(where(is.numeric))

rows_before <- nrow(covariates_all)

# Remove incomplete rows
good_idx    <- complete.cases(cov_df)
if (!any(good_idx)) stop("All rows have NA in covariates; nothing to analyze.")
dropped_rows <- rows_before - sum(good_idx)


## dropping from both info_df as well as cov_df 
covariates_all_clean <- dplyr::bind_cols(
  info_df[good_idx, , drop = FALSE],
  cov_df [good_idx, , drop = FALSE]
)

message(sprintf("Rows before NA filtering: %d | After: %d | Dropped: %d",
                rows_before, nrow(covariates_all_clean), dropped_rows))

# check!! - same number of rows and same ordering before using good_idx
stopifnot(nrow(data_spde1) == nrow(covariates_all))
stopifnot(identical(data_spde1$geo, covariates_all$geo))
stopifnot(identical(data_spde1$year_idx,  covariates_all$year_idx))


keep_keys <- covariates_all_clean %>% dplyr::select(geo, year_idx) %>% distinct()

data_spde <- data_spde1 %>% semi_join(keep_keys, by = c("geo","year_idx"))


# checks of new data_spde with clean covairates
stopifnot(nrow(data_spde) == nrow(covariates_all_clean))
stopifnot(setequal(paste(data_spde$geo, data_spde$year_idx),
                   paste(covariates_all_clean$geo, covariates_all_clean$year_idx)))



###########CORRELATION##############

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

info_df_cleaned_agg <- covariates_by_region %>% dplyr::select(all_of(info_cols2))
cov_df_cleaned_agg  <- covariates_by_region %>% dplyr::select(-all_of(info_cols2)) %>% dplyr::select(where(is.numeric))


# Assess covariate correlation
C_full <- cor(cov_df_cleaned_agg, use = "pairwise.complete.obs", method = "pearson")

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
png(".../Plots/A_final_plots/Correlation_plot/corr_heatmap_before_pruning_(cutoff).png", width = 1500, height = 1300, res = 180)
grid::grid.newpage(); grid::grid.draw(hm_full$gtable)
dev.off()


## saving helper (TIFF + SVG) for correlation plot
clean_cov_names <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

C_full   <- cor(cov_df_cleaned_agg, use = "pairwise.complete.obs", method = "pearson")
ord_full <- hclust(as.dist(1 - abs(C_full)))$order

mat <- abs(C_full)[ord_full, ord_full]

labs_clean <- clean_cov_names(rownames(mat))
rownames(mat) <- labs_clean
colnames(mat) <- labs_clean

hm_full <- pheatmap::pheatmap(
  mat,
  clustering_distance_rows = as.dist(1 - abs(C_full)),
  clustering_distance_cols = as.dist(1 - abs(C_full)),
  clustering_method = "complete",
  
  treeheight_col = 0,
  treeheight_row = 0,
  
  legend = TRUE,
  main = "",
  display_numbers = FALSE,
  silent = TRUE,
  
  fontsize      = 7,
  fontsize_row  = 7,
  fontsize_col  = 7,
  angle_col     = 90,
  border_color  = NA
)

out_dir_base <- ".../Plots/A_final_plots/Correlation_plot"

ext_tiff    <- "tiff"
ext_svg     <- "svg"
width_cm    <- 12.1
height_cm   <- 12.1
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

dir.create(file.path(out_dir_base, "Tiff"), recursive = TRUE, showWarnings = FALSE)
file_tiff <- file.path(out_dir_base, "Tiff", "corr_heatmap_before_pruning_(cutoff).tiff")

ragg::agg_tiff(
  filename = file_tiff,
  width = width_cm, height = height_cm, units = "cm",
  res = dpi,
  compression = compression,
  background = bg
)
grid::grid.newpage(); grid::grid.draw(hm_full$gtable)
dev.off()

dir.create(file.path(out_dir_base, "svg"), recursive = TRUE, showWarnings = FALSE)
file_svg <- file.path(out_dir_base, "svg", "corr_heatmap_before_pruning_(cutoff).svg")

w_cm <- 12.1
h_cm <- 12.1

svglite::svglite(
  file   = file_svg,
  width  = w_cm / 2.54,
  height = h_cm / 2.54,
  bg     = "transparent"
)
grid::grid.newpage()
grid::grid.draw(hm_full$gtable)
dev.off()


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
    
    # --- protect "keep" variables ---
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
pruned <- prune_cor_matrix(C_full, cutoff = 0.7, keep = protected_keep)   ## standard cutoff

cat("\n---- Pruning summary (cutoff = 0.70, with protection) ----\n")
cat(sprintf("Initial covariates: %d\n", ncol(C_full)))
cat(sprintf("Dropped: %d\n", length(pruned$dropped)))
if (length(pruned$dropped)) cat("Dropped vars:\n  ", paste(pruned$dropped, collapse = ", "), "\n")
cat(sprintf("Kept: %d\n", length(pruned$kept)))
cat("Kept vars:\n  ", paste(pruned$kept, collapse = ", "), "\n\n")


## heatmap with the protected set
if (ncol(pruned$C_kept) >= 2) {
  ord_kept <- hclust(as.dist(1 - abs(pruned$C_kept)))$order
  hm_kept <- pheatmap::pheatmap(
    abs(pruned$C_kept)[ord_kept, ord_kept],
    clustering_distance_rows = as.dist(1 - abs(pruned$C_kept)),
    clustering_distance_cols = as.dist(1 - abs(pruned$C_kept)),
    clustering_method = "complete",
    legend = TRUE,
    main = "Correlation heatmap (|r|) after pruning (proximity protected)",
    display_numbers = FALSE,
    silent = TRUE
  )
  png(".../Plots/A_final_plots/Correlation_plot/corr_heatmap_after_pruning_proximity_protected.png",
      width = 1500, height = 1300, res = 180)
  grid::grid.newpage(); grid::grid.draw(hm_kept$gtable)
  dev.off()
}



info_df_cleaned <- covariates_all_clean %>% dplyr::select(all_of(info_cols))
cov_df_cleaned  <- covariates_all_clean %>% dplyr::select(-all_of(info_cols))

prev_col <- ("presence_previous_year")     # we do not need this column although we add it for consistency
cols_to_keep <- union(pruned$kept, prev_col)

cov_df_pruned <- cov_df_cleaned %>%
  dplyr::select(dplyr::any_of(cols_to_keep))

covariates_all_pruned <- dplyr::bind_cols(info_df_cleaned, cov_df_pruned)

final_covariates <- covariates_all_pruned

colnames(final_covariates)



####################### NOW LOADING THE PREDICTION DATA also and preparing it ###############
out <- load(paste0(covariates_folder_data, "/Final_20CRV3-ERA5_counterfactual_climate_data_yearly.RData"))
df_climate_counter_2010_2021 <- climate_summaries_region

unique(df_climate_counter_2010_2021$year_idx)

### Since climate data in prediction is also  from year 2001 to 2021 (original year_idx is from 1 to 21), so taking only the data from year 2010 to 2021
df_climate_counter_2010_2021 <- df_climate_counter_2010_2021 %>%
  filter(year_idx > 9) %>%
  mutate(year_idx = year_idx - 9)

unique(df_climate_counter_2010_2021$year_idx)

# Since the climate table runs only to year_idx == 12, just trim the others to that same range.
max_year <- max(df_climate_counter_2010_2021$year_idx)


# only climates covariates are going to change in counterfactual scenario , we join every other socioeconomic covariates as it is in factual scenarios
df_population_pred <- dplyr::filter(df_population, year_idx <= max_year)
unique(df_population_pred$year_idx)
df_mobility_pred   <- dplyr::filter(df_mobility,   year_idx <= max_year)
df_landuse_pred    <- dplyr::filter(df_landuse,    year_idx <= max_year)
df_proximity_pred  <- dplyr::filter(df_proximity,  year_idx <= max_year)
unique(df_proximity_pred$year_idx)

coerce_keys <- function(x) {
  x %>%
    mutate(
      geo = as.character(geo),
      year_idx = as.integer(year_idx)
    )
}

df_population_pred  <- coerce_keys(df_population_pred)  %>% distinct(geo, year_idx, .keep_all = TRUE)
df_mobility_pred    <- coerce_keys(df_mobility_pred)    %>% distinct(geo, year_idx, .keep_all = TRUE)
df_climate_counter_mid     <- coerce_keys(df_climate_counter_2010_2021)     %>% distinct(geo, year_idx, .keep_all = TRUE)
df_landuse_pred     <- coerce_keys(df_landuse_pred)     %>% distinct(geo, year_idx, .keep_all = TRUE)
df_proximity_pred     <- coerce_keys(df_proximity_pred)     %>% distinct(geo, year_idx, .keep_all = TRUE)

df_population_sel <- df_population_pred %>%
  dplyr::select(geo, year_idx, year, LocationNa, population, area, population_density)

df_mobility_sel <- df_mobility_pred %>%
  dplyr::select(geo, year_idx, mobility_idx, net_flux_out, total_flux_in)

df_climate_counter_2010_2021_sel <- df_climate_counter_mid  %>%
  dplyr::select(-LocationNa) %>%
  dplyr::select(geo, year_idx, everything())

df_landuse_sel <- df_landuse_pred %>%
  dplyr::select(-LocationNa, -year) %>%
  dplyr::select(geo, year_idx, everything())

df_proximity_sel <- df_proximity_pred %>%
  dplyr::select(geo, year_idx, everything())

stopifnot(nrow(df_population_sel) == nrow(distinct(df_population_sel, geo, year_idx)))
stopifnot(nrow(df_mobility_sel)   == nrow(distinct(df_mobility_sel,   geo, year_idx)))
stopifnot(nrow(df_climate_counter_2010_2021_sel)    == nrow(distinct(df_climate_counter_2010_2021_sel,    geo, year_idx)))
stopifnot(nrow(df_landuse_sel)    == nrow(distinct(df_landuse_sel,    geo, year_idx)))
stopifnot(nrow(df_proximity_sel)  == nrow(distinct(df_proximity_sel,  geo, year_idx)))

counter_covariates_all_2010_2021 <- df_population_sel %>%
  left_join(df_mobility_sel,  by = c("geo", "year_idx")) %>%
  left_join(df_climate_counter_2010_2021_sel,   by = c("geo", "year_idx")) %>%
  left_join(df_landuse_sel,   by = c("geo", "year_idx")) %>%
  left_join(df_proximity_sel, by = c("geo", "year_idx")) %>%
  clean_names()


#  Quick diagnostics (i.e. number of pop rows are equal to all covariates joined)
n0 <- nrow(df_population_sel)
n_all <- nrow(counter_covariates_all_2010_2021)

message(sprintf("Rows after join: %s (base had %s)", n_all, n0))



## extra double check step
miss_mob   <- df_population_sel %>% anti_join(df_mobility_sel,  by = c("geo","year_idx")) %>% nrow()
miss_clim_counter  <- df_population_sel %>% anti_join(df_climate_counter_2010_2021_sel,   by = c("geo","year_idx")) %>% nrow()
miss_land  <- df_population_sel %>% anti_join(df_landuse_sel,   by = c("geo","year_idx")) %>% nrow()
miss_prox  <- df_population_sel %>% anti_join(df_proximity_sel, by = c("geo","year_idx")) %>% nrow()

message(sprintf("Unmatched rows vs mobility: %d | climate: %d | land-use: %d | proximity: %d",
                miss_mob, miss_clim_counter, miss_land, miss_prox))
## Now data loading of prediction is done


## check the number of rows
rows_before <- nrow(counter_covariates_all_2010_2021)
covariates_all_counter  <- counter_covariates_all_2010_2021 %>%
  filter(!(year_idx == 1 | year == 2010))

rows_after <- nrow(covariates_all_counter)
message(sprintf("Dropped %d rows (first year). Remaining: %d",
                rows_before - rows_after, rows_after))

# Shift global year_idx of prediction covariates down by 1 (so year_idx starts from 1)
covariates_all_counter  <- covariates_all_counter %>%
  mutate(year_idx = year_idx - 1L)


# Remove all the uninformative columns as in training (stored in cols_to_drop variable (check in training drop) and keep everything else
covariates_all_counter <- covariates_all_counter %>% dplyr::select(-dplyr::any_of(cols_to_drop))
colnames(covariates_all_counter)


#keeping only those covariates which are used in training after removing the correlation
desired_cols <- names(final_covariates)

extra_in_counter  <- setdiff(names(covariates_all_counter), desired_cols)
missing_in_counter <- setdiff(desired_cols, names(covariates_all_counter))

message("Extra columns in counter (will be dropped): ",
        paste(extra_in_counter, collapse = ", "))
message("Missing columns in counter (will be created as NA): ",
        paste(missing_in_counter, collapse = ", "))

## Error when something is present in training but not in prediction
if (length(missing_in_counter)) {
  stop("Counter dataset is missing required columns: ",
       paste(missing_in_counter, collapse = ", "))
}


# Build a type map if some column is in training and not in prediction, so we can create NAs of the right type for that column in prediction dataset
type_map <- vctrs::vec_ptype_abbr(covariates_all_pruned)
na_of_type <- function(t) {
  switch(t,
         int = NA_integer_,
         dbl = NA_real_,
         chr = NA_character_,
         lgl = NA,
         NA_real_)
}

covariates_all_counter_aligned <- covariates_all_counter %>%
  dplyr::select(any_of(intersect(desired_cols, names(.)))) %>%
  {
    if (length(missing_in_counter)) {
      for (col in missing_in_counter) {
        t <- type_map[[col]]
        .[[col]] <- na_of_type(t)
      }
    }
    .
  } %>%
  
  dplyr::select(all_of(desired_cols))

colnames(covariates_all_counter_aligned)

#strict checks
stopifnot(identical(names(covariates_all_counter_aligned), desired_cols))




####### check for any NA values of any covariates in prediction data
id_cols <- c("geo","year_idx","year","location_na")
cols_to_check <- setdiff(names(covariates_all_counter_aligned), id_cols)

good_idx <- complete.cases(covariates_all_counter_aligned[cols_to_check])

n_total <- nrow(covariates_all_counter_aligned)
n_good  <- sum(good_idx)
n_bad   <- n_total - n_good

if (n_bad == 0) {
  message("✅ No missing or NA values in selected columns (", length(cols_to_check), " columns; ", n_total, " rows).")
} else {
  message("⚠️ Found ", n_bad, " / ", n_total, " rows with at least one NA in selected columns.")
  message("First few problematic row indices: ", paste(head(which(!good_idx), 10), collapse = ", "))
}



## BINNING process for covariates
final_counter_covariates <- covariates_all_counter_aligned

binning_cols <- names(final_covariates)[
  grepl("temperature|precipitation|humidity|proximity|population|proximity", names(final_covariates))
]

binning_cols

final_covariates_binned <- final_covariates
final_counter_covariates_binned <- final_counter_covariates

# Choose the column one by one we want to bin by training rules and repeat this for all climate and proximity cols
#col <- binning_cols[1]

#col <- binning_cols[2]

#col <- binning_cols[3]

#col <- binning_cols[4]

#col <- binning_cols[5]

#col <- binning_cols[6]

#col <- binning_cols[7]

col <- binning_cols[8]

col



#Binning with NO merging first, inspect keys
out0 <- bin_and_merge_one(
  train_df   = final_covariates,
  pred_df    = final_counter_covariates,
  col        = col,
  n_bins     = 15,       # number of equal-width bins for the initial cut() 
  key_fmt    = "[%.1f,%.1f]",
  right = TRUE,
  include_lowest = TRUE,
  merge_spec = NULL        # if NULL then unmerged
  
)

#training and pred keys (pre-merge)
print(out0$keys_train)
print(out0$keys_pred)


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
  if (var == "population_density") 15 else 15    ## same initial bins for population too
}



merge_spec <- get_merge_spec(col)     
n_bins     <- get_n_bins(col)

out <- bin_and_merge_one(
  train_df   = final_covariates,
  pred_df    =final_counter_covariates,
  col        = col,
  n_bins     = n_bins,      # number of equal-width bins for the initial cut()
  right = TRUE,              #  right-closed intervals: (a, b] for bins
  include_lowest = TRUE,
  merge_spec = merge_spec,        # if NULL → unmerged, otherwise merging as per merge_spec
  key_fmt    = "[%.1f,%.1f]"
)

# Inspect merged keys (training + prediction)

print(out$keys_train)
print(out$keys_train_merged)
print(out$keys_pred_merged)


##### IMP - If the covariate is proximity covariate, save its binning spec for prediction proximity calculation in recurring manner for each year prediction
if (col == "y_proximity_exp300") {
  
  make_proximity_binning_spec <- function(out,
                                          col,
                                          merge_spec     = NULL,
                                          right          = TRUE,
                                          include_lowest = TRUE,
                                          key_fmt        = "[%.1f,%.1f]") {
    list(
      col              = col,
      
      breaks           = out$breaks,
      empty_bins       = out$empty_bins,
      remap_to_nearest = out$remap_to_nearest,
      compress_map     = out$compress_map,
      
      merge_spec       = merge_spec,
      
      keys_train        = out$keys_train,
      keys_train_merged = out$keys_train_merged,
      
      right            = right,
      include_lowest   = include_lowest,
      key_fmt          = key_fmt
    )
  }
  
  prox300_spec <- make_proximity_binning_spec(
    out           = out,
    col           = col,
    merge_spec    = merge_spec,
    right         = TRUE,
    include_lowest = TRUE,
    key_fmt       = "[%.1f,%.1f]"
  )
  
  out_dir <- ".../Output/A-Final_output/Proximity_binning_info_for_pred"
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  save(
    prox300_spec,
    file = file.path(out_dir, "prox300_binning_spec.RData")
  )
  
  message("Saved proximity binning spec to: ",
          file.path(out_dir, "prox300_binning_spec.RData"))
}




merged_bin_scenario <- out     ## also contain information about original bins
merged_bin_scenario



was_merged <- function(out) {
  if (!("bin_idx_merged" %in% names(out$train_bin))) return(FALSE)
  !isTRUE(all(out$train_bin$bin_idx_merged == out$train_bin$bin_idx))
}

was_merged(out)       ## TRUE if merging is done (FALSE if not)



#------------------------------------------------------------
#   -  attach a SINGLE canonical bin index + keys
#    - If merging happened: use merged indices and merged keys
#    - Else: use unmerged indices and unmerged keys
# ------------------------------------------------------------
attach_bins <- function(train_df, pred_df, out, col) {
  merged <- was_merged(out)
  
  idx_col  <- paste0("bin_idx_",   col)
  rng_col  <- paste0("bin_range_", col)
  mean_col <- paste0("bin_mean_",  col)
  
  if (merged) {
    train_df <- dplyr::mutate(train_df, !!idx_col := out$train_bin$bin_idx_merged)
    pred_df  <- dplyr::mutate(pred_df,  !!idx_col := out$pred_bin$bin_idx_merged)
    
    keys_tr <- out$keys_train_merged %>%
      dplyr::rename(!!rng_col := dplyr::starts_with("bin_range_"),
                    !!mean_col := dplyr::starts_with("bin_mean_"))
    keys_pr <- out$keys_pred_merged %>%
      dplyr::rename(!!rng_col := dplyr::starts_with("bin_range_"),
                    !!mean_col := dplyr::starts_with("bin_mean_"))
    
    msg <- "MERGED (since merging has been done)"
  } else {
    train_df <- dplyr::mutate(train_df, !!idx_col := out$train_bin$bin_idx)
    pred_df  <- dplyr::mutate(pred_df,  !!idx_col := out$pred_bin$bin_idx)
    
    keys_tr <- out$keys_train %>%
      dplyr::rename(!!rng_col := dplyr::starts_with("bin_range_"),
                    !!mean_col := dplyr::starts_with("bin_mean_"))
    keys_pr <- out$keys_pred %>%
      dplyr::rename(!!rng_col := dplyr::starts_with("bin_range_"),
                    !!mean_col := dplyr::starts_with("bin_mean_"))
    
    msg <- "UNMERGED (since no merging has been done)"
  }
  
  train_df <- train_df %>%
    dplyr::left_join(keys_tr %>% dplyr::select(bin_idx, !!rng_col, !!mean_col),
                     by = dplyr::join_by(!!rlang::sym(idx_col) == bin_idx))
  pred_df  <- pred_df %>%
    dplyr::left_join(keys_pr %>% dplyr::select(bin_idx, !!rng_col, !!mean_col),
                     by = dplyr::join_by(!!rlang::sym(idx_col) == bin_idx))
  
  message(sprintf("[%s] Using %s bins; wrote %s and joined %s/%s.",
                  col, msg, idx_col, rng_col, mean_col))
  list(train_df = train_df, pred_df = pred_df)
}


# --- Attach ONE canonical bin index + hreadable keys ---
res <- attach_bins(
  train_df = final_covariates_binned,
  pred_df  = final_counter_covariates_binned,
  out      = out,
  col      = col
)


unique(res$train_df[[paste0("bin_range_", col)]])
unique(res$train_df[[paste0("bin_range_", col)]])


# Extracting the bin details  columns from res and writng back the three new columns from res (created above)
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

# Overwrite the accumulating data frames with the updated versions
final_covariates_binned         <- write_back_by_keys(final_covariates_binned,         res$train_df, col)
final_counter_covariates_binned <- write_back_by_keys(final_counter_covariates_binned, res$pred_df,  col)


######!!!!!!!! SAVING both the final covariate and final counter covariate data!!!!!!!!!!!!!
out_dir <- ".../Data/climate/Final_covariates_data/Current_data/Binned_final_data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save(final_covariates_binned,
     file = file.path(out_dir, "final_binned_factual_data.RData"))

save(final_counter_covariates_binned,
     file = file.path(out_dir, "final_binned_counter_factual_data.RData"))

save(data_spde,
     file = file.path(out_dir, "final_mosq_data.RData"))

