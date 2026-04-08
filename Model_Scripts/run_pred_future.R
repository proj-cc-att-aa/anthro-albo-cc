# ================== Future prediction workflow ==================

# Setup
setwd(".../VectorNet_ 2_original")

# Load packages and functions
packages <- c("tidyverse", "stringr", "caret")
lapply(packages, library, character.only = TRUE)

library(dplyr)

source("functions_general.R")
source("functions_plotting.R")
source("functions_INLA.R")
library(INLA)
source("functions_N.R")
source("functions_matrix_operations.R")
source("functions_spatial_networks.R")
source("functions_risk_from_proximity.R")
source("functions_INLA_spde2.R")
source("bin_covariates.R")
source("COMPUTE_yearly_proximity_and_prox_bins.R") 



#climate_data_type_name = "pi_control"
proximity_exp = "050"

folder_data <- "Data"

binned_data_folder <- ".../Data/climate/Final_covariates_data/Current_data/Binned_final_data/"



# Loading the cleaned data one by one
# Mosquito cleaned data
out <- load(paste0(binned_data_folder, "final_mosq_data.RData"))
mosq_data <- data_spde

# Training cleaned data
out <- load(paste0(binned_data_folder, "final_binned_factual_data.RData"))
train_data <- final_covariates_binned


# Loading prediction cleaned data for future
binned_data_folder_future <- ".../Data/climate/Final_covariates_data/Future_data/Binned_final_data/"

out_pred <- load(paste0(binned_data_folder_future, "final_GFDL-ESM4_picontrol_climate_with_ssp370_socioeconomic_binned_data.RData"))

prediction_data <- final_GFDL_ESM4_picontrol_climate_with_ssp370_socioeconomic_binned_data



# Now join the mosquito data to training covariates data
# check: keys must be 1–1
stopifnot(all(c("geo","year_idx") %in% names(mosq_data)))
stopifnot(all(c("geo","year_idx") %in% names(train_data)))
stopifnot(nrow(train_data) == nrow(dplyr::distinct(train_data, geo, year_idx)))

# checks
stopifnot(
  nrow(train_data) ==
    nrow(distinct(train_data, geo, year_idx))
)


# checking that the two tables should cover the same (geo, year_idx) pairs
pairs_spde  <- mosq_data  %>% transmute(geo, year_idx)
pairs_cov   <- train_data %>% transmute(geo, year_idx)
stopifnot(anti_join(pairs_spde, pairs_cov, by=c("geo","year_idx")) %>% nrow() == 0)
stopifnot(anti_join(pairs_cov,  pairs_spde, by=c("geo","year_idx")) %>% nrow() == 0)


# Join: keep ALL columns from mosq_data; add year/location + all covariates
data_spde <- mosq_data %>%
  dplyr::left_join(
    train_data, # keep keys here
    by = c("geo","year_idx")
  )

# Training data is done ###########



# check the year column in future prediction data (and assign the year column)
if (!("year" %in% names(prediction_data))) {
  prediction_data <- prediction_data %>%
    mutate(year = 2022L + year_idx - 1L)%>% # year starting from 2022
    relocate(year, .after = year_idx)
}

pred_data <- prediction_data %>%
  mutate(
    data_idx   = dplyr::row_number(), # unique row id: 1..n
    region_id  = as.integer(factor(geo, levels = unique(geo))), # per-geo id, repeats across years
    date_start = as.integer(year), # mirror of `year`
    presence   = NA_real_, # prediction set → unknown presence
    entry_type = "predicted" # label rows as predicted (optional=)
  ) %>%
  dplyr::select( # put the key/meta cols first
    data_idx, region_id, year_idx, geo, date_start, entry_type, presence,
    everything()
  )

# Loading the map now
out <- attach("data_INLA.RData")
# For training
map_spde <- map_inla
# Prediction map
pred_map <- map_spde  ## same a straining map in our case



# Run INLA model
# Settings
proximity <- TRUE # use proximity covariate
#year1 <- FALSE # use year1 covariate
#presence_previous_year <- FALSE # use presence_previous_year covariate
model_id <- "vI"




# tuning parameters for the spde  (scale or precision, range or shape) but discarded in final version
# Used only for spatial effects assessment
Q_params <- c(log(1e-1),log(3))        

n_binom <-  1               # Bernoulli model: n_binom == 1, binomial model: n_binom > 1 

plot_larger <- TRUE
run_full_cv <- FALSE 

model_fits <- NULL
join_categories <- FALSE     ## FALSE for now




# Prepare SPDE inputs
# Spde mesh
out <- create_mesh(map_spde, mesh_precision="low")
spde_mesh <- out$spde_mesh
region_coords <- out$region_coords
n_spde <- out$n_spde

cat("Mesh, size: ", spde_mesh$n, "\n")
par(mfrow=c(1,1))
plot(spde_mesh)
plot_geometry(map_spde, add=TRUE)
points(region_coords, pch=20, col=2)

# Spde model 
spde_model <- inla.spde2.matern(spde_mesh, alpha=2)

# Tighten/identify the SPDE
# Use PC priors, mean-zero, and scaling: using bbox diagonal as a quick proxy
bb <- sf::st_bbox(sf::st_as_sf(map_spde))
study_diam <- sqrt((bb["xmax"]-bb["xmin"])^2 + (bb["ymax"]-bb["ymin"])^2)
r0 <- as.numeric(study_diam) / 3

## used in final version
spde_model <- inla.spde2.pcmatern(
  mesh = spde_mesh, alpha = 2,
  prior.range = c(r0, 0.3),
  prior.sigma = c(0.25, 0.01),
  constr = TRUE
)



# Precision matrix for model with manual tuning of the parameters, discarded in final version
Q_spde <- inla.spde2.precision(spde_model,theta=Q_params)         # Used only for spatial effects assessment




# INLA model function's and covariate's effects
response_column <- "presence"

data_spde__ <- data_spde # save original data 


parameters <- list(
  population_density            = TRUE,
  total_flux_in                 = TRUE,
  
  precipitation_winter_mean     = FALSE, 
  precipitation_spring_mean     = FALSE,
  precipitation_summer_mean     = FALSE,
  precipitation_autumn_mean     = FALSE,
  
  temperature_winter_mean       = TRUE,
  temperature_spring_mean       = FALSE,
  temperature_summer_mean       = TRUE,
  temperature_autumn_mean       = FALSE,
  
  relative_humidity_winter_mean = TRUE,
  relative_humidity_spring_mean = TRUE,
  relative_humidity_summer_mean = FALSE,
  relative_humidity_autumn_mean = FALSE,
  
  primary_forested_area         = TRUE,
  primary_non_forested_area     = TRUE,
  secondary_forested_area       = TRUE,
  croplands                     = FALSE,
  pasture_land                  = TRUE,
  rangelands                    = TRUE,
  urban_land                    = TRUE,
  
  proximity = TRUE, # Replace 'proximity' with TRUE or modify as required
  presence_previous_year = FALSE,
  year1 = FALSE # not used
)



# coerce logical (or 0/1) to integer 0/1 robustly
b <- function(x) as.integer(isTRUE(x))

# Creating a parameter string for selected/non selected model covariate's
parameters_str <- paste0(
  "param",
  b(parameters$population_density),
  b(parameters$total_flux_in),
  
  # Precipitation (seasonal)
  b(parameters$precipitation_winter_mean),
  b(parameters$precipitation_spring_mean),
  b(parameters$precipitation_summer_mean),
  b(parameters$precipitation_autumn_mean),
  
  # Temperature (seasonal)
  b(parameters$temperature_winter_mean),
  b(parameters$temperature_spring_mean),
  b(parameters$temperature_summer_mean),
  b(parameters$temperature_autumn_mean),
  
  # Relative humidity (seasonal)
  b(parameters$relative_humidity_winter_mean),
  b(parameters$relative_humidity_spring_mean),
  b(parameters$relative_humidity_summer_mean),
  b(parameters$relative_humidity_autumn_mean),
  
  # Land cover + others
  b(parameters$primary_forested_area),
  b(parameters$primary_non_forested_area),
  b(parameters$secondary_forested_area),
  b(parameters$croplands),
  b(parameters$pasture_land),
  b(parameters$rangelands),
  b(parameters$urban_land),
  
  # Proximity + fixed-effect flags
  b(parameters$proximity),
  b(parameters$presence_previous_year),
  b(parameters$year1)
)



tau_value <- 0.5     # spde Q parameter 



## fixed effects
fixed_effect_cols <- c() # create and empty list

if (parameters$presence_previous_year) {
  fixed_effect_cols <- c(fixed_effect_cols, "presence_previous_year")
}
if (parameters$year1) {
  fixed_effect_cols <- c(fixed_effect_cols, "year1")
}
if (parameters$primary_forested_area) {
  fixed_effect_cols <- c(fixed_effect_cols, "primary_forested_area")
}
if (parameters$primary_non_forested_area) {
  fixed_effect_cols <- c(fixed_effect_cols, "primary_non_forested_area")
}
if (parameters$secondary_forested_area) {
  fixed_effect_cols <- c(fixed_effect_cols, "secondary_forested_area")
}
if (parameters$croplands) {
  fixed_effect_cols <- c(fixed_effect_cols, "croplands")
}
if (parameters$pasture_land) {
  fixed_effect_cols <- c(fixed_effect_cols, "pasture_land")
}
if (parameters$rangelands) {
  fixed_effect_cols <- c(fixed_effect_cols, "rangelands")
}
if (parameters$urban_land) {
  fixed_effect_cols <- c(fixed_effect_cols, "urban_land")
}

if (parameters$total_flux_in) {
  fixed_effect_cols <- c(fixed_effect_cols, "total_flux_in")
}

## if fixed_effect_cols is empty then assign a default value
if (is.null(fixed_effect_cols) || length(fixed_effect_cols) == 0) {
  fixed_effect_cols <- ""
}

# Print(fixed_effect_cols)




### proximity random effect
proximity_type <- "besag" # or "rw1" (i.e. correlated to all neighboring bins (ICAR) or only previous (sequential) bins (rw1))

# Specify nb-graph (neighborhood graph) for proximity covariate 
if (parameters$proximity) {
  
  proximity_col <- "bin_idx_y_proximity_exp300"
  
  if (proximity_type %in% c("ICAR", "besag")) {
    
    # NOTE : ICAR / Besag-York model for proximity.
    # One connected graph over all bins
    nb_grid_proximity <-
      specify_network_by_nodes_and_links(
        n_nodes = max(data_spde__[[proximity_col]], na.rm = TRUE)
      )
    
    file_nb_proximity <- tempfile("graph_proximity_", fileext = ".adj") # # temporary directory
    
    nb2INLA(file_nb_proximity, nb_grid_proximity)
    
  } else if (identical(proximity_type, "rw1")) {           # define like all other covariates
    
    # F(proximity_bin_idx, model = "rw1", ...)
    
  } else {
    
    stop("proximity_type must be 'besag' (or 'ICAR') or 'rw1'")
  }
  
} else {
  proximity_col <- ""
}




## defining the spatio-temporal effects (precision matrix) + intercept 
is_free <- TRUE        ## set FALSE if we want to control and manually define tuning parameters of precision matrix for the spde (scale or precision, range or shape)

if (!is_free) {
  
  f_basic <- paste0("-1 + Intercept + f(spde_idx, model = \"generic0\", Cmatrix = Q_spde, ",      ## Q_spde is defined above with manual tuning parameters, not used in final version (used only for spatial effects assessment)
                    "group = spde_idx.group, control.group=list(model=\"iid\"))")
} else {
  f_basic <- paste0("-1 + Intercept + f(spde_idx, model = spde_model, ",
                    "group = spde_idx.group, control.group=list(model=\"iid\"))")   ## iid with respect to time (i.e. no correlation in time)
}



# Define full INLA model (in standard syntax)
f_list <- define_model( 
  response = response_column,
  f_basic = f_basic,                       ## spatio-temporal structure of model
  fixed_effect_cols = fixed_effect_cols, # all the fixed effect are here ()
  proximity_col = proximity_col,            ## besag or rw1
  proximity_type  = proximity_type,
  
  # Random effects
  # Seasonal precipitation
  precipitation_winter_mean = parameters$precipitation_winter_mean,        ## TRUE OR FALSE as defined in parameters
  precipitation_spring_mean = parameters$precipitation_spring_mean,
  precipitation_summer_mean = parameters$precipitation_summer_mean,
  precipitation_autumn_mean = parameters$precipitation_autumn_mean,
  
  # Seasonal temperature
  temperature_winter_mean   = parameters$temperature_winter_mean,
  temperature_spring_mean   = parameters$temperature_spring_mean,
  temperature_summer_mean   = parameters$temperature_summer_mean,
  temperature_autumn_mean   = parameters$temperature_autumn_mean,
  
  # Seasonal relative humidity
  relative_humidity_winter_mean = parameters$relative_humidity_winter_mean,
  relative_humidity_spring_mean = parameters$relative_humidity_spring_mean,
  relative_humidity_summer_mean = parameters$relative_humidity_summer_mean,
  relative_humidity_autumn_mean = parameters$relative_humidity_autumn_mean,
  
  # Other covariates exposed in define_model()
  population_density = parameters$population_density,
  spatial_idx_col = NULL
)




strategy <- "simplified.laplace"      ## can be set as auto too (for best and fast computation)




# Add model choices (or model id) to model_str
model_str <- model_id # # depends on numbers of parameters we are including

model_str <- paste0(model_str, parameters_str)


if (strategy=="laplace") {
  model_str <- paste0(model_str, "_LAPL")
} else if (strategy=="simplified.laplace") {
  model_str <- paste0(model_str, "_SIMP")     ## adding approximation strategy to model name string
}



# So, replicating the precision matrix Q for all years for spatio-temporal structure in final INLA run
if (grepl("spde_idx.group", toString(f_list$f[[3]]))) {                      ## if iid temporal effect is present  
  n_years_Q <- length(unique(data_spde__$year_idx))                       #length of years of traiaing data - 14 ( 2010 to 2023)
} else if (grepl("Q_spde|spde_model", toString(f_list$f[[3]]))) {
  n_years_Q <- 1
} else {
  stop("not implemented")
}



n_years_pred <- length(unique(pred_data$year_idx)) # prediction_years i.e. length of years of prediction data (future)




# Run inla model for all years prediction in loop

########### Important-------------------
##For prediction year 1 (2022) it uses the proximity which is created from training data (ECDC data) using the observed presence/absence for year 2021 (which is same as proximity for training data for 2022, used to initalize the proximity for prediction year 2022)

#and if there is a next prediction year:
###Uses compute_proximity_from_posterior() on the posterior samples for the current year,
###Extracts the next year’s y_proximity_exp300,
###Reorders to match the geo order of pred_data for that next year,
###Bins it with attach_proximity_bins_from_spec(),
###Keeps the resulting binned proximity to plug into the next year’s pred_data.
######## Important-----------------------


# Before running code we need two things -
# For binning of each year proximity created from posterior predictions, load the training bins information (saved during binning process of training and prediction data in code load_and_process_data.R)
out_prox <- load(".../Output/A-Final_output/Proximity_binning_info_for_pred/prox300_binning_spec.RData")
prox300_spec <- prox300_spec

#  check: prediction tibble must have a 'year' column
if (!"year" %in% names(pred_data)) {
  stop("`pred_data` must have a `year` column for prediction data (e.g. calendar year).")
}

## loading the 2022 year proximity created from observed presence/absence data to initialize the proximity for future prediction staring from 2022
proximity_2022 <- load(".../Output/A-Final_output/Proximity_for_some years_for_prediction/Proximity_for_year_2022_for_prediction.RData")
proximity_2022_initial <- proximity_for_year_2022_for_prediction


## ------------------------------------------------------------------------
## Basic setup: years, containers, coords
## ------------------------------------------------------------------------


# Prediction year_idx values (we assume they are 1,2,...,K in order)
years_pred <- sort(unique(pred_data$year_idx))

# Lists to store results over all prediction years
all_spde_region_pred <- list()
all_map_full_pred    <- list()
all_posterior_pred   <- list()
all_prox_binned      <- list()

prox_binned_by_year <- list()


# Loop over prediction years
for (i in seq_along(years_pred)) {
  
  # I = 1
  
  year_idx_curr <- years_pred[i]         ## each year_idx one by one (starting from year_idx ==1)
  
  # Current prediction data (one year_idx at a time)
  pred_data_curr <- pred_data %>%
    dplyr::filter(year_idx == year_idx_curr)
  
  ## Calendar year assignment
  year_curr <- unique(pred_data_curr$year)
  stopifnot(length(year_curr) == 1)
  
  if (i == 1) {
    
    ## We use the observed data proximity (2022) to initialize the proximity for 1st year future prediction (assigning the year idx)
    df_prox300_init <- proximity_2022_initial %>%
      dplyr::select(geo, y_proximity_exp300) %>%
      dplyr::mutate(year_idx = year_idx_curr)
    
    #aligning initial proximity (2022) geo order to the geo ordering of pred_data 
    geo_order_first <- pred_data %>%
      dplyr::filter(year_idx == year_idx_curr) %>%
      dplyr::pull(geo)
    
    df_prox300_init <- df_prox300_init %>%
      dplyr::mutate(geo = factor(geo, levels = geo_order_first)) %>%
      dplyr::arrange(geo) %>%
      dplyr::mutate(geo = as.character(geo)) # drop factor again
    
    ### binning of proximity using the information of bins created during training 
    df_prox300_init_binned <- attach_proximity_bins_from_spec(
      new_df   = df_prox300_init,
      bin_spec = prox300_spec
    )
    
    # Drop any existing proximity columns in current prediction data
    pred_data_curr_clean <- pred_data_curr %>%
      dplyr::select(
        -dplyr::any_of(c("y_proximity_exp300",
                         "bin_idx_y_proximity_exp300"))
      )
    
    # Join the initial proximity + bin index to prediction data by geo and year_idx
    pred_data_curr_effective <- pred_data_curr_clean %>%
      dplyr::left_join(
        df_prox300_init_binned %>%
          dplyr::select(geo, year_idx,
                        y_proximity_exp300, # join the proximity column
                        bin_idx_y_proximity_exp300), # then join the binning information of proximity column
        by = c("geo", "year_idx")
      )
    
  } else {
    # For all subsequent prediction years, we update proximity using the proximity computed from the posterior predictions of the previous year.
    prox_df <- prox_binned_by_year[[as.character(year_idx_curr)]]      ## "prox_binned_by_year" is created later in this code
    if (is.null(prox_df)) {
      stop("No binned proximity found for year_idx = ", year_idx_curr,
           ". Make sure the loop fills prox_binned_by_year correctly.")
    }
    
    # Drop any old proximity columns from current pred_data
    pred_data_curr_clean <- pred_data_curr %>%
      dplyr::select(
        -dplyr::any_of(c("y_proximity_exp300",
                         "bin_idx_y_proximity_exp300"))
      )
    
    # Join new proximity + bin index by geo and year_idx
    pred_data_curr_effective <- pred_data_curr_clean %>%
      dplyr::left_join(
        prox_df %>% dplyr::select(geo, year_idx,
                                  y_proximity_exp300,
                                  bin_idx_y_proximity_exp300),
        by = c("geo", "year_idx")
      )
  }
  
  # Run the SPDE model for this prediction year
  out <- run_spde_model2(
    parameters      = parameters,
    data_spde       = data_spde,
    map_spde        = map_spde,
    n_spde          = n_spde,
    n_years         = n_years_Q,
    n_years_pred    = length(unique(pred_data_curr_effective$year_idx)), # 1; mostly legacy as we go year by year prediction
    f_list          = f_list,
    response_column = response_column,
    Q_spde          = Q_spde,
    strategy        = strategy,
    n_binom         = n_binom,
    model_str       = model_str,
    pred_data       = pred_data_curr_effective # current year prediction data only
  )
  
  fitted_model <- out$fitted_model
  print(summary(fitted_model))
  
  
  A_pred <- out$A_pred # # testing data A_pred
  
  spde_indices <- out$spde_indices
  
  idx_linear_predictor_pred = out$idx_linear_predictor_pred # prediction idx_linear_predictor
  
  # Extract SPDE results for prediction (region level)
  # - pred_data (geo × year_idx layout for prediction)
  
  collected_results_spde_pred <- extract_spde_results_prediction(
    fitted_model              = fitted_model,
    response_column           = response_column,
    spde_mesh                 = spde_mesh, # mesh object used in the SPDE model
    map_spde                  = map_spde, # sf with columns geo, region_id, geometry
    data_spde                 = data_spde, # training data (used to infer K_train)
    pred_data                 = pred_data_curr_effective, # prediction data (rows in prediction stack)
    idx_linear_predictor_pred = idx_linear_predictor_pred, # indices of prediction rows
    A_pred                    = A_pred, # optional; we use projector inside (as in training)
    includes_time             = TRUE,
    spde_indices              = spde_indices, # from inla.spde.make.index, must contain $spde_idx.group
    order_of_year_of_prediction_data =
      sort(unique(pred_data_curr_effective$year_idx)), # here just current year_idx
    projection_locations      = region_coords # centroid locations
  )
  
  
  tb_region_results_pred <- collected_results_spde_pred$tb_region_results %>%
    dplyr::mutate(year = year_curr) %>% # add calendar year
    dplyr::relocate(year, .after = year_idx) # move year right after year_idx
  
  
  ##  Extract prediction-only fitted values (map_full_model_pred)
  collected_results_pred <- extract_main_results_prediction(
    model                     = fitted_model,
    data_pred                 = pred_data_curr_effective,
    y_column                  = response_column, # e.g. "presence" 
    idx_linear_predictor_pred = idx_linear_predictor_pred # will be linear_predictor pred for prediction indexes
  )
  
  # Important notes for collected_results_pred:
  #   - y_column in prediction data will typically be NA (no observed labels), but we keep it in the output for consistency.
  #   - We do NOT compute error, deviance, DIC, WAIC, or any summaries here.
  #   - linear predictor is computed manually as logit(p.mean), i.e. using the mean fitted probabilities, not model$summary.linear.predictor.
  
  
  map_full_model_pred <- collected_results_pred$map_full_model_pred %>%
    dplyr::mutate(year = year_curr) %>% # add calendar year
    dplyr::relocate(year, .after = year_idx) # move year right after year_idx
  
  
  
  ## Draw posterior samples for prediction year
  p_post_pred <- extract_posterior_fitted_samples_prediction(
    model                     = fitted_model,
    data_pred                 = pred_data_curr_effective,
    idx_linear_predictor_pred = idx_linear_predictor_pred,
    n_samples                 = 1000L,
    likelihood                = "binomial", # or "bernoulli"; with n_binom = 1 they are equivalent
    seed                      = 42
  ) %>%
    dplyr::mutate(year = year_curr) %>% # add calendar year
    dplyr::relocate(year, .after = year_idx) # move year right after year_idx
  
  # Store all the above three outputs for this year
  all_spde_region_pred[[i]] <- tb_region_results_pred
  all_map_full_pred[[i]]    <- map_full_model_pred
  all_posterior_pred[[i]]   <- p_post_pred
  
  
  if (i < length(years_pred)) {                     ## if there is a next prediction year, compute proximity for it
    
    year_idx_next <- years_pred[i + 1]
    
    # construct a “posterior-based invasion history” and its distance-based proximity for prediction runs 
    df_proximity_post <- compute_proximity_from_posterior(
      p_post_pred = p_post_pred, # posterior samples for current prediction year
      map_inla    = map_spde, # sf object with "geo" and geometry
      exponent_vec = c(0.5, 1, 2, 3)
    )
    
    
    df_prox300 <- df_proximity_post %>%
      dplyr::filter(year_idx == year_idx_next) %>% # keep only next year's proximity
      dplyr::select(geo, year_idx, y_proximity_exp300)
    
    geo_order_next <- pred_data %>%
      dplyr::filter(year_idx == year_idx_next) %>%
      dplyr::pull(geo)
    
    df_prox300 <- df_prox300 %>%
      dplyr::mutate(geo = factor(geo, levels = geo_order_next)) %>%
      dplyr::arrange(geo) %>%
      dplyr::mutate(geo = as.character(geo))
    
    ##  Apply binning to this next year proximity computed from predicted distribution using the training bin specification
    df_prox300_binned <- attach_proximity_bins_from_spec(
      new_df   = df_prox300,
      bin_spec = prox300_spec
    )
    
    # Add calendar year and move it after year_idx
    # the proximity is for next year so next year should be added (not the year_curr which is added in collected_results_spde or map_full_model)
    year_next <- pred_data %>%
      dplyr::filter(year_idx == year_idx_next) %>%
      dplyr::pull(year) %>% # # pull the year from pred data
      unique()
    
    if (length(year_next) != 1L) {
      stop("Prediction data for year_idx_next does not have a unique `year` value.")
    }
    
    df_prox300_binned <- df_prox300_binned %>%
      dplyr::mutate(year = year_next) %>%
      dplyr::relocate(year, .after = year_idx)
    
    prox_binned_by_year[[as.character(year_idx_next)]] <- df_prox300_binned     ## used to update the next year proximity in prediction data in beginning of loop for next year prediction
    
    all_prox_binned[[i + 1]] <- df_prox300_binned
  }
  
  # End of loop iteration for year_idx_curr
  cat("Finished prediction for year_idx =", year_idx_curr, # progress message for this loop iteration
      "(calendar year", year_curr, ")\n")
  flush.console()
  
}


# Bind results over all prediction years and return
tb_region_results_pred_all <- dplyr::bind_rows(all_spde_region_pred)
map_full_model_pred_all   <- dplyr::bind_rows(all_map_full_pred)
p_post_pred_all           <- dplyr::bind_rows(all_posterior_pred)

# Proximity binned list may have NULL for the first year (as we did not recompute proximity for year 1, as we used initial proximity from training data 2010 to do proximity for year 2011) so we only bind non-NULL elements.
prox_binned_non_null <- all_prox_binned[!vapply(all_prox_binned, is.null, logical(1))]

df_prox300_binned_all <- if (length(prox_binned_non_null) > 0) {
  dplyr::bind_rows(prox_binned_non_null)
} else {
  NULL
}
# df_prox300_binned_all starts with yer_idx = 2



# Save future prediction outputs (2022–2100)
full_model_dir <- ".../Output/A-Final_output/Main_model_output/"

# Make sure the directory exists
dir.create(full_model_dir, recursive = TRUE, showWarnings = FALSE)

pred_output_file <- file.path(
  full_model_dir,
  "future_prediction_picontrol_climate_ssp370_socioeconomic_2022_2080_output.RData"
)

save(
  tb_region_results_pred_all, # region-level summaries for all prediction years
  map_full_model_pred_all, # map-level summaries for all prediction years
  p_post_pred_all, # posterior samples for all predicted region-year combos
  df_prox300_binned_all, # all binned proximity for future years
  file = pred_output_file
)

message("Saved counterfactual prediction outputs to: ", pred_output_file)

