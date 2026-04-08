# ================== Factual model fitting and analysis of random and fixed effects ==================


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



climate_data_type_name = "pi_control"
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

# Prediction cleaned data
out <- load(paste0(binned_data_folder, "final_binned_counter_factual_data.RData"))
prediction_data <- final_counter_covariates_binned


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



### Now prediction data 
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
pred_map <- map_spde  ## same as training map in our case





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


## number of prediction years (if present)
n_years_pred <- length(unique(pred_data$year_idx))

model_str__ <- model_str

# NO Cross validation (RUN the entire Training data)





# Run INLA model (if prediction data is there)

# Parameters = parameters,
# Data_spde = data_spde,
# Map_spde = map_spde,
# N_spde = n_spde,
# N_years = n_years_Q,
# N_years_pred = n_years_pred,      ## number of years in prediction data
# F_list = f_list,
# Response_column = response_column,
# Q_spde = Q_spde,
# Strategy = strategy,
# N_binom = n_binom,
# Model_str = model_str,
# Pred_data = pred_data # non-NULL (supply the prediction data)
# )

## no prediction data (only training and model analysis)
out <- run_spde_model2(
  parameters   = parameters,
  data_spde    = data_spde,
  map_spde     = map_spde,
  n_spde       = n_spde,
  n_years      = n_years_Q,
  n_years_pred = 0, # not really used any more as prediction data is NULL
  f_list       = f_list,
  response_column = response_column,
  Q_spde       = Q_spde,
  strategy     = strategy,
  n_binom      = n_binom,
  model_str    = model_str,
  pred_data    = NULL          ## NULL if no prediction data needed
)

fitted_model <- out$fitted_model
summary(out$fitted_model)

## out (out of model) has a list  containing ----
# list(
#   fitted_model             = fitted_model,
#   A_train                  = A_train,
#   A_pred                   = A_pred,               ## will be NULL if pred_data  is NULL
#   idx_linear_predictor     = idx_linear_predictor,
#   idx_linear_predictor_pred = idx_linear_predictor_pred,    ## will be NULL if pred_data is NULL
#   spde_indices             = spde_indices,
#   computation_time         = computation_time


A <- out$A_train
# A_pred will be NULL if pred_datais NULL
A_pred <- out$A_pred 

idx_linear_predictor <- out$idx_linear_predictor
spde_indices <- out$spde_indices

idx_linear_predictor_pred = out$idx_linear_predictor_pred # prediction idx_linear_predictor  (if preiction data is supplied)




# Collected_results_spde collects the spatial effects over training
collected_results_spde <- extract_spde_results_training(fitted_model=fitted_model,
                                                        response_column=response_column,
                                                        spde_mesh=spde_mesh,
                                                        map_spde=map_spde, # map data
                                                        data_spde=data_spde,
                                                        projection_locations=region_coords,
                                                        idx_linear_predictor=idx_linear_predictor,
                                                        A=A,     #  fitted spde integrated to region level using the projection matrix A projection matrix A, it changes w.r.t space as well as time 
                                                        ##spde effect prior is matern correlated in space and iid in time
                                                        includes_time=TRUE,
                                                        spde_indices=spde_indices,
                                                        order_of_year_of_training_data = unique(data_spde$year_idx))
# Don't need to supply any prediction information here



# Collect relevant results
collected_results <- extract_main_results_training(model=fitted_model,
                                                   data=data_spde,
                                                   likelihood="binomial", n_binom=n_binom,      ## binomial with n_binom = 1
                                                   spatial_random_effect="spde_idx",           ## for spatial_random_effect = "spde_idx"
                                                   idx_linear_predictor=idx_linear_predictor,
                                                   random_effect_with_time_pattern="proximity|temperature|humidity|population_density",     ## assigned bin changes w.r.t for each of these covariates for all regions and for all years
                                                   data_pred=NULL,                             # no need for prediction data right now
                                                   idx_linear_predictor_pred= NULL,
                                                   y_column="presence",
                                                   cross_val=FALSE)                            # no cross validation for now



## adding the overall standard deviation (hyperparameter)  of spde effects for all regions/year to collected_results$results_summary$stdev_spatial
collected_results$results_summary$stdev_spatial <-
  sd(collected_results_spde$tb_region_results$spde_integrated) # sd of spde effects for all regions/year


# Collect statistics of fit (only needed if we manually control the hyperpriors of spatial effects (more precisely precision matrix Q)
if (!is_free) {
  rank_Q <- rankMatrix(Q_spde, method="qr.R")
  
  year_vec <- collected_results_spde$tb_spde_results$year_idx %>% unique()
  x <- collected_results_spde$tb_spde_results %>%
    filter(year_idx==year_vec[1]) %>% pull(spde_weights.mean)
  stats <-
    mvnorm_stats_from_Q(
      x=x,
      Q=Q_spde,
      rank_Q=rank_Q,
      tau=collected_results$results_summary$hyperpar_tau
    )
  collected_results$results_summary <- bind_cols(collected_results$results_summary,
                                                 log_density=stats$log_density,
                                                 log_density_tau1=stats$log_density_tau1,
                                                 norm_z=stats$norm_z)
} else {
  collected_results$results_summary <- bind_cols(collected_results$results_summary,
                                                 log_density=NA,
                                                 log_density_tau1=NA,
                                                 norm_z=NA)
}



tb_spde_results <- collected_results_spde$tb_spde_results
tb_region_results <- collected_results_spde$tb_region_results


# Taking out all the variables from collected_results
results_summary <- collected_results$results_summary
map_fixed_effects <- collected_results$map_fixed_effects
map_spatial_effect <- collected_results$map_spatial_effect
map_random_effects <- collected_results$map_random_effects
map_random_effects_with_time <- collected_results$map_random_effects_with_time
functions_random_effects <- collected_results$functions_random_effects
map_full_model <- collected_results$map_full_model

results_summary$corr_spde_with_data <- compute_correlation_with_data(tb_region_results, data_spde)


# Printing some results
print(results_summary)
print(results_summary$fixed_effects)
print(results_summary$hyperpar)



# this functions collects results used in plotting
map_complete_model_mapped <- compute_p_from_effects(map_full_model,
                                                    map_random_effects_with_time,
                                                    map_random_effects,
                                                    map_fixed_effects,
                                                    tb_region_results,
                                                    results_summary)


## This function is useful for plotting the covariate contribution of each covariate to prediction each region for all years 
map_all <- collect_all_model_effects(map_full_model, tb_region_results, map_fixed_effects,
                                     map_random_effects, map_random_effects_with_time)



#Drawing posteriors samples for each region/year 
# Factual model
p_samples_long <- extract_posterior_fitted_samples(
  model                = fitted_model,
  data                 = data_spde,
  idx_linear_predictor = idx_linear_predictor,
  n_samples            = 1000,
  likelihood           = "binomial",
  seed                 = 42
)

# p_samples_long returns a LONG tibble with one row per (region, year, draw):
#   columns:
#     - data_idx  : row index in the original data_spde
#     - geo       : NUTS3 (or region) ID
#     - region_id : numeric region index
#     - year_idx  : numeric year index (1..13 in our case)
#     - draw      : posterior sample index m = 1..n_samples
#     - p         : sampled fitted probability p_i^(fact)(m) for that
#                  region/year in that posterior draw


## SAVING the model results
# For main model
full_model_dir <- ".../Output/A-Final_output/Main_model_output/"

save(collected_results_spde, collected_results, map_full_model, map_all, p_samples_long,
     file=paste0(full_model_dir, "main_model", model_str, "_id", model_id, ".RData"))











# LOADING the saved results

m = load(".../Output/A-Final_output/Main_model_output/main_modelvIparam110000101011001110111100_SIMP_idvI.Rdata")
m



# ############### PLOTTING the results ############################################


# Plotting all the random effects curve in separate plots
plot_functions_random_effects_one_by_one <- function(fun_re,
                                                     data_spde,
                                                     out_dir = NULL,
                                                     width = 6,
                                                     height = 4,
                                                     dpi = 300,
                                                     ext = "png") {
  # Pretty x-axis label from effect name (bin_idx_*)
  pretty_xlabel <- function(idx_name) {
    base <- idx_name %>%
      str_remove("^bin_idx_") %>%
      str_remove("_mean$")
    
    if (str_detect(base, "proximity"))            return("Proximity")
    if (str_detect(base, "population_density"))   return(expression("Population density(people/km"^2*")"))
    
    # Season-first ordering for common climate variables
    if (str_detect(base, "^temperature_(winter|spring|summer|autumn)")) {
      season <- str_match(base, "^temperature_(winter|spring|summer|autumn)")[,2]
      return(paste(str_to_title(season), "temperature (°C)"))
    }
    if (str_detect(base, "^relative_humidity_(winter|spring|summer|autumn)")) {
      season <- str_match(base, "^relative_humidity_(winter|spring|summer|autumn)")[,2]
      return(paste(str_to_title(season), "relative humidity (%)"))
    }
    if (str_detect(base, "^precipitation_(winter|spring|summer|autumn)")) {
      season <- str_match(base, "^precipitation_(winter|spring|summer|autumn)")[,2]
      return(paste(str_to_title(season), "precipitation (mm)"))
    }
    
    # Fallback
    base %>% str_replace_all("_", " ") %>% str_to_sentence()
  }
  
  # Stack and standardize column names
  df_long <- imap_dfr(fun_re, function(tb, nm) {
    idx_col <- names(tb)[1] # first column is the bin-index col
    tb %>%
      transmute(
        effect  = nm,
        bin_idx = .data[[idx_col]],
        mean    = .data$mean,
        lower   = .data$lower_bound,
        upper   = .data$upper_bound
      )
  })
  
  # Build plots one-by-one with x = bin_mean_* 
  plots <- df_long %>%
    group_split(effect, .keep = TRUE) %>%
    set_names(map_chr(., ~ unique(.$effect))) %>%
    imap(function(df_eff, eff_idx_name) {
      
      # Find the matching bin_mean_* column in data_spde
      mean_col <- str_replace(eff_idx_name, "^bin_idx_", "bin_mean_")
      
      stopifnot(all(c(eff_idx_name, mean_col) %in% names(data_spde)))
      lut <- data_spde %>%
        distinct(
          bin_idx  = .data[[eff_idx_name]],
          bin_mean = .data[[mean_col]]
        ) %>%
        filter(!is.na(bin_idx), !is.na(bin_mean)) %>%
        distinct(bin_idx, .keep_all = TRUE) %>%
        arrange(bin_idx)
      
      # Join bin_mean onto df_eff
      df_eff2 <- df_eff %>%
        left_join(lut, by = "bin_idx") %>%
        arrange(bin_mean)
      
      # X/y labels
      x_lab <- pretty_xlabel(eff_idx_name)
      y_lab <- "Random effect"
      
      limit_y <- stringr::str_detect(eff_idx_name, "precipitation") |
        stringr::str_detect(eff_idx_name, "population")
      
      p <- ggplot2::ggplot(df_eff2, ggplot2::aes(x = bin_mean, y = mean)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2) +
        ggplot2::geom_line(linewidth = 0.3) +
        ggplot2::geom_point(size = 0.4) +
        ggplot2::labs(x = x_lab, y = y_lab) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          panel.grid   = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(color = "black", size = 7, family = "Arial"),
          axis.title.y = ggplot2::element_text(color = "black", size = 7, family = "Arial"),
          axis.text.x  = ggplot2::element_text(color = "black", size = 6, family = "Arial"),
          axis.text.y  = ggplot2::element_text(color = "black", size = 6, family = "Arial"),
          axis.ticks   = ggplot2::element_line(color = "black", linewidth = 0.2),
          axis.line    = ggplot2::element_line(color = "black", linewidth = 0.2)
        )
      
      if (limit_y) p <- p + ggplot2::coord_cartesian(ylim = c(-2, 2))   ## manual limit imposed on y axis for consistent y axis among all plots for visual comaprison
      
      p
    })
  
  # Optionally save each plot
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    walk2(plots, names(plots), ~ {
      f <- file.path(out_dir, paste0(.y, ".", ext))
      ggsave(filename = f, plot = .x, width = width, height = height, dpi = dpi)
    })
  }
  
  plots
}


plots_list <- plot_functions_random_effects_one_by_one(
  fun_re    = collected_results$functions_random_effects,
  data_spde = data_spde,
  out_dir   = NULL # or NULL to skip saving
)


# Specifying where to save
out_dir_save = ".../Output/A-Final_output/Programmatic_covariate_search/Full_model/Plots/random_effects_plot"


# Saving Params
ext         <- "tiff"
width_cm    <- 5.7
height_cm   <- 5.0
dpi         <- 600 # Higher DPI 
compression <- "lzw" # Required for TIFF to keep file size small
bg          <- "transparent"

dir.create(out_dir_save, recursive = TRUE, showWarnings = FALSE)
# Make filenames safe
sanitize_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

walk(names(plots_list), function(nm) {
  p <- plots_list[[nm]]
  file <- file.path(out_dir_save, paste0(sanitize_name(nm), ".", ext))
  
  ggplot2::ggsave(
    filename   = file,
    plot       = p,
    device     = ragg::agg_tiff,
    width      = width_cm,
    height     = height_cm,
    units      = "cm", # Switched to cm for precision
    dpi        = dpi,
    compression = compression,
    background = bg
  )
})

#  random effects plots in svg format
out_dir_save = ".../Output/A-Final_output/Programmatic_covariate_search/Full_model/Plots/random_effects_plot/svg"

ext      <- "svg"
width_cm <- 5.7
height_cm<- 5.0
bg       <- "transparent" # or "white" 

dir.create(out_dir_save, recursive = TRUE, showWarnings = FALSE)

# standardize filenames 
sanitize_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

purrr::walk(names(plots_list), function(nm) {
  p <- plots_list[[nm]]
  file <- file.path(out_dir_save, paste0(sanitize_name(nm), ".", ext))
  
  ggplot2::ggsave(
    filename   = file,
    plot       = p,
    device     = svglite::svglite, # vector SVG backend
    width      = width_cm,
    height     = height_cm,
    units      = "cm",
    bg         = bg
  )
})






# Plotting the fixed effects (waterfall plot for all fixed effects)
fixed_effects <- collected_results$fixed_effects

# What we do:
# We drop the intercept (we only want covariate effects).
# clean names: "_" -> " ", Title Case.
# rename "total_flux_in" to "Mobility" for the plot.
# and draw a forest plot with a 0-reference line.

plot_and_save_fixed_effects_forest <- function(fixed_eff_tbl,
                                               out_dir = ".../Plots/A_final_plots/Fixed_effects_plots",
                                               ext = "svg", # "svg" or "tiff"
                                               width_cm = 5.7,
                                               height_cm = 5.0,
                                               dpi = 600, # used for TIFF only
                                               compression = "lzw", # used for TIFF only
                                               bg = "transparent", # or "white" 
                                               sort_by = c("mean", "absmean", "none")) {
  
  sort_by <- match.arg(sort_by)
  ext <- tolower(ext)
  
  df <- as.data.frame(fixed_eff_tbl)
  df$term <- rownames(df)
  
  df <- df %>%
    dplyr::filter(term != "Intercept") %>%
    dplyr::transmute(
      term,
      mean = mean,
      lwr  = `0.025quant`,
      upr  = `0.975quant`
    ) %>%
    dplyr::mutate(
      term  = dplyr::if_else(term == "total_flux_in", "Mobility", term),
      label = stringr::str_replace_all(term, "_", " "),
      label = stringr::str_to_title(label)
    )
  
  # We build the forest plot (all effects in forestgreen)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean, y = label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = lwr, xmax = upr),
      height = 0.18, linewidth = 0.4, color = "darkorange4"
    ) +
    ggplot2::geom_point(size = 1.0, color = "darkorange4") +
    ggplot2::labs(
      x = "Scaled covariate fixed effect (log odds)",
      y = "Driver"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(color = "black", size = 7, family = "Arial"),
      axis.title.y = ggplot2::element_text(color = "black", size = 7, family = "Arial"),
      axis.text.x  = ggplot2::element_text(color = "black", size = 6, family = "Arial"),
      axis.text.y  = ggplot2::element_text(color = "black", size = 6, family = "Arial"),
      axis.ticks   = ggplot2::element_line(color = "black", linewidth = 0.2),
      axis.line    = ggplot2::element_line(color = "black", linewidth = 0.2)
    )
  
  # save the plot (svg or tiff)
  if (!is.null(out_dir)) {
    # keep outputs organized by format
    out_dir2 <- file.path(out_dir, ifelse(ext %in% c("tiff", "tif"), "Tiff", "svg"))
    dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)
    
    file_out <- file.path(out_dir2, paste0("fixed_effects_forest.", ext))
    
    if (ext == "svg") {
      ggplot2::ggsave(
        filename = file_out,
        plot     = p,
        device   = svglite::svglite, # vector SVG backend
        width    = width_cm,
        height   = height_cm,
        units    = "cm",
        bg       = bg
      )
    } else if (ext %in% c("tiff", "tif")) {
      ggplot2::ggsave(
        filename    = file_out,
        plot        = p,
        device      = ragg::agg_tiff, # high-quality TIFF backend
        width       = width_cm,
        height      = height_cm,
        units       = "cm",
        dpi         = dpi,
        compression = compression,
        background  = bg
      )
    } else {
      stop("ext must be 'svg' or 'tiff' (or 'tif').")
    }
  }
  
  return(p)
}


# Save SVG to .../Fixed_effects_plots/svg/
p_fx <- plot_and_save_fixed_effects_forest(
  fixed_eff_tbl = fixed_effects[[1]],
  out_dir = ".../Plots/A_final_plots/Fixed_effects_plots",
  ext = "svg",
  width_cm = 9.0,
  height_cm = 8.0,
  bg = "transparent",
  sort_by = "none"
)

p_fx # # execute the saving code

# Save Tiff to .../Fixed_effects_plots/svg/
p_fx <- plot_and_save_fixed_effects_forest(
  fixed_eff_tbl = fixed_effects[[1]],
  out_dir = ".../Plots/A_final_plots/Fixed_effects_plots",
  ext = "tiff",
  width_cm = 9.0,
  height_cm = 8.0,
  bg = "transparent",
  sort_by = "none"
)

p_fx # # execute the saving code

