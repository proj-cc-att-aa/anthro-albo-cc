# ================== Complete cross-validation ==================
# Cross-validation workflow

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

# Now both train and pred data is ready


# Loading the map now
out <- attach("data_INLA.RData")
# For training
map_spde <- map_inla
# Prediction map
pred_map <- map_spde  ## same as training map in our case



proximity <- TRUE # use proximity covariate
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



## ############### CROSS VALIDATIONS ###############################################

## CONFUSION MATRIX - function for creating confusion matrix 
draw_confusion_matrix_present <- function(y_true, p_hat, thr = 0.5, main = "") {
  stopifnot(length(y_true) == length(p_hat))
  stopifnot(all(y_true %in% c(0,1)))
  stopifnot(all(is.finite(p_hat)))
  
  # Labels (force "present" as event/positive)
  truth_lab <- ifelse(y_true == 1, "present", "absent")
  pred_lab  <- ifelse(p_hat >= thr,  "present", "absent")
  
  # 2x2 counts with rows = Actual, cols = Predicted
  lv <- c("present", "absent") # put 'present' first (positive)
  tbl <- table(Actual = factor(truth_lab, levels = lv),
               Pred   = factor(pred_lab,  levels = lv))
  # Extract counts
  TP <- as.integer(tbl["present","present"])
  FN <- as.integer(tbl["present","absent"])
  FP <- as.integer(tbl["absent","present"])
  TN <- as.integer(tbl["absent","absent"])
  N  <- TP + TN + FP + FN
  
  # Metrics (all w.r.t. 'present')
  acc  <- (TP + TN) / N
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_ # recall
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_ # PPV
  f1   <- if (is.finite(prec + sens) && (prec + sens) > 0) 2*prec*sens/(prec+sens) else NA_real_
  
  # Simple Cohen's kappa (optional)
  p_obs <- acc
  p_exp <- ((TP+FN)/N)*((TP+FP)/N) + ((FP+TN)/N)*((FN+TN)/N)
  kappa <- if (p_exp < 1) (p_obs - p_exp)/(1 - p_exp) else NA_real_
  
  # Draw (base graphics)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(matrix(c(1,1,2))) # big matrix + details
  par(mar = c(2,2,2,2))
  
  # Main confusion matrix panel
  plot(c(100, 345), c(300, 450), type = "n", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  title(main, cex.main = 2)
  
  # Header row (Predicted classes)
  rect(150, 430, 240, 370, col = '#3F97D0') ; text(195, 435, 'present', cex = 1.2)
  rect(250, 430, 340, 370, col = '#F7AD50') ; text(295, 435, 'absent',  cex = 1.2)
  text(125, 370, 'Predicted', cex = 1.3, srt = 90, font = 2)
  text(245, 450, 'Actual',    cex = 1.3, font = 2)
  
  # Actual rows
  rect(150, 305, 240, 365, col = '#F7AD50') # (present, predicted present) cell background
  rect(250, 305, 340, 365, col = '#3F97D0') # (present, predicted absent) cell background
  text(140, 400, 'present', cex = 1.2, srt = 90)
  text(140, 335, 'absent',  cex = 1.2, srt = 90)
  
  # Counts in cells (white text, bold)
  text(195, 400, TP, cex = 1.6, font = 2, col = 'white') # True Positive
  text(295, 400, FP, cex = 1.6, font = 2, col = 'white') # False Positive
  text(195, 335, FN, cex = 1.6, font = 2, col = 'white') # False Negative
  text(295, 335, TN, cex = 1.6, font = 2, col = 'white') # True Negative
  
  # Details panel
  plot(c(100, 0), c(100, 0), type = "n", xlab = "", ylab = "", main = "DETAILS", xaxt = 'n', yaxt = 'n')
  # Column 1–5 labels & values
  txt_lab <- c("Sensitivity", "Specificity", "Precision", "F1", "Accuracy")
  txt_val <- c(sens, spec, prec, f1, acc)
  xs <- c(10, 30, 50, 70, 90)
  for (i in seq_along(xs)) {
    text(xs[i], 85, txt_lab[i],                   cex = 1.2, font = 2)
    text(xs[i], 70, round(as.numeric(txt_val[i]), 3), cex = 1.2)
  }
  # Bottom: Kappa + Threshold
  text(30, 35, "Kappa", cex = 1.3, font = 2)
  text(30, 20, round(kappa, 3), cex = 1.3)
  text(70, 35, "Threshold", cex = 1.3, font = 2)
  text(70, 20, thr, cex = 1.3)
  
  invisible(list(
    table  = tbl,
    counts = c(TP = TP, FP = FP, FN = FN, TN = TN),
    metrics = c(accuracy = acc, sensitivity = sens, specificity = spec,
                precision = prec, F1 = f1, kappa = kappa)
  ))
}



## functions of important metrics for cross validation 
bin_metrics <- function(y_true, p_hat, thr = 0.5) {
  stopifnot(length(y_true) == length(p_hat))
  # Threshold
  y_pred <- as.integer(p_hat >= thr)
  
  TP <- sum(y_true == 1 & y_pred == 1, na.rm = TRUE)
  TN <- sum(y_true == 0 & y_pred == 0, na.rm = TRUE)
  FP <- sum(y_true == 0 & y_pred == 1, na.rm = TRUE)
  FN <- sum(y_true == 1 & y_pred == 0, na.rm = TRUE)
  
  acc <- (TP + TN) / (TP + TN + FP + FN)
  sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_) # recall
  spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  prec <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  f1   <- ifelse(is.na(prec) | is.na(sens) | (prec + sens) == 0,
                 NA_real_, 2 * prec * sens / (prec + sens))
  tibble(TP, TN, FP, FN, accuracy = acc, sensitivity = sens, specificity = spec,
         precision = prec, F1 = f1)
}




# NOW RUNNING the cross validation
################################################################################
## temporal CV
################################################################################

# Cross-validation
cv_type <- "cv_years"

run_full_cv <- TRUE # !

n_cv <- 13L # # numbers of years of cross validation
# If n_cv = 1 then it hold out last 1 years (n_CV is 2L then years of cross validation is last two years and it goes one by one in ascending order of years i.e. first for year 12 and then for  year13)

# Choose CV years
if (!run_full_cv) {
  n_cv <- 0 # no CV
} else {
  n_cv <- n_cv
}

unique_years <- sort(unique(data_spde__$year_idx)) # # unique years in data
i_cv         <- seq.int(from = length(unique_years) - n_cv + 1L, to = length(unique_years)) # ##set of years to be cross validated
cv_set       <- unique_years[i_cv]
cv_set

# Safety: nothing to do?
if (length(cv_set) == 0L) {
  stop("No CV years selected (set run_full_cv=TRUE or increase n_cv).")
}


# Run CV
cv_results <- vector("list", length(cv_set))                   ## vector holding all results for each year cross validation

for (j_cv in seq_along(cv_set)) {
  year_holdout <- cv_set[[j_cv]]
  
  message(sprintf(">>> CV %d/%d — holding out year_idx = %s", j_cv, length(cv_set), year_holdout))
  
  # copy; keep original response in y_orig
  data_spde_train_test <- data_spde__ %>%
    mutate(
      training_set = as.integer(year_idx != year_holdout & !is.na(.data[[response_column]])),
      test_set     = as.integer(year_idx == year_holdout & !is.na(.data[[response_column]])),
      y_orig       = .data[[response_column]] # Original data column
    )
  
  # Mask the hold-out response
  data_spde_train_test[[response_column]][data_spde_train_test$test_set == 1] <- NA_real_
  
  # Tag string (optional)
  model_str <- sprintf("cv_year%s_%s", year_holdout, model_str__)
  
  # Fit once via wrapper model training
  out <- run_spde_model2(
    parameters     = parameters,
    data_spde      = data_spde_train_test, #train-test data set created
    map_spde       = map_spde, # map of training dataset
    n_spde         = n_spde,
    n_years        = n_years_Q,
    n_years_pred   = n_years_pred, ## numbers of years of prediction   (optional and can be set to 0)
    f_list         = f_list, ## model definition
    response_column= response_column, ## response coluumn - presence (set above)
    Q_spde         = Q_spde, 
    strategy       = strategy, 
    n_binom        = n_binom, 
    model_str      = model_str, 
    pred_data      = pred_data       ## prediction data is not needed here (just for consistency, can be set to NULL)
  )
  
  fitted_model <- out$fitted_model
  summary(fitted_model)
  
  
  A <- out$A_train
  A_pred <- out$A_pred # # testing data A_pred
  idx_linear_predictor <- out$idx_linear_predictor
  spde_indices <- out$spde_indices_train
  
  idx_eta <- idx_linear_predictor
  
  # INLA returns fitted probabilities for binomial/logit in summary.fitted.values
  p_all <- fitted_model$summary.fitted.values$mean[idx_eta] ## probabilities
  
  # pulling held-out rows only 
  idx_test <- which(data_spde_train_test$test_set == 1L)
  y_true   <- data_spde_train_test$y_orig[idx_test]
  p_hat    <- p_all[idx_test] # predicted value of test set
  
  
  # Setting STANDARD THRESHOLD
  # Metrics @ 0.5 threshold
  # Computing the metrics 
  m <- bin_metrics(y_true, p_hat, thr = 0.5) %>%
    mutate(year_idx = year_holdout)
  
  
  # ROC / AUC
  roc_obj <- tryCatch(pROC::roc(response = y_true, predictor = p_hat, quiet = TRUE),
                      error = function(e) NULL)
  auc_val <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
  
  # Draw + save confusion matrix 
  fig_width_px  <- 3000
  fig_height_px <- 2800
  fig_dpi       <- 600
  thr_std       <- 0.5
  
  
  # Directory
  out_dir_conf <- ".../Plots/A_final_plots/Cross_validation/Confusion_matrix_temporal"
  dir.create(out_dir_conf, recursive = TRUE, showWarnings = FALSE)
  
  fig_name <- sprintf(
    "confus_matrix_holdoutyear%03d_threshold%0.2f.png",
    as.integer(year_holdout), thr_std
  )
  
  fig_path <- file.path(out_dir_conf, fig_name)
  
  png(filename = fig_path,
      width = fig_width_px, height = fig_height_px, res = fig_dpi)
  
  draw_confusion_matrix_present(
    y_true, p_hat, thr = thr_std,
    main = "" # No title
  )
  dev.off()
  
  # Stash per-fold
  cv_results[[j_cv]] <- list(
    year_idx  = year_holdout,
    metrics   = m %>% mutate(AUC = auc_val),
    roc       = roc_obj,
    n_test    = length(idx_test)
  )
}

# Saving the results of Cross validation to .Rdata file

out_dir_cv <- ".../Output/A-Final_output/Cross_Validation_Output_temporal"
dir.create(out_dir_cv, recursive = TRUE, showWarnings = FALSE)

span_tag <- sprintf("%03d-%03d",
                    as.integer(min(cv_set)),
                    as.integer(max(cv_set)))

# Use span_tag in names
fig_name <- sprintf(
  "cv_results_temporal_holdoutyear_%s.RData",
  span_tag
)

# Saving .RData (contains the object named `cv_results`)
rdata_path <- file.path(out_dir_cv, fig_name)
save(cv_results, file = rdata_path)



# Now loading the CV metrics stored
out <- load(paste0(out_dir_cv, "/cv_results_temporal_holdoutyear_001-013.Rdata"))
cv_results <- cv_results

# Further processing of cross validation results

# Saving each year CV
year_map <- data_spde %>%
  dplyr::select(year_idx, year) %>%
  dplyr::distinct() %>%
  dplyr::arrange(year_idx)

# Collect + add year
cv_metrics_by_year <- dplyr::bind_rows(lapply(cv_results, `[[`, "metrics")) %>%
  dplyr::relocate(year_idx)

cv_metrics_by_year <- cv_metrics_by_year %>%
  dplyr::left_join(year_map, by = "year_idx") %>%
  dplyr::relocate(year, .after = year_idx)

print(cv_metrics_by_year)

# Save 
out_dir_cv_csv <- ".../Output/A-Final_output/Cross_Validation_Output_temporal"

csv_path <- file.path(out_dir_cv_csv, sprintf("cv_metrics_temporal_holdoutyears_%s.csv", span_tag))
readr::write_csv(cv_metrics_by_year, csv_path)

# Overall mean of CV 
cv_metrics_overall <- cv_metrics_by_year %>%
  summarise(
    across(c(accuracy, sensitivity, specificity, precision, F1, AUC), ~ mean(.x, na.rm = TRUE))
  )
cv_metrics_overall


out_dir_cv_csv_overall <- ".../Output/A-Final_output/Cross_Validation_Output_temporal"

csv_path_overall <- file.path(out_dir_cv_csv_overall, sprintf("cv_metrics_temporal_%s_overall.csv", span_tag))
readr::write_csv(cv_metrics_overall, csv_path_overall)

# Plot ROC curves for all held-out years 
roc_df <- map_dfr(cv_results, function(x) {
  if (is.null(x$roc)) return(NULL)
  tibble(
    year_idx = x$year_idx,
    tpr = rev(x$roc$sensitivities),
    fpr = rev(1 - x$roc$specificities)
  )
})

# ROC Curve
if (nrow(roc_df) > 0) {
  p <- ggplot2::ggplot(
    roc_df,
    ggplot2::aes(x = fpr, y = tpr, color = factor(year_idx), group = year_idx)
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::labs(
      x = "False Positive rate",
      y = "True Positve rate",
      color = "Test year",
      title = "ROC - AUC"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.title.x = element_text(color = "black", size = 20),
      axis.title.y = element_text(color = "black", size = 20),
      axis.line  = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.text  = ggplot2::element_text(color = "black"),
      panel.background = ggplot2::element_blank(),
      plot.background  = ggplot2::element_blank(),
      # Legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      # Legend sizes
      legend.title      = element_text(size = 18, face = "bold", colour = "black"),
      legend.text       = element_text(size = 16, colour = "black")
    )
  
  print(p) 
} else {
  message("No ROC data to plot.")
}

# Save block for ROC plot ARTICLE description 
ext         <- "tiff" # "tiff", "png", "pdf",
width_in    <- 8
height_in   <- 7
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

out_dir_ROC_AUC <- ".../Plots/A_final_plots/Cross_validation/ROC-AUC-Temporal"
dir.create(out_dir_ROC_AUC, recursive = TRUE, showWarnings = FALSE)

# File name 
out_file <- file.path(out_dir_ROC_AUC, paste0("roc_by_all_year.", ext))

# Save (TIFF, transparent) using ragg 
# Install.packages("ragg") 
ggplot2::ggsave(
  filename   = out_file,
  plot       = p, #  ggplot object
  device     = ragg::agg_tiff, #  TIFF device
  width      = width_in,
  height     = height_in,
  units      = "in",
  dpi        = dpi,
  compression = compression,
  background  = bg
)

message("Saved: ", out_file)


## Saving with standard plotting style
if (nrow(roc_df) > 0) {
  p <- ggplot2::ggplot(
    roc_df,
    ggplot2::aes(x = fpr, y = tpr, color = factor(year_idx), group = year_idx)
  ) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.2) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::labs(
      x = "False positive rate",
      y = "True positive rate",
      color = "Test year" # Shortened title to save space
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = element_text(family = "Arial"),
      # standard sizes
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6),
      axis.line  = element_line(linewidth = 0.2), 
      axis.ticks = element_line(linewidth = 0.2),
      
      # COMPACT LEGEND SETTING
      legend.title = element_text(
        size = 6,
        face = "bold",
        margin = margin(b = 2.5, unit = "pt")
      ),
      legend.text  = element_text(size = 5), # # legend text
      # Shrink the colored lines in the legend
      legend.key.size = unit(0.15, "cm"),
      legend.key.height = unit(0.1, "cm"), # # reduce the gap between legend
      # Legend position
      legend.position = c(0.92, 0.38),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    ) +
    # Force legend into one column to keep it narrow
    ggplot2::guides(color = guide_legend(ncol = 1))
  
  print(p)
}

out_dir_ROC_AUC_sci <- ".../Plots/A_final_plots/Cross_validation/ROC-AUC-Temporal/Science_journal"
dir.create(out_dir_ROC_AUC_sci, recursive = TRUE, showWarnings = FALSE)

# Save as TIFF
out_file_tiff_temp <- file.path(out_dir_ROC_AUC_sci, "roc_by_all_year.tiff")

ggplot2::ggsave(
  filename   = out_file_tiff_temp,
  plot       = p,
  device     = ragg::agg_tiff,
  width      = 5.7, # 1-column width
  height     = 5.7,
  units      = "cm",
  dpi        = 600,
  compression = "lzw",
  background  = "white"
)

# Save as SVG 
out_file_svg_temp <- file.path(out_dir_ROC_AUC_sci, "roc_by_all_year.svg")

ggplot2::ggsave(
  filename   = out_file_svg_temp,
  plot       = p,
  device     = "svg", # or svglite::svglite
  width      = 5.7,
  height     = 5.7,
  units      = "cm",
  bg = "white"
)

message("Saved TIFF and SVG to: ", out_dir_ROC_AUC_sci)

# TEMPORAL FINISHED




################################################################################
## spatial CV
################################################################################


set.seed(42)
n_folds_spat <- 5L # folds of Cross validation (~ 80/20 split each)

regions <- sort(unique(data_spde__$region_id)) # # Unique region IDs

regions_shuf <- sample(regions, length(regions), replace = FALSE)

fold_id <- rep(seq_len(n_folds_spat), length.out = length(regions_shuf))

fold_map <- tibble(region_id = regions_shuf, fold = fold_id)
fold_map

# Quick check for number of regions in each fold
table(fold_map$fold)


# Run Spatial CV
cv_results_spatial <- vector("list", n_folds_spat)

for (k in seq_len(n_folds_spat)) {
  # k =1
  
  # Regions held out in this fold
  regions_test_k <- fold_map %>% filter(fold == k) %>% pull(region_id)
  
  message(sprintf(">>> Spatial CV fold %d/%d — holding out %d regions",
                  k, n_folds_spat, length(regions_test_k)))
  
  # Build a lookup table 
  lkp_region <- data_spde__ %>%
    dplyr::select(region_id, geo) %>%
    dplyr::distinct()
  
  # Join polygons to region_id via 'geo'
  regions_sf <- map_spde %>% # joining the data_spde__ to map
    dplyr::left_join(lkp_region, by = "geo")
  
  region_status <- tibble::tibble(region_id = unique(na.omit(regions_sf$region_id))) %>%
    dplyr::mutate(set = dplyr::if_else(region_id %in% regions_test_k,
                                       "Test (held-out)", "Train"))
  
  plot_sf <- regions_sf %>%
    dplyr::left_join(region_status, by = "region_id") %>%
    dplyr::mutate(set = dplyr::if_else(is.na(set), "Train", set)) %>%
    sf::st_as_sf(crs = sf::st_crs(regions_sf))
  
  # plotting the train test map for each fold
  # Simple legend + clean map (no grids), no titles
  p_map <- ggplot2::ggplot(plot_sf) +
    ggplot2::geom_sf(
      ggplot2::aes(fill = set, color = set),
      linewidth = 0.4
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::scale_fill_manual(
      name = NULL, # no legend title
      values = c("Train" = "grey50", "Test (held-out)" = "white")
    ) +
    ggplot2::scale_color_manual(
      name = NULL, # no legend title
      values = c("Train" = "grey20", "Test (held-out)" = "red")
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "right",
                   legend.title = element_blank(),
                   legend.key = ggplot2::element_blank(),
                   # Legend sizes
                   legend.text       = element_text(size = 16, colour = "black"))
  # Legend.position = "right")
  
  # Save TIFF (8x7 in, 600 dpi, transparent)
  ext         <- "tiff"
  width_in    <- 8
  height_in   <- 7
  dpi         <- 600
  compression <- "lzw"
  bg          <- "transparent"
  
  out_dir_map_cv <- ".../Plots/A_final_plots/Cross_validation/Hold_out_regions_spatial"
  dir.create(out_dir_map_cv, recursive = TRUE, showWarnings = FALSE)
  
  out_file_map <- file.path(out_dir_map_cv, sprintf("spatial_cv_fold_%02d_map.%s", k, ext))
  
  ggplot2::ggsave(
    filename    = out_file_map,
    plot        = p_map,
    device      = ragg::agg_tiff,
    width       = width_in,
    height      = height_in,
    units       = "in",
    dpi         = dpi,
    compression = compression,
    background  = bg
  )
  
  message("Saved map: ", out_file_map)
  
  # Saving the plot with standard style guidelines
  p_map <- ggplot2::ggplot(plot_sf) +
    ggplot2::geom_sf(
      ggplot2::aes(fill = set, color = set),
      linewidth = 0.08 # Thinner borders 
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c("Train" = "grey50", "Test (held-out)" = "white")
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = c("Train" = "grey20", "Test (held-out)" = "red")
    ) +
    ggplot2::theme_void() + # Start with void
    ggplot2::theme(
      text = element_text(family = "Arial"), # Required sans-serif font
      legend.position = "right",
      legend.title = element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.text = element_text(
        size = 6,
        margin = margin(l = 0.9, unit = "pt") # # sppce between legend symbol and legend text
      ),
      legend.key.size = unit(0.2, "cm"),
      # REMOVE GAP BETWEEN KEY AND TEXT
      # Set spacing to a very small or zero value
      legend.spacing.x = unit(0.01, "cm"),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      legend.margin = margin(l = 2, r = 0, unit = "pt") # Tiny gap for legend only
    )
  
  # 2. Saving Setup
  out_dir_base <- ".../Plots/A_final_plots/Cross_validation/Hold_out_regions_spatial"
  out_dir_sci  <- file.path(out_dir_base, "Science_journal")
  dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)
  
  width_cm  <- 5.7
  height_cm <- 3.0
  dpi       <- 600
  
  # 3. Save as TIFF (High-res production raster)
  tiff_path <- file.path(out_dir_sci, sprintf("spatial_cv_fold_%02d_map.tiff", k))
  ggplot2::ggsave(
    filename    = tiff_path,
    plot        = p_map,
    device      = ragg::agg_tiff,
    width       = width_cm,
    height      = height_cm,
    units       = "cm",
    dpi         = dpi,
    compression = "lzw",
    background  = "white"
  )
  
  # 4. Save as SVG 
  svg_path <- file.path(out_dir_sci, sprintf("spatial_cv_fold_%02d_map.svg", k))
  ggplot2::ggsave(
    filename = svg_path,
    plot     = p_map,
    device   = "svg",
    width    = width_cm,
    height   = height_cm,
    units    = "cm",
    bg       = "white"
  )
  
  message("Saved Science-ready TIFF and SVG to: ", out_dir_sci)
  
  # Tag train/test rows across all years based on region membership in each fold
  data_spde_train_test <- data_spde__ %>%
    mutate(
      # See the details in temporal validation
      training_set = as.integer(!region_id %in% regions_test_k & !is.na(.data[[response_column]])),
      test_set     = as.integer( region_id %in% regions_test_k  & !is.na(.data[[response_column]])),
      y_orig       = .data[[response_column]] # keep the untouched original response
    )
  
  data_spde_train_test[[response_column]][data_spde_train_test$test_set == 1L] <- NA_real_
  
  
  # A descriptive model tag 
  model_str <- sprintf("cv_spatial_fold%02d_%s", k, model_str__)
  
  out <- run_spde_model2(
    parameters      = parameters,
    data_spde       = data_spde_train_test,
    map_spde        = map_spde,
    n_spde          = n_spde,
    n_years         = n_years_Q,
    n_years_pred    = n_years_pred,
    f_list          = f_list,
    response_column = response_column,
    Q_spde          = Q_spde,
    strategy        = strategy,
    n_binom         = n_binom,
    model_str       = model_str,
    pred_data       = pred_data
  )
  
  fitted_model <- out$fitted_model
  summary(fitted_model)
  
  
  A <- out$A_train
  A_pred <- out$A_pred 
  idx_linear_predictor <- out$idx_linear_predictor
  spde_indices <- out$spde_indices_train
  
  p_all <- fitted_model$summary.fitted.values$mean[idx_linear_predictor] #  probabilities of all regions fitted in INLA model
  
  idx_test <- which(data_spde_train_test$test_set == 1L)
  y_true   <- data_spde_train_test$y_orig[idx_test]
  p_hat    <- p_all[idx_test]
  
  # Metrics at threshold 0.5
  m <- bin_metrics(y_true, p_hat, thr = 0.5) %>%
    mutate(fold = k)
  
  # ROC / AUC
  roc_obj <- tryCatch(
    pROC::roc(response = y_true, predictor = p_hat, quiet = TRUE),
    error = function(e) NULL
  )
  auc_val <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
  
  # Confusion matrix per fold
  fig_width_px  <- 3000
  fig_height_px <- 2800
  fig_dpi       <- 600
  thr_std       <- 0.5
  
  out_dir_conf_sp <- ".../Plots/A_final_plots/Cross_validation/Confusion_matrix_spatial"
  dir.create(out_dir_conf_sp, recursive = TRUE, showWarnings = FALSE)
  
  fig_name <- sprintf("confus_matrix_spatial_fold%02d_threshold%0.2f.png", k, thr_std)
  fig_path <- file.path(out_dir_conf_sp, fig_name)
  
  png(filename = fig_path, width = fig_width_px, height = fig_height_px, res = fig_dpi)
  draw_confusion_matrix_present(y_true, p_hat, thr = thr_std, main = "")
  dev.off()
  
  # Store results for this fold
  cv_results_spatial[[k]] <- list(
    fold             = k,
    held_out_regions = regions_test_k, # keep which regions were held out
    metrics          = m %>% mutate(AUC = auc_val),
    roc              = roc_obj,
    n_test           = length(idx_test)
  )
}

# Save spatial CV results object
out_dir_cv_sp <- ".../Output/A-Final_output/Cross_Validation_Output_spatial"
dir.create(out_dir_cv_sp, recursive = TRUE, showWarnings = FALSE)

save(cv_results_spatial,
     file = file.path(out_dir_cv_sp, "cv_results_spatial_region5fold.RData"))

# Now loading the CV metrics stored
out <- load(paste0(out_dir_cv_sp, "/cv_results_spatial_region5fold.RData"))
cv_results_spatial <- cv_results_spatial

# Building tidy metrics table 
cv_metrics_spatial <- dplyr::bind_rows(lapply(cv_results_spatial, `[[`, "metrics")) %>%
  arrange(fold)

cv_metrics_spatial <- cv_metrics_spatial %>%
  dplyr::relocate(fold) %>% # move fold to first
  dplyr::relocate(dplyr::any_of(c("TP","TN","FP","FN")),
                  .after = dplyr::last_col()) # move counts to the end

readr::write_csv(cv_metrics_spatial,
                 file.path(out_dir_cv_sp, "cv_metrics_spatial_region5fold.csv"))

# Overall mean across folds (ignore NA)
cv_metrics_spatial_overall <- cv_metrics_spatial %>%
  summarise(across(c(accuracy, sensitivity, specificity, precision, F1, AUC),
                   ~ mean(.x, na.rm = TRUE)))

readr::write_csv(cv_metrics_spatial_overall,
                 file.path(out_dir_cv_sp, "cv_metrics_spatial_region5fold_overall.csv"))

# ROC curves by spatial fold (for all folds together)
roc_df_sp <- purrr::map_dfr(cv_results_spatial, function(x) {
  if (is.null(x$roc)) return(NULL)
  tibble(
    fold = x$fold,
    tpr  = rev(x$roc$sensitivities),
    fpr  = rev(1 - x$roc$specificities)
  )
})

if (nrow(roc_df_sp) > 0) {
  p_sp <- ggplot2::ggplot(
    roc_df_sp,
    ggplot2::aes(x = fpr, y = tpr, color = factor(fold), group = fold)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::labs(
      x = "False Positive rate",
      y = "True Positve rate",
      color = "Spatial fold",
      title = "ROC - AUC"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.title.x = element_text(color = "black", size = 20),
      axis.title.y = element_text(color = "black", size = 20),
      axis.line  = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.text  = ggplot2::element_text(color = "black"),
      panel.background = ggplot2::element_blank(),
      plot.background  = ggplot2::element_blank(),
      # Legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      # Legend sizes
      legend.title      = element_text(size = 18, face = "bold", colour = "black"),
      legend.text       = element_text(size = 16, colour = "black")
    )
  
  print(p_sp)
} else {
  message("No ROC data to plot (spatial).")
}

# Save ROC figure (spatial)
ext         <- "tiff" # "tiff", "png", "pdf",
width_in    <- 8
height_in   <- 7
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

out_dir_ROC_sp <- ".../Plots/A_final_plots/Cross_validation/ROC-AUC-Spatial"
dir.create(out_dir_ROC_sp, recursive = TRUE, showWarnings = FALSE)

out_file_sp <- file.path(out_dir_ROC_sp, paste0("roc_spatial_region5fold.", ext))

# Install.packages("ragg") 
ggplot2::ggsave(
  filename    = out_file_sp,
  plot        = p_sp,
  device      = ragg::agg_tiff, # use ragg::agg_png if ext == "png"
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  background  = bg
)

message("Saved: ", out_file_sp)

## Final plot description
if (nrow(roc_df_sp) > 0) {
  p <- ggplot2::ggplot(
    roc_df_sp,
    ggplot2::aes(x = fpr, y = tpr, color = factor(fold), group = fold)
  ) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.2) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::labs(
      x = "False positive rate",
      y = "True positive rate",
      color = "Spatial fold" # Shortened title to save space
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = element_text(family = "Arial"),
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6),
      axis.line  = element_line(linewidth = 0.2), # and axis.line also to 0.2 standard
      axis.ticks = element_line(linewidth = 0.2),
      
      # COMPACT LEGEND SETTING
      legend.title = element_text(
        size = 6,
        face = "bold",
        margin = margin(b = 2.5, unit = "pt")
      ),
      legend.text  = element_text(size = 5), # # legend text
      # Shrink the colored lines in the legend
      legend.key.size = unit(0.15, "cm"),
      legend.key.height = unit(0.1, "cm"), #reduce the gap between legend
      # Legend position
      legend.position = c(0.90, 0.38),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    ) +
    # Force legend into one column to keep it narrow
    ggplot2::guides(color = guide_legend(ncol = 1))
  
  print(p)
}

out_dir_ROC_sp_sci <- ".../Plots/A_final_plots/Cross_validation/ROC-AUC-Spatial/Science_journal"
dir.create(out_dir_ROC_sp_sci, recursive = TRUE, showWarnings = FALSE)

# Save as TIFF 
out_file_tiff_sp_sci <- file.path(out_dir_ROC_sp_sci, "roc_spatial_region5fold.tiff")

ggplot2::ggsave(
  filename   = out_file_tiff_sp_sci,
  plot       = p,
  device     = ragg::agg_tiff,
  width      = 5.7, #
  height     = 5.7,
  units      = "cm",
  dpi        = 600,
  compression = "lzw",
  background  = "white"
)

# Save as SVG (Vector format for maximum quality and edit
out_file_svg_sci <- file.path(out_dir_ROC_sp_sci, "roc_spatial_region5fold.svg")

ggplot2::ggsave(
  filename   = out_file_svg_sci,
  plot       = p,
  device     = "svg", # or svglite::svglite
  width      = 5.7,
  height     = 5.7,
  units      = "cm",
  bg = "white"
)

message("Saved TIFF and SVG to: ", out_dir_ROC_sp_sci)

# End spatial CV




################################################################################
## Spatially-confined  CV 
################################################################################

#STEPS- 
## 1 -- Contiguous folds built from NUTS3 adjacency (sf::st_touches) (compact contiguous cluster (blob) per fold (no scattered clusters) as far as possible with limitation that each region was used in testing only once)
## 2 -- Optional buffer ring around held-out block (default 35 km with trivial case as 0 km)
## 3 -- Holds out full regions across *all years*
## 4 -- Saves per-fold maps, confusion matrices, held-out region lists,
## 5--  and CV results objects + metrics CSV.
###############################################################################

# Settings
set.seed(42) # reproducibility
n_folds_spat <- 5L # ~80/20 each fold (5)
buffer_km    <- 35


#assign unique ID to each region (starting from 1)
lkp_region <- data_spde__ %>%
  dplyr::select(region_id, geo) %>% # # picking up all distinct regions
  distinct()

# Joining the unique id column from lkp_region with map with geo as joint indentifier
regions_sf <- map_spde %>%
  left_join(lkp_region, by = "geo") %>%
  st_as_sf() %>%
  filter(!st_is_empty(geometry))

# Checking that region i
stopifnot("region_id" %in% names(regions_sf))

N <- nrow(regions_sf) # # numbers of regions


# Build adjacency (contiguity graph)
nb_list <- sf::st_touches(regions_sf)

# Working in a metric CRS for compactness (for seed picking and compact growth we need distances; LAEA Europe (EPSG:3035) )
regions_sf_etrs <- sf::st_transform(regions_sf, 3035)

# Centroid of each region (x,y) in meters; matrix with N rows, 2 columns
cent_xy <- sf::st_coordinates(sf::st_centroid(regions_sf_etrs))[, 1:2, drop = FALSE]


# Seed pickers
# anchor each blob in a distinct “macro-area”.
pick_corners_plus_center <- function(coords, k) {
  
  # K =5L ## for example 5 fold
  
  if (k != 5L) return(NULL) # only defined for 5-way split
  x <- coords[,1]; y <- coords[,2] # # picking up coordinates
  
  # Scale function just scales each x coordinate in LAEA (EPSG:3035) across the whole set of coordinates
  idx_SW <- which.min(scale(x) + scale(y))
  idx_NE <- which.max(scale(x) + scale(y)) # northeast-ish (large x and large y)
  idx_NW <- which.max(scale(-x) + scale(y)) # northwest-ish (small x, large y)
  idx_SE <- which.max(scale(x) + scale(-y)) # southeast-ish (large x, small y)
  
  # Center: nearest to median point
  cx <- median(x)
  cy <- median(y)
  
  d2 <- (x - cx)^2 + (y - cy)^2
  idx_C  <- which.min(d2)
  seeds <- c(idx_SW, idx_NW, idx_NE, idx_SE, idx_C) # # # ensure uniqueness
  
  seeds <- unique(seeds)
  if (length(seeds) < 5) {
    # Farthest-point fill-in
    n <- nrow(C)
    while (length(seeds) < 5) {
      # Distance to nearest chosen seed
      d2 <- rep(Inf, n)
      for (s in seeds) d2 <- pmin(d2, rowSums((C - matrix(C[s,], nrow=n, ncol=2, byrow=TRUE))^2))
      d2[seeds] <- -Inf
      seeds <- c(seeds, which.max(d2))
    }
  }
  seeds
}


#Farthest-point seeding (fallback for k != 5 or corner failure) 
pick_seeds_farthest_point <- function(coords, k) {
  
  # K =5L ## for example 5 fold
  
  n <- nrow(coords) # # total number of coordinates or regions
  first <- sample.int(n, 1) # # sample any region from it (random sample)
  seeds <- integer(k)
  seeds[1] <- first # # store the first seed's row index
  for (i in 2:k) {
    # I = 2
    d2 <- rep(Inf, n)
    
    #The inner loop computes squared distances from every region to seed 650 (650vis example, any id can be used for initialization), then takes a pmin with the current d2
    for (s in seeds[1:(i-1)]) { # # loop over the indices of already-chosen seeds
      d2 <- pmin(d2, rowSums((coords - matrix(coords[s,], nrow=n, ncol=2, byrow=TRUE))^2))
    }
    d2[seeds[1:(i-1)]] <- -Inf      # avoid re-picking existing seeds ((d2[650] <- -Inf, which prevents us from re-selecting region 650)
    seeds[i] <- which.max(d2)       # pick the region that is farthest from ALL existing seeds
    
    # Now for i = 3 the existing seeds are {650, r2}:
    #   Reinitialize d2 <- rep(Inf, 1335).
    # Loop over s ∈ {650, r2}:
    # compute squared distances to 650 → update d2 with pmin.
    # 
    # compute squared distances to r2 → update d2 again with pmin.
    # 
    # Now d2[j] equals the squared distance from j to its nearest seed among {650, r2}.
    # Set d2[c(650, r2)] <- -Inf.
    # Pick seeds[3] <- which.max(d2) — the region farthest from both seeds {650, r2} (i.e., farthest from its nearest seed).
  }
  # Repeat this until we have k seeds:
  seeds                         # the seeds will be at farthest distance from each other (initialized by a region picked by a random sample)
}


# Single-blob fold builder
# Build k folds; each fold is one contiguous blob, grown  from its seed.
build_single_blob_folds <- function(nb_list, coords, k, target_size = NULL, seeds = NULL) {
  
  # Target_size = NULL
  # #seeds = NULL
  
  n <- nrow(coords)
  if (is.null(target_size)) target_size <- ceiling(n / k)       ## average number of regions in each fold,  if no particular target number of fold is given 
  
  if (is.null(seeds)) seeds <- pick_corners_plus_center(coords, k)
  
  if (is.null(seeds) || length(seeds) < k) {
    seeds <- pick_seeds_farthest_point(coords, k)
  }
  
  # Assignment vectors
  assign <- rep(NA_integer_, n) # fold id for each node
  taken  <- rep(FALSE, n) # TRUE once a node is assigned
  
  # Pre compute squared distances from each seed to all nodes 
  dist2_seed <- lapply(seeds, function(si) {
    rowSums((coords - matrix(coords[si,], nrow=n, ncol=2, byrow=TRUE))^2)
  })
  
  for (i in seq_len(k)) {
    # I = 1
    si <- seeds[i] # # select each farthest region (or corner region)
    
    if (taken[si]) {
      candidates <- which(!taken)
      si <- candidates[which.min(dist2_seed[[i]][candidates])]
    }
    
    # Initialize blob with the seed
    assign[si] <- i
    taken[si]  <- TRUE      #also mark taken as TRUE
    
    frontier <- setdiff(nb_list[[si]], si)
    frontier <- frontier[!taken[frontier]] #check that this ID of neighbours must not be taken
    
    # Grow until the blob reaches target_size or frontier is exhausted
    while (sum(assign == i, na.rm = TRUE) < target_size && length(frontier) > 0) {     #The loop condition : Keep growing while: The blob size (how many regions already assigned to fold i) is less than target_size, and there are still frontier candidates to add
      
      # Choose among frontier the node closest to this fold's seed (compactness)
      pick_idx <- frontier[which.min(dist2_seed[[i]][frontier])]
      assign[pick_idx] <- i
      taken[pick_idx]  <- TRUE
      
      # Update frontier: add neighbors of the picked node; remove taken nodes
      new_neigh <- nb_list[[pick_idx]] # Get all neighbors of the region we just added
      frontier  <- unique(c(setdiff(frontier, pick_idx), new_neigh))
      frontier  <- frontier[!taken[frontier]]
    }
  }
  
  # leftover regions which are not picked in any fold
  leftovers <- which(!taken)                # rarely needed (blob growth typically covers the graph fairly)
  if (length(leftovers) > 0) {
    for (v in leftovers) {
      dists <- sapply(seq_len(k), function(i) dist2_seed[[i]][v])
      fi <- which.min(dists)
      assign[v] <- fi
      taken[v]  <- TRUE
    }
  }
  
  # Return mapping to region_id 
  tibble(
    region_row = seq_len(nrow(regions_sf)),
    region_id  = regions_sf$region_id,
    fold       = assign
  )
}

# Build the fold map 
fold_map_idx <- build_single_blob_folds(nb_list, cent_xy, n_folds_spat)

fold_map <- dplyr::select(fold_map_idx, region_id, fold)
print(table(fold_map$fold)) # check for near-balanced fold sizes



# Buffer helper to prevent spatial leakage (optional)
to_etrs   <- function(x) sf::st_transform(x, 3035)         # conversion to LAEA system for distance computation between regions 
from_etrs <- function(x, crs_target) sf::st_transform(x, crs_target)   # from LAEA to original CRS

make_buffer_ids <- function(test_ids, regions_sf, buffer_km = 0) {
  # test_ids <- regions_test_k 
  # regions_sf <- regions_sf   
  # buffer_km <- 35         
  
  if (buffer_km <= 0) return(integer(0))
  crs0 <- sf::st_crs(regions_sf)        # assigning the standard WGS84 crs system 
  
  test_sf <- regions_sf %>% dplyr::filter(region_id %in% test_ids)
  
  # Buffer regions
  buf <- test_sf %>%
    to_etrs() %>%         # conversion to ETRS89-LAEA (EPSG:3035) for distance computation
    
    sf::st_buffer(dist = buffer_km * 1000) %>%      # meters (convert the distance to meter as we define buffer distance in km)
    
    from_etrs(crs0)         #bring the buffered polygons back to original CRS 
  
  hits <- sf::st_intersects(regions_sf, sf::st_union(buf), sparse = TRUE)
  
  buf_ids <- regions_sf$region_id[lengths(hits) > 0] # # Pick regions from that intersect the buffer
  
  setdiff(buf_ids, test_ids) # buffer excludes the test set itself (test regions are already held-out, the buffers are meant to remove only nearby train regions)
}

# Run the spatially-confined CV
cv_results_spatial_conf <- vector("list", n_folds_spat)

buffer_km <- 35


# Run spatially-confined CV
cv_results_spatial <- vector("list", n_folds_spat)

for (k in seq_len(n_folds_spat)) { # # running through each fold
  # k = 1
  
  regions_test_k <- fold_map %>% filter(fold == k) %>% pull(region_id)
  regions_buff_k <- make_buffer_ids(regions_test_k, regions_sf, buffer_km)
  
  message(sprintf(">>> Spatial CV fold %d/%d — Test: %d regions  Buffer: %d regions",
                  k, n_folds_spat, length(regions_test_k), length(regions_buff_k)))
  
  #  Train/Buffer/Test map
  region_status <- tibble(region_id = regions_sf$region_id) %>%
    mutate(set = case_when(
      region_id %in% regions_test_k ~ "Test (held-out)",
      region_id %in% regions_buff_k ~ "Buffer zone (excluded)",
      TRUE                          ~ "Train"
    ))
  
  plot_sf <- regions_sf %>%
    left_join(region_status, by = "region_id") %>%
    mutate(set = factor(set, levels = c("Train","Buffer zone (excluded)","Test (held-out)")))
  
  p_map <- ggplot(plot_sf) +
    geom_sf(aes(fill = set, color = set), linewidth = 0.4) +
    coord_sf(expand = FALSE) +
    scale_fill_manual(
      name = NULL,
      values = c("Train" = "grey50", "Buffer zone (excluded)" = "grey80", "Test (held-out)" = "#ffefdd")
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Train" = "grey20", "Buffer zone (excluded)" = "blue", "Test (held-out)" = "red")
    ) +
    theme_void(base_size = 12) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.key = ggplot2::element_blank(),
          # Legend sizes
          legend.text       = element_text(size = 16, colour = "black"))
  # Legend.position = "right")
  
  # Save the fold map as TIFF (8x7", 600 dpi, transparent)
  out_dir_map_cv_conf_spat <- ".../Plots/A_final_plots/Cross_validation/Hold_out_region_confined_sptial"
  dir.create(out_dir_map_cv_conf_spat, recursive = TRUE, showWarnings = FALSE)
  
  ggplot2::ggsave(
    filename    = file.path(out_dir_map_cv_conf_spat, sprintf("spatial_cv_fold_%02d_map_conf.tiff", k)),
    plot        = p_map,
    device      = ragg::agg_tiff, # use ragg for crisp text/lines
    width       = 8, height = 7, units = "in", dpi = 600,
    compression = "lzw", background = "transparent"
  )
  
  # Save lists of held-out and buffer regions (helps in reproducibility) with map attributes
  cols_prefer <- c("region_id","geo","LocationNa","LocationPa","LocationTy",
                   "CountryNam","CountryISO","isEU","isEEA","isWHOEurop")
  cols_keep   <- intersect(cols_prefer, names(regions_sf))
  
  # Save a tabular metadata (no geometry) for each fold
  # Heldout test regions
  heldout_meta <- regions_sf %>%
    dplyr::filter(region_id %in% regions_test_k) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(cols_keep)) %>%
    dplyr::mutate(fold = k, set = "Test (held-out)")
  
  # Heldout buffer region
  buffer_meta <- regions_sf %>%
    dplyr::filter(region_id %in% regions_buff_k) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(cols_keep)) %>%
    dplyr::mutate(fold = k, set = "Buffer (excluded)")
  
  out_dir_cv_conf_spat_reg    <- ".../Output/A-Final_output/Cross_Validation_Output_confined_spatial/Fold's_test_and_buffer_region"
  dir.create(out_dir_cv_conf_spat_reg,   recursive = TRUE, showWarnings = FALSE)
  
  # Write CSVs (held-out, buffer, and combined)
  readr::write_csv(heldout_meta,
                   file.path(out_dir_cv_conf_spat_reg, sprintf("heldout_region_ids_fold%02d_with_map.csv", k)))
  
  if (nrow(buffer_meta) > 0) {
    readr::write_csv(buffer_meta,
                     file.path(out_dir_cv_conf_spat_reg, sprintf("buffer_region_ids_fold%02d_with_map.csv", k)))
  }
  
  # Saving the plot clean and as per standaed requirements
  p_map <- ggplot2::ggplot(plot_sf) +
    ggplot2::geom_sf(
      ggplot2::aes(fill = set, color = set),
      linewidth = 0.08 # Thinner borders 
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    scale_fill_manual(
      name = NULL,
      values = c("Train" = "grey50", "Buffer zone (excluded)" = "grey80", "Test (held-out)" = "#ffefdd")
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Train" = "grey20", "Buffer zone (excluded)" = "blue", "Test (held-out)" = "red")
    ) +
    ggplot2::theme_void() + 
    ggplot2::theme(
      text = element_text(family = "Arial"), # Required sans-serif font
      legend.position = "right",
      legend.title = element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.text = element_text(
        size = 6,
        margin = margin(l = 0.9, unit = "pt") # # space between legend symbol and legend text
      ),
      legend.key.size = unit(0.2, "cm"),
      # REMOVE GAP BETWEEN KEY AND TEXT
      # Set spacing to a very small or zero value
      legend.spacing.x = unit(0.01, "cm"),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      legend.margin = margin(l = 2, r = 0, unit = "pt") # Tiny gap for legend only
    )
  
  # 2. Saving Setup
  out_dir_map_cv_conf_spat_sci <- ".../Plots/A_final_plots/Cross_validation/Hold_out_region_confined_sptial"
  out_dir_sci_conf_spat_sci  <- file.path(out_dir_map_cv_conf_spat_sci, "Science_journal")
  dir.create(out_dir_sci_conf_spat_sci, recursive = TRUE, showWarnings = FALSE)
  
  width_cm  <- 5.7
  height_cm <- 3.0
  dpi       <- 600
  
  # 3. Save as TIFF 
  tiff_path_sci <- file.path(out_dir_sci_conf_spat_sci, sprintf("spatial_cv_fold_%02d_map_conf.tiff", k))
  ggplot2::ggsave(
    filename    = tiff_path_sci,
    plot        = p_map,
    device      = ragg::agg_tiff,
    width       = width_cm,
    height      = height_cm,
    units       = "cm",
    dpi         = dpi,
    compression = "lzw",
    background  = "white"
  )
  
  # 4. Save as SVG
  svg_path_sci <- file.path(out_dir_sci_conf_spat_sci, sprintf("spatial_cv_fold_%02d_map_conf.svg", k))
  ggplot2::ggsave(
    filename = svg_path_sci,
    plot     = p_map,
    device   = "svg",
    width    = width_cm,
    height   = height_cm,
    units    = "cm",
    bg       = "white"
  )
  
  message("Saved Science-ready TIFF and SVG to: ", out_dir_sci)
  
  # Build train/test masks across *all years* of panel
  data_spde_train_test <- data_spde__ %>%
    dplyr::mutate(
      in_test   = region_id %in% regions_test_k,
      in_buffer = region_id %in% regions_buff_k,
      # Training_set = 1 only if: in neither test nor buffer AND response is non-NA
      #held out the same training region for all years
      training_set = as.integer(!in_test & !in_buffer & !is.na(.data[[response_column]])),
      test_set     = as.integer( in_test               & !is.na(.data[[response_column]])),
      y_orig       = .data[[response_column]]        # keep a copy of the original response
    )
  
  # Mask the *response* for test rows 
  data_spde_train_test[[response_column]][data_spde_train_test$test_set == 1L] <- NA_real_
  
  model_str <- sprintf("cv_spatial_fold%02d_%s", k, model_str__)
  
  out <- run_spde_model2(
    parameters      = parameters,
    data_spde       = data_spde_train_test,
    map_spde        = map_spde,
    n_spde          = n_spde,
    n_years         = n_years_Q,
    n_years_pred    = n_years_pred,
    f_list          = f_list,
    response_column = response_column,
    Q_spde          = Q_spde,
    strategy        = strategy,
    n_binom         = n_binom,
    model_str       = model_str,
    pred_data       = pred_data
  )
  
  fitted_model <- out$fitted_model
  summary(fitted_model)
  
  
  A <- out$A_train
  A_pred <- out$A_pred
  idx_linear_predictor <- out$idx_linear_predictor
  spde_indices <- out$spde_indices_train
  
  p_all <- fitted_model$summary.fitted.values$mean[idx_linear_predictor]
  
  # Extract held-out rows and their truth/predictions
  idx_test <- which(data_spde_train_test$test_set == 1L)
  y_true   <- data_spde_train_test$y_orig[idx_test]
  p_hat    <- p_all[idx_test]
  
  #  Metrics, ROC, Confusion matrix plot
  m <- bin_metrics(y_true, p_hat, thr = 0.5) %>% dplyr::mutate(fold = k)
  
  roc_obj <- tryCatch(pROC::roc(response = y_true, predictor = p_hat, quiet = TRUE),
                      error = function(e) NULL)
  auc_val <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
  
  out_dir_map_cv_conf_spat_confus <- ".../Plots/A_final_plots/Cross_validation/Confusion_matrix_confined_spatial"
  dir.create(out_dir_map_cv_conf_spat_confus, recursive = TRUE, showWarnings = FALSE)
  
  fig_name <- sprintf("confus_matrix_spatial_fold%02d_threshold%0.2f.png", k, 0.5)
  png(filename = file.path(out_dir_map_cv_conf_spat_confus, fig_name),
      width = 3000, height = 2800, res = 600)
  draw_confusion_matrix_present(y_true, p_hat, thr = 0.5, main = "")
  dev.off()
  
  # Stash fold results
  cv_results_spatial_conf[[k]] <- list(
    fold             = k,
    held_out_regions = regions_test_k,
    buffer_regions   = regions_buff_k,
    metrics          = m %>% dplyr::mutate(AUC = auc_val),
    roc              = roc_obj,
    n_test           = length(idx_test)
  )
}

# Save confined spatial results & metrics table
out_dir_cv_spat_conf_data    <- ".../Output/A-Final_output/Cross_Validation_Output_confined_spatial"

save(cv_results_spatial_conf,
     file = file.path(out_dir_cv_spat_conf_data , "cv_results_confined_spatial_region5fold.RData"))

message("Done. Built confined folds, ran CV, and saved artifacts.")

# Now loading the CV metrics stored
out <- load(paste0(out_dir_cv_spat_conf_data, "/cv_results_confined_spatial_region5fold.RData"))
cv_results_spatial_conf <- cv_results_spatial_conf

cv_metrics_spatial_conf <- dplyr::bind_rows(lapply(cv_results_spatial_conf, `[[`, "metrics")) %>%
  arrange(fold)

cv_metrics_spatial_conf <- cv_metrics_spatial_conf %>%
  dplyr::relocate(fold) %>% # move fold to first
  dplyr::relocate(dplyr::any_of(c("TP","TN","FP","FN")),
                  .after = dplyr::last_col()) # move counts to the end

readr::write_csv(cv_metrics_spatial_conf,
                 file.path(out_dir_cv_spat_conf_data, "cv_metrics_confined_spatial_region5fold.csv"))

# Overall mean across folds (all 5) (ignore NA)
cv_metrics_spatial_overall_conf <- cv_metrics_spatial_conf %>%
  summarise(across(c(accuracy, sensitivity, specificity, precision, F1, AUC),
                   ~ mean(.x, na.rm = TRUE)))

readr::write_csv(cv_metrics_spatial_overall_conf,
                 file.path(out_dir_cv_spat_conf_data, "cv_metrics_confined_spatial_region5fold_overall.csv"))

# ROC curves by spatial fold (for all folds together)
roc_df_sp_conf <- purrr::map_dfr(cv_results_spatial_conf, function(x) {
  if (is.null(x$roc)) return(NULL)
  tibble(
    fold = x$fold,
    tpr  = rev(x$roc$sensitivities),
    fpr  = rev(1 - x$roc$specificities)
  )
})


# Filter Data: Select only Folds 1, 3, and 5
# EXCLUSION REASONING:
# - Fold 2 is excluded: It includes the UK and Ireland, which contains mostly absence data, resulting in very imbalanced validation metrics (see supplementary for detials)
# - Fold 4 is excluded: It includes the South of Italy, which contains mostly presence data, also resulting in very imbalanced validation metrics.
########################################

roc_df_filtered <- roc_df_sp_conf[roc_df_sp_conf$fold %in% c(1, 3, 5), ]

# Plotting

if (nrow(roc_df_filtered) > 0) {
  p_sp <- ggplot2::ggplot(
    roc_df_filtered,
    ggplot2::aes(x = fpr, y = tpr, color = factor(fold), group = fold)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::labs(
      x = "False Positive rate",
      y = "True Positive rate", # Fixed typo 'Positve' to 'Positive'
      color = "Spatial fold",
      title = "ROC - AUC"
    ) +
    # Clean look: dark axes, ticks; no grid
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(color = "black", size = 20),
      axis.title.y = ggplot2::element_text(color = "black", size = 20),
      axis.line  = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.text  = ggplot2::element_text(color = "black"),
      panel.background = ggplot2::element_blank(),
      plot.background  = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      # Legend sizes
      legend.title       = ggplot2::element_text(size = 18, face = "bold", colour = "black"),
      legend.text        = ggplot2::element_text(size = 16, colour = "black")
    )
  
  print(p_sp)
  
} else {
  message("No ROC data to plot (spatial) after filtering.")
}

# Save ROC figure 
ext         <- "tiff" # "tiff", "png", "pdf",
width_in    <- 8
height_in   <- 7
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

out_dir_ROC_sp_con <- ".../Plots/A_final_plots/Cross_validation/ROC-AUC-Confined_spatial"
dir.create(out_dir_ROC_sp_con, recursive = TRUE, showWarnings = FALSE)

out_file_sp <- file.path(out_dir_ROC_sp_con, paste0("roc_confined_spatial_region5fold.", ext))


ggplot2::ggsave(
  filename    = out_file_sp,
  plot        = p_sp,
  device      = ragg::agg_tiff, # use ragg::agg_png if ext == "png"
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  background  = bg
)

message("Saved: ", out_file_sp)


## saving for clean standard figure
if (nrow(roc_df_filtered) > 0) {
  p_sp <- ggplot2::ggplot(
    roc_df_filtered,
    ggplot2::aes(x = fpr, y = tpr, color = factor(fold), group = fold)
  ) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.2) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::labs(
      x = "False positive rate",
      y = "True positive rate",
      color = "Spatial fold" # Shortened title to save space
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = element_text(family = "Arial"),
      
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6),
      axis.line  = element_line(linewidth = 0.2), # and axis.line also to 0.2 standard
      axis.ticks = element_line(linewidth = 0.2),
      
      # COMPACT LEGEND SETTING
      legend.title = element_text(
        size = 6,
        face = "bold",
        margin = margin(b = 2.5, unit = "pt")
      ),
      legend.text  = element_text(size = 5), # # legend text
      # Shrink the colored lines in the legend
      legend.key.size = unit(0.15, "cm"),
      legend.key.height = unit(0.1, "cm"), # # reduce the gap between legend
      # Legend position
      legend.position = c(0.90, 0.38),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    ) +
    # Force legend into one column to keep it narrow
    ggplot2::guides(color = guide_legend(ncol = 1))
  
  print(p_sp)
} else {
  message("No ROC data to plot (spatial) after filtering.")
}


out_dir_ROC_sp_con_sci <- ".../Plots/A_final_plots/Cross_validation/ROC-AUC-Confined_spatial/Science_journal"
dir.create(out_dir_ROC_sp_con_sci, recursive = TRUE, showWarnings = FALSE)

# Save as TIFF 
out_file_tiff_conf <- file.path(out_dir_ROC_sp_con_sci, "roc_confined_spatial_region5fold.tiff")

ggplot2::ggsave(
  filename   = out_file_tiff_conf,
  plot       = p_sp,
  device     = ragg::agg_tiff,
  width      = 5.7, 
  height     = 5.7,
  units      = "cm",
  dpi        = 600,
  compression = "lzw",
  background  = "white"
)

# Save as SVG 
out_file_svg_conf <- file.path(out_dir_ROC_sp_con_sci, "roc_confined_spatial_region5fold.svg")

ggplot2::ggsave(
  filename   = out_file_svg_conf,
  plot       = p_sp,
  device     = "svg", # or svglite::svglite
  width      = 5.7,
  height     = 5.7,
  units      = "cm",
  bg = "white"
)

message("Saved TIFF and SVG to: ", out_dir_ROC_sp_con_sci)





################################################################################
## spatio-temporal CV  (random 80% - 20% train test picking)
################################################################################

set.seed(42) # reproducible fold assignment

n_all <- nrow(data_spde__)
eligible <- which(!is.na(data_spde__[[response_column]]))

n_folds <- 5L
fold_id_vec <- sample(rep(seq_len(n_folds), length.out = length(eligible)))
fold_map_df <- tibble(row_id = eligible, fold = fold_id_vec)

table(fold_map_df$fold) # size of each fold

# vector for holding results
cv_results_st <- vector("list", n_folds)

# Loop over folds
for (k in seq_len(n_folds)) {
  
  # K = 1
  
  message(sprintf(">>> Spatio-temporal CV (random) — fold %d/%d", k, n_folds))
  
  # Indices of rows assigned to this fold (test set)
  idx_test_rows <- fold_map_df$row_id[fold_map_df$fold == k] 
  
  # Tag train/test & keep y_orig for evaluation
  data_spde_train_test <- data_spde__ %>%
    mutate(
      training_set = as.integer(!(row_number() %in% idx_test_rows) &
                                  !is.na(.data[[response_column]])),
      test_set     = as.integer( (row_number() %in% idx_test_rows) &
                                   !is.na(.data[[response_column]]) ),
      y_orig       = .data[[response_column]]
    )
  
  data_spde_train_test[[response_column]][data_spde_train_test$test_set == 1L] <- NA_real_
  
  # Label for logs/artifacts
  model_str <- sprintf("cv_st_rand5fold_fold%02d_%s", k, model_str__)
  
  # Fit model 
  out <- run_spde_model2(
    parameters      = parameters,
    data_spde       = data_spde_train_test,
    map_spde        = map_spde,
    n_spde          = n_spde,
    n_years         = n_years_Q,
    n_years_pred    = n_years_pred,
    f_list          = f_list,
    response_column = response_column,
    Q_spde          = Q_spde,
    strategy        = strategy,
    n_binom         = n_binom,
    model_str       = model_str,
    pred_data       = pred_data
  )
  
  fitted_model         <- out$fitted_model
  summary(out$fitted_model)
  idx_linear_predictor <- out$idx_linear_predictor
  
  p_all <- fitted_model$summary.fitted.values$mean[idx_linear_predictor]
  
  # Evaluation on test rows only
  idx_test <- which(data_spde_train_test$test_set == 1L)
  y_true   <- data_spde_train_test$y_orig[idx_test]
  p_hat    <- p_all[idx_test]
  
  # Metrics @ 0.5 threshold
  m <- bin_metrics(y_true, p_hat, thr = 0.5) %>%
    mutate(fold = k)
  
  # ROC / AUC
  roc_obj <- tryCatch(pROC::roc(response = y_true, predictor = p_hat, quiet = TRUE),
                      error = function(e) NULL)
  auc_val <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
  
  
  # Directory for confusion matrix plot
  out_dir_conf_st <- ".../Plots/A_final_plots/Cross_validation/Confusion_matrix_spatio-temporal"
  dir.create(out_dir_conf_st, recursive = TRUE, showWarnings = FALSE)
  
  fig_width_px  <- 3000
  fig_height_px <- 2800
  fig_dpi       <- 600
  thr_std       <- 0.5
  
  fig_path <- file.path(out_dir_conf_st,
                        sprintf("confus_matrix_spatiotemporal_fold%02d_threshold%0.2f.png",
                                k, thr_std))
  png(filename = fig_path, width = fig_width_px, height = fig_height_px, res = fig_dpi)
  draw_confusion_matrix_present(y_true, p_hat, thr = thr_std, main = "")
  dev.off()
  
  # Store this fold’s result
  cv_results_st[[k]] <- list(
    fold    = k,
    metrics = m %>% mutate(AUC = auc_val),
    roc     = roc_obj,
    n_test  = length(idx_test)
  )
}

# Saving results
out_dir_cv_st <- ".../Output/A-Final_output/Cross_validation_spatio-temporal"
dir.create(out_dir_cv_st, recursive = TRUE, showWarnings = FALSE)

save(cv_results_st,
     file = file.path(out_dir_cv_st, "cv_results_spatiotemporal_5fold.RData"))

# Now loading the CV metrics spatio-temporal stored
out <- load(paste0(out_dir_cv_st , "/cv_results_spatiotemporal_5fold.RData"))
cv_results_st<- cv_results_st

# Per-fold metrics table + overall mean
cv_metrics_st <- dplyr::bind_rows(lapply(cv_results_st, `[[`, "metrics")) %>%
  arrange(fold) %>%
  dplyr::relocate(fold) %>%
  dplyr::relocate(dplyr::any_of(c("TP","TN","FP","FN")), .after = dplyr::last_col())

readr::write_csv(cv_metrics_st,
                 file.path(out_dir_cv_st, "cv_metrics_spatiotemporal_k5.csv"))

cv_metrics_st_overall <- cv_metrics_st %>%
  summarise(across(c(accuracy, sensitivity, specificity, precision, F1, AUC),
                   ~ mean(.x, na.rm = TRUE)))

readr::write_csv(cv_metrics_st_overall,
                 file.path(out_dir_cv_st, "cv_metrics_spatiotemporal_k5_overall.csv"))

# Save ROC curves across folds (combined figure)
roc_df_st <- purrr::map_dfr(cv_results_st, function(x) {
  if (is.null(x$roc)) return(NULL)
  tibble(
    fold = x$fold,
    tpr  = rev(x$roc$sensitivities),
    fpr  = rev(1 - x$roc$specificities)
  )
})

out_dir_roc_st <- ".../Plots/A_final_plots/Cross_validation/ROC-AUC-Spatio-temporal"
dir.create(out_dir_roc_st, recursive = TRUE, showWarnings = FALSE)

if (nrow(roc_df_st) > 0) {
  p_st <- ggplot2::ggplot(
    roc_df_st,
    ggplot2::aes(x = fpr, y = tpr, color = factor(fold), group = fold)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::labs(x = "False Positive rate", y = "True Positive rate", color = "Fold") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.title.x = element_text(color = "black", size = 20),
      axis.title.y = element_text(color = "black", size = 20),
      axis.line    = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.ticks   = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.text    = ggplot2::element_text(color = "black"),
      panel.background = ggplot2::element_blank(),
      plot.background  = ggplot2::element_blank(),
      legend.key       = ggplot2::element_blank(),
      legend.title     = element_text(size = 18, face = "bold", colour = "black"),
      legend.text      = element_text(size = 16, colour = "black")
    )
  
  ggplot2::ggsave(
    filename    = file.path(out_dir_roc_st, "roc_spatiotemporal_k5.tiff"),
    plot        = p_st,
    device      = ragg::agg_tiff,
    width       = 8, height = 7, units = "in", dpi = 600,
    compression = "lzw", background = "transparent"
  )
} else {
  message("No ROC data to plot (spatio-temporal k-fold).")
}

message("Spatio-temporal (random 5-fold, ~80/20) CV — done.")


## plot with standard guidelines
if (nrow(roc_df_st) > 0) {
  p_st <- ggplot2::ggplot(
    roc_df_st,
    ggplot2::aes(x = fpr, y = tpr, color = factor(fold), group = fold)
  ) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.2) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::labs(
      x = "False positive rate",
      y = "True positive rate",
      color = "Fold" # Shortened title to save space
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = element_text(family = "Arial"),
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6),
      axis.line  = element_line(linewidth = 0.2), # and axis.line also to 0.2 standard
      axis.ticks = element_line(linewidth = 0.2),
      
      # COMPACT LEGEND SETTING
      legend.title = element_text(
        size = 6,
        face = "bold",
        margin = margin(b = 2.5, unit = "pt")
      ),
      legend.text  = element_text(size = 5), # # legend text
      # Shrink the colored lines in the legend
      legend.key.size = unit(0.15, "cm"),
      legend.key.height = unit(0.1, "cm"), # # reduce the gap between legend
      # Legend position
      legend.position = c(0.90, 0.38),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    ) +
    # Force legend into one column to keep it narrow
    ggplot2::guides(color = guide_legend(ncol = 1))
  
  print(p_st)
} else {
  message("No ROC data to plot (spatial) after filtering.")
}

# For main model
out_dir_roc_st_sci <- ".../Plots/A_final_plots/Cross_validation/ROC-AUC-Spatio-temporal/Science_journal"
dir.create(out_dir_roc_st_sci, recursive = TRUE, showWarnings = FALSE)

# Save as TIFF 
out_file_tiff_st_sci <- file.path(out_dir_roc_st_sci, "roc_spatiotemporal_k5.tiff")

ggplot2::ggsave(
  filename   = out_file_tiff_st_sci,
  plot       = p_st,
  device     = ragg::agg_tiff,
  width      = 5.7, 
  height     = 5.7,
  units      = "cm",
  dpi        = 600,
  compression = "lzw",
  background  = "white"
)

# Save as SVG 
out_file_svg_st_sci <- file.path(out_dir_roc_st_sci, "roc_spatiotemporal_k5.svg")

ggplot2::ggsave(
  filename   = out_file_svg_st_sci,
  plot       = p_st,
  device     = "svg", # or svglite::svglite
  width      = 5.7,
  height     = 5.7,
  units      = "cm",
  bg = "white"
)

message("Saved TIFF and SVG to: ", out_dir_roc_st_sci)


# CV completed





#############################################################
####### PLOTS FOR VALIDATION METRIC ############
#############################################################

# Plots of accuracy, sensitivity, specificity, precision,
# FIRST TEMPORAL
metrics <- c("accuracy","precision","sensitivity","specificity","F1", "AUC")

# Load the file containing temporal validation metrics
cv_metrics_by_year <- read.csv(".../Output/A-Final_output/Cross_Validation_Output_temporal/cv_metrics_temporal_holdoutyears_001-013.csv")
print(cv_metrics_by_year)

# Converting this to long format
temporal_long <- cv_metrics_by_year %>% 
  dplyr::select(year, all_of(metrics)) %>%
  pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("accuracy","precision","sensitivity","specificity","F1", "AUC"),
                         labels = c("Accuracy","Precision","Sensitivity","Specificity","F1", "AUC")))

p_temporal_all <- ggplot(temporal_long, aes(x = year, y = Value, color = Metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_y_continuous(limits = c(0.6,1), breaks = seq(0,1,0.1)) +
  scale_x_continuous(breaks = pretty(temporal_long$year)) +
  labs(x = "Year", y = "Value", color = NULL) + # no legend title
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20),
    axis.line  = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.text  = ggplot2::element_text(color = "black"),
    panel.background = ggplot2::element_blank(),
    plot.background  = ggplot2::element_blank(),
    # Legend.title = ggplot2::element_text(face = "bold"),
    legend.key = ggplot2::element_blank(),
    # Legend sizes
    legend.title      = element_text(size = 18, face = "bold", colour = "black"),
    legend.text       = element_text(size = 16, colour = "black"),
    legend.position = "right"
  )

p_temporal_all

## SPATIAL
cv_metrics_by_spatial <- read.csv(".../Output/A-Final_output/Cross_Validation_Output_spatial/cv_metrics_spatial_region5fold.csv")
print(cv_metrics_by_spatial)

spatial_long <- cv_metrics_by_spatial %>%
  dplyr::select(fold, all_of(metrics)) %>%
  pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("accuracy","precision","sensitivity","specificity","F1", "AUC"),
                         labels = c("Accuracy","Precision","Sensitivity","Specificity","F1", "AUC")
  ))

p_spatial_all <- ggplot(spatial_long, aes(x = fold, y = Value, color = Metric, group = Metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = sort(unique(spatial_long$fold))) +
  scale_y_continuous(limits = c(0.6,1), breaks = seq(0,1,0.1)) +
  labs(x = "Spatial Fold", y = "Value", color = NULL) + # no legend title
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20),
    axis.line  = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.text  = ggplot2::element_text(color = "black"),
    panel.background = ggplot2::element_blank(),
    plot.background  = ggplot2::element_blank(),
    # Legend.title = ggplot2::element_text(face = "bold"),
    legend.key = ggplot2::element_blank(),
    # Legend sizes
    legend.title      = element_text(size = 18, face = "bold", colour = "black"),
    legend.text       = element_text(size = 16, colour = "black"),
    legend.position = "right"
  )

print(p_spatial_all)


## CONFINED SPATIAL
cv_metrics_by_confined_spatial <- read.csv(".../Output/A-Final_output/Cross_Validation_Output_confined_spatial/cv_metrics_confined_spatial_region5fold.csv")
print(cv_metrics_by_confined_spatial)

# Filter Data: Select only Folds 1, 3, and 5 (see supplementary)
cv_metrics_by_confined_spatial_filtered <- cv_metrics_by_confined_spatial[cv_metrics_by_confined_spatial$fold %in% c(1, 3, 5), ]

spatial_conf_long <- cv_metrics_by_confined_spatial_filtered %>%
  dplyr::select(fold, all_of(metrics)) %>%
  pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("accuracy","precision","sensitivity","specificity","F1", "AUC"),
                         labels = c("Accuracy","Precision","Sensitivity","Specificity","F1", "AUC")
  ))

p_spatial_conf_all <- ggplot(spatial_conf_long, aes(x = fold, y = Value, color = Metric, group = Metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = sort(unique(spatial_conf_long$fold))) +
  scale_y_continuous(limits = c(0.5,1), breaks = seq(0,1,0.1)) +
  labs(x = "Spatial Fold", y = "Value", color = NULL) + # no legend title
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20),
    axis.line  = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.text  = ggplot2::element_text(color = "black"),
    panel.background = ggplot2::element_blank(),
    plot.background  = ggplot2::element_blank(),
    # Legend.title = ggplot2::element_text(face = "bold"),
    legend.key = ggplot2::element_blank(),
    # Legend sizes
    legend.title      = element_text(size = 18, face = "bold", colour = "black"),
    legend.text       = element_text(size = 16, colour = "black"),
    legend.position = "right"
  )

print(p_spatial_conf_all)


## SPATIO-TEMPORAL
cv_metrics_by_spatio_temp <- read.csv(".../Output/A-Final_output/Cross_validation_spatio-temporal/cv_metrics_spatiotemporal_k5.csv")
print(cv_metrics_by_spatio_temp)

spatio_temp_long <- cv_metrics_by_spatio_temp %>%
  dplyr::select(fold, all_of(metrics)) %>%
  pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("accuracy","precision","sensitivity","specificity","F1", "AUC"),
                         labels = c("Accuracy","Precision","Sensitivity","Specificity","F1", "AUC")
  ))

p_spatio_temp_all <- ggplot(spatio_temp_long, aes(x = fold, y = Value, color = Metric, group = Metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = sort(unique(spatial_long$fold))) +
  scale_y_continuous(limits = c(0.6,1), breaks = seq(0,1,0.1)) +
  labs(x = "Spatial Fold", y = "Value", color = NULL) + # no legend title
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20),
    axis.line  = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.text  = ggplot2::element_text(color = "black"),
    panel.background = ggplot2::element_blank(),
    plot.background  = ggplot2::element_blank(),
    # Legend.title = ggplot2::element_text(face = "bold"),
    legend.key = ggplot2::element_blank(),
    # Legend sizes
    legend.title      = element_text(size = 18, face = "bold", colour = "black"),
    legend.text       = element_text(size = 16, colour = "black"),
    legend.position = "right"
  )

print(p_spatio_temp_all)



# Save combined-metrics plots (spatial & temporal)
ext         <- "tiff" # "tiff", "png", "pdf",
width_in    <- 8
height_in   <- 7
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

out_dir_metrics <- ".../Plots/A_final_plots/Cross_validation/Metrics_Combined_plot"
dir.create(out_dir_metrics, recursive = TRUE, showWarnings = FALSE)

# Filenames
file_spatial   <- file.path(out_dir_metrics, paste0("spatial.",  ext))
file_temporal  <- file.path(out_dir_metrics, paste0("temporal.", ext))
file_conf_spatial  <- file.path(out_dir_metrics, paste0("confined_spatial.", ext))
file_conf_spatio_temp  <- file.path(out_dir_metrics, paste0("spatio-temporal.", ext))


# Save all 4 validation combined metrics in separate plot
ggplot2::ggsave(
  filename    = file_spatial,
  plot        = p_spatial_all,
  device      = ragg::agg_tiff,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  background  = bg
)

ggplot2::ggsave(
  filename    = file_temporal,
  plot        = p_temporal_all,
  device      = ragg::agg_tiff,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  background  = bg
)

ggplot2::ggsave(
  filename    = file_conf_spatial,
  plot        = p_spatial_conf_all,
  device      = ragg::agg_tiff,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  background  = bg
)

ggplot2::ggsave(
  filename    = file_conf_spatio_temp,
  plot        = p_spatio_temp_all,
  device      = ragg::agg_tiff,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  background  = bg
)

message("Saved:\n- ", file_spatial, "\n- ", file_temporal, "\n- ", file_conf_spatial, "\n- ", file_conf_spatio_temp)




# plotting with standard style
metrics <- c("accuracy","precision","sensitivity","specificity","F1", "AUC")

theme_science_small <- function(legend_pos = c(0.98, 0.30)) {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text = element_text(family = "Arial"),
      
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6),
      
      axis.line  = element_line(linewidth = 0.2), # and axis.line also to 0.2 standard
      axis.ticks = element_line(linewidth = 0.2),
      
      # Compact legend settings 
      legend.title = element_text(
        size = 6,
        face = "bold",
        margin = margin(b = 2.5, unit = "pt")
      ),
      legend.text = element_text(size = 4.5, margin = margin(l = 1.2, unit = "pt")),
      
      # Shrink the colored lines in the legend
      legend.key.size   = unit(0.15, "cm"),
      legend.key.height = unit(0.10, "cm"), # # reduce the gap between legend
      
      legend.position = legend_pos, # # legend position
      legend.justification = c(1, 0),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      
      panel.grid = element_blank()
    )
}

## fuction - useful for saving the plots
save_science_tiff_svg <- function(p, out_dir, stem,
                                  w_cm = 5.7, h_cm = 5.7, dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  f_tiff <- file.path(out_dir, paste0(stem, ".tiff"))
  f_svg  <- file.path(out_dir, paste0(stem,  ".svg"))
  
  ggplot2::ggsave(
    filename    = f_tiff,
    plot        = p,
    device      = ragg::agg_tiff,
    width       = w_cm,
    height      = h_cm,
    units       = "cm",
    dpi         = dpi,
    compression = "lzw",
    background  = "white"
  )
  
  ggplot2::ggsave(
    filename  = f_svg,
    plot      = p,
    device    = "svg",
    width     = w_cm,
    height    = h_cm,
    units     = "cm",
    bg        = "white"
  )
  
  message("Saved: ", f_tiff, " and ", f_svg)
}



# TEMPORAL metrics plot
cv_metrics_by_year <- read.csv(
  ".../Output/A-Final_output/Cross_Validation_Output_temporal/cv_metrics_temporal_holdoutyears_001-013.csv"
)

temporal_long <- cv_metrics_by_year %>%
  dplyr::select(year, all_of(metrics)) %>%
  tidyr::pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  dplyr::mutate(
    Metric = factor(
      Metric,
      levels = c("accuracy","precision","sensitivity","specificity","F1","AUC"),
      labels = c("Accuracy","Precision","Sensitivity","Specificity","F1","AUC")
    )
  )

p_temporal_sci <- ggplot2::ggplot(
  temporal_long,
  ggplot2::aes(x = year, y = Value, color = Metric, group = Metric)
) +
  ggplot2::geom_line(linewidth = 0.3) +
  ggplot2::geom_point(size = 0.3) + # # metrics point size
  ggplot2::scale_y_continuous(
    limits = c(0.6, 1),
    breaks = seq(0.6, 1.0, 0.1),
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5)
  ) +
  ggplot2::labs(
    x = "Year",
    y = "Value",
    color = NULL
  ) +
  theme_science_small(legend_pos = c(1.07, 0.80)) + # # top right
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1))

print(p_temporal_sci)



out_dir_metrics_sci <- ".../Plots/A_final_plots/Cross_validation/Metrics_Combined_plot/Science_journal"
save_science_tiff_svg(p_temporal_sci, out_dir_metrics_sci, "metrics_temporal")


# SPATIAL metrics plot
cv_metrics_by_spatial <- read.csv(
  ".../Output/A-Final_output/Cross_Validation_Output_spatial/cv_metrics_spatial_region5fold.csv"
)

spatial_long <- cv_metrics_by_spatial %>%
  dplyr::select(fold, all_of(metrics)) %>%
  tidyr::pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  dplyr::mutate(
    Metric = factor(
      Metric,
      levels = c("accuracy","precision","sensitivity","specificity","F1","AUC"),
      labels = c("Accuracy","Precision","Sensitivity","Specificity","F1","AUC")
    )
  )

p_spatial_sci <- ggplot2::ggplot(
  spatial_long,
  ggplot2::aes(x = fold, y = Value, color = Metric, group = Metric)
) +
  ggplot2::geom_line(linewidth = 0.3) +
  ggplot2::geom_point(size = 0.3) +
  ggplot2::scale_x_continuous(breaks = sort(unique(spatial_long$fold))) +
  ggplot2::scale_y_continuous(
    limits = c(0.6, 1),
    breaks = seq(0.6, 1.0, 0.1)
  ) +
  ggplot2::labs(
    x = "Spatial fold",
    y = "Value",
    color = NULL
  ) +
  theme_science_small(legend_pos = c(1.07, 0.05)) + # ## bottom right
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1))

print(p_spatial_sci)

out_dir_metrics_sci <- ".../Plots/A_final_plots/Cross_validation/Metrics_Combined_plot/Science_journal"

save_science_tiff_svg(p_spatial_sci, out_dir_metrics_sci, "metrics_spatial")



# CONFINED SPATIAL metrics plot
cv_metrics_by_confined_spatial <- read.csv(
  ".../Output/A-Final_output/Cross_Validation_Output_confined_spatial/cv_metrics_confined_spatial_region5fold.csv"
)

cv_metrics_by_confined_spatial_filtered <-
  cv_metrics_by_confined_spatial[cv_metrics_by_confined_spatial$fold %in% c(1, 3, 5), ]

spatial_conf_long <- cv_metrics_by_confined_spatial_filtered %>%
  dplyr::select(fold, all_of(metrics)) %>%
  tidyr::pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  dplyr::mutate(
    Metric = factor(
      Metric,
      levels = c("accuracy","precision","sensitivity","specificity","F1","AUC"),
      labels = c("Accuracy","Precision","Sensitivity","Specificity","F1","AUC")
    )
  )

p_spatial_conf_sci <- ggplot2::ggplot(
  spatial_conf_long,
  ggplot2::aes(x = fold, y = Value, color = Metric, group = Metric)
) +
  ggplot2::geom_line(linewidth = 0.3) +
  ggplot2::geom_point(size = 0.3) +
  ggplot2::scale_x_continuous(breaks = sort(unique(spatial_conf_long$fold))) +
  ggplot2::scale_y_continuous(
    limits = c(0.5, 1),
    breaks = seq(0.5, 1.0, 0.1)
  ) +
  ggplot2::labs(
    x = "Spatial fold",
    y = "Value",
    color = NULL
  ) +
  theme_science_small(legend_pos = c(1.07, 0.80)) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1))

print(p_spatial_conf_sci)

# Saving the file
out_dir_metrics_sci <- ".../Plots/A_final_plots/Cross_validation/Metrics_Combined_plot/Science_journal"

save_science_tiff_svg(p_spatial_conf_sci, out_dir_metrics_sci, "metrics_confined_spatial")

# SPATIO-TEMPORAL metrics plot
cv_metrics_by_spatio_temp <- read.csv(
  ".../Output/A-Final_output/Cross_validation_spatio-temporal/cv_metrics_spatiotemporal_k5.csv"
)

spatio_temp_long <- cv_metrics_by_spatio_temp %>%
  dplyr::select(fold, all_of(metrics)) %>%
  tidyr::pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  dplyr::mutate(
    Metric = factor(
      Metric,
      levels = c("accuracy","precision","sensitivity","specificity","F1","AUC"),
      labels = c("Accuracy","Precision","Sensitivity","Specificity","F1","AUC")
    )
  )

p_spatio_temp_sci <- ggplot2::ggplot(
  spatio_temp_long,
  ggplot2::aes(x = fold, y = Value, color = Metric, group = Metric)
) +
  ggplot2::geom_line(linewidth = 0.3) +
  ggplot2::geom_point(size = 0.3) +
  ggplot2::scale_x_continuous(breaks = sort(unique(spatio_temp_long$fold))) + # (fixed)
  ggplot2::scale_y_continuous(
    limits = c(0.6, 1),
    breaks = seq(0.6, 1.0, 0.1)
  ) +
  ggplot2::labs(
    x = "Spatio-temporal fold",
    y = "Value",
    color = NULL
  ) +
  theme_science_small(legend_pos = c(1.07, 0.05)) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1))

print(p_spatio_temp_sci)

# SAVE 
out_dir_metrics_sci <- ".../Plots/A_final_plots/Cross_validation/Metrics_Combined_plot/Science_journal"
save_science_tiff_svg(p_spatio_temp_sci, out_dir_metrics_sci, "metrics_spatiotemporal_k5")

message("All Science-style metric plots saved to: ", out_dir_metrics_sci)

