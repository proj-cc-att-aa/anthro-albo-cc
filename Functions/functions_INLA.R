# ================== Inla model and output functions ==================

packages <- c("spdep")

lapply(packages, library, character.only = TRUE)

library("INLA")


# Function: logit.
logit <- function(p) log(p / (1-p))


# Function: inverse_logit.
inverse_logit <- function(x) 1 / (1 + exp(-x))


#  update_f (for preparation of standard INLA model structure)
update_f <- function(f, str_to_add) {
  f_updated <- f %>% as.character()
  
  f_updated[[3]] <- paste0(f_updated[[3]], " + ", str_to_add)
  
  f_updated <- sprintf("%s ~ %s",
                       as.character(f_updated[[2]]),
                       f_updated[[3]]) %>%
    as.formula()
}

# define_model- model formula
define_model <- function(
    response = "presence",
    f_basic  = "1",
    
    precipitation_winter_mean     = FALSE,
    precipitation_spring_mean     = FALSE,
    precipitation_summer_mean     = FALSE,
    precipitation_autumn_mean     = FALSE,
    
    temperature_winter_mean       = FALSE,
    temperature_spring_mean       = FALSE,
    temperature_summer_mean       = FALSE,
    temperature_autumn_mean       = FALSE,
    
    relative_humidity_winter_mean = FALSE,
    relative_humidity_spring_mean = FALSE,
    relative_humidity_summer_mean = FALSE,
    relative_humidity_autumn_mean = FALSE,
    
    population_density = FALSE,
    
    spatial_idx_col = NULL,
    spatial_model   = "",
    proximity_col   = "",
    fixed_effect_cols = "",
    proximity_type     = ""
) {
  
  extraconstr <- list()
  
  f <- sprintf("%s ~ %s", response, f_basic) %>% as.formula()
  
  if (!identical(fixed_effect_cols, "")) {
    for (fixed_effect in fixed_effect_cols) {
      f <- update_f(f, fixed_effect)
    }
  }
  
  if (precipitation_winter_mean) {
    precipitation_winter_mean_str <- paste0(
      "f(bin_idx_precipitation_winter_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, precipitation_winter_mean_str)
  }
  if (precipitation_spring_mean) {
    precipitation_spring_mean_str <- paste0(
      "f(bin_idx_precipitation_spring_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, precipitation_spring_mean_str)
  }
  if (precipitation_summer_mean) {
    precipitation_summer_mean_str <- paste0(
      "f(bin_idx_precipitation_summer_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, precipitation_summer_mean_str)
  }
  if (precipitation_autumn_mean) {
    precipitation_autumn_mean_str <- paste0(
      "f(bin_idx_precipitation_autumn_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, precipitation_autumn_mean_str)
  }
  
  if (temperature_winter_mean) {
    temperature_winter_mean_str <- paste0(
      "f(bin_idx_temperature_winter_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, temperature_winter_mean_str)
  }
  if (temperature_spring_mean) {
    temperature_spring_mean_str <- paste0(
      "f(bin_idx_temperature_spring_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, temperature_spring_mean_str)
  }
  if (temperature_summer_mean) {
    temperature_summer_mean_str <- paste0(
      "f(bin_idx_temperature_summer_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, temperature_summer_mean_str)
  }
  if (temperature_autumn_mean) {
    temperature_autumn_mean_str <- paste0(
      "f(bin_idx_temperature_autumn_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, temperature_autumn_mean_str)
  }
  
  if (relative_humidity_winter_mean) {
    relative_humidity_winter_mean_str <- paste0(
      "f(bin_idx_relative_humidity_winter_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, relative_humidity_winter_mean_str)
  }
  if (relative_humidity_spring_mean) {
    relative_humidity_spring_mean_str <- paste0(
      "f(bin_idx_relative_humidity_spring_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, relative_humidity_spring_mean_str)
  }
  if (relative_humidity_summer_mean) {
    relative_humidity_summer_mean_str <- paste0(
      "f(bin_idx_relative_humidity_summer_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, relative_humidity_summer_mean_str)
  }
  if (relative_humidity_autumn_mean) {
    relative_humidity_autumn_mean_str <- paste0(
      "f(bin_idx_relative_humidity_autumn_mean, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, relative_humidity_autumn_mean_str)
  }
  
  if (population_density) {
    population_density_str <- paste0(
      "f(bin_idx_population_density, model=\"rw1\", scale.model=TRUE, ",
      "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))))"
    )
    f <- update_f(f, population_density_str)
  }
  
  if (!is.null(spatial_idx_col)) {
    
    if (!spatial_model %in% c("iid", "besagproper"))
      stop("spatial_model error")
    
    if (spatial_idx_col=="spatial_idx") {
      
      spatial_model_str <- paste0("f(spatial_idx, model = \"", spatial_model,
                                  "\", graph = file_ginla, constr = FALSE)")
    } else {
      
      stop("not implemented")
    }
    
    f <- update_f(f, spatial_model_str)
  }
  
  if (proximity_col != "") {
    
    if (identical(proximity_type, "rw1")) {
      
      proximity_str <- paste0(
        "f(", proximity_col, ", model=\"rw1\", scale.model=TRUE, ",
        "hyper = list(theta = list(prior=\"pc.prec\",param = c(0.3, 0.05))), ",
        "constr = TRUE)"
      )
      
    } else if (proximity_type %in% c("ICAR", "besag")) {
      
      proximity_str <- paste0(
        "f(", proximity_col, ", model=\"besag\", ",
        "graph=file_nb_proximity, scale.model=TRUE, ",
        "hyper = list(theta = list(prior=\"pc.prec\", param=c(0.5,0.01))), ",
        "constr = TRUE)"
      )
      
    } else {
      
      stop("proximity_type must be either 'rw1' or 'ICAR'/'besag'")
    }
    
    f <- update_f(f, proximity_str)
  }
  
  return(list(f=f, extraconstr=extraconstr))
}


#  prepare_inla_input.
prepare_inla_input <- function(map_inla, tb_inla, file_ginla=NULL) {
  
  if (!all(tb_inla$geo %in% map_inla$geo)) {
    stop("some regions in tb_inla are missing in map_inla")
  }
  
  geo_neighs_inla <- geo_to_neighs(map_inla)
  nb2INLA(file_ginla, geo_neighs_inla$nb_map)
  
  map_inla$spatial_idx <- 1:nrow(map_inla)
  
  data_inla <- full_join(tb_inla %>% dplyr::select(-any_of("spatial_idx")), map_inla %>% dplyr::select(geo, spatial_idx)) %>%
    arrange(year_idx, spatial_idx) %>%
    relocate(geo, year_idx, spatial_idx)
  
  n_years <- length(unique(data_inla$year_idx))
  if (!identical(data_inla$spatial_idx, rep(1:nrow(map_inla), times=n_years)))
    stop("data_inla should have one row for each combination (year_idx, spatial_idx)")
  
  if (!identical(map_inla$spatial_idx, 1:nrow(map_inla)))
    stop("map_inla should be sorted by spatial_idx")
  
  data_inla$data_idx <- 1:nrow(data_inla)
  
  data_inla <- data_inla %>% relocate(data_idx, spatial_idx, year_idx)
  map_inla <- map_inla %>% relocate(spatial_idx)
  
  return(list(map_inla=map_inla, data_inla = data_inla))
}



# link_for_NA_handling.
link_for_NA_handling <- function(df, y_column="presence") {
  NA_handling <- rep(NA, nrow(df))
  ind_NA <- which(is.na(df[[y_column]]))
  NA_handling[ind_NA] = 1
  
  return(NA_handling)
}




# extract_main_results_training.
extract_main_results_training <- function(model, data, remove_NA=FALSE,
                                          y_column="presence", likelihood="bernoulli",
                                          n_binom=NULL, n_time=1,
                                          spatial_random_effect="spatial_idx",
                                          random_effect_with_time_pattern="",
                                          idx_linear_predictor=NULL,
                                          data_pred=NULL,
                                          idx_linear_predictor_pred=NULL,
                                          cross_val=FALSE) {
  
  if (!likelihood %in% c("bernoulli", "binomial"))
    stop("only implemented for bernoulli and binomial likelihoods")
  
  if (likelihood=="binomial" & !is.numeric(n_binom)) {
    
    if (!is.numeric(n_binom))
      stop("n_binom must be supplied for the binomial likelihood")
    
    if (!(x>0 && identical(round(x), x)))
      stop("n_binom should be a positive integer")
  }
  
  if (spatial_random_effect=="spde_idx") {
    if (is.null(idx_linear_predictor)) {
      stop("idx_linear_predictor needs to be suplied for an spde model")
    } else if (length(idx_linear_predictor) != nrow(data)) {
      stop("idx_linear_predictor does not match data")
    }
  }
  
  if (is.null(idx_linear_predictor))
    idx_linear_predictor <- 1:nrow(data)
  
  intercept_str <- rownames(model$summary.fixed) %>% str_subset("Intercept")
  tau_str <- rownames(model$summary.hyperpar) %>% str_subset("spatial_idx|spde_idx") %>%
    str_subset("Precision")
  
  random_effect_col_names <- names(model$summary.random)
  spatial_effect <- intersect(random_effect_col_names,
                              c("spatial_idx", "spde_idx"))
  
  random_effect_col_names <- setdiff(random_effect_col_names, spatial_effect)
  
  random_effects_with_time <- str_subset(random_effect_col_names,
                                         random_effect_with_time_pattern)
  
  random_effects_no_time <- setdiff(random_effect_col_names, random_effects_with_time)
  
  results_summary <- tibble(DIC = round(model$dic$dic,2),
                            n_eff_parameters.DIC = model$dic$p.eff,
                            logscore = round(-mean(log(model$cpo$cpo), na.rm = T), 3),
                            WAIC = round(model$waic$waic, 2),
                            n_eff_parameters.WAIC = model$waic$p.eff,
                            fixed_effects = list(model$summary.fixed),
                            hyperpar = list(model$summary.hyperpar),
                            intercept = model$summary.fixed[intercept_str,"mean"],
                            hyperpar_tau = ifelse(length(model$summary.hyperpar[tau_str, "mean"])==1,
                                                  model$summary.hyperpar[tau_str, "mean"], NA),
                            fixed_effect_col_names =
                              list(model$model.matrix@Dimnames[[2]] %>% setdiff(intercept_str)),
                            random_effect_col_names = list(list(no_time=random_effects_no_time,
                                                                with_time=random_effects_with_time)),
                            spatial_random_effect = list(spatial_effect))
  
  ## Hyperparameters
  rnames <- rownames(results_summary$hyperpar[[1]])
  results_summary$hyperpar_processed <- list(list())
  
  i <- str_detect(rnames, "Precision")
  if (any(i)) {
    h <- results_summary$hyperpar[[1]][which(i),]
    
    v1 <- as.matrix(t(h$mean))
    col_str <- rownames(h) %>% str_extract("(?<=Precision for ).*") %>% paste0(".pr.mean")
    colnames(v1) <- col_str
    
    v2 <- as.matrix(t(h$mean))
    col_str <- rownames(h) %>% str_extract("(?<=Precision for ).*") %>% paste0(".pr.m05")
    colnames(v2) <- col_str
    
    results_summary$hyperpar_processed[[1]]$precision <- as_tibble(cbind(v2, v1))
  }
  
  i <- str_detect(rnames, "Theta")
  if (any(i)) {
    h <- results_summary$hyperpar[[1]][which(i),]
    
    v1 <- as.matrix(t(h$mean))
    col_str <- rownames(h) %>% str_extract("Theta.") %>% paste0(".meanlog")
    colnames(v1) <- col_str
    
    v2 <- as.matrix(t(h$mean))
    col_str <- rownames(h) %>% str_extract("Theta.") %>% paste0(".m05log")
    colnames(v2) <- col_str
    
    v3 <- exp(v2)
    col_str <- rownames(h) %>% str_extract("Theta.") %>% paste0(".m05exp")
    colnames(v3) <- col_str
    
    results_summary$hyperpar_processed[[1]]$theta <- as_tibble(cbind(v3, v2, v1))
    
    names(results_summary$hyperpar_processed[[1]]$theta) <-
      names(results_summary$hyperpar_processed[[1]]$theta) %>%
      str_replace("Theta1", "tau") %>% str_replace("Theta2", "kappa")
  }
  
  n_fixed_effects <- length(results_summary$fixed_effect_col_names[[1]])
  if (n_fixed_effects > 0) {
    map_fixed_effects <- tibble(data_idx = data$data_idx,
                                geo = data$geo,
                                region_id = data$region_id,
                                year_idx = data$year_idx)
    
    for (j in 1:n_fixed_effects) {
      fixed_effect_name <- results_summary$fixed_effect_col_names[[1]][j]
      fitted_coef <- results_summary$fixed_effects[[1]][fixed_effect_name, "mean"]
      
      map_fixed_effects[[paste0("feffect_", fixed_effect_name, ".mean")]] <-
        data[[fixed_effect_name]] * fitted_coef
    }
  } else {
    map_fixed_effects <- NULL
  }
  
  n_random_effects_no_time <- length(results_summary$random_effect_col_names[[1]]$no_time)
  if (n_random_effects_no_time > 0) {
    year_idx_ <- sort(unique(data$year_idx)) %>% `[`(1)
    map_random_effects <- data %>% filter(year_idx==year_idx_) %>% dplyr::select(data_idx, geo, region_id)
    functions_random_effects <- list()
    
    for (random_effect_col in results_summary$random_effect_col_names[[1]]$no_time) {
      
      functions_random_effects[[random_effect_col]] <-
        tibble(!!random_effect_col := model$summary.random[[random_effect_col]]$ID,
               mean = model$summary.random[[random_effect_col]]$mean,
               sd = model$summary.random[[random_effect_col]]$sd,
               lower_bound = model$summary.random[[random_effect_col]]$`0.025quant`,
               upper_bound = model$summary.random[[random_effect_col]]$`0.975quant`)
      
      random_effect_str <- paste0("reffect_", random_effect_col)
      
      map_random_effects_curr <-
        left_join(data %>% dplyr::select(data_idx, region_id, !!sym(random_effect_col)),
                  functions_random_effects[[random_effect_col]] %>%
                    rename(!!paste0(random_effect_str, ".mean") := mean,
                           !!paste0(random_effect_str, ".sd") := sd,
                           !!paste0(random_effect_str, ".lower_bound") := lower_bound,
                           !!paste0(random_effect_str, ".upper_bound") := upper_bound))
      map_random_effects <- left_join(map_random_effects, map_random_effects_curr)
      
    }
  } else {
    map_random_effects <- NULL
    functions_random_effects <- NULL
  }
  
  n_random_effects_with_time <- length(results_summary$random_effect_col_names[[1]]$with_time)
  if (n_random_effects_with_time > 0) {
    map_random_effects_with_time <- data %>% dplyr::select(data_idx, geo, region_id, year_idx)
    
    if (!exists("functions_random_effects"))
      functions_random_effects <- list()
    
    for (random_effect_col in results_summary$random_effect_col_names[[1]]$with_time) {
      
      functions_random_effects[[random_effect_col]] <-
        tibble(!!random_effect_col := model$summary.random[[random_effect_col]]$ID,
               mean = model$summary.random[[random_effect_col]]$mean,
               sd = model$summary.random[[random_effect_col]]$sd,
               lower_bound = model$summary.random[[random_effect_col]]$`0.025quant`,
               upper_bound = model$summary.random[[random_effect_col]]$`0.975quant`)
      
      random_effect_str <- paste0("reffect_", random_effect_col)
      
      map_random_effects_curr <-
        left_join(data %>% dplyr::select(data_idx, region_id, year_idx, !!sym(random_effect_col)),
                  functions_random_effects[[random_effect_col]] %>%
                    rename(!!paste0(random_effect_str, ".mean") := mean,
                           !!paste0(random_effect_str, ".sd") := sd,
                           !!paste0(random_effect_str, ".lower_bound") := lower_bound,
                           !!paste0(random_effect_str, ".upper_bound") := upper_bound))
      map_random_effects_with_time <- left_join(map_random_effects_with_time, map_random_effects_curr)
      
    }
  } else {
    
    map_random_effects_with_time <- NULL
  }
  
  if (!is.null(model$summary.random$spatial_idx$ID)) {
    
    map_spatial_effect <- tibble(spatial_idx = model$summary.random$spatial_idx$ID,
                                 spatial_field.mean = model$summary.random$spatial_idx$mean,
                                 spatial_field.sd = model$summary.random$spatial_idx$sd)
    
    if (nrow(map_spatial_effect) < nrow(data)) {
      map_spatial_effect <- left_join(map_spatial_effect,
                                      data %>% dplyr::select(geo, spatial_idx) %>% distinct()) %>%
        relocate(geo, spatial_idx)
    } else if (nrow(map_spatial_effect) == nrow(data)) {
      map_spatial_effect <- bind_cols(data %>% dplyr::select(geo, spatial_idx, year_idx), map_spatial_effect)
    } else {
      stop("nrow error")
    }
    
    results_summary$stdev_spatial <- sd(map_spatial_effect$spatial_field.mean)
  } else {
    
    map_spatial_effect <- NULL
    results_summary$stdev_spatial <- NA
  }
  
  map_full_model <- tibble(data_idx = data$data_idx,
                           geo = data$geo,
                           region_id = data$region_id,
                           year_idx = data$year_idx,
                           y = data[[y_column]],
                           p.mean = model$summary.fitted.values$mean[idx_linear_predictor],
                           p.sd = model$summary.fitted.values$sd[idx_linear_predictor],
                           fitted_log_odds.mean = logit(p.mean))
  
  if (cross_val) {
    map_full_model$y_orig <- data$y_orig
    map_full_model$training_set <- data$training_set
    map_full_model$test_set <- data$test_set
  }
  
  if (!is.null(idx_linear_predictor_pred)) {
    
    map_full_model_pred <- data_pred %>% dplyr::select(geo, year_idx)
    map_full_model_pred$p_pred.mean <-
      model$summary.fitted.values$mean[idx_linear_predictor_pred]
    map_full_model_pred$p_pred.sd <-
      model$summary.fitted.values$sd[idx_linear_predictor_pred]
    map_full_model_pred$fitted_log_odds.median <-
      logit(model$summary.fitted.values$mean[idx_linear_predictor_pred])
  } else {
    
    map_full_model_pred <- NULL
  }
  
  if (remove_NA) {
    map_full_model <- map_full_model %>%
      mutate(p.mean = ifelse(is.na(data[[y_column]]), NA, p.mean),
             p.sd = ifelse(is.na(data[[y_column]]), NA, p.sd),
             fitted_log_odds.median = ifelse(is.na(data[[y_column]]), NA, fitted_log_odds.median))
  }
  
  if (likelihood=="bernoulli") {
    map_full_model$error <- compute_error(map_full_model$p.mean, map_full_model$y)
  } else if (likelihood=="binomial") {
    map_full_model$error <- compute_error(map_full_model$p.mean, map_full_model$y, n=n_binom)
  }
  results_summary$total_abs_error <- sum(map_full_model$error, na.rm=TRUE)
  results_summary$mean_abs_error <- mean(map_full_model$error, na.rm=TRUE)
  
  if (likelihood=="bernoulli") {
    map_full_model$deviance <- compute_bernoulli_deviance(map_full_model$p.mean, map_full_model$y)
  } else if (likelihood=="binomial") {
    map_full_model$deviance <- compute_binomial_deviance(map_full_model$p.mean, map_full_model$y, n=n_binom)
  }
  results_summary$total_deviance <- sum(map_full_model$deviance, na.rm=TRUE)
  results_summary$mean_deviance <- mean(map_full_model$deviance, na.rm=TRUE)
  
  if (cross_val) {
    
    if (sum(map_full_model$test_set)==0) {
      
      warning("test set is empty")
      
      map_full_model$error_cv <- NA
      map_full_model$deviance_cv <- NA
      results_summary$total_abs_error_cv <- NA
      result_summary_mean_abs_error_cv <- NA
      results_summary$total_deviance_cv <- NA
      results_summary$mean_deviance_cv <- NA
    } else {
      
      if (likelihood=="bernoulli") {
        map_full_model$error_cv <-
          compute_error(map_full_model$p.mean, map_full_model$y_orig)
      } else if (likelihood=="binomial") {
        map_full_model$error_cv <-
          compute_error(map_full_model$p.mean, map_full_model$y_orig, n=n_binom)
      }
      results_summary$total_abs_error_cv <-
        sum(map_full_model %>% filter(test_set==TRUE) %>% pull(error_cv), na.rm=TRUE)
      results_summary$mean_abs_error_cv <-
        mean(map_full_model %>% filter(test_set==TRUE) %>% pull(error_cv), na.rm=TRUE)
      
      if (likelihood=="bernoulli") {
        map_full_model$deviance_cv <-
          compute_bernoulli_deviance(map_full_model$p.mean, map_full_model$y_orig)
      } else if (likelihood=="binomial") {
        map_full_model$deviance_cv <-
          compute_binomial_deviance(map_full_model$p.mean, map_full_model$y_orig, n=n_binom)
      }
      results_summary$total_deviance_cv <-
        sum(map_full_model %>% filter(test_set==TRUE) %>% pull(deviance_cv), na.rm=TRUE)
      results_summary$mean_deviance_cv <-
        mean(map_full_model %>% filter(test_set==TRUE) %>% pull(deviance_cv), na.rm=TRUE)
    }
  }
  
  results_summary <- results_summary %>%
    relocate(DIC, n_eff_parameters.DIC, matches("mean_deviance"), matches("mean_abs_error"))
  
  return(list(results_summary=results_summary,
              fixed_effects = list(model$summary.fixed),
              map_fixed_effects=map_fixed_effects,
              map_random_effects=map_random_effects,
              map_random_effects_with_time=map_random_effects_with_time,
              functions_random_effects=functions_random_effects,
              map_spatial_effect=map_spatial_effect,
              map_full_model=map_full_model,
              map_full_model_pred=map_full_model_pred))
}




# extract_main_results_prediction.
extract_main_results_prediction <- function(model,
                                            data_pred,
                                            y_column = "presence",
                                            idx_linear_predictor_pred = NULL) {
  
  if (is.null(idx_linear_predictor_pred)) {
    stop("idx_linear_predictor_pred must be supplied for prediction output.")
  }
  
  if (length(idx_linear_predictor_pred) != nrow(data_pred)) {
    stop("idx_linear_predictor_pred length does not match nrow(data_pred).")
  }
  
  if (!y_column %in% names(data_pred)) {
    stop(paste0("y_column '", y_column, "' not found in data_pred."))
  }
  
  p_mean_pred <- model$summary.fitted.values$mean[idx_linear_predictor_pred]
  p_sd_pred   <- model$summary.fitted.values$sd[idx_linear_predictor_pred]
  
  map_full_model_pred <- tibble(
    data_idx = data_pred$data_idx,
    geo      = data_pred$geo,
    region_id = data_pred$region_id,
    year_idx  = data_pred$year_idx,
    y         = data_pred[[y_column]],
    p.mean    = p_mean_pred,
    p.sd      = p_sd_pred,
    fitted_log_odds.mean = logit(p_mean_pred)
  )
  
  return(list(
    map_full_model_pred = map_full_model_pred
  ))
}




# compute_error.
compute_error <- function(p, y, n=1, abs=TRUE) {
  err <- rep(NA, length(p))
  
  idx_non_na <- which(!is.na(y))
  if (any(is.na(p[idx_non_na])))
    stop("na error")
  
  p_num <- p[idx_non_na]
  if (any(p_num[p_num] > 1-1e-10) || any(p_num[p_num] < 1e-10))
    warning("very high/low values of p")
  p_num[p_num > 1-1e-10] <- 1-1e-10
  p_num[p_num < 1e-10] <- 1e-10
  
  if (abs) {
    err[idx_non_na] <- abs(p_num*n - y[idx_non_na])
  } else {
    err[idx_non_na] <- p_num*n - y[idx_non_na]
  }
  
  return(err)
}




# Function: compute_bernoulli_deviance.
compute_bernoulli_deviance <- function(p, y) {
  deviance <- rep(NA, length(y))
  
  idx_non_na <- which(!is.na(y))
  
  if (any(is.na(p[idx_non_na])))
    stop("na error")
  
  p_num <- p[idx_non_na]
  y_num <- y[idx_non_na]
  
  if (any(p_num > 1-1e-10) || any(p_num < 1e-10))
    warning("very high/low values of p")
  p_num[p_num > 1-1e-10] <- 1-1e-10
  p_num[p_num < 1e-10] <- 1e-10
  
  deviance[idx_non_na] <- -2 * ((y_num * log(p_num) + (1-y_num) * log(1-p_num)))
  
  return(deviance)
}




# compute_binomial_deviance.
compute_binomial_deviance <- function(p, y, n) {
  
  if (length(n) != 1)
    stop("n length error")
  
  deviance <- rep(NA, length(y))
  
  idx_non_na <- which(!is.na(y))
  
  if (any(is.na(p[idx_non_na])))
    stop("na error")
  
  p_num <- p[idx_non_na]
  y_num <- y[idx_non_na]
  
  if (any(p_num > 1-1e-10) || any(p_num < 1e-10))
    warning("very high/low values of p")
  p_num[p_num > 1-1e-10] <- 1-1e-10
  p_num[p_num < 1e-10] <- 1e-10
  
  deviance[idx_non_na] <- -2 * (log(choose(n, y_num)) + y_num * log(p_num) + (n-y_num) * log(1-p_num))
  
  return(deviance)
}



# extract_posterior_fitted_samples
extract_posterior_fitted_samples <- function(model,
                                             data,
                                             idx_linear_predictor,
                                             n_samples = 1000L,
                                             likelihood = "bernoulli",
                                             seed = NULL) {
  
  if (!likelihood %in% c("bernoulli", "binomial")) {
    stop("Only implemented for bernoulli / binomial (logit link).")
  }
  
  if (is.null(idx_linear_predictor)) {
    stop("idx_linear_predictor must be supplied (indices of Predictor in latent field).")
  }
  if (length(idx_linear_predictor) != nrow(data)) {
    stop("length(idx_linear_predictor) must equal nrow(data).")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  samp <- INLA::inla.posterior.sample(
    n      = n_samples,
    result = model,
    seed   = seed
  )
  
  eta_mat <- sapply(
    samp,
    function(x) as.numeric(x$latent[idx_linear_predictor])
  )
  
  p_mat <- plogis(eta_mat)
  
  n_data <- nrow(data)
  stopifnot(n_data == nrow(p_mat))
  
  p_df <- as_tibble(p_mat, .name_repair = "minimal")
  names(p_df) <- paste0("draw_", seq_len(n_samples))
  
  p_long <- p_df %>%
    mutate(
      data_idx  = data$data_idx,
      geo       = data$geo,
      region_id = data$region_id,
      year_idx  = data$year_idx
    ) %>%
    pivot_longer(
      cols      = starts_with("draw_"),
      names_to  = "draw",
      values_to = "p"
    ) %>%
    mutate(draw = as.integer(sub("^draw_", "", draw))) %>%
    arrange(draw, data_idx)
  
  return(p_long)
}




#  extract_posterior_fitted_samples_prediction.
extract_posterior_fitted_samples_prediction <- function(model,
                                                        data_pred,
                                                        idx_linear_predictor_pred,
                                                        n_samples  = 1000L,
                                                        likelihood = "bernoulli",
                                                        seed       = NULL) {
  
  if (!likelihood %in% c("bernoulli", "binomial")) {
    stop("Only implemented for bernoulli / binomial (logit link).")
  }
  
  if (is.null(idx_linear_predictor_pred)) {
    stop("idx_linear_predictor_pred must be supplied (prediction indices of Predictor in latent field).")
  }
  
  if (length(idx_linear_predictor_pred) != nrow(data_pred)) {
    stop("length(idx_linear_predictor_pred) must equal nrow(data_pred).")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  samp <- INLA::inla.posterior.sample(
    n      = n_samples,
    result = model,
    seed   = seed
  )
  
  eta_mat_pred <- sapply(
    samp,
    function(x) as.numeric(x$latent[idx_linear_predictor_pred])
  )
  
  p_mat_pred <- plogis(eta_mat_pred)
  
  n_pred <- nrow(data_pred)
  stopifnot(n_pred == nrow(p_mat_pred))
  
  p_df_pred <- as_tibble(p_mat_pred, .name_repair = "minimal")
  
  names(p_df_pred) <- paste0("draw_", seq_len(n_samples))
  
  p_long_pred <- p_df_pred %>%
    mutate(
      data_idx  = data_pred$data_idx,
      geo       = data_pred$geo,
      region_id = data_pred$region_id,
      year_idx  = data_pred$year_idx
    ) %>%
    tidyr::pivot_longer(
      cols      = dplyr::starts_with("draw_"),
      names_to  = "draw",
      values_to = "p"
    ) %>%
    mutate(
      draw = as.integer(sub("^draw_", "", draw))
    ) %>%
    arrange(draw, data_idx)
  
  return(p_long_pred)
}




# compute_bin_keys.
compute_bin_keys <- function(data, bin_column, original_column,
                             print_format="[%.1f,%.1f]") {
  
  bin_range_column <- paste0("bin_range_", original_column)
  bin_mean_column <- paste0("bin_mean_", original_column)
  
  keys_to_bins <- tibble(bin_idx = sort(unique(data[[bin_column]])),
                         !!sym(bin_range_column) := NA,
                         bin_min = NA, bin_max = NA,
                         !!sym(bin_mean_column) := NA,
                         n_in_bin = NA)
  
  for (j in 1:nrow(keys_to_bins)) {
    
    curr_bin_idx <- keys_to_bins$bin_idx[j]
    data_curr_bin <- data %>% filter(!!sym(bin_column) == curr_bin_idx)
    
    bin_range <- range(data_curr_bin[[original_column]])
    bin_range <- round(bin_range, digits=13)
    keys_to_bins[[bin_range_column]][j] <- sprintf(print_format, bin_range[1], bin_range[2])
    keys_to_bins$bin_min[j] <- bin_range[1]
    keys_to_bins$bin_max[j] <- bin_range[2]
    
    keys_to_bins[[bin_mean_column]][j] <- mean(data_curr_bin[[original_column]])
    
    keys_to_bins$n_in_bin[j] <- nrow(data_curr_bin)
    
  }
  
  keys_to_bins <- keys_to_bins %>% mutate(!!sym(bin_range_column) :=
                                            factor(!!sym(bin_range_column),
                                                   levels=!!sym(bin_range_column), ordered=TRUE))
  
  return(keys_to_bins)
}





# linked_random_effect_and_histogram_plot.
linked_random_effect_and_histogram_plot <- function(df_random_effect,
                                                    x_column_type="bin_mean",
                                                    interaction_x_column=NULL,
                                                    random_effect_title="",
                                                    bin_width_constant=TRUE,
                                                    vector_to_plot_with_rug=NULL,
                                                    plot_larger=FALSE,
                                                    ylim=ylim,
                                                    reduce_window_size=TRUE) {
  
  n_column <- names(df_random_effect) %>% str_subset("n_in_bin")
  
  if (x_column_type=="bin_mean") {
    
    x_column <- names(df_random_effect) %>% str_subset("bin_mean")
  } else if (x_column_type=="bin_idx") {
    
    x_column <- names(df_random_effect) %>% str_subset("bin_idx")
  } else if (x_column_type=="interaction_bin_idx") {
    
    x_column <- names(df_random_effect) %>% str_subset("bin_idx.*inter")
    
    if (!bin_width_constant) {
      warning("setting bin_width_constant=TRUE, which is required when x_column_type=\"interaction_bin_idx\"")
      bin_width_constant <- TRUE
    }
  } else {
    
    stop("x_column_type should equal either \"bin_mean\", \"bin_idx\", or \"interaction_bin_idx\"")
  }
  
  if (length(x_column) != 1 || length(n_column) != 1)
    stop("interaction_x_column does not match any of the interaction_bin_mean-columns")
  
  df_random_effect_tall <- bind_rows(tibble(!!sym(x_column) := df_random_effect[[x_column]],
                                            random_effect = df_random_effect$lower_bound,
                                            type = "lower_bound"),
                                     tibble(!!sym(x_column) := df_random_effect[[x_column]],
                                            random_effect = df_random_effect$mean,
                                            type = "mean"),
                                     tibble(!!sym(x_column) := df_random_effect[[x_column]],
                                            random_effect = df_random_effect$upper_bound,
                                            type = "upper_bound"))
  
  if (plot_larger) {
    palette_idx <- NA
    col_rgb <- col2rgb("gray60")
    col_hex <- rgb(col_rgb[1], col_rgb[2], col_rgb[3], maxColorValue=255)
    manual_palette <- c(col_hex, NA, col_hex)
    col_rgb <- col2rgb("gray5")
    col_hex <- rgb(col_rgb[1], col_rgb[2], col_rgb[3], maxColorValue=255)
    manual_palette[2] <- col_hex
  } else {
    palette_idx <- 4
    manual_palette <- NULL
  }
  
  if (bin_width_constant) {
    
    p_random_effect <-
      plot_scatter_lines(df_random_effect_tall, x_column = x_column, y_column = "random_effect",
                         color_column = "type", plot_constant_lines_y = 0,
                         palette_idx=palette_idx, manual_palette = manual_palette,
                         title=random_effect_title,
                         plot_larger=plot_larger,
                         ylim=ylim)
    
    xlimits <- layer_scales(p_random_effect)$x$range$range
    p_hist <- ggplot(df_random_effect, aes(x = !!sym(x_column), y = !!sym(n_column))) +
      geom_col() +
      coord_cartesian(xlim = xlimits) +
      theme_classic() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank()
      )
  } else if (!bin_width_constant && reduce_window_size) {
    
    p_random_effect <-
      plot_scatter_lines(df_random_effect_tall, x_column = x_column, y_column = "random_effect",
                         color_column = "type", plot_constant_lines_y = 0,
                         palette_idx=palette_idx, manual_palette = manual_palette,
                         title=random_effect_title,
                         plot_larger=plot_larger,
                         ylim=ylim)
    
    xlimits <- layer_scales(p_random_effect)$x$range$range
    
    first_bin <- list()
    first_bin <- df_random_effect %>% head(1) %>% dplyr::select(bin_min, bin_max) %>%
      rename(original_bin_min=bin_min)
    first_bin$bin_min <- xlimits[1]
    first_bin <- first_bin %>% mutate(original_length=bin_max-original_bin_min,
                                      new_length=bin_max-bin_min)
    
    df_random_effect[1,"bin_min"] <- first_bin$bin_min
    
    last_bin <- list()
    last_bin <- df_random_effect %>% tail(1) %>% dplyr::select(bin_min, bin_max) %>%
      rename(original_bin_max=bin_max)
    last_bin$bin_max <- xlimits[2]
    last_bin <- last_bin %>% mutate(
      original_length=original_bin_max-bin_min,
      new_length=bin_max-bin_min)
    
    df_random_effect[nrow(df_random_effect),"bin_max"] <- last_bin$bin_max
    
    print(first_bin)
    print(last_bin)
    
    df_random_effect <- df_random_effect %>%
      mutate(bin_height=n_in_bin/(bin_max-bin_min))
    
    p_hist <- ggplot(df_random_effect, aes(x = !!sym(x_column), y = bin_height)) +
      geom_rect(aes(xmin = bin_min, xmax = bin_max,
                    ymin = 0, ymax = bin_height)) +
      theme_classic() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank()
      )
  } else if (!bin_width_constant && !reduce_window_size) {
    
    df_random_effect <- df_random_effect %>%
      mutate(bin_height=n_in_bin/(bin_max-bin_min))
    
    p_hist <- ggplot(df_random_effect, aes(x = !!sym(x_column), y = bin_height)) +
      geom_rect(aes(xmin = bin_min, xmax = bin_max,
                    ymin = 0, ymax = bin_height)) +
      theme_classic() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank()
      )
    
    xlimits <- layer_scales(p_hist)$x$range$range
    p_random_effect <-
      plot_scatter_lines(df_random_effect_tall, x_column = x_column, y_column = "random_effect",
                         color_column = "type", plot_constant_lines_y = 0,
                         palette_idx=palette_idx, manual_palette = manual_palette,
                         title=random_effect_title,
                         xlim=xlimits,
                         plot_larger=plot_larger,
                         ylim=ylim)
  }
  
  if (!is.null(vector_to_plot_with_rug)) {
    
    p_random_effect <- p_random_effect +
      geom_rug(data=data.frame(x=vector_to_plot_with_rug), aes(x=x, y=0),
               inherit.aes = FALSE, alpha = 1/2, position = "jitter", sides="b")
  }
  
  g <- align_plots(p_random_effect, p_hist)
  grid.draw(g)
  
  return(list(g=g, p_random_effect=p_random_effect, p_hist=p_hist))
}




# calibration_curve.
calibration_curve <- function(p, y, y_orig=NULL, n_binom=1,
                              n_intervals = 100, by_quantiles=TRUE,
                              n_min=NULL, min_p_range=NULL) {
  
  if (by_quantiles) {
    interval <- 0.20
    if (is.null(n_min)) {
      n_min <- 100
    }
    if (is.null(min_p_range)) {
      min_p_range <- 0.05
    }
  } else {
    interval <- 0.05
    if (is.null(n_min)) {
      n_min <- 500
    }
    if (is.null(min_p_range)) {
      min_p_range <- 0.05
    }
  }
  
  breaks <- seq(from=0, to=1, length.out=(n_intervals+1))
  if (breaks[length(breaks)] != 1)
    stop("breaks error")
  
  n_per_interval <- which(breaks >= interval)[1]
  
  breaks_lower <- breaks[1:(length(breaks)-n_per_interval+1)]
  breaks_upper <- breaks[n_per_interval:length(breaks)]
  
  if (by_quantiles) {
    level_values <- tibble(n = NA, q_lower = breaks_lower, q_upper = breaks_upper,
                           p_lower = NA, p_upper = NA, p_mean = NA,
                           empirical_p = NA, empirical_p_w_CV = NA) %>% mutate(level = 1:n()) %>%
      relocate(level)
  } else {
    level_values <- tibble(n = NA, q_lower = NA, q_upper = NA,
                           p_lower = breaks_lower, p_upper = breaks_upper, p_mean = NA,
                           empirical_p = NA, empirical_p_w_CV = NA) %>% mutate(level = 1:n()) %>%
      relocate(level)
  }
  
  ind_groups <- rep(NA, nrow(level_values))
  for (j in 1:nrow(level_values)) {
    if (by_quantiles) {
      q_lower <- breaks_lower[j]
      q_upper <- breaks_upper[j]
      p_quantiles <- quantile(p, c(q_lower, q_upper))
      level_values$p_lower[j] <- p_quantiles[1]
      level_values$p_upper[j] <- p_quantiles[2]
    } else {
      p_ecdf <- ecdf(p)
      p_lower <- breaks_lower[j]
      p_upper <- breaks_upper[j]
      qs <- p_ecdf(c(p_lower, p_upper))
      level_values$q_lower[j] <- qs[1]
      level_values$q_upper[j] <- qs[2]
    }
    
    ind_level <- which(p <= level_values$p_upper[j] & p >= level_values$p_lower[j])
    
    level_values$n[j] <- length(ind_level)
    
    while (level_values$n[j] < n_min) {
      level_values$p_upper[j] <- level_values$p_upper[j]*1.05
      level_values$p_lower[j] <- level_values$p_lower[j]/1.05
      level_values$p_upper[j] <- min(1, level_values$p_upper[j])
      level_values$p_lower[j] <- max(0, level_values$p_lower[j])
      
      ind_level <- which(p <= level_values$p_upper[j] & p >= level_values$p_lower[j])
      
      level_values$n[j] <- length(ind_level)
    }
    p_range <- level_values$p_upper[j] - level_values$p_lower[j]
    while (p_range < min_p_range) {
      level_values$p_upper[j] <- level_values$p_upper[j]*1.05
      level_values$p_lower[j] <- level_values$p_lower[j]/1.05
      level_values$p_upper[j] <- min(1, level_values$p_upper[j])
      level_values$p_lower[j] <- max(0, level_values$p_lower[j])
      
      ind_level <- which(p <= level_values$p_upper[j] & p >= level_values$p_lower[j])
      
      p_range <- level_values$p_upper[j] - level_values$p_lower[j]
      level_values$n[j] <- length(ind_level)
    }
    
    level_values$p_mean[j] <- mean(p[ind_level])
    ind_groups[[j]] <- list(ind_level)
    
    n_tot <- sum(n_binom * !is.na(y[ind_level]))
    
    frac_1s <- sum(y[ind_level], na.rm=TRUE) / n_tot
    level_values$empirical_p[j] <- frac_1s
    if (!is.null(y_orig)) {
      n_tot_orig <- sum(n_binom * !is.na(y_orig[ind_level]))
      
      frac_1s_orig <- sum(y_orig[ind_level], na.rm=TRUE) / n_tot_orig
      level_values$empirical_p_w_CV[j] <- frac_1s_orig
    }
    
  }
  
  level_vaues <- level_values %>% dplyr::select(-level) %>% distinct() %>% mutate(level=1:n()) %>%
    arrange(level, empirical_p)
  
  return(list(level_values=level_values, ind_groups=ind_groups))
}




#  plot_calibration_curve.
plot_calibration_curve <- function(level_values, title="") {
  
  ps <- list()
  
  ps[[1]] <- plot_scatter_lines(level_values,
                                "empirical_p", "p_mean", plot_diagonal=TRUE)
  
  logit_transformation <- trans_new(
    "logit",
    trans=function(x) {shift = 1e-5; z = x * (1-2*shift) + shift;
    log( z / (1-z) )},
    inverse=function(y) {shift = 1e-5; t = (y - shift) / (1-2*shift) ;
    1 / (1 + exp(-t))},
    domain = c(0,1)
  )
  
  level_values <- level_values %>% mutate(n_size = n-min(n),
                                          n_size = n_size/max(n_size)*2 + 0.5)
  
  level_values_long <- level_values %>%
    select(p_mean, p_lower) %>%
    rename(y = p_lower) %>%
    mutate(type = "min. p",
           n_size = 0.5) %>%
    bind_rows(level_values %>%
                select(p_mean, empirical_p, n_size) %>%
                rename(y = empirical_p) %>%
                mutate(type = "empirical p")) %>%
    bind_rows(level_values %>%
                select(p_mean, p_upper) %>%
                rename(y = p_upper) %>%
                mutate(type = "max. p",
                       n_size = 0.5)) %>%
    mutate(type = factor(type))
  
  ps[[2]] <-
    plot_scatter_lines(level_values_long, "p_mean", "y", color_column = "type", size_column="n_size",
                       plot_diagonal=TRUE, plot_lines=TRUE, plot_smooth=FALSE,
                       xlim=c(2.5e-3, 1-1e-3),
                       ylim=c(2.5e-3, 1-1e-3),
                       axis_transformation_x = logit_transformation,
                       axis_transformation_y = logit_transformation,
                       palette_idx=4,
                       title=title)
  
  return(ps)
}




#  plot_random_effect.
plot_random_effect <- function(functions_random_effects, keys,
                               col_subscript, vector_to_plot_with_rug=NULL,
                               filename=NULL,
                               plot_larger=FALSE,
                               ylim=NULL,
                               resolution_scaling_width=1,
                               resolution_scaling_height=1) {
  
  bin_idx_col <- paste0("bin_idx_", col_subscript)
  bin_range_col <- paste0("bin_range_", col_subscript)
  bin_mean_col <- paste0("bin_mean_", col_subscript)
  
  tb_reffect <- functions_random_effects[[bin_idx_col]]
  tb_reffect <- left_join(tb_reffect, keys,
                          by = join_by(!!sym(bin_idx_col) == bin_idx)) %>%
    relocate(sym(bin_idx_col), sym(bin_range_col), sym(bin_mean_col))
  
  re_plots <-
    linked_random_effect_and_histogram_plot(tb_reffect,
                                            x_column_type="bin_mean",
                                            bin_width_constant=FALSE,
                                            vector_to_plot_with_rug=vector_to_plot_with_rug,
                                            plot_larger=plot_larger,
                                            ylim=ylim)
  print(tb_reffect)
  if (!is.null(filename)) {
    save_to_png(re_plots$g, filename=filename,
                resolution_scaling_width=resolution_scaling_width,
                resolution_scaling_height=resolution_scaling_height)
  }
  
  return(re_plots$g)
}




#  plot_proximity_random_effect.
plot_proximity_random_effect <- function(functions_random_effects, keys, data,
                                         filename=NULL,
                                         plot_rug_vector=FALSE,
                                         plot_larger=FALSE,
                                         ylim=NULL,
                                         resolution_scaling_width=1) {
  
  bin_idx_proximity_col <- names(functions_random_effects) %>% str_subset("proximity")
  proximity_exp <- str_extract(bin_idx_proximity_col, "(?<=proximity)[0-9].*")
  keys_col <- str_subset(names(keys), paste0("exp", proximity_exp))
  
  reffect_proximity <- functions_random_effects[[bin_idx_proximity_col]]
  if (plot_larger) {
    reffect_proximity <- reffect_proximity %>% filter(bin_idx_proximity050>=3)
  }
  reffect_proximity <- left_join(reffect_proximity, keys[[keys_col]],
                                 by = join_by(!!sym(bin_idx_proximity_col) == bin_idx))
  
  if (plot_rug_vector) {
    rug_vector <- data$log_population
  } else {
    rug_vector <- NULL
  }
  
  re_plots_proximity <-
    linked_random_effect_and_histogram_plot(reffect_proximity,
                                            x_column_type="bin_mean",
                                            bin_width_constant=FALSE,
                                            vector_to_plot_with_rug=rug_vector,
                                            plot_larger=plot_larger,
                                            ylim=ylim)
  print(reffect_proximity)
  
  save_to_png(re_plots_proximity$g, filename=filename, resolution_scaling_width=resolution_scaling_width)
  
  return(re_plots_proximity$g)
}




#  plot_proximity_random_effect_on_map.
plot_proximity_random_effect_on_map <- function(map_random_effects,
                                                map_inla,
                                                special_group_bin_idx,
                                                year_to_plot=9,
                                                plot_larger=FALSE) {
  proximity_col <- names(map_random_effects) %>% str_subset("reffect_bin_idx_proximity.*mean")
  
  ps <- list()
  ps[[1]] <- plot_sf(left_join(map_inla, map_random_effects %>% filter(year_idx==year_to_plot)),
                     proximity_col, plot_larger=plot_larger)
  
  bin_idx_col <- names(map_random_effects) %>% str_subset("bin_idx_proximity[0-9]*$")
  row_idx_special_groups <- map_random_effects[[bin_idx_col]] %in% special_group_bin_idx
  
  map_random_effects$proximity_excl_special_groups <-
    map_random_effects[[proximity_col]]
  map_random_effects$proximity_excl_special_groups[row_idx_special_groups] <- NA
  
  ps[[2]] <- plot_sf(left_join(map_inla, map_random_effects %>% filter(year_idx==year_to_plot)),
                     "proximity_excl_special_groups",
                     title=paste0("year ", year_to_plot),
                     plot_larger=plot_larger)
  
  if (!plot_larger) {
    p <- do.call(ggarrange, c(ps, nrow=2, ncol=1))
  } else {
    p <- ps[[2]]
  }
  
  return(p)
}




#  plot_random_effects_and_full_model_on_map.
plot_random_effects_and_full_model_on_map <- function(map_random_effects,
                                                      map_full_model,
                                                      map_spatial_effect=NULL,
                                                      map_geometry,
                                                      year_to_plot=9,
                                                      plot_larger=FALSE) {
  ps_map <- list()
  
  if (!is.null(map_random_effects) && nrow(map_random_effects)>0) {
    
    if ("reffect_bin_idx_max_temp.mean" %in% names(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_maxtemp = reffect_bin_idx_max_temp.mean),
                "reffect_maxtemp",
                title="random effect, max temp",
                plot_larger=plot_larger)
    }
    
    if ("reffect_bin_idx_median_temp.mean" %in% names(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_mediantemp = reffect_bin_idx_median_temp.mean),
                "reffect_mediantemp",
                title="random effect, median temp",
                plot_larger=plot_larger)
    }
    
    if ("reffect_bin_idx_min_temp.mean" %in% names(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_mintemp =
                           reffect_bin_idx_min_temp.mean),
                "reffect_mintemp",
                title="random effect, min temp",
                plot_larger=plot_larger)
    }
    
    if ("reffect_bin_idx_min_temp_merged.mean" %in% names(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_mintemp_merged =
                           reffect_bin_idx_min_temp_merged.mean),
                "reffect_mintemp_merged",
                title="random effect, min temp, merged",
                plot_larger=plot_larger)
    }
    
    if ("reffect_bin_idx_min_temp_CROSS.mean" %in% names(map_random_effects)) {
      map_random_effects <- map_random_effects %>%
        mutate(reffect_bin_idx_min_temp_CROSS.meanNA=
                 if_else(bin_idx_min_temp_CROSS==1, NA,
                         reffect_bin_idx_min_temp_CROSS.mean))
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_mintemp_CROSS =
                           reffect_bin_idx_min_temp_CROSS.meanNA),
                "reffect_mintemp_CROSS",
                title="random effect, min temp, CROSS",
                plot_larger=plot_larger)
    }
    
    if ("reffect_bin_idx_median_precip.mean" %in% names(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_medianprecip =
                           reffect_bin_idx_median_precip.mean),
                "reffect_medianprecip",
                title="random effect, median precip",
                plot_larger=plot_larger)
    }
    
    if ("reffect_bin_idx_median_relative_humidity.mean" %in% names(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_medianhumidity =
                           reffect_bin_idx_median_relative_humidity.mean),
                "reffect_medianhumidity",
                title="random effect, median relative humidity",
                plot_larger=plot_larger)
    }
    
    if ("reffect_bin_idx_log_population.mean" %in% names(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_logpop =
                           reffect_bin_idx_log_population.mean),
                "reffect_logpop",
                title="random effect, log population",
                plot_larger=plot_larger)
    }
    
    if ("reffect_bin_idx_log_pop_merged.mean" %in% names(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_logpop_merged =
                           reffect_bin_idx_log_pop_merged.mean),
                "reffect_logpop_merged",
                title="random effect, log population, merged",
                plot_larger=plot_larger)
    }
    
    if ("reffect_bin_idx_mean_adults_pow_merged.mean" %in% names(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_random_effects %>% left_join(map_geometry) %>%
                  rename(reffect_meanODE_merged =
                           reffect_bin_idx_mean_adults_pow_merged.mean),
                "reffect_meanODE_merged",
                title="random effect, mean ODE, merged",
                plot_larger=plot_larger)
    }
    
  }
  
  n_time <- length(unique(map_full_model$year_idx))
  if (n_time==1) {
    year_to_plot_ <- unique(map_full_model$year_idx)
  } else if (n_time>1) {
    year_to_plot_ <- year_to_plot
    if (!year_to_plot %in% unique(map_full_model$year_idx)) {
      stop("year does not exist")
    }
  } else {
    stop("not implemented")
  }
  
  if (!is.null(map_spatial_effect)) {
    spatial_column <- str_subset(names(map_spatial_effect), "spde_integrated|spatial_field.mean")
    
    if (nrow(map_spatial_effect) == nrow(map_random_effects)) {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_spatial_effect %>% left_join(map_geometry), spatial_column,
                title="spatial effect",
                plot_larger=plot_larger)
    } else {
      ps_map[[length(ps_map)+1]] <-
        plot_sf(map_spatial_effect %>% left_join(map_geometry) %>%
                  filter(year_idx==year_to_plot_), spatial_column,
                title="spatial effect",
                plot_larger=plot_larger)
    }
  }
  
  ps_map[[length(ps_map)+1]] <-
    plot_sf(map_full_model %>% filter(year_idx==year_to_plot_) %>%
              left_join(map_geometry),
            "p.mean",
            title="probability",
            plot_larger=plot_larger)
  ps_map[[length(ps_map)+1]] <-
    plot_sf(map_full_model %>% filter(year_idx==year_to_plot_) %>%
              mutate(p_from_sumNA = ifelse(presence_previous_year==1,NA,p_from_sum)) %>%
              left_join(map_geometry),
            "p_from_sumNA",
            title="probability",
            plot_larger=plot_larger)
  ps_map[[length(ps_map)+1]] <- plot_sf(map_full_model %>% filter(year_idx==year_to_plot_) %>%
                                          left_join(map_geometry),
                                        "y",
                                        title="data",
                                        plot_larger=plot_larger)
  
  n_figures <- length(ps_map)
  ncol <- floor(sqrt(n_figures))
  nrow <- ceiling(n_figures / ncol)
  p <- do.call(ggarrange, c(ps_map, nrow=nrow, ncol=ncol))
  
  return(list(p=p, ps_map=ps_map))
}



# plot_random_effects_as_functions.
plot_random_effects_as_functions <- function(functions_random_effects, keys,
                                             data_inla, model_str,
                                             plot_rug_vector=FALSE,
                                             plot_larger=FALSE,
                                             ylim=NULL,
                                             resolution_scaling_width=1) {
  
  random_effects_col <- names(functions_random_effects)
  
  ps <- list()
  if ("bin_idx_max_temp" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$max_temp
    } else {
      rug_vector <- NULL
    }
    
    p_max_temp <- plot_random_effect(functions_random_effects, keys=keys$max_temp,
                                     col_subscript="max_temp",
                                     vector_to_plot_with_rug=rug_vector,
                                     filename=paste0("Plots/INLA_MODEL", model_str,
                                                     "_max_temp_reffect.png"),
                                     plot_larger=plot_larger,
                                     ylim=ylim,
                                     resolution_scaling_width=resolution_scaling_width)
    
    ps$p_max_temp <- p_max_temp
  }
  
  if ("bin_idx_median_temp" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$median_temp
    } else {
      rug_vector <- NULL
    }
    
    p_median_temp <- plot_random_effect(functions_random_effects, keys=keys$median_temp,
                                        col_subscript="median_temp",
                                        vector_to_plot_with_rug=rug_vector,
                                        filename=paste0("Plots/INLA_MODEL", model_str,
                                                        "_median_temp_reffect.png"),
                                        plot_larger=plot_larger,
                                        ylim=ylim,
                                        resolution_scaling_width=resolution_scaling_width)
    
    ps$p_median_temp <- p_median_temp
  }
  
  if ("bin_idx_min_temp" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$min_temp
    } else {
      rug_vector <- NULL
    }
    
    p_min_temp <- plot_random_effect(functions_random_effects, keys=keys$min_temp,
                                     col_subscript="min_temp",
                                     vector_to_plot_with_rug=rug_vector,
                                     filename=paste0("Plots/INLA_MODEL", model_str,
                                                     "_min_temp_reffect.png"),
                                     plot_larger=plot_larger,
                                     ylim=ylim,
                                     resolution_scaling_width=resolution_scaling_width)
    
    ps$p_min_temp <- p_min_temp
  }
  
  if ("bin_idx_min_temp_merged" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$min_temp
    } else {
      rug_vector <- NULL
    }
    
    p_min_temp_merged <- plot_random_effect(functions_random_effects, keys=keys$min_temp_merged,
                                            col_subscript="min_temp_merged",
                                            vector_to_plot_with_rug=rug_vector,
                                            filename=paste0("Plots/INLA_MODEL", model_str,
                                                            "_min_temp_merged_reffect.png"),
                                            plot_larger=plot_larger,
                                            ylim=ylim,
                                            resolution_scaling_width=resolution_scaling_width)
    
    ps$p_min_temp_merged <- p_min_temp_merged
  }
  
  if ("bin_idx_min_temp_CROSS" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$min_temp
    } else {
      rug_vector <- NULL
    }
    
    functions_random_effects[["bin_idx_min_temp_CROSS"]] <-
      functions_random_effects[["bin_idx_min_temp_CROSS"]] %>%
      filter(bin_idx_min_temp_CROSS >= 2)
    
    p_min_temp_CROSS <- plot_random_effect(functions_random_effects, keys=keys$min_temp_CROSS,
                                           col_subscript="min_temp_CROSS",
                                           vector_to_plot_with_rug=rug_vector,
                                           filename=paste0("Plots/INLA_MODEL", model_str,
                                                           "_min_temp_CROSS_reffect.png"),
                                           plot_larger=plot_larger,
                                           ylim=ylim,
                                           resolution_scaling_width=resolution_scaling_width)
    
    ps$p_min_temp_CROSS <- p_min_temp_CROSS
  }
  
  if ("bin_idx_median_precip" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$median_precip
    } else {
      rug_vector <- NULL
    }
    
    p_median_precip <- plot_random_effect(functions_random_effects, keys=keys$median_precip,
                                          col_subscript="median_precip",
                                          vector_to_plot_with_rug=rug_vector,
                                          filename=paste0("Plots/INLA_MODEL", model_str,
                                                          "_median_precip_reffect.png"),
                                          plot_larger=plot_larger,
                                          ylim=ylim,
                                          resolution_scaling_width=resolution_scaling_width)
    
    ps$p_median_precip <- p_median_precip
  }
  
  if ("bin_idx_median_relative_humidity" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$median_relative_humidity
    } else {
      rug_vector <- NULL
    }
    
    p_median_relative_humidity <-
      plot_random_effect(functions_random_effects, keys=keys$median_relative_humidity,
                         col_subscript="median_relative_humidity",
                         vector_to_plot_with_rug=rug_vector,
                         filename=paste0("Plots/INLA_MODEL", model_str,
                                         "_median_relative_humidity_reffect.png"),
                         plot_larger=plot_larger,
                         ylim=ylim,
                         resolution_scaling_width=resolution_scaling_width)
    
    ps$p_median_relative_humidity <- p_median_relative_humidity
  }
  
  if ("bin_idx_log_population" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$log_population
    } else {
      rug_vector <- NULL
    }
    
    p_log_pop <- plot_random_effect(functions_random_effects, keys=keys$log_pop,
                                    col_subscript="log_population",
                                    vector_to_plot_with_rug=rug_vector,
                                    filename=paste0("Plots/INLA_MODEL", model_str,
                                                    "_log_pop.png"),
                                    plot_larger=plot_larger,
                                    ylim=ylim,
                                    resolution_scaling_width=resolution_scaling_width)
    
    ps$p_log_pop <- p_log_pop
  }
  
  if ("bin_idx_log_pop_merged" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$log_population
    } else {
      rug_vector <- NULL
    }
    
    p_log_pop_merged <- plot_random_effect(functions_random_effects, keys=keys$log_pop_merged,
                                           col_subscript="log_pop_merged",
                                           vector_to_plot_with_rug=rug_vector,
                                           filename=paste0("Plots/INLA_MODEL", model_str,
                                                           "_log_pop_merged.png"),
                                           plot_larger=plot_larger,
                                           ylim=ylim,
                                           resolution_scaling_width=resolution_scaling_width)
    
    ps$p_log_pop_merged <- p_log_pop_merged
  }
  
  if ("bin_idx_mean_adults_pow_merged" %in% random_effects_col) {
    if (plot_rug_vector) {
      rug_vector <- data_inla$mean_adults_pow_merged
    } else {
      rug_vector <- NULL
    }
    
    p_mean_ODE <- plot_random_effect(functions_random_effects, keys=keys$mean_adults_pow_merged,
                                     col_subscript="mean_adults_pow_merged",
                                     vector_to_plot_with_rug=rug_vector,
                                     filename=paste0("Plots/INLA_MODEL", model_str,
                                                     "_mean_ODE_merged.png"),
                                     plot_larger=plot_larger,
                                     ylim=ylim,
                                     resolution_scaling_width=resolution_scaling_width)
    
    ps$p_mean_ODE <- p_mean_ODE
  }
  
  if (any(str_detect(random_effects_col, "proximity"))) {
    p_proximity <- plot_proximity_random_effect(functions_random_effects, keys=keys, data=data_inla,
                                                filename=paste0("Plots/INLA_MODEL", model_str,
                                                                "_proximity.png"),
                                                plot_rug_vector=plot_rug_vector,
                                                plot_larger=plot_larger,
                                                ylim=ylim,
                                                resolution_scaling_width=resolution_scaling_width)
    
    ps$p_proximity <- p_proximity
  }
  
  return(list(ps))
}








# bin_covariate_for_random_effect__CROSS.
bin_covariate_for_random_effect__CROSS <- function(data_spde, keys,
                                                   join_categories=FALSE) {
  
  df_to_cross <- data_spde %>% dplyr::select(data_idx, year_idx,
                                             bin_idx_median_temp,
                                             bin_idx_min_temp_merged,
                                             min_temp,
                                             presence) %>%
    mutate(spatial_idx=1:n())
  
  df_to_cross <- df_to_cross %>% mutate(
    median_temp_class = case_when(
      bin_idx_median_temp <= 2 ~ 1,
      TRUE ~ 2
    )
  )
  
  out <- cross_variables(df_to_cross, "median_temp_class", "bin_idx_min_temp_merged",
                         "medianXmin_temp")
  
  df_cross <- left_join(df_to_cross, out$data_cross)
  
  coords_altered <- out$coords
  coords_altered <- coords_altered %>% mutate(value_str = as.character(value),
                                              value_new = ifelse(
                                                x==2, value_str, "1x-"
                                              )) %>%
    rename(medianXmin_temp_CROSS=value,
           min_temp_CROSS=value_new) %>%
    mutate(bin_idx_min_temp_CROSS = match(min_temp_CROSS,
                                          sort(unique(min_temp_CROSS))))
  
  df_cross <- left_join(df_cross, coords_altered %>% dplyr::select(matches("CROSS")))
  
  if (join_categories) {
    df_cross <- df_cross %>% mutate(
      bin_idx_min_temp_CROSS = case_when(
        bin_idx_min_temp_CROSS %in% 5:6 ~ 5,
        TRUE ~ bin_idx_min_temp_CROSS
      ),
      bin_idx_min_temp_CROSS = match(bin_idx_min_temp_CROSS,
                                     sort(unique(bin_idx_min_temp_CROSS)))
    )
  }
  
  df_cross %>% group_by(min_temp_CROSS) %>%
    summarise(presence = mean(presence, na.rm=TRUE), n=n())
  
  keys$min_temp_CROSS <- compute_bin_keys(df_cross, bin_column="bin_idx_min_temp_CROSS",
                                          original_column="min_temp")
  keys$min_temp_CROSS <-
    keys$min_temp_CROSS %>% rename_with(~ sub('min_temp', 'min_temp_CROSS', .x))
  print(keys$min_temp_CROSS)
  
  return(list(keys=keys, df_cross=df_cross))
}



#  plot_precision_prior_and_posterior.
plot_precision_prior_and_posterior <- function(fitted_model, covariate_name="spde_idx") {
  
  covariate_name_vec <- lapply(fitted_model$all.hyper$random, function(x) {x$hyperid}) %>% unlist()
  i <- which(covariate_name_vec==covariate_name)
  prior <- fitted_model$all.hyper$random[[i]]$hyper
  
  if (prior$theta$prior != "loggamma")
    stop("prior error")
  
  prior_params <- prior$theta$param[1:2]
  
  precision_column <- paste0("Precision for ", covariate_name)
  
  i <- which(names(fitted_model$marginals.hyperpar) == precision_column)
  x_posterior <- fitted_model$marginals.hyperpar[[i]][,1]
  dens_posterior <- fitted_model$marginals.hyperpar[[i]][,2]
  if (is.null(x_posterior))
    stop("marginals error")
  
  tau_posterior <- tibble(x = x_posterior,
                          y = dens_posterior,
                          type = "posterior")
  
  dens_prior <- dgamma(x_posterior, prior_params[1], rate = prior_params[2], log = FALSE)
  tau_prior <- tibble(x = x_posterior,
                      y = dens_prior,
                      type = "prior")
  
  post_tau_mean <- fitted_model$summary.hyperpar[precision_column, "mean"]
  tau_title <- paste0(precision_column, ", post. mean: ", sprintf("%.2f", post_tau_mean))
  
  ps <- list()
  ps[["tau"]] <- plot_scatter_lines(tau_prior,
                                    "x", "y", color_column = "type",
                                    logx = TRUE,
                                    title = tau_title)
  
  ps[["tau_prior"]] <- plot_scatter_lines(tau_posterior,
                                          "x", "y", color_column = "type",
                                          logx = TRUE,
                                          title = tau_title)
  
  p <- ggarrange(ps[["tau"]], ps[["tau_prior"]], ncol=1, nrow=2)
  
  return(list(p=p, tau_prior=tau_prior, tau_posterior=tau_posterior))
}



# evaluate_precision_prior.
evaluate_precision_prior <- function(fitted_model, x, covariate_name="spde_idx") {
  
  covariate_name_vec <- lapply(fitted_model$all.hyper$random, function(x) {x$hyperid}) %>% unlist()
  i <- which(covariate_name_vec==covariate_name)
  prior <- fitted_model$all.hyper$random[[i]]$hyper
  
  if (prior$theta$prior != "loggamma")
    stop("prior error")
  
  prior_params <- prior$theta$param[1:2]
  
  prior_at_x <- dgamma(x, prior_params[1], rate = prior_params[2], log = FALSE)
  
  return(prior_at_x)
}



# compute_p_from_effects.
compute_p_from_effects <- function(map_full_model,
                                   map_random_effects_with_time,
                                   map_random_effects,
                                   map_fixed_effects,
                                   tb_region_results,
                                   results_summary) {
  
  if (!is.null(map_random_effects_with_time)) {
    map_all <- map_full_model %>%
      left_join(
        map_random_effects_with_time
      )
  } else {
    map_all <- map_full_model
  }
  
  if (!is.null(map_random_effects)) {
    map_all <- map_all %>%
      left_join(
        map_random_effects %>% dplyr::select(-data_idx)
      )
  }
  
  if (!is.null(map_fixed_effects)) {
    map_all <- map_all %>%
      left_join(
        map_fixed_effects %>% dplyr::select(-data_idx)
      )
  }
  
  if (nrow(tb_region_results)==nrow(map_all)) {
    map_all <- map_all %>%
      left_join(
        tb_region_results %>% dplyr::select(geo, year_idx, spde_integrated) %>%
          rename(spde.mean=spde_integrated)
      )
  } else {
    stop("not implemented")
  }
  
  effect_columns <- names(map_all) %>% str_subset("mean$")
  
  map_full_model$effect_sum <- rowSums(map_all[,effect_columns]) + results_summary$intercept
  map_full_model$p_from_sum <- exp(map_full_model$effect_sum) / (1 + exp(map_full_model$effect_sum))
  
  return(map_all)
}






# collect_all_model_effects.
collect_all_model_effects <- function(map_full_model, tb_region_results, map_fixed_effects,
                                      map_random_effects, map_random_effects_with_time) {
  map_all <- map_full_model
  
  if (nrow(tb_region_results)>0) {
    contribution_curr <- tb_region_results %>% dplyr::select(region_id, year_idx, spde_integrated)
    if (nrow(contribution_curr) == nrow(map_all)) {
      map_all <- left_join(map_all, contribution_curr)
    } else {
      stop("nrow do not match")
    }
  }
  
  if (!is.null(map_fixed_effects) && nrow(map_fixed_effects)>0) {
    contribution_curr <- map_fixed_effects %>% dplyr::select(region_id, year_idx, matches("feffect.*mean$"))
    if (nrow(contribution_curr) == nrow(map_all)) {
      map_all <- left_join(map_all, contribution_curr)
    } else {
      stop("nrow do not match")
    }
  }
  
  if (!is.null(map_random_effects) && nrow(map_random_effects)>0) {
    contribution_curr <- map_random_effects %>% dplyr::select(region_id, matches("reffect.*mean$"))
    if (nrow(contribution_curr) == nrow(map_all %>% filter(year_idx==min(map_all$year_idx)))) {
      map_all <- left_join(map_all, contribution_curr)
    } else {
      stop("nrow do not match")
    }
  }
  
  if (!is.null(map_random_effects_with_time) && nrow(map_random_effects_with_time)>0) {
    contribution_curr <- map_random_effects_with_time %>%
      dplyr::select(region_id, year_idx, matches("reffect.*mean$"))
    if (nrow(contribution_curr) == nrow(map_all)) {
      map_all <- left_join(map_all, contribution_curr)
    } else {
      stop("nrow do not match")
    }
  }
  
  return(map_all)
}
