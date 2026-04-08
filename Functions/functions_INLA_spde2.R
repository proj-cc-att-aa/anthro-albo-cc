# ================== Inla Spde function used in the model ==================

library(INLA)


packages <- c("sf", "ggpubr")

lapply(packages, library, character.only = TRUE)



# Function: create_grid_covering_map.
create_grid_covering_map <- function(map, n_x, n_y) {
  lims <- compute_boundary_of_sf(map)
  
  x_proj <- seq(lims$xlim[1], lims$xlim[2], length = n_x)
  y_proj <- seq(lims$ylim[1], lims$ylim[2], length = n_y)
  
  mesh_XY <- meshgrid(x_proj, y_proj)
  loc_XY <- cbind(c(mesh_XY$X), c(mesh_XY$Y))
}



# Function: meshgrid.
meshgrid <- function(x,y) {
  m <- length(x)
  n <- length(y)
  X <- matrix(rep(x, each = n), nrow = n)
  Y <- matrix(rep(y, m), nrow = n)
  return(list(X = X, Y = Y))
}



# Function: match_location_to_region.
match_location_to_region <- function(loc, map_idx, do_plot = FALSE) {
  
  if (!"region_id" %in% map_idx)
    map_idx$region_id <- 1:nrow(map_idx)
  
  region_geometries_planar <- st_transform(map_idx$geometry, 2163)
  
  tens <- 10
  n_matches_per_region = rep(0, nrow(map_idx))
  n_projections_to_region = rep(0, nrow(map_idx))
  loc_to_region <- NULL
  if (do_plot)
    has_plotted <- FALSE
  for (j in 1:nrow(loc)) {
    
    point_planar <- st_transform(st_sfc(st_point(loc[j,]), crs = 4326), 2163)
    
    ind <- which(st_intersects(point_planar, region_geometries_planar, sparse = FALSE))
    if (length(ind) == 1) {
      
      n_matches_per_region[map_idx$region_id[ind]] <- n_matches_per_region[map_idx$region_id[ind]] + 1
      
      loc_to_region <- bind_rows(loc_to_region,
                                 tibble(longitude=loc[j,1], latitude=loc[j,2],
                                        loc_id = j, region_id_match = map_idx$region_id[ind],
                                        type = "intersects_region"))
      
      if (do_plot && !has_plotted) {
        map_match <- map_idx %>% mutate(match=region_id==map_idx$region_id[ind])
        plot_sf(map_match, "match")
        
        all_regions <- as_Spatial(map_idx$geometry)
        matching_region <- as_Spatial(map_idx$geometry[ind])
        location <- SpatialPoints(as.matrix(t(loc[j,])))
        
        plot(all_regions)
        plot(matching_region, add=TRUE, col="yellow")
        plot(location, add=TRUE, col = "red", pch = 20, cex = 1.5)
        
        has_plotted <- TRUE
      }
    } else if (length(ind) > 1) {
      
      loc_to_region <- bind_rows(loc_to_region,
                                 tibble(longitude=loc[j,1], latitude=loc[j,2],
                                        loc_id = j, region_id_match = map_idx$region_id[ind],
                                        type = "intersects_multiple_regions"))
      n_matches_per_region[map_idx$region_id[ind]] <- n_matches_per_region[map_idx$region_id[ind]] + 1
    }
    
    progress <- (j/nrow(loc) * 100)
    if (progress >= tens) {
      cat(sprintf("Progress: %.0f%%\n", progress))
      tens <- tens + 10
    }
  }
  
  region_match_stats <- tibble(region_id=map_idx$region_id,
                               n_matches=n_matches_per_region)
  
  return(list(loc_to_region=loc_to_region,
              region_match_stats=region_match_stats))
}




# Function: add_extra_locations_in_regions_with_no_matches.
add_extra_locations_in_regions_with_no_matches <- function(map_idx, loc_to_region, region_match_stats,
                                                           do_plot=FALSE) {
  loc_to_region$location_type <- "original"
  region_match_stats$n_extra_locations <- 0
  if (any(region_match_stats$n_matches == 0)) {
    
    id_no_matches <- region_match_stats %>% filter(n_matches == 0) %>% pull(region_id)
    
    coords <- map_idx %>% filter(region_id %in% id_no_matches) %>% pull(geometry) %>%
      as_Spatial() %>% coordinates()
    
    coords_sampled <- NULL
    for (j in 1:length(id_no_matches)) {
      geometry_projected <- st_transform(map_idx$geometry[id_no_matches[j]], crs = 3578)
      coords_j_projected <- st_sample(geometry_projected,5,type="regular")
      coords_j_lonlat <- st_transform(coords_j_projected, crs = 4326)
      coords_sampled <- rbind(coords_sampled, as.matrix(st_coordinates(coords_j_lonlat)))
      
      if (do_plot && j==1) {
        map_no_match <- map_idx %>% filter(region_id==id_no_matches[j])
        lon_range <- c(map_no_match$longitude - 0.5, map_no_match$longitude + 0.5)
        lat_range <- c(map_no_match$latitude - 0.5, map_no_match$latitude + 0.5)
        
        plot_geometry(map_idx %>% filter(longitude > lon_range[1] & longitude < lon_range[2],
                                         latitude > lat_range[1] & latitude < lat_range[2]))
        plot_geometry(map_no_match, add=TRUE, color="red", opacity = 1)
        
        coords_midpoint <- map_idx %>% filter(region_id==id_no_matches[j]) %>% pull(geometry) %>%
          as_Spatial() %>% coordinates()
        plot(SpatialPoints(coords_midpoint), add=TRUE)
        
        plot(coords_j_lonlat, add=TRUE, pch = 20, cex = 0.95)
      }
    }
    coords <- rbind(coords, coords_sampled)
    
    out <- match_location_to_region(coords, map_idx)
    
    region_match_stats_extra <- out$region_match_stats %>% rename(n_extra_locations=n_matches)
    region_match_stats <- left_join(region_match_stats %>% select(-any_of("n_extra_locations")),
                                    region_match_stats_extra)
    loc_to_region_extra <- out$loc_to_region
    rm("out")
    
    loc_to_region_extra$location_type <- "extra"
    loc_to_region_extra$loc_id <- loc_to_region_extra$loc_id +
      max(loc_to_region$loc_id)
    
    loc_to_region <- bind_rows(loc_to_region, loc_to_region_extra)
    
    n_points_in_region <- region_match_stats$n_matches + region_match_stats$n_extra_locations
    if (any(n_points_in_region == 0)) {
      stop("some regions are still empty")
      
    }
  }
  
  return(list(loc_to_region=loc_to_region,
              region_match_stats=region_match_stats))
}




# Function: extract_spde_results_training.
extract_spde_results_training <- function(fitted_model, response_column, spde_mesh, map_spde, data_spde,
                                          idx_linear_predictor,
                                          A,
                                          includes_time=FALSE, spde_indices=NULL,
                                          order_of_year_of_training_data=NULL,
                                          loc_to_region_spde=NULL,
                                          projection_locations=NULL) {
  
  if (includes_time && is.null(spde_indices))
    stop("if the spde includes time you need to supply spde_indices")
  if (!includes_time) {
    order_of_year_of_training_data <- unique(data_spde$year_idx)
    if (length(order_of_year_of_training_data) != 1)
      stop("if the spde does not include time then length(unique(year_idx)) should equal one")
  } else {
    year_to_plot_vec <- intersect(order_of_year_of_training_data, unique(data_spde$year_idx))
  }
  
  if (includes_time) {
    
    tb_spde_results <- NULL
    
    for (each_year in order_of_year_of_training_data) {
      
      spde_mean <- fitted_model$summary.random$spde_idx$mean[spde_indices$spde_idx.group==each_year]
      spde_sd <- fitted_model$summary.random$spde_idx$sd[spde_indices$spde_idx.group==each_year]
      
      tb_spde_results <- bind_rows(tb_spde_results,
                                   tibble(longitude=spde_mesh$loc[,1],
                                          latitude=spde_mesh$loc[,2],
                                          spde_weights.mean=spde_mean,
                                          spde_weights.sd=spde_sd,
                                          year_idx=each_year))
    }
    
  } else {
    
    spde_mean <- fitted_model$summary.random$spde_idx$mean
    spde_sd <- fitted_model$summary.random$spde_idx$sd
    
    tb_spde_results <- tibble(longitude=spde_mesh$loc[,1],
                              latitude=spde_mesh$loc[,2],
                              spde_weights.mean=spde_mean,
                              spde_weights.sd=spde_sd,
                              year_idx=each_year)
  }
  
  if (includes_time) {
    
    tb_regions <- data_spde %>% dplyr::select(geo, region_id) %>% distinct()
    
    data_spde_with_na <-
      tidyr::crossing(tb_regions, tibble(year_idx=unique(data_spde$year_idx))) %>%
      arrange(year_idx, region_id) %>% left_join(data_spde)
    
    tb_region_results <- NULL
    proj <- inla.mesh.projector(spde_mesh, loc=projection_locations)
    
    for (each_year in order_of_year_of_training_data) {
      
      spde_integrated <-
        inla.mesh.project(
          proj,
          fitted_model$summary.random$spde_idx$mean[spde_indices$spde_idx.group==each_year]
        )
      
      tb_region_results <- bind_rows(tb_region_results,
                                     bind_cols(data_spde_with_na %>% filter(year_idx == each_year) %>%
                                                 dplyr::select(!!response_column, geo, region_id, year_idx),
                                               spde_integrated=spde_integrated))
    }
    
    tb_region_results <- rename(tb_region_results, y = !!response_column)
    
  } else {
    
    tb_region_results <- tibble(y = data_spde[[response_column]])
    proj <- inla.mesh.projector(spde_mesh, loc=projection_locations)
    tb_region_results$spde_integrated <- inla.mesh.project(proj, tb_spde_results$spde_weights.mean)
    tb_region_results <- bind_cols(data_spde %>% dplyr::select(region_id, geo, year_idx), tb_region_results)
  }
  
  if (includes_time) {
    
    idx_subset <- idx_linear_predictor[data_spde$year_idx %in% order_of_year_of_training_data]
  } else {
    
    idx_subset <- idx_linear_predictor
  }
  
  if (includes_time) {
    
    data_sub <- data_spde %>% filter(year_idx %in% order_of_year_of_training_data)
    
    data_sub$val_subset <- fitted_model$summary.linear.predictor[idx_subset, "mean"]
    data_sub_with_na <- data_spde_with_na %>% filter(year_idx %in% order_of_year_of_training_data) %>%
      left_join(data_sub %>% dplyr::select(geo, year_idx, val_subset))
    tb_region_results$linear_predictor.mean <- data_sub_with_na$val_subset
    
    data_sub$val_subset <- fitted_model$summary.linear.predictor[idx_subset, "sd"]
    data_sub_with_na <- data_spde_with_na %>% filter(year_idx %in% order_of_year_of_training_data) %>%
      left_join(data_sub %>% dplyr::select(geo, year_idx, val_subset))
    tb_region_results$linear_predictor.sd <- data_sub_with_na$val_subset
    
    data_sub$val_subset <- fitted_model$summary.fitted.values[idx_subset, "mean"]
    data_sub_with_na <- data_spde_with_na %>% filter(year_idx %in% order_of_year_of_training_data) %>%
      left_join(data_sub %>% dplyr::select(geo, year_idx, val_subset))
    tb_region_results$p.mean <- data_sub_with_na$val_subset
    
    data_sub$val_subset <- fitted_model$summary.fitted.values[idx_subset, "sd"]
    data_sub_with_na <- data_spde_with_na %>% filter(year_idx %in% order_of_year_of_training_data) %>%
      left_join(data_sub %>% dplyr::select(geo, year_idx, val_subset))
    tb_region_results$p.sd <- data_sub_with_na$val_subset
  } else {
    
    tb_region_results$linear_predictor.mean <- fitted_model$summary.linear.predictor[idx_subset, "mean"]
    tb_region_results$linear_predictor.sd <- fitted_model$summary.linear.predictor[idx_subset, "sd"]
    tb_region_results$p.mean <- fitted_model$summary.fitted.values[idx_subset, "mean"]
    tb_region_results$p.sd <- fitted_model$summary.fitted.values[idx_subset, "sd"]
  }
  
  return(list(tb_spde_results=tb_spde_results,
              tb_region_results=tb_region_results))
}





# Function: extract_spde_results_prediction.
extract_spde_results_prediction <- function(fitted_model, response_column,
                                            spde_mesh, map_spde,
                                            data_spde,
                                            pred_data,
                                            idx_linear_predictor_pred,
                                            A_pred = NULL,
                                            includes_time = TRUE,
                                            spde_indices = NULL,
                                            order_of_year_of_prediction_data = NULL,
                                            projection_locations = NULL) {
  
  if (includes_time && is.null(spde_indices))
    stop("if the spde includes time you need to supply spde_indices")
  
  if (!includes_time) {
    order_of_year_of_prediction_data <- unique(pred_data$year_idx)
    if (length(order_of_year_of_prediction_data) != 1)
      stop("if the spde does not include time then length(unique(year_idx)) should equal one")
  }
  
  if (is.null(order_of_year_of_prediction_data)) {
    order_of_year_of_prediction_data <- sort(unique(pred_data$year_idx))
  }
  
  years_train <- sort(unique(data_spde$year_idx))
  K_train     <- length(years_train)
  
  years_pred <- sort(unique(pred_data$year_idx))
  K_pred     <- length(years_pred)
  
  if (includes_time) {
    
    tb_spde_results <- NULL
    
    for (each_year in order_of_year_of_prediction_data) {
      
      k_local <- match(each_year, years_pred)
      if (is.na(k_local)) {
        stop("year_idx = ", each_year,
             " in order_of_year_of_prediction_data not found in pred_data$year_idx")
      }
      
      group_idx <- K_train + k_local
      
      spde_mean <- fitted_model$summary.random$spde_idx$mean[
        spde_indices$spde_idx.group == group_idx
      ]
      spde_sd   <- fitted_model$summary.random$spde_idx$sd[
        spde_indices$spde_idx.group == group_idx
      ]
      
      tb_spde_results <- dplyr::bind_rows(
        tb_spde_results,
        tibble::tibble(
          longitude         = spde_mesh$loc[, 1],
          latitude          = spde_mesh$loc[, 2],
          spde_weights.mean = spde_mean,
          spde_weights.sd   = spde_sd,
          year_idx          = each_year
        )
      )
    }
    
  } else {
    
    spde_mean <- fitted_model$summary.random$spde_idx$mean
    spde_sd   <- fitted_model$summary.random$spde_idx$sd
    
    tb_spde_results <- tibble::tibble(
      longitude         = spde_mesh$loc[, 1],
      latitude          = spde_mesh$loc[, 2],
      spde_weights.mean = spde_mean,
      spde_weights.sd   = spde_sd,
      year_idx          = order_of_year_of_prediction_data[1]
    )
  }
  
  if (is.null(projection_locations)) {
    stop("projection_locations (e.g., region centroids) must be supplied")
  }
  
  tb_regions <- pred_data %>%
    dplyr::select(geo, region_id) %>%
    dplyr::distinct()
  
  data_pred_with_na <- tidyr::crossing(
    tb_regions,
    tibble::tibble(year_idx = order_of_year_of_prediction_data)
  ) %>%
    dplyr::arrange(year_idx, region_id)
  
  proj <- inla.mesh.projector(spde_mesh, loc = projection_locations)
  
  tb_region_results <- NULL
  
  if (includes_time) {
    
    for (each_year in order_of_year_of_prediction_data) {
      
      k_local <- match(each_year, years_pred)
      if (is.na(k_local)) {
        stop("year_idx = ", each_year,
             " in order_of_year_of_prediction_data not found in pred_data$year_idx")
      }
      group_idx <- K_train + k_local
      
      spde_integrated <- inla.mesh.project(
        proj,
        fitted_model$summary.random$spde_idx$mean[
          spde_indices$spde_idx.group == group_idx
        ]
      )
      
      tb_region_results <- dplyr::bind_rows(
        tb_region_results,
        dplyr::bind_cols(
          data_pred_with_na %>%
            dplyr::filter(year_idx == each_year) %>%
            dplyr::select(geo, region_id, year_idx),
          tibble::tibble(
            y              = NA_real_,
            spde_integrated = spde_integrated
          )
        )
      )
    }
    
  } else {
    
    spde_integrated <- inla.mesh.project(
      proj,
      tb_spde_results$spde_weights.mean
    )
    
    tb_region_results <- dplyr::bind_cols(
      data_pred_with_na %>% dplyr::select(geo, region_id, year_idx),
      tibble::tibble(
        y               = NA_real_,
        spde_integrated = spde_integrated
      )
    )
  }
  
  if (includes_time) {
    
    idx_subset <- idx_linear_predictor_pred[
      pred_data$year_idx %in% order_of_year_of_prediction_data
    ]
    
    data_sub <- pred_data %>%
      dplyr::filter(year_idx %in% order_of_year_of_prediction_data)
    
    data_sub$val_subset <- fitted_model$summary.linear.predictor[idx_subset, "mean"]
    data_pred_with_na_lp <- data_pred_with_na %>%
      dplyr::left_join(
        data_sub %>% dplyr::select(geo, year_idx, val_subset),
        by = c("geo", "year_idx")
      )
    tb_region_results$linear_predictor.mean <- data_pred_with_na_lp$val_subset
    
    data_sub$val_subset <- fitted_model$summary.linear.predictor[idx_subset, "sd"]
    data_pred_with_na_lp <- data_pred_with_na %>%
      dplyr::left_join(
        data_sub %>% dplyr::select(geo, year_idx, val_subset),
        by = c("geo", "year_idx")
      )
    tb_region_results$linear_predictor.sd <- data_pred_with_na_lp$val_subset
    
    data_sub$val_subset <- fitted_model$summary.fitted.values[idx_subset, "mean"]
    data_pred_with_na_lp <- data_pred_with_na %>%
      dplyr::left_join(
        data_sub %>% dplyr::select(geo, year_idx, val_subset),
        by = c("geo", "year_idx")
      )
    tb_region_results$p.mean <- data_pred_with_na_lp$val_subset
    
    data_sub$val_subset <- fitted_model$summary.fitted.values[idx_subset, "sd"]
    data_pred_with_na_lp <- data_pred_with_na %>%
      dplyr::left_join(
        data_sub %>% dplyr::select(geo, year_idx, val_subset),
        by = c("geo", "year_idx")
      )
    tb_region_results$p.sd <- data_pred_with_na_lp$val_subset
    
  } else {
    
    idx_subset <- idx_linear_predictor_pred
    
    tb_region_results$linear_predictor.mean <-
      fitted_model$summary.linear.predictor[idx_subset, "mean"]
    tb_region_results$linear_predictor.sd   <-
      fitted_model$summary.linear.predictor[idx_subset, "sd"]
    tb_region_results$p.mean                <-
      fitted_model$summary.fitted.values[idx_subset, "mean"]
    tb_region_results$p.sd                  <-
      fitted_model$summary.fitted.values[idx_subset, "sd"]
  }
  
  return(list(
    tb_spde_results   = tb_spde_results,
    tb_region_results = tb_region_results
  ))
}




# Function: plot_spde_results.
plot_spde_results <- function(tb_spde_results, tb_grid_results, tb_region_results,
                              tb_region_results_pred = NULL) {
  ps_spde <- list()
  ps_spde$weights.mean <- plot_scatter_lines(tb_spde_results,
                                             "longitude", "latitude", "spde_weights.mean",
                                             plot_lines=FALSE,
                                             title="spde, weights, mean")
  ps_spde$weights.sd <- plot_scatter_lines(tb_spde_results, "longitude",
                                           "latitude", "spde_weights.sd",
                                           plot_lines=FALSE,
                                           title="spde, weights, sd")
  
  ps_spde$interp.mean <- plot_scatter_lines(tb_grid_results, "longitude", "latitude", "spde_interp.mean",
                                            plot_lines=FALSE,
                                            title="spde, interpolated to grid, mean")
  ps_spde$interp.sd <- plot_scatter_lines(tb_grid_results, "longitude", "latitude", "spde_interp.sd",
                                          plot_lines=FALSE,
                                          title="spde, interpolated to grid, st.dev.")
  
  ps_spde$integrated.mean <- plot_sf(tb_region_results, "spde_integrated",
                                     title="spde, integrated from grid, mean")
  
  ps_region <- list()
  
  ps_region$spde_integrated.mean <- ps_spde$integrated.mean
  
  ps_region$data <- plot_sf(tb_region_results, "y", title="data")
  
  ps_region$linear_predictor.mean <- plot_sf(tb_region_results, "linear_predictor.mean",
                                             title="full model, linear predictor, mean")
  ps_region$linear_predictor.sd <- plot_sf(tb_region_results, "linear_predictor.sd",
                                           title="full model, linear predictor, st.dev.")
  
  ps_region$p.mean <- plot_sf(tb_region_results, "p.mean",
                              title="full model, probability, mean")
  ps_region$p.sd <- plot_sf(tb_region_results, "p.sd",
                            title="full model, probability, st.dev")
  
  if (!is.null(tb_region_results_pred)) {
    
    tb_region_with_pred <- bind_rows(tb_region_results %>% rename(spde_joined=spde_integrated),
                                     tb_region_results_pred %>% rename(spde_joined=spde_pred))
    
    ps_region_with_pred <- list()
    
    ps_region_with_pred$spde_joined.mean <-
      plot_sf(tb_region_with_pred, "spde_joined", title="spde, joined, mean")
    
    ps_region_with_pred$data <- plot_sf(tb_region_with_pred, "y", title="data")
    
    ps_region_with_pred$linear_predictor.mean <-
      plot_sf(tb_region_with_pred, "linear_predictor.mean", title="joined model, linear predictor, mean")
    ps_region_with_pred$linear_predictor.sd <-
      plot_sf(tb_region_with_pred, "linear_predictor.sd", title="joined model, linear predictor, st.dev.")
    
    ps_region_with_pred$p.mean <-
      plot_sf(tb_region_with_pred, "p.mean", title="joined model, probability, mean")
    ps_region_with_pred$p.sd <-
      plot_sf(tb_region_with_pred, "p.sd", title="joined model, probability, st.dev")
    
  } else {
    
    ps_region_with_pred <- NULL
  }
  
  return(list(ps_spde=ps_spde,
              ps_region=ps_region,
              ps_region_with_pred=ps_region_with_pred))
}



# Function: mvnorm_stats_from_Q.
mvnorm_stats_from_Q <- function(x, Q, rank_Q, tau) {
  
  Q_chol <- chol(Q)
  log_determinant <- sum(log(diag(Q_chol)))*2
  
  n <- length(x)
  log_density_tau1 <- - n/2 * log(2*pi) +
    log_determinant / 2 +
    log(tau) * rank_Q / 2 -
    tau * as.numeric(t(x) %*% Q %*% x) / 2
  
  log_density <- - n/2 * log(2*pi) +
    log_determinant / 2 -
    as.numeric(t(x) %*% Q %*% x) / 2
  
  z <- Q_chol %*% x
  
  norm_z <- norm(z, type="2")
  norm_x <- norm(x, type="2")
  
  return(list(log_density=log_density, log_density_tau1=log_density_tau1,
              norm_z=norm_z, norm_x=norm_x))
}



# Function: prepare_for_spde_model.
prepare_for_spde_model <- function(parameters, response_column = "y", n_spde, data_spde, map_spde, A_train = NULL,
                                   data_na = NULL, A_pred = NULL, n_group = 1, n_group_pred = 1) {
  
  if (n_group == 1) {
    spde_indices_train <- inla.spde.make.index("spde_idx", n.spde = n_spde$full_mesh)
  } else {
    spde_indices_train <- inla.spde.make.index("spde_idx", n.spde = n_spde$full_mesh,
                                               n.group = length(unique(data_spde$year_idx)))
  }
  
  effects_list <- list(Intercept = 1, year_idx = data_spde$year_idx)
  
  if (parameters$precipitation_mean) {
    effects_list$bin_idx_precipitation_mean <- data_spde$bin_idx_precipitation_mean
  }
  if (parameters$precipitation_warm_mean) {
    effects_list$bin_idx_precipitation_warm_mean <- data_spde$bin_idx_precipitation_warm_mean
  }
  if (parameters$precipitation_summer_mean) {
    effects_list$bin_idx_precipitation_summer_mean <- data_spde$bin_idx_precipitation_summer_mean
  }
  if (parameters$precipitation_cold_mean) {
    effects_list$bin_idx_precipitation_cold_mean <- data_spde$bin_idx_precipitation_cold_mean
  }
  if (parameters$temperature_mean) {
    effects_list$bin_idx_temperature_mean <- data_spde$bin_idx_temperature_mean
  }
  if (parameters$temperature_warm_mean) {
    effects_list$bin_idx_temperature_warm_mean <- data_spde$bin_idx_temperature_warm_mean
  }
  if (parameters$temperature_summer_mean) {
    effects_list$bin_idx_temperature_summer_mean <- data_spde$bin_idx_temperature_summer_mean
  }
  if (parameters$temperature_cold_mean) {
    effects_list$bin_idx_temperature_cold_mean <- data_spde$bin_idx_temperature_cold_mean
  }
  if (parameters$relative_humidity_mean) {
    effects_list$bin_idx_relative_humidity_mean <- data_spde$bin_idx_relative_humidity_mean
  }
  if (parameters$relative_humidity_warm_mean) {
    effects_list$bin_idx_relative_humidity_warm_mean <- data_spde$bin_idx_relative_humidity_warm_mean
  }
  if (parameters$relative_humidity_summer_mean) {
    effects_list$bin_idx_relative_humidity_summer_mean <- data_spde$bin_idx_relative_humidity_summer_mean
  }
  if (parameters$relative_humidity_cold_mean) {
    effects_list$bin_idx_relative_humidity_cold_mean <- data_spde$bin_idx_relative_humidity_cold_mean
  }
  if (parameters$population_density) {
    effects_list$bin_idx_population_density <- data_spde$bin_idx_population_density
  }
  if (parameters$proximity) {
    effects_list$bin_idx_proximity050 <- data_spde$bin_idx_proximity050
  }
  if (parameters$presence_previous_year) {
    effects_list$presence_previous_year <- data_spde$presence_previous_year
  }
  if (parameters$year1) {
    effects_list$year1 <- data_spde$year1
  }
  if (parameters$primary_forested_area) {
    effects_list$primary_forested_area <- data_spde$primary_forested_area
  }
  if (parameters$primary_non_forested_area) {
    effects_list$primary_non_forested_area <- data_spde$primary_non_forested_area
  }
  if (parameters$secondary_forested_area) {
    effects_list$secondary_forested_area <- data_spde$secondary_forested_area
  }
  if (parameters$croplands) {
    effects_list$croplands <- data_spde$croplands
  }
  if (parameters$pasture_land) {
    effects_list$pasture_land <- data_spde$pasture_land
  }
  if (parameters$rangelands) {
    effects_list$rangelands <- data_spde$rangelands
  }
  if (parameters$urban_land) {
    effects_list$urban_land <- data_spde$urban_land
  }
  
  stack <- inla.stack(
    data = data_spde[response_column],
    A = list(A_train, 1),
    effects = list(spde_indices_train, effects_list),
    tag = "spde_data"
  )
  
  if (!is.null(data_na) && !is.null(A_pred)) {
    if (n_group_pred == 1) {
      spde_indices_pred <- inla.spde.make.index("spde_idx", n.spde = n_spde$full_mesh)
    } else {
      spde_indices_pred <- inla.spde.make.index("spde_idx", n.spde = n_spde$full_mesh,
                                                n.group = length(unique(data_na$year_idx)))
    }
    
    effects_list_pred <- list(Intercept = 1, year_idx = data_na$year_idx)
    
    if (parameters$precipitation_mean) {
      effects_list_pred$bin_idx_precipitation_mean <- data_na$bin_idx_precipitation_mean
    }
    if (parameters$precipitation_warm_mean) {
      effects_list_pred$bin_idx_precipitation_warm_mean <- data_na$bin_idx_precipitation_warm_mean
    }
    if (parameters$precipitation_summer_mean) {
      effects_list_pred$bin_idx_precipitation_summer_mean <- data_na$bin_idx_precipitation_summer_mean
    }
    if (parameters$precipitation_cold_mean) {
      effects_list_pred$bin_idx_precipitation_cold_mean <- data_na$bin_idx_precipitation_cold_mean
    }
    if (parameters$temperature_mean) {
      effects_list_pred$bin_idx_temperature_mean <- data_na$bin_idx_temperature_mean
    }
    if (parameters$temperature_warm_mean) {
      effects_list_pred$bin_idx_temperature_warm_mean <- data_na$bin_idx_temperature_warm_mean
    }
    if (parameters$temperature_summer_mean) {
      effects_list_pred$bin_idx_temperature_summer_mean <- data_na$bin_idx_temperature_summer_mean
    }
    if (parameters$temperature_cold_mean) {
      effects_list_pred$bin_idx_temperature_cold_mean <- data_na$bin_idx_temperature_cold_mean
    }
    if (parameters$relative_humidity_mean) {
      effects_list_pred$bin_idx_relative_humidity_mean <- data_na$bin_idx_relative_humidity_mean
    }
    if (parameters$relative_humidity_warm_mean) {
      effects_list_pred$bin_idx_relative_humidity_warm_mean <- data_na$bin_idx_relative_humidity_warm_mean
    }
    if (parameters$relative_humidity_summer_mean) {
      effects_list_pred$bin_idx_relative_humidity_summer_mean <- data_na$bin_idx_relative_humidity_summer_mean
    }
    if (parameters$relative_humidity_cold_mean) {
      effects_list_pred$bin_idx_relative_humidity_cold_mean <- data_na$bin_idx_relative_humidity_cold_mean
    }
    if (parameters$population_density) {
      effects_list_pred$bin_idx_population_density <- data_na$bin_idx_population_density
    }
    if (parameters$proximity) {
      effects_list_pred$bin_idx_proximity050 <- data_na$bin_idx_proximity050
    }
    if (parameters$presence_previous_year) {
      effects_list_pred$presence_previous_year <- data_na$presence_previous_year
    }
    if (parameters$year1) {
      effects_list_pred$year1 <- data_na$year1
    }
    if (parameters$primary_forested_area) {
      effects_list_pred$primary_forested_area <- data_na$primary_forested_area
    }
    if (parameters$primary_non_forested_area) {
      effects_list_pred$primary_non_forested_area <- data_na$primary_non_forested_area
    }
    if (parameters$secondary_forested_area) {
      effects_list_pred$secondary_forested_area <- data_na$secondary_forested_area
    }
    if (parameters$croplands) {
      effects_list_pred$croplands <- data_na$croplands
    }
    if (parameters$pasture_land) {
      effects_list_pred$pasture_land <- data_na$pasture_land
    }
    if (parameters$rangelands) {
      effects_list_pred$rangelands <- data_na$rangelands
    }
    if (parameters$urban_land) {
      effects_list_pred$urban_land <- data_na$urban_land
    }
    
    stack_pred <- inla.stack(
      data = list(presence_pred = NA),
      A = list(A_pred, 1),
      effects = list(spde_indices_pred, effects_list_pred),
      tag = "pred_data"
    )
    
    if (nrow(A_pred) != nrow(data_na)) {
      stop("Mismatch between A_pred and pred_data dimensions.")
    }
    
    stack_joined <- inla.stack(stack, stack_pred)
    idx_linear_predictor_pred <- inla.stack.index(stack_joined, "pred_data")$data
  } else {
    stack_joined <- stack
    idx_linear_predictor_pred <- NULL
  }
  
  idx_linear_predictor <- inla.stack.index(stack_joined, "spde_data")$data
  
  link <- rep(NA, length(stack_joined$data$data[[response_column]]))
  ind_NA <- which(is.na(stack_joined$data$data[[response_column]]))
  link[ind_NA] <- 1
  
  return(list(stack_joined = stack_joined,
              idx_linear_predictor = idx_linear_predictor,
              idx_linear_predictor_pred = idx_linear_predictor_pred,
              link = link,
              spde_indices_train = spde_indices_train,
              spde_indices_pred = spde_indices_pred))
}




# Function: create_mesh.
create_mesh <- function(map_spde, mesh_precision="mid") {
  
  n_spde <- list()
  n_spde$geo <- nrow(map_spde)
  
  shapes <- as_Spatial(map_spde$geometry)
  region_coords <- coordinates(shapes)
  
  map_spde$region_id <- 1:nrow(map_spde)
  
  boundary <- inla.nonconvex.hull(region_coords, convex = 1.2)
  lines(boundary, add = FALSE)
  points(region_coords)
  
  if (mesh_precision=="high") {
    spde_mesh <- inla.mesh.2d(loc=region_coords, cutoff=0.3, max.edge=c(0.35,1), boundary=boundary)
  } else if (mesh_precision=="mid") {
    spde_mesh <- inla.mesh.2d(loc=region_coords, cutoff=0.45, max.edge=c(0.45,2), boundary=boundary)
  } else if (mesh_precision=="low") {
    spde_mesh <- inla.mesh.2d(loc=region_coords, cutoff=0.6, max.edge=c(0.6,3), boundary=boundary)
  }
  
  n_spde$full_mesh <- nrow(spde_mesh$loc)
  
  return(list(spde_mesh=spde_mesh, region_coords=region_coords, n_spde=n_spde))
}




# prepare_for_spde_model_new.
prepare_for_spde_model_new <- function(parameters, response_column = "presence",
                                       n_spde, data_spde, map_spde,
                                       A_train = NULL, data_na = NULL, A_pred = NULL,
                                       spde_indices) {
  
  effects_list <- list(Intercept = 1, year_idx = data_spde$year_idx)
  
  # Function: add_if.
  add_if <- function(flag, name, src, lst) {
    if (isTRUE(flag)) lst[[name]] <- src[[name]]
    lst
  }
  
  # Function: add_block.
  add_block <- function(lst, flags, src) {
    for (nm in names(flags)) lst <- add_if(flags[[nm]], nm, src, lst)
    lst
  }
  
  flags <- list(
    bin_idx_precipitation_winter_mean  = parameters$precipitation_winter_mean,
    bin_idx_precipitation_spring_mean  = parameters$precipitation_spring_mean,
    bin_idx_precipitation_summer_mean  = parameters$precipitation_summer_mean,
    bin_idx_precipitation_autumn_mean  = parameters$precipitation_autumn_mean,
    
    bin_idx_temperature_winter_mean    = parameters$temperature_winter_mean,
    bin_idx_temperature_spring_mean    = parameters$temperature_spring_mean,
    bin_idx_temperature_summer_mean    = parameters$temperature_summer_mean,
    bin_idx_temperature_autumn_mean    = parameters$temperature_autumn_mean,
    
    bin_idx_relative_humidity_winter_mean = parameters$relative_humidity_winter_mean,
    bin_idx_relative_humidity_spring_mean = parameters$relative_humidity_spring_mean,
    bin_idx_relative_humidity_summer_mean = parameters$relative_humidity_summer_mean,
    bin_idx_relative_humidity_autumn_mean = parameters$relative_humidity_autumn_mean,
    
    bin_idx_population_density = parameters$population_density,
    
    bin_idx_y_proximity_exp300       = parameters$proximity,
    presence_previous_year     = parameters$presence_previous_year,
    year1                      = parameters$year1,
    
    primary_forested_area      = parameters$primary_forested_area,
    primary_non_forested_area  = parameters$primary_non_forested_area,
    secondary_forested_area    = parameters$secondary_forested_area,
    croplands                  = parameters$croplands,
    pasture_land               = parameters$pasture_land,
    rangelands                 = parameters$rangelands,
    urban_land                 = parameters$urban_land,
    total_flux_in = parameters$total_flux_in
  )
  
  effects_list <- add_block(effects_list, flags, data_spde)
  
  stack <- inla.stack(
    data    = list(presence = data_spde[[response_column]]),
    A       = list(A_train, 1),
    effects = list(spde_indices, effects_list),
    tag     = "spde_data"
  )
  
  stack_pred <- NULL
  stack_joined <- stack
  idx_linear_predictor_pred <- NULL
  
  if (!is.null(data_na) && !is.null(A_pred)) {
    # Function: align_levels.
    align_levels <- function(src_train, src_pred, cols) {
      for (cn in cols) {
        if (cn %in% names(src_train) && cn %in% names(src_pred)) {
          if (is.factor(src_train[[cn]])) {
            src_pred[[cn]] <- factor(src_pred[[cn]], levels = levels(src_train[[cn]]))
          }
        }
      }
      src_pred
    }
    
    bin_cols <- names(flags)[unlist(flags)]
    data_na  <- align_levels(data_spde, data_na, bin_cols)
    
    effects_list_pred <- list(Intercept = 1, year_idx = data_na$year_idx)
    effects_list_pred <- add_block(effects_list_pred, flags, data_na)
    
    stack_pred <- inla.stack(
      data    = list(presence = rep(NA_real_, nrow(data_na))),
      A       = list(A_pred, 1),
      effects = list(spde_indices, effects_list_pred),
      tag     = "pred_data"
    )
    
    if (nrow(A_pred) != nrow(data_na)) stop("Mismatch between A_pred and pred_data dimensions.")
    stack_joined <- inla.stack(stack, stack_pred)
    idx_linear_predictor_pred <- inla.stack.index(stack_joined, "pred_data")$data
    idx_linear_predictor <- inla.stack.index(stack_joined, "spde_data")$data
  } else {
    stack_joined <- stack
    idx_linear_predictor_pred <- NULL
    idx_linear_predictor      <- inla.stack.index(stack_joined, "spde_data")$data
  }
  
  y_all <- inla.stack.data(stack_joined)$presence
  link  <- rep(1, length(y_all))
  
  list(
    stack_joined             = stack_joined,
    idx_linear_predictor     = idx_linear_predictor,
    idx_linear_predictor_pred = idx_linear_predictor_pred,
    link                     = link,
    spde_indices             = spde_indices
  )
}





# Function: run_spde_model2.
run_spde_model2 <- function(parameters, data_spde, map_spde, n_spde,
                            n_years, n_years_pred,
                            f_list, response_column, Q_spde,
                            strategy = "simplified.laplace", n_binom = 1,
                            model_str = "",
                            pred_data = NULL) {
  
  Sys.setenv(OPENBLAS_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1",
             MKL_NUM_THREADS = "1")
  INLA::inla.setOption(num.threads = "8:1")
  
  if (!identical(Q_spde, 0)) {
    x <- as.vector(inla.qsample(n = 1, Q_spde))
    
    if (length(x) == n_spde$full_mesh) {
      tb_sim <- tibble(longitude = spde_mesh$loc[, 1],
                       latitude = spde_mesh$loc[, 2],
                       x = x)
      p <- plot_scatter_lines(tb_sim, "longitude", "latitude", "x",
                              plot_lines = FALSE)
    } else if (length(x) == n_spde$geo) {
      tb_sim <- bind_cols(map_spde %>% dplyr::select(region_id, geo, geometry), x = x)
      p <- plot_sf(tb_sim, "x")
    }
    
  }
  
  shapes <- as_Spatial(map_spde$geometry)
  region_coords <- sp::coordinates(shapes)
  
  coords_df <- tibble::tibble(
    geo = map_spde$geo,
    lon = region_coords[, 1],
    lat = region_coords[, 2]
  )
  
  tb_coords <- data_spde %>%
    dplyr::select(data_idx, geo, year_idx) %>%
    dplyr::left_join(coords_df, by = "geo")
  
  stopifnot(!anyNA(tb_coords$lon), !anyNA(tb_coords$lat))
  loc_train <- as.matrix(tb_coords[, c("lon","lat")])
  
  if (is.null(pred_data)) {
    
    years_train <- sort(unique(data_spde$year_idx))
    K_train     <- length(years_train)
    
    g_train <- if (K_train == 1L) {
      rep(1L, nrow(tb_coords))
    } else {
      as.integer(factor(tb_coords$year_idx, levels = years_train))
    }
    
    if (K_train == 1L) {
      A_train <- inla.spde.make.A(mesh = spde_mesh, loc = loc_train)
    } else {
      A_train <- inla.spde.make.A(
        mesh    = spde_mesh,
        loc     = loc_train,
        group   = g_train,
        n.group = K_train
      )
    }
    
    print(dim(A_train))
    
    spde_indices <- inla.spde.make.index(
      "spde_idx",
      n.spde  = n_spde$full_mesh,
      n.group = K_train
    )
    
    A_pred  <- NULL
    pred_df <- NULL
    
  } else {
    
    years_train <- sort(unique(data_spde$year_idx))
    K_train     <- length(years_train)
    
    g_train <- if (K_train == 1L) {
      rep(1L, nrow(tb_coords))
    } else {
      as.integer(factor(tb_coords$year_idx, levels = years_train))
    }
    
    tb_coords_pred <- pred_data %>%
      dplyr::select(data_idx, geo, year_idx) %>%
      dplyr::left_join(coords_df, by = "geo")
    
    stopifnot(!anyNA(tb_coords_pred$lon), !anyNA(tb_coords_pred$lat))
    loc_pred <- as.matrix(tb_coords_pred[, c("lon","lat")])
    
    years_pred <- sort(unique(pred_data$year_idx))
    K_pred     <- length(years_pred)
    
    if (K_pred == 1L) {
      g_pred <- rep(K_train + 1L, nrow(tb_coords_pred))
    } else {
      g_pred_local <- as.integer(factor(tb_coords_pred$year_idx, levels = years_pred))
      g_pred       <- K_train + g_pred_local
    }
    
    K_total <- K_train + K_pred
    
    spde_indices <- inla.spde.make.index(
      "spde_idx",
      n.spde  = n_spde$full_mesh,
      n.group = K_total
    )
    
    if (K_total == 1L) {
      A_train <- inla.spde.make.A(mesh = spde_mesh, loc = loc_train)
      A_pred  <- inla.spde.make.A(mesh = spde_mesh, loc = loc_pred)
    } else {
      A_train <- inla.spde.make.A(
        mesh    = spde_mesh,
        loc     = loc_train,
        group   = g_train,
        n.group = K_total
      )
      
      A_pred <- inla.spde.make.A(
        mesh    = spde_mesh,
        loc     = loc_pred,
        group   = g_pred,
        n.group = K_total
      )
    }
    
    print(dim(A_train))
    print(dim(A_pred))
    
    pred_df <- pred_data
  }
  
  stopifnot(nrow(A_train) == nrow(tb_coords))
  stopifnot(Matrix::rowSums(abs(A_train)) |> min() > 0)
  
  out <- prepare_for_spde_model_new(
    parameters      = parameters,
    response_column = response_column,
    n_spde          = n_spde,
    data_spde       = data_spde,
    map_spde        = map_spde,
    A_train         = A_train,
    A_pred          = A_pred,
    data_na         = pred_df,
    spde_indices    = spde_indices
  )
  
  stack_joined             <- out$stack_joined
  idx_linear_predictor     <- out$idx_linear_predictor
  idx_linear_predictor_pred <- out$idx_linear_predictor_pred
  link                     <- out$link
  
  start_time <- Sys.time()
  
  extraconst <- f_list$extraconstr
  
  parallel::detectCores(logical = FALSE)
  
  Sys.setenv(OPENBLAS_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1",
             MKL_NUM_THREADS = "1")
  
  INLA::inla.setOption(num.threads = "8:1")
  
  fitted_model <- tryCatch(
    {
      inla(formula = f_list$f, data = inla.stack.data(stack_joined),
           family = "binomial", Ntrials = n_binom,
           control.family = list(link = "logit"),
           control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE),
           
           control.predictor = list(A = inla.stack.A(stack_joined),
                                    link = link,
                                    compute = TRUE),
           
           control.inla = list(strategy = strategy, int.strategy = "auto", tolerance = 1e-4),
           
           verbose = F)
    },
    error = function(e) {
      message('inla call failed.')
      print(e)
      return(NULL)
    }
  )
  
  end_time <- Sys.time()
  computation_time <- end_time - start_time
  
  summary(fitted_model)
  
  return(list(
    fitted_model             = fitted_model,
    A_train                  = A_train,
    A_pred                   = A_pred,
    idx_linear_predictor     = idx_linear_predictor,
    idx_linear_predictor_pred = idx_linear_predictor_pred,
    spde_indices             = spde_indices,
    computation_time         = computation_time
  ))
}




# Function: compute_correlation_with_data.
compute_correlation_with_data <- function(tb_region_results, data_spde) {
  tb_all <- left_join(tb_region_results, data_spde)
  tb_select <- tb_all %>% filter(!is.na(presence)) %>% dplyr::select(spde_integrated, presence)
  tb_select$spde_integrated <- tb_select$spde_integrated - mean(tb_select$spde_integrated)
  tb_select$presence <- tb_select$presence - mean(tb_select$presence)
  corr_spde_with_data <- cor(tb_select$spde_integrated, tb_select$presence)
  return(corr_spde_with_data)
}
