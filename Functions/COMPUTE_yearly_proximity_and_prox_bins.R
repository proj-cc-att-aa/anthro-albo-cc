# ================== Yearly proximity covariates and bin assignment ==================

source("functions_general.R")
source("functions_plotting.R")
source("functions_risk_from_proximity.R")
source("functions_spatial_networks.R")
source("functions_matrix_operations.R")



# Compute proximity from posterior generated from each year prediction
compute_proximity_from_posterior <- function(p_post_pred,
                                             map_inla,
                                             exponent_vec = c(0.5, 1, 2, 3)) {
  
  tb_data <- p_post_pred %>%
    dplyr::group_by(geo, region_id, year_idx) %>%
    dplyr::summarise(
      p_mean = mean(p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      presence = as.numeric(p_mean >= 0.5),         ## threshold
      y_to_convolve = dplyr::if_else(is.na(presence), 0, presence)
    )
  
  geo_neighs <- geo_to_neighs(map_inla)
  
  Q_besag1 <- compute_Q_besag_centroid_dist_weights(
    geo_neighs$nb_map,
    geo_neighs$coords,
    geo_neighs$n_geo,
    proper_version = FALSE,
    lambda = NULL
  )
  
  Q_with_eig <- collect_Q_and_eig(Q_besag1)
  
  y_column <- "y_to_convolve"
  
  df_proximity <- proximity_from_Q(
    Q_with_eig,
    tb_data,
    y_column,
    map_inla,
    exponent_vec
  )
  
  df_proximity <- df_proximity %>%
    dplyr::rename(presence_previous_year = y)
  
  drop_years <- sort(unique(tb_data$year_idx))
  
  df_proximity <- df_proximity %>%
    dplyr::filter(!year_idx %in% drop_years)
  
  df_surrounded <- get_surrounded_regions()
  
  proximity_col_vec <- names(df_proximity) %>% stringr::str_subset("y_proximity")
  
  for (j in seq_len(nrow(df_proximity))) {
    df_curr <- df_proximity[j, ]
    if (df_curr$geo %in% df_surrounded$inner) {
      geo_outer <- df_surrounded %>%
        dplyr::filter(inner == df_curr$geo) %>%
        dplyr::pull(outer)
      df_curr_match <- df_proximity %>%
        dplyr::filter(geo == geo_outer, year_idx == df_curr$year_idx)
      if (df_curr$presence_previous_year != 1 &&
          df_curr_match$presence_previous_year != 1) {
        for (proximity_col in proximity_col_vec) {
          df_proximity[[proximity_col]][j] <- df_curr_match[[proximity_col]]
        }
      }
    }
  }
  
  return(df_proximity)
}



# Attach proximity bins from spec (yearly assignment to proximity bins)
attach_proximity_bins_from_spec <- function(new_df, bin_spec) {
  
  if (is.null(bin_spec$merge_spec)) {
    stop("This helper assumes merged bins (merge_spec must not be NULL).")
  }
  
  col <- bin_spec$col
  stopifnot(col %in% names(new_df))
  vals <- new_df[[col]]
  stopifnot(is.numeric(vals))
  
  vmin <- min(bin_spec$breaks, na.rm = TRUE)
  vmax <- max(bin_spec$breaks, na.rm = TRUE)
  vals_clamped <- pmin(pmax(vals, vmin), vmax)
  
  idx0 <- as.numeric(cut(
    vals_clamped,
    breaks        = bin_spec$breaks,
    right         = bin_spec$right,
    include.lowest = bin_spec$include_lowest
  ))
  
  idx_nearest <- bin_spec$remap_to_nearest[idx0]
  
  idx_compact <- bin_spec$compress_map[idx_nearest]
  
  K_full <- length(bin_spec$breaks) - 1L
  idx_merged <- remap_monotone(
    idx      = idx_compact,
    spec     = bin_spec$merge_spec,
    K        = K_full,
    compress = FALSE
  )
  
  idx_col <- paste0("bin_idx_", col)
  
  new_df2 <- new_df %>%
    dplyr::mutate(!!idx_col := idx_merged)
  
  message(sprintf(
    "[%s] MERGED: wrote %s (no range/mean columns).",
    col, idx_col
  ))
  
  new_df2
}



# Remap monotone
# Take original bin indices (1..K) and a set of "merge rules" (spec), and return new bin indices where some of the original bins are merge into larger groups, while:
##     - preserving the left-to-right (monotone) order of bins,
##     - producing consecutive integer IDs (1,2,3,...) for the merged bins.

remap_monotone <- function(idx, spec, K, compress = FALSE) {
  
  if (is.null(spec)) return(idx)
  
  stopifnot(is.list(spec))
  
  # To int groups.
  to_int_groups <- function(s) {
    lapply(s, function(x) {
      if (is.character(x) && length(x) == 1L &&
          grepl("^\\d+\\s*:\\s*\\d+$", x)) {
        a <- as.integer(sub(":.*", "", x))
        b <- as.integer(sub(".*:", "", x))
        a:b
      } else {
        as.integer(x)
      }
    })
  }
  
  groups <- to_int_groups(spec)
  
  lab <- integer(K)
  lab[] <- 0L
  for (g in seq_along(groups)) {
    lab[groups[[g]]] <- -g
  }
  
  new_id   <- integer(K)
  next_id  <- 1L
  seen_grp <- logical(length(groups))
  
  for (b in seq_len(K)) {
    if (lab[b] < 0L) {
      g <- -lab[b]
      if (!seen_grp[g]) {
        new_id[groups[[g]]] <- next_id
        next_id <- next_id + 1L
        seen_grp[g] <- TRUE
      }
    } else if (lab[b] == 0L) {
      new_id[b] <- next_id
      next_id   <- next_id + 1L
    }
  }
  
  out <- new_id[idx]
  
  if (!compress) return(out)
  uniq <- sort(na.omit(unique(out)))
  match(out, uniq)
}
