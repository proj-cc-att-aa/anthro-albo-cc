# ================== Bin Covariates ==================

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



# compute_bin_keys_pred.
compute_bin_keys_pred <- function(data, data_pred, bin_column, original_column,
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
    data_pred_curr_bin <- data_pred %>% filter(!!sym(bin_column) == curr_bin_idx)
    
    bin_range <- range(data_curr_bin[[original_column]])
    bin_range <- round(bin_range, digits=13)
    keys_to_bins[[bin_range_column]][j] <- sprintf(print_format, bin_range[1], bin_range[2])
    keys_to_bins$bin_min[j] <- bin_range[1]
    keys_to_bins$bin_max[j] <- bin_range[2]
    
    keys_to_bins[[bin_mean_column]][j] <- mean(data_curr_bin[[original_column]])
    
    keys_to_bins$n_in_bin[j] <- nrow(data_pred_curr_bin)
    
  }
  
  keys_to_bins <- keys_to_bins %>% mutate(!!sym(bin_range_column) :=
                                            factor(!!sym(bin_range_column),
                                                   levels=!!sym(bin_range_column), ordered=TRUE))
  
  return(keys_to_bins)
}






## GENRIC BINNING FUNCTIONs

# .compute_keys.
# -- helper: compute keys (range string, min, max, mean, n)
.compute_keys <- function(df, bin_col, value_col, fmt="[%.2f,%.2f]") {
  rng_col  <- paste0("bin_range_", value_col)
  mean_col <- paste0("bin_mean_",  value_col)
  
  tibble(bin_idx = sort(unique(df[[bin_col]]))) %>%
    rowwise() %>%
    mutate(
      vals     = list(df[[value_col]][df[[bin_col]] == bin_idx]),
      bin_min  = suppressWarnings(min(vals, na.rm=TRUE)),
      bin_max  = suppressWarnings(max(vals, na.rm=TRUE)),
      !!rng_col  := if (is.finite(bin_min) && is.finite(bin_max))
        sprintf(fmt, bin_min, bin_max) else NA_character_,
      !!mean_col := if (length(vals) && any(is.finite(vals)))
        mean(vals, na.rm=TRUE) else NA_real_,
      n_in_bin = sum(df[[bin_col]] == bin_idx, na.rm=TRUE)
    ) %>%
    ungroup() %>%
    mutate(!!rng_col := factor(!!sym(rng_col), levels = !!sym(rng_col), ordered = TRUE)) %>%
    dplyr::select(bin_idx, all_of(rng_col), bin_min, bin_max, all_of(mean_col), n_in_bin)
}



#  bin_and_merge_one.
bin_and_merge_one <- function(train_df, pred_df, col, n_bins = 15,
                              right = TRUE, include_lowest = TRUE,
                              merge_spec = NULL, key_fmt = "[%.2f,%.2f]") {
  
  stopifnot(col %in% names(train_df), col %in% names(pred_df))
  stopifnot(is.numeric(train_df[[col]]), is.numeric(pred_df[[col]]))
  
  train_bins <- cut(train_df[[col]], breaks = n_bins, right = right, include.lowest = include_lowest)
  cut_levels <- levels(train_bins)
  
  # get_bounds.
  get_bounds <- function(lbls) {
    lbls <- cut_levels
    cleaned <- gsub("[\\[\\]()]", "", cut_levels, perl = TRUE)
    
    parts <- strsplit(cleaned, ",")
    
    ends <- sort(unique(unlist(parts)))
    
    sort(as.numeric(ends))
  }
  
  breaks <- get_bounds(cut_levels)
  
  if (length(breaks) < 2) {
    rng <- range(train_df[[col]], na.rm = TRUE)
    breaks <- seq(rng[1], rng[2], length.out = n_bins + 1)
  }
  
  train_idx0 <- as.numeric(cut(train_df[[col]], breaks = breaks, right = right, include_lowest = include_lowest))
  
  train_min <- min(breaks, na.rm = TRUE)
  train_max <- max(breaks, na.rm = TRUE)
  
  tmp <- pmax(pred_df[[col]], train_min)
  
  pred_vals <- pmin(tmp, train_max)
  
  pred_idx0  <- as.numeric(cut(pred_vals, breaks = breaks, right = right, include_lowest = include_lowest))
  
  counts <- as.integer(tabulate(train_idx0, nbins = n_bins))
  nonempty <- which(counts > 0)
  empty    <- which(counts == 0)
  
  # Function: rng_lab.
  rng_lab <- function(i) sprintf(key_fmt, breaks[i], breaks[i + 1])
  
  if (length(empty)) {
    cat("\nFollowing uniform division of bins, these are then empty training bins (bins with no data points):\n")
    empties_tbl <- tibble(
      bin_idx_orig = empty,
      range = map_chr(empty, rng_lab)
    )
    print(empties_tbl, n = nrow(empties_tbl))
  } else {
    cat("\nNo empty training bins under uniform division.\n")
  }
  
  remap_to_nearest <- seq_len(n_bins)
  if (length(empty)) {
    # Function: nearest_idx.
    nearest_idx <- function(k) nonempty[which.min(abs(nonempty - k))]
    remap_to_nearest[empty] <- vapply(empty, nearest_idx, integer(1))
  }
  
  pred_idx_nearest <- remap_to_nearest[pred_idx0]
  
  if (length(empty)) {
    moved_counts <- tibble(orig_pred_bin = pred_idx0) |>
      mutate(is_empty = orig_pred_bin %in% empty) |>
      filter(is_empty) |>
      count(orig_pred_bin, name = "n_moved") |>
      mutate(
        from_range = map_chr(orig_pred_bin, rng_lab),
        to_bin     = remap_to_nearest[orig_pred_bin],
        to_range   = map_chr(to_bin, rng_lab)
      ) |>
      dplyr::select(orig_pred_bin, from_range, to_bin, to_range, n_moved)
    
    if (nrow(moved_counts)) {
      cat("\nPrediction data in empty bins will be moved to the nearest non-empty bin:\n")
      print(moved_counts, n = nrow(moved_counts))
    } else {
      cat("\nNo prediction rows fell into empty training bins.\n")
    }
  }
  
  if (!length(nonempty)) stop("All training bins are empty; cannot reindex.")
  compress_map <- integer(length(breaks) - 1)
  compress_map[nonempty] <- seq_along(nonempty)
  
  train_idx <- compress_map[train_idx0]
  pred_idx  <- compress_map[pred_idx_nearest]
  
  
  
  
  # Function: remap_monotone.
  remap_monotone <- function(idx,
                             spec,
                             K,
                             compress = FALSE) {
    
    if (is.null(spec)) return(idx)
    
    # Function: to_int_groups.
    to_int_groups <- function(s) {
      stopifnot(is.list(s))
      lapply(s, function(x) {
        if (is.character(x) && length(x) == 1 && grepl("^\\d+\\s*:\\s*\\d+$", x)) {
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
    lab[] <- 0
    for (g in seq_along(groups)) {
      lab[groups[[g]]] <- -g
    }
    
    new_id   <- integer(K)
    next_id  <- 1
    seen_grp <- logical(length(groups))
    
    for (b in seq_len(K)) {
      if (lab[b] < 0) {
        g <- -lab[b]
        if (!seen_grp[g]) {
          new_id[groups[[g]]] <- next_id
          next_id <- next_id + 1
          seen_grp[g] <- TRUE
        }
      } else if (lab[b] == 0) {
        new_id[b] <- next_id
        next_id <- next_id + 1
      }
    }
    
    out <- new_id[idx]
    
    if (!compress) return(out)
    
    uniq <- sort(na.omit(unique(out)))
    match(out, uniq)
  }
  
  K_full <- length(breaks) - 1
  train_idx_merged <- remap_monotone(train_idx, merge_spec, K = K_full, compress = FALSE)
  pred_idx_merged  <- remap_monotone(pred_idx,  merge_spec, K = K_full, compress = FALSE)
  
  train_tmp <- train_df %>% mutate(!!paste0("bin_idx_", col) := train_idx)
  pred_tmp  <- pred_df  %>% mutate(!!paste0("bin_idx_", col) := pred_idx)
  keys_train <- .compute_keys(train_tmp, bin_col = paste0("bin_idx_", col), value_col = col, fmt = key_fmt)
  keys_pred  <- .compute_keys(pred_tmp,  bin_col = paste0("bin_idx_", col), value_col = col, fmt = key_fmt)
  
  train_tmp_m <- train_df %>% mutate(!!paste0("bin_idx_", col, "_merged") := train_idx_merged)
  pred_tmp_m  <- pred_df  %>% mutate(!!paste0("bin_idx_", col, "_merged") := pred_idx_merged)
  keys_train_merged <- .compute_keys(train_tmp_m, bin_col = paste0("bin_idx_", col, "_merged"), value_col = col, fmt = key_fmt)
  keys_pred_merged  <- .compute_keys(pred_tmp_m,  bin_col = paste0("bin_idx_", col, "_merged"), value_col = col, fmt = key_fmt)
  
  list(
    breaks             = breaks,
    empty_bins         = empty,
    remap_to_nearest   = remap_to_nearest,
    compress_map       = compress_map,
    train_bin          = tibble(bin_idx = train_idx, bin_idx_merged = train_idx_merged),
    pred_bin           = tibble(bin_idx = pred_idx,  bin_idx_merged = pred_idx_merged),
    keys_train         = keys_train,
    keys_pred          = keys_pred,
    keys_train_merged  = keys_train_merged,
    keys_pred_merged   = keys_pred_merged
  )
}




# Function: proximity_binning_aligned.
proximity_binning_aligned <- function(train_df,
                                      pred_df,
                                      proximity_col,
                                      presence_prev_col,
                                      n_bins = 12,
                                      right = TRUE,
                                      include_lowest = TRUE,
                                      key_fmt = "[%.2f,%.2f]",
                                      label_prev = "albo_prev_year",
                                      compress_nonempty = TRUE) {
  
  stopifnot(proximity_col %in% names(train_df), proximity_col %in% names(pred_df))
  stopifnot(presence_prev_col %in% names(train_df), presence_prev_col %in% names(pred_df))
  stopifnot(is.numeric(train_df[[proximity_col]]), is.numeric(pred_df[[proximity_col]]))
  
  # Function: make_bin_col_name.
  make_bin_col_name <- function(x) paste0("bin_idx_", gsub("[^A-Za-z0-9]+", "_", x))
  bin_col        <- make_bin_col_name(proximity_col)
  bin_col_merged <- bin_col
  
  train_core <- train_df %>% filter(.data[[presence_prev_col]] != 1)
  rng        <- range(train_core[[proximity_col]], na.rm = TRUE)
  breaks     <- seq(rng[1], rng[2], length.out = n_bins + 1)
  
  # Function: cut_index.
  cut_index <- function(x) as.integer(cut(x, breaks = breaks, right = right, include_lowest = include_lowest))
  
  train_idx_core <- cut_index(train_df[[proximity_col]])
  pred_idx_core  <- cut_index(pred_df[[proximity_col]])
  
  train_idx <- ifelse(train_df[[presence_prev_col]] == 1, 1L, train_idx_core + 1L)
  pred_idx  <- ifelse(pred_df[[presence_prev_col]]   == 1, 1L, pred_idx_core  + 1L)
  
  compress_map <- NULL
  if (compress_nonempty) {
    used      <- sort(unique(na.omit(train_idx)))
    used_core <- used[used >= 2L]
    if (length(used_core)) {
      new_ids   <- seq(2L, 1L + length(used_core))
      map_core  <- setNames(new_ids, used_core)
      # Function: remap_vec.
      remap_vec <- function(v) {
        v2 <- v
        mask <- !is.na(v2) & v2 >= 2L
        v2[mask] <- unname(map_core[as.character(v2[mask])])
        v2
      }
      train_idx <- remap_vec(train_idx)
      pred_idx  <- remap_vec(pred_idx)
      compress_map <- list(core_map = map_core, reserved_bin = 1L)
    } else {
      compress_map <- list(core_map = integer(0), reserved_bin = 1L)
    }
  } else {
    compress_map <- list(core_map = integer(0), reserved_bin = 1L)
  }
  
  train_out <- train_df %>% mutate(!!bin_col := train_idx)
  pred_out  <- pred_df  %>% mutate(!!bin_col := pred_idx)
  
  keys_train <- .compute_keys(train_out, bin_col = bin_col, value_col = proximity_col, fmt = key_fmt)
  keys_pred  <- .compute_keys(pred_out,  bin_col = bin_col, value_col = proximity_col, fmt = key_fmt)
  
  rng_col <- paste0("bin_range_", proximity_col)
  if (nrow(keys_train)) {
    lv <- levels(keys_train[[rng_col]])
    if (!is.null(lv) && length(lv) >= 1) { lv[1] <- label_prev; levels(keys_train[[rng_col]]) <- lv }
  }
  if (nrow(keys_pred)) {
    lv <- levels(keys_pred[[rng_col]])
    if (!is.null(lv) && length(lv) >= 1) { lv[1] <- label_prev; levels(keys_pred[[rng_col]]) <- lv }
  }
  
  join_cols <- keys_train %>% dplyr::select(bin_idx, matches("bin_range|bin_mean"))
  train_out <- train_out %>% left_join(join_cols, by = join_by(!!sym(bin_col) == bin_idx))
  pred_out  <- pred_out  %>% left_join(join_cols, by = join_by(!!sym(bin_col) == bin_idx))
  
  empty_bins       <- integer(0)
  remap_to_nearest <- NULL
  keys_train_merged <- keys_train
  keys_pred_merged  <- keys_pred
  
  cat("\nProximity binning (aligned) complete.\n")
  cat("Reserved bin 1 for previous-year presence rows.\n")
  cat("Bin column:        ", bin_col, "\n", sep = "")
  cat("Merged bin column: ", bin_col_merged, " (identical)\n", sep = "")
  cat("Train counts:\n"); print(table(train_out[[bin_col]], useNA = "ifany"))
  cat("Pred  counts:\n"); print(table(pred_out[[bin_col]],  useNA = "ifany"))
  
  train_bin <- tibble(
    bin_idx        = train_idx,
    bin_idx_merged = train_idx
  )
  pred_bin <- tibble(
    bin_idx        = pred_idx,
    bin_idx_merged = pred_idx
  )
  
  cat("\nProximity binning (aligned) complete.\n")
  cat("Reserved bin 1 for previous-year presence rows.\n")
  cat("Bin column:        ", bin_col, "\n", sep = "")
  cat("Merged bin column: ", bin_col_merged, " (identical)\n", sep = "")
  cat("Train counts:\n"); print(table(train_out[[bin_col]], useNA = "ifany"))
  cat("Pred  counts:\n"); print(table(pred_out[[bin_col]],  useNA = "ifany"))
  
  list(
    breaks            = breaks,
    empty_bins        = integer(0),
    remap_to_nearest  = NULL,
    compress_map      = compress_map,
    
    bin_col           = bin_col,
    bin_col_merged    = bin_col_merged,
    
    train_bin         = train_bin,
    pred_bin          = pred_bin,
    
    keys_train        = keys_train,
    keys_pred         = keys_pred,
    keys_train_merged = keys_train_merged,
    keys_pred_merged  = keys_pred_merged
  )
}

