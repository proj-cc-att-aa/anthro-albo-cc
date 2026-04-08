# ================== Spatial network construction functions ==================

# Geo to neighs
geo_to_neighs <- function(map, use_snap=TRUE, use_queen=TRUE) {
  
  if ("sfc" %in% class(map)) {
    
    warning("creating data.frame from sfc input")
    
    map_geometry <- map
    map <- tibble(geo_idx=1:length(map_geometry), geometry=map_geometry)
  } else if (!"sfc" %in% class(map$geometry)) {
    
    stop("no valid geometry column")
  }
  
  if (!"geo_idx" %in% colnames(map)) {
    
    map <- map %>% dplyr::select(geometry) %>% mutate(geo_idx=1:nrow(map))
  }
  
  n_geo <- map$geometry %>% unique() %>% length()
  if (n_geo != nrow(map)) {
    stop('duplicate geometries in map')
  }
  
  shapes <- as_Spatial(map$geometry)
  
  if (use_snap) {
    
    nb_map <- poly2nb(shapes, snap=sqrt(.Machine$double.eps)*5e5,
                      queen = use_queen)
  } else {
    
    nb_map <- poly2nb(shapes, queen = use_queen)
  }
  
  coords <- coordinates(shapes)
  
  return(list(geo_idx = map$geo_idx,
              shapes = shapes,
              nb_map = nb_map,
              n_geo = n_geo,
              coords = coords))
}



# Cross variables (optionals)
cross_variables <- function(df, x_col, y_col, var_name) {
  
  df <- df %>% as_tibble() %>% dplyr::select(data_idx, spatial_idx, year_idx, all_of(c(x_col, y_col)))
  
  cross_var_str <- paste0(var_name, "_CROSS")
  df <- df %>% mutate(!!sym(cross_var_str) := paste0(!!sym(x_col), "x", !!sym(y_col)))
  df %>% pull(sym(cross_var_str)) %>% table() %>% print()
  
  cross_var_num_str <- paste0(cross_var_str, "Num")
  combinations <- sapply(unique(df[[x_col]]) %>% sort(),
                         function(x) {paste0(x, "x",
                                             unique(df[[y_col]]) %>% sort())}) %>%
    unlist() %>% as.vector()
  data_cross_combinations <- tibble(!!sym(cross_var_str) := combinations) %>%
    add_column(cross_join(tibble(x=sort(unique(df[[x_col]]))),
                          tibble(y=sort(unique(df[[y_col]])))))
  
  df <- left_join(df, data_cross_combinations)
  data_cross <- df %>% dplyr::select(data_idx, spatial_idx, year_idx, !!sym(x_col), !!sym(y_col),
                                     !!sym(cross_var_str))
  
  n1 <- length(unique(data_cross[[x_col]]))
  n2 <- length(unique(data_cross[[y_col]]))
  nb_grid <- cell2nb(nrow=n1,
                     ncol=n2)
  coords <- tibble(x = as.numeric(as.factor(data_cross_combinations[["x"]])),
                   y = as.numeric(as.factor(data_cross_combinations[["y"]])),
                   value = data_cross_combinations[[cross_var_str]],
                   x_str = x_col,
                   y_str = y_col,
                   x_orig = data_cross_combinations[["x"]],
                   y_orig = data_cross_combinations[["y"]]) %>%
    distinct() %>% arrange(x, y)
  coords$value <- factor(coords$value, level=coords$value, ordered=TRUE)
  
  data_cross[[cross_var_str]] <- factor(data_cross[[cross_var_str]], level=coords$value, ordered=TRUE)
  data_cross[[cross_var_num_str]] <- match(data_cross[[cross_var_str]], sort(unique(data_cross[[cross_var_str]])))
  
  tb_n_key <- data_cross %>% group_by(!!sym(cross_var_str)) %>% summarise(n=n())
  coords <- left_join(coords, tb_n_key, by=c("value" = cross_var_str))
  
  plot(nb_grid, coords %>% dplyr::select(x,y) %>% as.matrix())
  
  return(list(data_cross=data_cross, coords=coords, nb_grid=nb_grid))
}




# # Control network connections (useful for crossed variables) --------------------------
specify_network_by_nodes_and_links <- function(n_nodes,
                                               nodes_to_disconnect=NULL,
                                               links_to_remove=NULL,
                                               links_to_add=NULL,
                                               save_to_file=NULL) {
  
  nb_grid <- cell2nb(nrow=1, ncol=n_nodes)
  if (!is.null(nodes_to_disconnect)) {
    for (j in 1:n_nodes) {
      if (j %in% nodes_to_disconnect) {
        nb_grid[[j]] <- as.integer(0)
      } else {
        nb_grid[[j]] <- setdiff(nb_grid[[j]], nodes_to_disconnect)
      }
    }
  }
  
  if (!is.null(links_to_remove)) {
    for (j in 1:length(links_to_remove)) {
      link_to_disconnect <- links_to_remove[[j]]
      
      node1 <- link_to_disconnect[1]
      node2 <- link_to_disconnect[2]
      nb_grid[[node1]] <- setdiff(nb_grid[[node1]], node2)
      nb_grid[[node2]] <- setdiff(nb_grid[[node2]], node1)
    }
  }
  
  if (!is.null(links_to_add)) {
    for (j in 1:length(links_to_add)) {
      nodes_j <- links_to_add[[j]]
      
      for (k in nodes_j) {
        nb_grid[[k]] <- c(nb_grid[[k]], as.integer(setdiff(nodes_j, k)))
      }
    }
  }
  
  coords <- tibble(x=1, y=1:n_nodes)
  coords$x[coords$y %in% nodes_to_disconnect] <- 2
  for (j in 1:length(links_to_add)) {
    nodes_j <- links_to_add[[j]]
    coords$x[coords$y %in% nodes_j] <- coords$x[coords$y %in% nodes_j] - j*0.25
  }
  plot(nb_grid, coords %>% as.matrix())
  
  if (!is.null(save_to_file))
    nb2INLA(save_to_file, nb_grid)
  
  return(nb_grid)
}
