# ================== functions to compute risk From Proximity ==================

packages <- c("igraph", "geosphere")

lapply(packages, library, character.only = TRUE)


# mobility from population (radiation model)
mobility_from_population <- function(map, population_column="population") {
  
  geo_neighs <- geo_to_neighs(map)
  
  mobility <-
    list(tb_regions = tibble(geo = map$geo,
                             mobility_idx = 1:nrow(map),
                             population = map[[population_column]],
                             coords = geo_neighs$coords,
                             flux_in = NULL,
                             flux_out = NULL),
         
         flux_matrix = matrix(nrow = geo_neighs$n_geo,
                              ncol = geo_neighs$n_geo,
                              dimnames =
                                list(sapply(1:geo_neighs$n_geo,
                                            function(i) sprintf("fr_i%d", i)),
                                     sapply(1:geo_neighs$n_geo,
                                            function(i) sprintf("to_i%d", i)))))
  
  geo_distances <- distm(geo_neighs$coords, fun = distGeo)
  
  flux_prob <- mobility$flux_matrix
  flux_prob_rescaled <- mobility$flux_matrix
  tot_pop <- sum(map[[population_column]])
  
  for (i_curr_region in 1:geo_neighs$n_geo) {
    
    dist_curr <- geo_distances[i_curr_region,]
    
    idx_order <- order(dist_curr)
    
    pop_ordered <- map[[population_column]][idx_order]
    
    if (dist_curr[idx_order[2]] == 0)
      stop("dist error, duplicate geos")
    
    if (any(is.na(dist_curr)))
      stop("dist error, na")
    
    s_pop_in_circle <- rep(0, length.out=geo_neighs$n_geo)
    
    s_pop_in_circle[3:geo_neighs$n_geo] <-
      cumsum(pop_ordered[2:(geo_neighs$n_geo-1)])
    
    s_pop_in_circle_original_order <- numeric(length(s_pop_in_circle))
    
    s_pop_in_circle_original_order[idx_order] <- s_pop_in_circle
    
    for (j in 1:geo_neighs$n_geo) {
      if (i_curr_region == j)
        next
      m = map[[population_column]][i_curr_region]
      n = map[[population_column]][j]
      s = s_pop_in_circle_original_order[j]
      if ((m + n + s) == 0 || (m + s) == 0) {
        flux_prob[i_curr_region, j] <- 0
        flux_prob_rescaled[i_curr_region, j] <- 0
      } else {
        flux_prob[i_curr_region, j] <- m * n / ((m + s) * (m + n + s))
        flux_prob_rescaled[i_curr_region, j] <- 1 / (1 - m / tot_pop) *
          flux_prob[i_curr_region, j]
      }
      mobility$flux_matrix[i_curr_region, j] <-
        0.9 * m * flux_prob[i_curr_region, j]
      if (is.nan(mobility$flux_matrix[i_curr_region, j]))
        stop("nan error")
    }
  }
  
  mobility$tb_regions$flux_out <-
    rowSums(mobility$flux_matrix, na.rm = TRUE)
  
  mobility$tb_regions$flux_in <-
    colSums(mobility$flux_matrix, na.rm = TRUE)
  
  return(mobility)
}



# flux matrix to tb
flux_matrix_to_tb <- function(flux_matrix, type="directed") {
  
  if (!type %in% c("directed", "undirected"))
    stop("type should equal either 'directed' or 'undirected'")
  
  n_nodes <- nrow(flux_matrix)
  
  adj_matrix <- matrix(1, nrow = n_nodes, ncol = n_nodes)
  diag(adj_matrix) <- 0
  
  if (type=="undirected") {
    adj_matrix[lower.tri(adj_matrix)] <- 0
    
    flux_matrix <- pmax(flux_matrix, t(flux_matrix))
  }
  
  flux_tb <- tibble(link = which(adj_matrix != 0, arr.ind = T))
  
  flux_tb$flux <- flux_matrix[flux_tb$link]
  
  return(flux_tb)
}



# plot mobility graph
plot_mobility_graph <- function(mobility,
                                nodes = "",
                                link_color="tomato",
                                sc_link_width=10, sc_node_size=6,
                                filtering_p=0.99,
                                title="") {
  
  n_nodes <- nrow(mobility$tb_regions)
  
  flux_tb <- flux_matrix_to_tb(mobility$flux_matrix, type="directed")
  
  netw <- tibble(links = flux_tb$link,
                 link_weights = flux_tb$flux,
                 link_color = link_color)
  
  quant <- quantile(netw$link_weights, filtering_p)
  netw_subset <- netw %>% filter(link_weights > quant)
  cat(sprintf("filtering links, only showing links with weights > %.10g%%-percentile\n",
              filtering_p))
  cat(sprintf("-> showing %d links (total %d links)\n", nrow(netw_subset), nrow(netw)))
  
  if (title != "") {
    title <- paste0(title, ": ")
  }
  title <- paste0(title, "showing links with size > ",
                  sprintf("%.10g", filtering_p*100), "%-percentile")
  
  node_weights <- NULL
  if (nodes == "net_flux_abs") {
    mobility$tb_regions <- mobility$tb_regions %>%
      mutate(net_flux = abs(flux_out - flux_in))
    
    node_weights <- mobility$tb_regions$net_flux
    
  } else if (nodes == "population") {
    node_weights <- mobility$tb_regions$population
  }
  
  plot_network_mobility(
    links = netw_subset$links,
    link_weights = netw_subset$link_weights,
    n_nodes = n_nodes,
    node_coords = mobility$tb_regions$coords,
    node_weights = node_weights,
    sc_link_width = sc_link_width,
    sc_node_size = sc_node_size,
    link_color = link_color,
    title = title)
  
}




# plot Q network
plot_Q_network <- function(Q, coords,
                           link_color="tomato",
                           sc_link_width=2, sc_node_size=6,
                           title="",
                           filtering_p=NULL) {
  
  if ((nrow(Q) != nrow(coords)) || (nrow(Q) != ncol(Q)))
    stop("dimension error, Q should be an n x n-matrix and coords should contain the corresponding n locations")
  
  n_nodes <- nrow(Q)
  
  node_weights <- diag(Q)
  
  diag(Q) <- 0
  Q[lower.tri(Q)] <- 0
  netw <- tibble(links=which(Q != 0, arr.ind = T),
                 link_weights=Q[links])
  
  if (!is.null(filtering_p)) {
    quant <- quantile(netw$link_weights, filtering_p)
    n_links_tot <- nrow(netw)
    netw <- netw %>% filter(link_weights >= quant)
    cat(sprintf("filtering links, only showing links with weights > %.10g%%-quantile\n",
                filtering_p))
    cat(sprintf("-> showing %d links (total %d links)\n", nrow(netw), n_links_tot))
  }
  
  plot_network(links=netw$links,
               link_weights=netw$link_weights,
               n_nodes=n_nodes,
               node_coords=coords,
               node_weights=node_weights,
               sc_link_width = sc_link_width, sc_node_size = sc_node_size,
               link_color=link_color,
               title = title)
  
  return(netw)
}




# plot Q network1
plot_Q_network1 <- function(Q, coords,
                            link_color="tomato",
                            sc_link_width=2, sc_node_size=6,
                            title="",
                            filtering_p=NULL, edge_colour_by_link_weights = FALSE) {
  
  if ((nrow(Q) != nrow(coords)) || (nrow(Q) != ncol(Q)))
    stop("dimension error, Q should be an n x n-matrix and coords should contain the corresponding n locations")
  
  n_nodes <- nrow(Q)
  
  node_weights <- diag(Q)
  
  diag(Q) <- 0
  Q[lower.tri(Q)] <- 0
  netw <- tibble(links=which(Q != 0, arr.ind = T),
                 link_weights=Q[links])
  
  if (!is.null(filtering_p)) {
    quant <- quantile(netw$link_weights, filtering_p)
    n_links_tot <- nrow(netw)
    netw <- netw %>% filter(link_weights >= quant)
    cat(sprintf("filtering links, only showing links with weights > %.10g%%-quantile\n",
                filtering_p))
    cat(sprintf("-> showing %d links (total %d links)\n", nrow(netw), n_links_tot))
  }
  
  plot_network1(links=netw$links,
                link_weights=netw$link_weights,
                n_nodes=n_nodes,
                node_coords=coords,
                node_weights=node_weights,
                sc_link_width = sc_link_width, sc_node_size = sc_node_size,
                link_color=link_color,
                title = title,
                edge_colour_by_link_weights = edge_colour_by_link_weights)
  
  return(netw)
}




# subset links by map
subset_links_by_map <- function(df_links, map) {
  df_links_subset <- df_links %>% filter(from_geo %in% map$geo,
                                         to_geo %in% map$geo)
  
  return(df_links_subset)
}



# links to Q (optional)
links_to_Q <- function(df_links, n_nodes) {
  
  df_links$node_start <- df_links$link[,1]
  df_links$node_end <- df_links$link[,2]
  df_links <- df_links %>% dplyr::select(-link)
  
  Q <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  for (j in 1:nrow(df_links)) {
    Q[df_links$node_start[j], df_links$node_end[j]] <- df_links$flux[j]
    Q[df_links$node_end[j], df_links$node_start[j]] <- df_links$flux[j]
  }
  
  diag(Q) <- -rowSums(Q)
  
  return(Q)
}




# compute Q for ICAR based neighbourhood structure
compute_Q_besag <- function(nb_map, n_geo, proper_version = FALSE,
                            lambda = NULL) {
  
  if (proper_version) {
    if (is.null(lambda))
      stop("lambda parameter should be supplied for proper version of Besag")
    if (!(lambda > 0 && lambda < 1))
      stop("invalid lambda-value")
  }
  if (!proper_version && !is.null(lambda))
    stop("lambda parameter is not required when proper_version == FALSE")
  
  if (n_geo != length(nb_map))
    stop("n_geo error")
  neighs <- vector(mode = "list", length = n_geo)
  n_neighs <- vector(length = n_geo)
  
  Q = matrix(0, nrow = n_geo, ncol = n_geo)
  
  for (j in seq_along(nb_map)) {
    ns <- nb_map[[j]]
    if (!all(ns == 0)) {
      
      neighs[[j]] <- ns
      n_neighs[j] <- length(ns)
      
      Q[neighs[[j]],j] = -1
      Q[j,neighs[[j]]] = -1
      Q[j,j] = n_neighs[j]
    }
  }
  if (proper_version) {
    Q <- lambda * Q
    for (i in 1:n_geo)
      Q[i,i] <- Q[i,i] + (1-lambda)
  }
  
  return(Q)
}



# compute Q with centroid dist weights neighbourhood structure
compute_Q_besag_centroid_dist_weights <- function(nb_map, coords, n_geo, proper_version = FALSE, lambda = NULL) {
  
  if (proper_version) {
    if (is.null(lambda))
      stop("lambda parameter should be supplied for proper version of Besag")
    if (!(lambda > 0 && lambda < 1))
      stop("invalid lambda-value")
  }
  if (!proper_version && !is.null(lambda))
    stop("lambda parameter is not required when proper_version == FALSE")
  
  if (n_geo != length(nb_map))
    stop("n_geo error")
  
  geo_distances <- distm(coords, fun = distGeo)
  
  all_nb_dists <- unlist(lapply(seq_along(nb_map), function(i) geo_distances[i, nb_map[[i]]]))
  summary(all_nb_dists)
  
  L <- max(all_nb_dists)
  
  geo_distances <- geo_distances/L
  
  kappa <- 1
  
  tau <- 0.5
  
  Q <- matrix(0, nrow = n_geo, ncol = n_geo)
  
  for (j in seq_along(nb_map)) {
    neighbors <- nb_map[[j]]
    deg_j <- length(neighbors)
    diag_sum <- ifelse(deg_j == 1, 0.5, 0)
    
    if (deg_j > 0 && !all(neighbors == 0)) {
      for (i in neighbors) {
        ell <- geo_distances[j, i]
        off_diag_weight <- -2 * kappa * tau * (exp(-kappa * ell) / (1 - exp(-2 * kappa * ell)))
        
        Q[i, j] <- off_diag_weight
        Q[j, i] <- off_diag_weight
        
        diag_sum <- diag_sum + ((0.5) + (exp(-2 * kappa * ell) / (1 - exp(-2 * kappa * ell))))
      }
      Q[j, j] <- 2 * kappa * tau * diag_sum
    }
  }
  
  if (proper_version) {
    
    Q <- lambda * Q
    
    diag(Q) <- diag(Q) + (1 - lambda)
  }
  
  return(Q)
}


# collect Q and eiginvalues
collect_Q_and_eig <- function(Q) {
  
  rank <- rankMatrix(Q)
  E <- eigen(Q)
  determinant_log <- sum(log(E$values[1:rank]))
  
  return(list(rank = rank,
              determinant_log = determinant_log,
              E = E,
              Q = Q))
}



# Q convolution (solving the expected Dirichlet's condition for Q)
Q_convolution <- function(Q_with_eig, exponent, x, NAval = NA, use_noise = FALSE) {
  
  n <- nrow(Q_with_eig$Q)
  
  if (!is.na(NAval)) {
    ind <- which(x==NAval)
  } else {
    ind <- which(is.na(x))
  }
  
  Q_rank <- Q_with_eig$rank
  
  if (Q_rank < nrow(Q_with_eig$Q)) {
    Q_with_eig$E$values[(Q_with_eig$rank[1]+1):n] <- 0
  }
  if (use_noise) {
    d <- median(Q_with_eig$E$values) / 5
  } else {
    d <- 0
  }
  
  Q <- Q_with_eig$E$vectors %*% diag(Q_with_eig$E$values^exponent + d) %*% t(Q_with_eig$E$vectors)
  
  Q_AA <- Q[ind, ind]
  Q_AB <- Q[ind, -ind]
  
  if (rankMatrix(Q_AA) != nrow(Q_AA))
    stop("rank error, Q_AA")
  
  if (length(x[-ind]) > 1) {
    x[ind] <- - solve(Q_AA, Q_AB %*% x[-ind])
  } else {
    x[ind] <- - solve(Q_AA, Q_AB * x[-ind])
  }
  
  return(x)
}




# proximity computation form Q
proximity_from_Q <- function(Q_with_eig, tb, y_column, map, exponent_vec=1) {
  
  year_vec <- sort(unique(tb$year_idx))
  df_proximity <- tb %>% dplyr::select(geo, year_idx)
  
  df_proximity_last_year <- tibble(geo=unique(df_proximity$geo), year_idx=max(tb$year_idx)+1)
  
  df_proximity <- bind_rows(df_proximity, df_proximity_last_year)
  
  for (exponent in exponent_vec) {
    
    cat(sprintf("Computing convolution of y and Q^exponent, with exponent: %.2g\n", exponent))
    
    y_proximity_str <- sprintf("y_proximity_exp%03d", exponent*100)
    
    df_proximity_curr <- NULL
    for (year in year_vec) {
      
      map_year <- left_join(map, tb %>% filter(year_idx == year))
      
      y <- map_year[[y_column]]
      
      y_convolved <- Q_convolution(Q_with_eig, exponent=exponent, y,
                                   NAval=0, use_noise=TRUE)
      
      df_proximity_curr <-  df_proximity_curr %>%
        bind_rows(tibble(!!sym(y_proximity_str) := y_convolved,
                         y = y,
                         year_idx = year + 1,
                         geo = map_year$geo))
      
    }
    
    df_proximity_curr <- df_proximity_curr %>%
      bind_rows(tibble(!!sym(y_proximity_str) := 0,
                       y = 0,
                       year_idx = 1,
                       geo = map_year$geo))
    
    df_proximity <- left_join(df_proximity, df_proximity_curr)
  }
  
  return(df_proximity)
}




# plotting the proximity
plot_proximity <- function(df_proximity, map, filename_start=NULL) {
  
  col_vec <- str_subset(names(df_proximity), "y_proximity")
  year_vec <- sort(unique(df_proximity$year_idx))
  
  ps <- list()
  ps_collected <- list()
  for (col in col_vec) {
    ps[[col]] <- list()
    
    i_plot <- 0
    for (year in year_vec) {
      i_plot <- i_plot + 1
      
      ps[[col]][[i_plot]] <- plot_sf(left_join(map, df_proximity %>% filter(year_idx==year)),
                                     col, title=paste0(col, ", year ", year))
      
    }
    
    if (!is.null(filename_start)) {
      n_plots <- length(ps[[col]])
      ncol <- floor(sqrt(n_plots))
      nrow <- ceiling(n_plots / ncol)
      ps_collected[[col]] <- ggarrange(plotlist=ps[[col]], nrow=nrow, ncol=ncol)
      
      save_to_png(ps_collected[[col]], filename=paste0(filename_start, "_proximity_", col, ".png"))
    }
  }
  
  return(ps_collected)
}



# get surrounded regions.
## To handle the problem of one region (especially cities) lying inside another region is because NUTS3 region is divided into population
get_surrounded_regions <- function() {
  
  tb_surrounded3 <- tibble(outer=c("CZ020",
                                   "HU120",
                                   "RO322",
                                   "EL306",
                                   "EL305",
                                   "EL305",
                                   "EL305",
                                   "BG412",
                                   "AT126"),
                           inner=c("CZ010",
                                   "HU110",
                                   "RO321",
                                   "EL302",
                                   "EL304",
                                   "EL303",
                                   "EL301",
                                   "BG411",
                                   "AT130"))
  
  tb_surrounded <- tibble(outer=c("DEG0L",
                                  "DEE0B",
                                  "DE40G",
                                  "DEA2C",
                                  "DE734",
                                  "DE716",
                                  "DE94E",
                                  "DE238",
                                  "DE145",
                                  "DED2F",
                                  "DE276",
                                  "DE219",
                                  "DE40E",
                                  "DE40E",
                                  "DEE07",
                                  "DEE05"),
                          inner=c("DEG02",
                                  "DEE02",
                                  "DE402",
                                  "DEA22",
                                  "DE731",
                                  "DE711",
                                  "DE944",
                                  "DE232",
                                  "DE144",
                                  "DED21",
                                  "DE271",
                                  "DE211",
                                  "DE404",
                                  "DE401",
                                  "DEE03",
                                  "DEE01"))
  
  tb_surrounded2 <-
    tibble(outer=c("DEB3K",
                   "DEB3H",
                   "DEB3F",
                   "DE124",
                   "DE12B",
                   "DE118",
                   "DEB25",
                   "DE124",
                   "DE27E",
                   "DE27B",
                   "DE256",
                   "DE26C",
                   "DE26B",
                   "DE245",
                   "DE247",
                   "DE246",
                   "DE249",
                   "DE237",
                   "DE234",
                   "DE21K",
                   "DE227",
                   "DE22B",
                   "DEG0G",
                   "DEG0P",
                   "DE734",
                   "DEG0J"),
           inner=c("DEB37",
                   "DEB33",
                   "DEB32",
                   "DE121",
                   "DE129",
                   "DE117",
                   "DEB21",
                   "DE121",
                   "DE273",
                   "DE272",
                   "DE251",
                   "DE263",
                   "DE262",
                   "DE241",
                   "DE243",
                   "DE242",
                   "DE244",
                   "DE233",
                   "DE231",
                   "DE213",
                   "DE221",
                   "DE223",
                   "DEG05",
                   "DEG0N",
                   "DE731",
                   "DEG03"))
  
  tb_surrounded <- bind_rows(tb_surrounded, tb_surrounded2) %>%
    bind_rows(tb_surrounded3) %>% distinct()
  
  return(tb_surrounded)
}
