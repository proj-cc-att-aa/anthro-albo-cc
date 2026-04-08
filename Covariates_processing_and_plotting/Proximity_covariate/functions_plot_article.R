# ================== Plotting functions==================

# plot geo neighs article
plot_geo_neighs_article <- function(geo_neighs,
                                    filename = "",
                                    selection_nodes = NULL,
                                    node_col = "red",
                                    width_in = 6,
                                    height_in = 4,
                                    dpi = 600,
                                    compression = "lzw",
                                    resolution_scaling = 1) {
  
  coords_sp <- sp::SpatialPoints(geo_neighs$coords)
  
  if (filename != "") {
    cat(sprintf("saving plot to file: %s\n", filename))
    
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    ext <- tolower(tools::file_ext(filename))
    
    if (ext %in% c("tif", "tiff")) {
      grDevices::tiff(filename = filename,
                      width = width_in * resolution_scaling,
                      height = height_in * resolution_scaling,
                      units = "in",
                      res = dpi,
                      compression = compression)
    } else if (ext == "svg") {
      grDevices::svg(filename = filename,
                     width = width_in * resolution_scaling,
                     height = height_in * resolution_scaling,
                     pointsize = 12,
                     bg = "transparent"
      )
    }
  }
  
  if ("shapes" %in% names(geo_neighs)) {
    plot(geo_neighs$shapes, border = gray(.5), lwd = 0.4)
    add <- TRUE
  } else {
    add <- FALSE
  }
  
  plot(geo_neighs$nb_map, geo_neighs$coords, add = add,
       col = "black", lwd = 0.4, pch = ".", cex = 1)
  
  plot(coords_sp, add = TRUE,
       pch = 21,
       bg  = node_col,
       col = "black",
       lwd = 0.4,
       cex = 0.3)
  
  if (!is.null(selection_nodes)) {
    plot(coords_sp[selection_nodes], add = TRUE, col = "green", pch = 20, cex = 0.2)
  }
  
  if (filename != "") dev.off()
}



# plot mobility graph article.
plot_mobility_graph_article <- function(
    mobility,
    nodes = "",
    link_color = "tomato",
    sc_link_width = 10, sc_node_size = 6,
    filtering_p = 0.99,
    title = "",
    filename = "",
    width_in = 6, height_in = 4, dpi = 600,
    compression = "lzw",
    resolution_scaling = 1
) {
  
  n_nodes <- nrow(mobility$tb_regions)
  
  flux_tb <- flux_matrix_to_tb(mobility$flux_matrix, type = "directed")
  netw <- tibble::tibble(
    links = flux_tb$link,
    link_weights = flux_tb$flux,
    link_color = link_color
  )
  
  quant <- stats::quantile(netw$link_weights, filtering_p, na.rm = TRUE)
  netw_subset <- dplyr::filter(netw, link_weights > quant)
  cat(sprintf(
    "filtering links, only showing links with weights > %.10g%%-percentile\n",
    filtering_p * 100
  ))
  cat(sprintf("-> showing %d links (total %d links)\n", nrow(netw_subset), nrow(netw)))
  
  if (title != "") {
    title <- paste0(title, "showing links with size > ",
                    sprintf("%.10g", filtering_p*100), "%-percentile")
  }
  
  node_weights <- NULL
  if (nodes == "net_flux_abs") {
    mobility$tb_regions <- mobility$tb_regions |>
      dplyr::mutate(net_flux = abs(flux_out - flux_in))
    node_weights <- mobility$tb_regions$net_flux
  } else if (nodes == "population") {
    node_weights <- mobility$tb_regions$population
  }
  
  plot_network_mobility_article(
    links = netw_subset$links,
    link_weights = netw_subset$link_weights,
    n_nodes = n_nodes,
    node_coords = mobility$tb_regions$coords,
    node_weights = node_weights,
    sc_link_width = sc_link_width,
    sc_node_size = sc_node_size,
    link_color = link_color,
    title = title,
    filename = filename,
    width_in = width_in,
    height_in = height_in,
    dpi = dpi,
    compression = compression,
    resolution_scaling = resolution_scaling
  )
}



# plot network mobility article
plot_network_mobility_article <- function(
    links, n_nodes, node_coords,
    node_weights = NULL, link_weights = NULL,
    link_color = "tomato",
    sc_link_width = 10, sc_node_size = 6,
    title = "",
    filename = "",
    width_in = 6, height_in = 4, dpi = 600,
    compression = "lzw",
    resolution_scaling = 1
) {
  stopifnot(is.matrix(node_coords) || is.data.frame(node_coords))
  
  if (filename != "") {
    cat(sprintf("saving plot to file: %s\n", filename))
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    ext <- tolower(tools::file_ext(filename))
    
    if (ext %in% c("tif","tiff")) {
      grDevices::tiff(
        filename = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        units = "in", res = dpi, compression = compression
      )
    } else if (ext == "png") {
      grDevices::png(
        filename = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        units = "in", res = dpi
      )
    } else if (ext == "svg") {
      grDevices::svg(
        filename = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        pointsize = 12, bg = "transparent"
      )
    } else if (ext == "pdf") {
      grDevices::pdf(
        file = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        onefile = FALSE
      )
    } else {
      stop("Unsupported file extension: use .tiff/.tif, .png, .svg, or .pdf")
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  nodes <- 1:n_nodes
  net <- igraph::graph.data.frame(links, vertices = nodes, directed = TRUE)
  
  if (!is.null(node_weights)) {
    if (length(node_weights) != n_nodes) stop("node_weights length must equal n_nodes")
    igraph::V(net)$size <- (node_weights / max(node_weights, na.rm = TRUE)) * sc_node_size
  } else {
    igraph::V(net)$size <- sc_node_size
  }
  
  if (!is.null(link_weights)) {
    if (length(link_weights) != igraph::ecount(net)) {
      stop("link_weights length must equal number of edges")
    }
    igraph::E(net)$width <- (link_weights / max(link_weights, na.rm = TRUE)) * sc_link_width + 0.2
  } else {
    igraph::E(net)$width <- sc_link_width
  }
  
  graphics::plot(
    net, layout = as.matrix(node_coords),
    edge.arrow.size = 0.05, edge.curved = 0.1,
    edge.color = link_color,
    vertex.label = NA,
    main = title
  )
}



# plot Q network - for neighbourhood based proximity network
plot_Q_network_article <- function(
    Q, coords,
    link_color = "tomato",
    sc_link_width = 0.8, sc_node_size = 2,
    title = "",
    filtering_p = NULL,
    filename = "",
    width_in = 6, height_in = 4, dpi = 600,
    compression = "lzw",
    resolution_scaling = 1
) {
  
  if ((nrow(Q) != nrow(coords)) || (nrow(Q) != ncol(Q)))
    stop("dimension error: Q must be n x n and coords must have n rows")
  
  n_nodes <- nrow(Q)
  
  node_weights <- diag(Q)
  
  diag(Q) <- 0
  Q[lower.tri(Q)] <- 0
  
  links_idx <- which(Q != 0, arr.ind = TRUE)
  link_vals  <- Q[links_idx]
  
  netw <- tibble::tibble(
    links = as.data.frame(links_idx),
    link_weights = link_vals
  )
  
  if (!is.null(filtering_p)) {
    n_links_tot <- nrow(netw)
    thr <- stats::quantile(netw$link_weights, filtering_p, na.rm = TRUE)
    netw <- dplyr::filter(netw, link_weights >= thr)
    cat(sprintf(
      "filtering links, only showing links with weights > %.10g%%-quantile\n",
      filtering_p * 100
    ))
    cat(sprintf("-> showing %d links (total %d links)\n", nrow(netw), n_links_tot))
  }
  
  # Delegate plotting (handles device + saving)
  plot_network_article(
    links = netw$links,
    link_weights = netw$link_weights,
    n_nodes = n_nodes,
    node_coords = coords,
    node_weights = node_weights,
    sc_link_width = sc_link_width,
    sc_node_size = sc_node_size,
    link_color = link_color,
    title = title,
    filename = filename,
    width_in = width_in,
    height_in = height_in,
    dpi = dpi,
    compression = compression,
    resolution_scaling = resolution_scaling
  )
  
  invisible(netw)
}



# plot network for article
plot_network_article <- function(
    links, n_nodes, node_coords,
    node_weights = NULL, link_weights = NULL,
    link_color = "tomato",
    sc_link_width = 10, sc_node_size = 6,
    title = "",
    filename = "",
    width_in = 6, height_in = 4, dpi = 600,
    compression = "lzw",
    resolution_scaling = 1
) {
  
  if (filename != "") {
    cat(sprintf("saving plot to file: %s\n", filename))
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    ext <- tolower(tools::file_ext(filename))
    if (ext %in% c("tif","tiff")) {
      grDevices::tiff(
        filename = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        units = "in", res = dpi, compression = compression
      )
    } else if (ext == "png") {
      grDevices::png(
        filename = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        units = "in", res = dpi
      )
    } else if (ext == "svg") {
      grDevices::svg(
        filename = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        pointsize = 12, bg = "transparent"
      )
    } else if (ext == "pdf") {
      grDevices::pdf(
        file = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        onefile = FALSE
      )
    } else stop("Unsupported file extension. Use .tiff/.tif, .png, .svg, or .pdf")
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  nodes <- 1:n_nodes
  net <- igraph::graph.data.frame(links, vertices = nodes, directed = F)
  
  if (!is.null(node_weights)) {
    if (length(node_weights) != n_nodes) stop("node_weights length must equal n_nodes")
    igraph::V(net)$size <- (node_weights / max(node_weights, na.rm = TRUE)) * sc_node_size
  } else {
    igraph::V(net)$size <- sc_node_size
  }
  
  if (!is.null(link_weights)) {
    if (length(link_weights) != igraph::ecount(net))
      stop("link_weights length must equal number of edges")
    igraph::E(net)$width <- (link_weights / max(link_weights, na.rm = TRUE)) * sc_link_width
  } else {
    igraph::E(net)$width <- sc_link_width
  }
  
  if (!is.null(node_weights)) {
    if (length(node_weights) != n_nodes) stop("node_weights length must equal n_nodes")
    node_sizes_px <- (node_weights / max(node_weights, na.rm = TRUE)) * sc_node_size
  } else {
    node_sizes_px <- rep(sc_node_size, n_nodes)
  }
  
  node_fill   <- "orange"
  node_border <- "black"
  node_lwd    <- 0.4
  
  cex_nodes <- pmax(0.001, node_sizes_px) * 0.45
  
  lay <- as.matrix(node_coords)
  
  lay_s <- igraph::norm_coords(lay, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  
  graphics::plot(
    net, layout = lay_s,
    edge.arrow.size = 0.05, edge.curved = 0.1,
    edge.color = link_color,
    vertex.size = 0, vertex.color = NA, vertex.frame.color = NA,
    vertex.label = NA,
    main = title, axes = FALSE
  )
  
  points(
    lay_s[,1], lay_s[,2],
    pch = 21,
    bg  = node_fill,
    col = node_border,
    lwd = node_lwd,
    cex = cex_nodes
  )
  
}



# plot Q network centroid (centroid based proximity network)
plot_Q_network_centroid_article <- function(
    Q, coords,
    link_color = "tomato",
    sc_link_width = 0.8, sc_node_size = 2,
    title = "",
    filtering_p = NULL,
    filename = "",
    width_in = 6, height_in = 4, dpi = 600,
    compression = "lzw",
    resolution_scaling = 1,
    edge_colour_by_link_weights = FALSE,
    legend_title = "Link strength (percentile)"
) {
  
  if ((nrow(Q) != nrow(coords)) || (nrow(Q) != ncol(Q)))
    stop("dimension error: Q must be n x n and coords must have n rows")
  
  n_nodes <- nrow(Q)
  
  node_weights <- diag(Q)
  
  diag(Q) <- 0
  Q[lower.tri(Q)] <- 0
  
  links_idx <- which(Q != 0, arr.ind = TRUE)
  link_vals <- Q[links_idx]
  
  netw <- tibble::tibble(
    links        = as.data.frame(links_idx),
    link_weights = link_vals
  )
  
  if (!is.null(filtering_p)) {
    n_links_tot <- nrow(netw)
    thr <- stats::quantile(netw$link_weights, filtering_p, na.rm = TRUE)
    netw <- dplyr::filter(netw, link_weights >= thr)
    cat(sprintf(
      "filtering links, only showing links with weights > %.10g%%-quantile\n",
      filtering_p * 100
    ))
    cat(sprintf("-> showing %d links (total %d links)\n", nrow(netw), n_links_tot))
  }
  
  plot_network_article_centroid(
    links        = netw$links,
    link_weights = netw$link_weights,
    n_nodes      = n_nodes,
    node_coords  = coords,
    node_weights = node_weights,
    sc_link_width = sc_link_width,
    sc_node_size  = sc_node_size,
    link_color    = link_color,
    title         = title,
    filename      = filename,
    width_in      = width_in,
    height_in     = height_in,
    dpi           = dpi,
    compression   = compression,
    resolution_scaling = resolution_scaling,
    edge_colour_by_link_weights = edge_colour_by_link_weights,
    legend_title  = legend_title
  )
  
  invisible(netw)
}



# plot network article centroid
plot_network_article_centroid <- function(
    links, n_nodes, node_coords,
    node_weights = NULL, link_weights = NULL,
    link_color = "tomato",
    sc_link_width = 2, sc_node_size = 5,
    title = "",
    filename = "",
    width_in = 6, height_in = 4, dpi = 600,
    compression = "lzw",
    resolution_scaling = 1,
    edge_colour_by_link_weights = FALSE,
    legend_title = "Link strength",
    legend_n = 80
) {
  
  if (filename != "") {
    cat(sprintf("saving plot to file: %s\n", filename))
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    ext <- tolower(tools::file_ext(filename))
    
    if (ext %in% c("tif", "tiff")) {
      grDevices::tiff(
        filename = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        units = "in", res = dpi, compression = compression
      )
    } else if (ext == "png") {
      grDevices::png(
        filename = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        units = "in", res = dpi
      )
    } else if (ext == "svg") {
      grDevices::svg(
        filename = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        pointsize = 12, bg = "transparent"
      )
    } else if (ext == "pdf") {
      grDevices::pdf(
        file = filename,
        width = width_in * resolution_scaling,
        height = height_in * resolution_scaling,
        onefile = FALSE
      )
    } else stop("Unsupported file extension. Use .tiff/.tif, .png, .svg, or .pdf")
    
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  
  if (op$mar[4] < 5) {
    graphics::par(mar = c(op$mar[1], op$mar[2], op$mar[3], 5))
  }
  
  nodes <- 1:n_nodes
  net   <- igraph::graph.data.frame(links, vertices = nodes, directed = FALSE)
  
  if (!is.null(node_weights)) {
    if (length(node_weights) != n_nodes)
      stop("node_weights length must equal n_nodes")
    node_sizes_px <- (node_weights / max(node_weights, na.rm = TRUE)) * sc_node_size
  } else {
    node_sizes_px <- rep(sc_node_size, n_nodes)
  }
  
  cex_nodes <- pmax(0.001, node_sizes_px) * 0.45
  
  if (!is.null(link_weights)) {
    if (length(link_weights) != igraph::ecount(net))
      stop("link_weights length must equal number of edges")
    
    w_abs <- abs(link_weights)
    ecdf_w   <- stats::ecdf(w_abs)
    w_scaled <- ecdf_w(w_abs)
    w_scaled[!is.finite(w_scaled)] <- 0
    
    gamma <- 0.7
    w_width_scaled <- w_scaled ^ gamma
    
    igraph::E(net)$width <- pmax(0.4, w_width_scaled * sc_link_width)
    
    if (edge_colour_by_link_weights) {
      pal <- grDevices::colorRampPalette(c("darkblue", "orange", "darkred"))(legend_n)
      idx <- pmax(1, pmin(legend_n, floor(w_scaled * (legend_n - 1)) + 1))
      edge_col <- pal[idx]
    } else {
      edge_col <- link_color
    }
    
  } else {
    igraph::E(net)$width <- sc_link_width
    edge_col <- link_color
  }
  
  lay   <- as.matrix(node_coords)
  lay_s <- igraph::norm_coords(lay, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  
  #draw EDGES ONLY with igraph
  graphics::plot(
    net, layout = lay_s,
    edge.arrow.size = 0.05, edge.curved = 0.1,
    edge.color = edge_col,
    vertex.size = 0.0, vertex.color = NA, vertex.frame.color = NA,
    vertex.label = NA,
    main = title, axes = FALSE
  )
  
  #overlay NODES with base points() 
  node_fill   <- "orange"
  node_border <- "black"
  node_lwd    <- 0.4
  
  graphics::points(
    lay_s[, 1], lay_s[, 2],
    pch = 21,
    bg  = node_fill,
    col = node_border,
    lwd = node_lwd,
    cex = cex_nodes
  )
  
  if (edge_colour_by_link_weights) {
    
    # draw colorbar.
    draw_colorbar <- function(
    zlim = c(0, 1), pal, n = length(pal),
    label = "Link strength (percentile)",
    n_ticks = 5,
    cex_axis = 0.9,
    cex_label = 1.1,
    height_frac = 0.7,    # <- fraction of plot height used by legend (0–1)
    width_frac  = 0.04
    ) {
      
      usr     <- graphics::par("usr")
      x_range <- usr[2] - usr[1]
      y_range <- usr[4] - usr[3]
      
      x_left  <- usr[2] + 0.02 * x_range
      x_right <- x_left + width_frac * x_range
      
      y_mid   <- (usr[3] + usr[4]) / 2
      half_h  <- (y_range * height_frac) / 2
      y_bottom <- y_mid - half_h
      y_top    <- y_mid + half_h
      
      y_vals <- seq(y_bottom, y_top, length.out = n + 1)
      
      graphics::par(xpd = TRUE)
      for (i in seq_len(n)) {
        graphics::rect(
          xleft  = x_left,
          ybottom = y_vals[i],
          xright = x_right,
          ytop   = y_vals[i + 1],
          col    = pal[i],
          border = NA
        )
      }
      
      tick_probs <- seq(0, 1, length.out = n_ticks)
      tick_pos   <- y_bottom + tick_probs * (y_top - y_bottom)
      tick_labs  <- paste0(round(tick_probs * 100), "%")
      
      graphics::axis(
        side = 4,
        at   = tick_pos,
        labels = tick_labs,
        las  = 1,
        cex.axis = cex_axis
      )
      
      graphics::mtext(
        label,
        side = 4,
        line = 1.8,         # distance of legend title from axis (smaller = closer) 
        cex  = cex_label
      )
      
      graphics::par(xpd = FALSE)
    }
    
    pal <- grDevices::colorRampPalette(c("darkblue", "orange", "darkred"))(legend_n)
    draw_colorbar(
      zlim  = c(0, 1),
      pal   = pal,
      label = legend_title,
      cex_axis  = 0.5,
      cex_label = 0.6,
      height_frac = 0.5,
      width_frac  = 0.03
    )
  }
}



# save plot device (evice helper, reused everywhere)
save_plot_device <- function(p, filename,
                             width_in = 6, height_in = 4, dpi = 600,
                             compression = "lzw", bg = "transparent") {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(filename))
  
  if (ext %in% c("tif", "tiff")) {
    grDevices::tiff(
      filename = filename,
      width = width_in, height = height_in, units = "in",
      res = dpi, compression = compression, bg = bg
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
    
  } else if (ext == "png") {
    grDevices::png(
      filename = filename,
      width = width_in, height = height_in, units = "in",
      res = dpi, bg = bg
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
    
  } else if (ext == "svg") {
    grDevices::svg(
      filename = filename,
      width = width_in, height = height_in,
      pointsize = 12, bg = bg
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
    
  } else if (ext == "pdf") {
    grDevices::pdf(
      file = filename, width = width_in, height = height_in,
      onefile = FALSE, paper = "special", bg = bg
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
    
  } else {
    stop("Unsupported file extension. Use .tiff/.tif, .png, .svg, or .pdf")
  }
}



# plot save proximities for year for multple strength (/single strength ) 
plot_save_proximities_for_year <- function(
    df_proximity, map_inla, year_idx,
    value_cols = c("y_proximity_exp050","y_proximity_exp100","y_proximity_exp200","y_proximity_exp300"),
    filename_stub,
    ext = "tiff",
    width_in = 3.5, height_in = 4, dpi = 600, compression = "lzw",
    legend_grad_filename     = NULL,
    legend_presence_filename = NULL,
    legend_grad_width_in     = 1.4,
    legend_grad_height_in    = 4,
    legend_presence_width_in  = 1.4,
    legend_presence_height_in = 0.9,
    legend_dpi               = NULL,
    legend_ext               = NULL,
    presence_fill = "grey50",
    na_fill = "grey95",
    border_col = "black", border_size = 0.2,
    
    outer_boundary_col  = "black",
    outer_boundary_size = 0.3,
    
    palette_name = "Zissou1", palette_n = 100, palette_reverse = FALSE
) {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr")
  if (!requireNamespace("sf", quietly = TRUE))     stop("Please install sf")
  if (!requireNamespace("ggplot2", quietly = TRUE))stop("Please install ggplot2")
  if (!requireNamespace("cowplot", quietly = TRUE))stop("Please install cowplot")
  if (!requireNamespace("wesanderson", quietly = TRUE)) stop("Please install wesanderson")
  
  library(dplyr); library(sf); library(ggplot2)
  
  pal <- wesanderson::wes_palette(palette_name, palette_n, type = "continuous")
  if (isTRUE(palette_reverse)) pal <- rev(pal)
  
  yr_df <- df_proximity %>%
    dplyr::filter(year_idx == !!year_idx) %>%
    dplyr::select(geo, presence_previous_year, dplyr::all_of(value_cols))
  
  sf_year <- map_inla %>%
    dplyr::select(geo, geometry) %>%
    dplyr::left_join(yr_df, by = "geo")
  
  sf_year <- sf::st_as_sf(sf_year, sf_column_name = "geometry")
  
  vals_for_scale <- sf_year %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(presence_previous_year == 0) %>%
    dplyr::select(dplyr::all_of(value_cols)) %>%
    unlist(use.names = FALSE)
  
  lims  <- range(vals_for_scale, na.rm = TRUE, finite = TRUE)
  
  step  <- 0.2
  min_b <- floor(lims[1] / step) * step
  max_b <- ceiling(lims[2] / step) * step
  breaks <- seq(min_b, max_b, by = step)
  
  map_inla <- sf::st_as_sf(map_inla, sf_column_name = "geometry")
  
  geom <- sf::st_geometry(map_inla)
  outline <- sf::st_as_sf(sf::st_union(geom))
  
  for (vc in value_cols) {
    
    sf_year <- sf_year %>%
      dplyr::mutate(value_masked = ifelse(presence_previous_year == 1, NA_real_, .data[[vc]]))
    
    p <- ggplot() +
      geom_sf(
        data = sf_year,
        aes(fill = value_masked),
        color = border_col, linewidth = border_size
      ) +
      scale_fill_gradientn(
        colors = pal,
        limits = lims,
        breaks = breaks,
        labels = scales::label_number(accuracy = 0.1),
        oob = scales::squish,
        na.value = NA,
        name = "Proximity"
      ) +
      geom_sf(
        data = sf_year %>% dplyr::filter(presence_previous_year == 1),
        fill = presence_fill, color = border_col, linewidth = border_size,
        inherit.aes = FALSE
      ) +
      
      geom_sf(
        data = outline,
        fill = NA,
        color = outer_boundary_col,
        linewidth = outer_boundary_size,
        inherit.aes = FALSE
      ) +
      
      coord_sf(datum = NA) +
      labs(title = "") +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none"
      )
    
    out_file <- sprintf("%s_%s.%s", filename_stub, vc, ext)
    save_plot_device(p, out_file, width_in, height_in, dpi, compression)
    cat("saved:", out_file, "\n")
  }
  
  # gradient legend (use the last value_masked column present in sf_year)
  p_grad_leg <- ggplot(sf_year) +
    geom_sf(aes(fill = value_masked), color = NA) +
    scale_fill_gradientn(
      colors = pal, limits = lims, breaks = breaks,
      oob = scales::squish,
      labels = scales::label_number(accuracy = 0.1),
      na.value = NA, name = NULL
    ) +
    theme_void(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text  = element_text(size = 8)
    )
  grad_grob <- cowplot::get_legend(p_grad_leg)
  
  if (is.null(legend_ext)) legend_ext <- ext
  if (is.null(legend_dpi)) legend_dpi <- dpi
  if (is.null(legend_grad_filename)) {
    legend_grad_filename <- sprintf("%s__LEGEND_GRAD.%s", filename_stub, legend_ext)
  }
  save_plot_device(
    cowplot::ggdraw(grad_grob),
    legend_grad_filename,
    legend_grad_width_in, legend_grad_height_in, legend_dpi, compression
  )
  cat("saved legend:", legend_grad_filename, "\n")
  
  presence_key <- data.frame(lbl = "Presence previous year", x = 0, y = 0)
  p_grey_leg <- ggplot(presence_key, aes(x, y, fill = lbl)) +
    geom_tile() +
    scale_fill_manual(
      values = c("Presence previous year" = presence_fill),
      name   = NULL,
      breaks = "Presence previous year",
      labels = ""
    ) +
    theme_void(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title    = element_blank(),
      legend.text     = element_blank()
    )
  grey_grob <- cowplot::get_legend(p_grey_leg)
  
  if (is.null(legend_presence_filename)) {
    legend_presence_filename <- sprintf("%s__LEGEND_PRESENCE.%s", filename_stub, legend_ext)
  }
  save_plot_device(
    cowplot::ggdraw(grey_grob),
    legend_presence_filename,
    legend_presence_width_in, legend_presence_height_in, legend_dpi, compression
  )
  cat("saved legend:", legend_presence_filename, "\n")
  
  invisible(TRUE)
}



# plot save proximity all years (for single strength (/single strength ) for all years with common scale legend)
plot_save_proximity_all_years <- function(
    df_proximity, map_inla,
    value_col = "y_proximity_exp100",
    years = NULL,
    filename_stub,
    ext = "tiff",
    width_cm  = 5.7,
    height_cm = 8.5,
    dpi = 600, compression = "lzw",
    legend_grad_filename     = NULL,
    legend_presence_filename = NULL,
    legend_grad_width_cm     = 2.8,
    legend_grad_height_cm    = 8.5,
    legend_presence_width_cm  = 2.8,
    legend_presence_height_cm = 2.3,
    legend_dpi               = NULL,
    legend_ext               = NULL,
    presence_fill = "grey50",
    na_fill = "grey95",
    border_col = "black", border_size = 0.2,
    
    outer_boundary_col  = "black",
    outer_boundary_size = 0.3,
    palette_name = "Zissou1", palette_n = 100, palette_reverse = FALSE
) {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr")
  if (!requireNamespace("wesanderson", quietly = TRUE)) stop("Please install wesanderson")
  
  library(dplyr); library(sf); library(ggplot2)
  
  pal <- wesanderson::wes_palette(palette_name, palette_n, type = "continuous")
  if (isTRUE(palette_reverse)) pal <- rev(pal)
  
  if (is.null(years)) years <- sort(unique(df_proximity$year_idx))
  
  vals_for_scale <- df_proximity %>%
    filter(year_idx %in% years, presence_previous_year == 0) %>%
    pull(all_of(value_col))
  
  lims  <- range(vals_for_scale, na.rm = TRUE, finite = TRUE)
  
  step  <- 0.2
  min_b <- floor(lims[1] / step) * step
  max_b <- ceiling(lims[2] / step) * step
  breaks <- seq(min_b, max_b, by = step)
  
  map_inla <- sf::st_as_sf(map_inla, sf_column_name = "geometry")
  
  geom <- sf::st_geometry(map_inla)
  outline <- sf::st_as_sf(sf::st_union(geom))
  
  # cm to in.
  cm_to_in <- function(x) x / 2.54
  
  width_in  <- cm_to_in(width_cm)
  height_in <- cm_to_in(height_cm)
  
  legend_grad_width_in      <- cm_to_in(legend_grad_width_cm)
  legend_grad_height_in     <- cm_to_in(legend_grad_height_cm)
  legend_presence_width_in  <- cm_to_in(legend_presence_width_cm)
  legend_presence_height_in <- cm_to_in(legend_presence_height_cm)
  
  for (yy in years) {
    yr_df <- df_proximity %>%
      filter(year_idx == !!yy) %>%
      select(geo, presence_previous_year, !!sym(value_col))
    
    sf_year <- map_inla %>%
      select(geo, geometry) %>%
      left_join(yr_df, by = "geo") %>%
      st_as_sf(sf_column_name = "geometry") %>%
      mutate(value_masked = ifelse(presence_previous_year == 1, NA_real_, .data[[value_col]]))
    
    p <- ggplot() +
      geom_sf(data = sf_year, aes(fill = value_masked), color = border_col, linewidth = border_size) +
      scale_fill_gradientn(
        colors = pal,
        limits = lims,
        breaks = breaks,
        labels = scales::label_number(accuracy = 0.1),
        oob = scales::squish,
        na.value = NA,
        name = "Proximity"
      ) +
      
      geom_sf(
        data = sf_year %>% filter(presence_previous_year == 1),
        fill = presence_fill, color = border_col, linewidth = border_size, inherit.aes = FALSE
      ) +
      
      geom_sf(
        data = outline,
        fill = NA,
        color = outer_boundary_col,
        linewidth = outer_boundary_size,
        inherit.aes = FALSE
      ) +
      coord_sf(datum = NA) +
      labs(title = "") +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none"
      )
    
    out_file <- sprintf("%s_%s_year%02d.%s", filename_stub, value_col, yy, ext)
    save_plot_device(p, out_file, width_in, height_in, dpi, compression)
    cat("saved:", out_file, "\n")
    
    if (yy == years[1]) {
      p_grad_leg <- ggplot(sf_year %>% mutate(value_unmasked = .data[[value_col]])) +
        geom_sf(aes(fill = value_unmasked), color = NA) +
        scale_fill_gradientn(
          colors = pal, limits = lims, breaks = breaks,
          labels = scales::label_number(accuracy = 0.1),
          oob = scales::squish, na.value = NA, name = NULL
        ) +
        theme_void(base_size = 10) +
        theme(
          legend.position = "right",
          legend.title = element_blank(),
          legend.text  = element_text(size = 8)
        )
      grad_grob <- cowplot::get_legend(p_grad_leg)
      
      if (is.null(legend_ext)) legend_ext <- ext
      if (is.null(legend_dpi)) legend_dpi <- dpi
      if (is.null(legend_grad_filename)) {
        legend_grad_filename <- sprintf("%s_%s__LEGEND_GRAD.%s", filename_stub, value_col, legend_ext)
      }
      save_plot_device(cowplot::ggdraw(grad_grob),
                       legend_grad_filename,
                       legend_grad_width_in, legend_grad_height_in, legend_dpi, compression)
      cat("saved legend:", legend_grad_filename, "\n")
      
      presence_key <- data.frame(lbl = "Presence previous year", x = 0, y = 0)
      p_grey_leg <- ggplot(presence_key, aes(x, y, fill = lbl)) +
        geom_tile() +
        scale_fill_manual(
          values = c("Presence previous year" = presence_fill),
          name   = NULL,
          breaks = "Presence previous year",
          labels = ""
        ) +
        theme_void(base_size = 10) +
        theme(
          legend.position = "right",
          legend.title    = element_blank(),
          legend.text     = element_blank()
        )
      grey_grob <- cowplot::get_legend(p_grey_leg)
      
      if (is.null(legend_presence_filename)) {
        legend_presence_filename <- sprintf("%s_%s__LEGEND_PRESENCE.%s", filename_stub, value_col, legend_ext)
      }
      save_plot_device(cowplot::ggdraw(grey_grob),
                       legend_presence_filename,
                       legend_presence_width_in, legend_presence_height_in, legend_dpi, compression)
      cat("saved legend:", legend_presence_filename, "\n")
    }
  }
  
  invisible(TRUE)
}




# create aggregated plot (land use  scaled version (line plot))
create_aggregated_plot_article <- function(
    adjusted_data,
    map_inla,
    variable,
    y_axis_label = NULL,
    break_year = 2015,
    save_filename = NULL,
    width_in = 6, height_in = 4, dpi = 600,
    compression = "lzw", bg = "transparent",
    legend_filename = NULL,
    legend_width_in = 2.0, legend_height_in = 0.8, legend_dpi = 600,
    legend_compression = "lzw", legend_bg = "transparent"
) {
  
  # clean variable name.
  clean_variable_name <- function(var_name) {
    pretty <- gsub("_mean$", "", var_name)
    toTitleCase(gsub("_", " ", pretty))
  }
  
  var_pretty      <- clean_variable_name(variable)
  scaled_variable <- paste0("scaled_", variable)
  if (!scaled_variable %in% names(adjusted_data)) {
    stop(paste0("Column '", scaled_variable, "' not found in adjusted_data."))
  }
  
  y_lab <- if (is.null(y_axis_label)) var_pretty else y_axis_label
  
  base <- adjusted_data %>%
    mutate(scenario_lower = trimws(tolower(scenario))) %>%
    left_join(map_inla %>% select(geo, area), by = "geo")
  
  if (any(is.na(base$area))) {
    base$area[is.na(base$area)] <- 1
  }
  
  yr_min <- suppressWarnings(min(base$year, na.rm = TRUE))
  yr_max <- suppressWarnings(max(base$year, na.rm = TRUE))
  if (!is.finite(yr_min) || !is.finite(yr_max)) stop("No valid years in adjusted_data.")
  if (break_year <= yr_min) warning("break_year <= min(year); historical segment may be empty.")
  if (break_year > yr_max)  warning("break_year > max(year); scaled/SSP370 segments may be empty.")
  
  hist_unscaled <- base %>%
    filter(scenario == "historical", year <= (break_year - 1)) %>%
    group_by(year) %>%
    summarise(
      mean_aw = sum(.data[[variable]] * area, na.rm = TRUE) / sum(area, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(series = paste(var_pretty, "Historical"))
  
  hist_scaled <- base %>%
    filter(scenario == "SSP370", year >= break_year) %>%
    group_by(year) %>%
    summarise(
      mean_aw = sum(.data[[scaled_variable]] * area, na.rm = TRUE) / sum(area, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(series = paste(var_pretty, "Historical (scaled)"))
  
  ssp370_unscaled <- base %>%
    filter(scenario == "SSP370", year >= break_year) %>%
    group_by(year) %>%
    summarise(
      mean_aw = sum(.data[[variable]] * area, na.rm = TRUE) / sum(area, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(series = paste(var_pretty, "SSP370"))
  
  plot_df <- bind_rows(hist_unscaled, hist_scaled, ssp370_unscaled)
  
  series_levels <- c(
    paste(var_pretty, "Historical"),
    paste(var_pretty, "Historical (scaled)"),
    paste(var_pretty, "SSP370")
  )
  plot_df$series <- factor(plot_df$series, levels = series_levels)
  
  color_map <- setNames(
    c("#1f77b4", "#2ca02c", "#d62728"),
    series_levels
  )
  
  p_with_legend <- ggplot(plot_df, aes(x = year, y = mean_aw, color = series, group = series)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = color_map, breaks = series_levels, drop = FALSE) +
    scale_x_continuous(
      limits = c(min(plot_df$year, na.rm = TRUE), max(plot_df$year, na.rm = TRUE) + 1),
      expand = c(0, 0)
    ) +
    labs(
      title = "",
      x = "Year",
      y = y_lab,
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.9),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(5, "pt"),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10, color = "black"),
      legend.position = "bottom",
      legend.box = "outside"
    )
  
  series_levels_leg <- c("Historical", "Historical (scaled)", "SSP370")
  
  legend_df <- data.frame(
    series_leg = factor(series_levels_leg, levels = series_levels_leg),
    x = 1, y = 1
  )
  
  color_map_leg <- setNames(
    as.character(color_map[series_levels]),
    series_levels_leg
  )
  
  p_legend_source <- ggplot(legend_df, aes(x = x, y = y, color = series_leg)) +
    geom_point(size = 3, show.legend = TRUE) +
    scale_color_manual(values = color_map_leg, breaks = series_levels_leg, drop = FALSE) +
    guides(color = guide_legend(
      ncol = 1, byrow = TRUE,
      keyheight = grid::unit(15, "pt")
    )) +
    theme_void(base_size = 12) +
    theme(
      legend.position   = "right",
      legend.direction  = "vertical",
      legend.title      = element_blank(),
      legend.key.height = grid::unit(12, "pt"),
      legend.spacing.y  = grid::unit(10,  "pt"),
      legend.text       = element_text(size = 10)
    )
  
  ## Extract legend via gtable (robust)
  gt <- ggplotGrob(p_legend_source)
  legend_grob  <- gtable::gtable_filter(gt, "guide-box", trim = TRUE)
  legend_plot  <- cowplot::ggdraw(legend_grob)
  
  p <- p_with_legend + theme(legend.position = "none")
  
  y_hist_prev <- hist_unscaled %>% filter(year == (break_year - 1)) %>% pull(mean_aw)
  y_scal_curr <- hist_scaled   %>% filter(year ==  break_year     ) %>% pull(mean_aw)
  if (length(y_hist_prev) == 1 && length(y_scal_curr) == 1) {
    bridge_df <- data.frame(x0 = break_year - 1, y0 = y_hist_prev,
                            x1 = break_year,     y1 = y_scal_curr)
    p <- p + geom_segment(
      data = bridge_df,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      inherit.aes = FALSE,
      linewidth = 1,
      color = color_map[paste(var_pretty, "Historical (scaled)")]
    )
  }
  
  # Save main plot (optional)
  if (!is.null(save_filename)) {
    save_plot_device(
      p, save_filename,
      width_in = width_in, height_in = height_in, dpi = dpi,
      compression = compression, bg = bg
    )
    cat("Saved plot:", save_filename, "\n")
  }
  
  if (!is.null(legend_filename)) {
    save_plot_device(
      legend_plot, legend_filename,
      width_in = legend_width_in, height_in = legend_height_in, dpi = legend_dpi,
      compression = legend_compression, bg = legend_bg
    )
    cat("Saved legend:", legend_filename, "\n")
  }
  
  return(p)
}



# create aggregated plot article population
##  CAUTION -- line plot for population for OLD ISMIP DATA PLOT NOT USED IN FINAL VERSION
create_aggregated_plot_article_pop <- function(
    adjusted_data,
    map_inla,
    variable,
    y_axis_label = NULL,
    break_year = 2015,
    save_filename = NULL,
    width_in = 6, height_in = 4, dpi = 600,
    compression = "lzw", bg = "transparent",
    legend_filename = NULL,
    legend_width_in = 2.0, legend_height_in = 0.8, legend_dpi = 600,
    legend_compression = "lzw", legend_bg = "transparent"
) {
  
  # clean variable name.
  clean_variable_name <- function(var_name) {
    pretty <- gsub("_mean$", "", var_name)
    toTitleCase(gsub("_", " ", pretty))
  }
  
  var_pretty      <- clean_variable_name(variable)
  scaled_variable <- paste0("scaled_", variable)
  if (!scaled_variable %in% names(adjusted_data)) {
    stop(paste0("Column '", scaled_variable, "' not found in adjusted_data."))
  }
  
  y_lab <- if (is.null(y_axis_label)) var_pretty else y_axis_label
  
  base <- adjusted_data %>%
    mutate(scenario_lower = trimws(tolower(scenario))) %>%
    left_join(map_inla %>% select(geo, area), by = "geo")
  
  if (any(is.na(base$area))) {
    base$area[is.na(base$area)] <- 1
  }
  
  yr_min <- suppressWarnings(min(base$year, na.rm = TRUE))
  yr_max <- suppressWarnings(max(base$year, na.rm = TRUE))
  if (!is.finite(yr_min) || !is.finite(yr_max)) stop("No valid years in adjusted_data.")
  if (break_year <= yr_min) warning("break_year <= min(year); historical segment may be empty.")
  if (break_year > yr_max)  warning("break_year > max(year); scaled/SSP370 segments may be empty.")
  
  hist_unscaled <- base %>%
    filter(scenario == "historical", year <= (break_year - 1)) %>%
    group_by(year) %>%
    summarise(
      mean_aw = sum(.data[[variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(series = paste(var_pretty, "Historical"))
  
  hist_scaled <- base %>%
    filter(scenario == "SSP370", year >= break_year) %>%
    group_by(year) %>%
    summarise(
      mean_aw = sum(.data[[scaled_variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(series = paste(var_pretty, "Historical (scaled)"))
  
  ssp370_unscaled <- base %>%
    filter(scenario == "SSP370", year >= break_year) %>%
    group_by(year) %>%
    summarise(
      mean_aw = sum(.data[[variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(series = paste(var_pretty, "SSP370"))
  
  plot_df <- bind_rows(hist_unscaled, hist_scaled, ssp370_unscaled)
  
  series_levels <- c(
    paste(var_pretty, "Historical"),
    paste(var_pretty, "Historical (scaled)"),
    paste(var_pretty, "SSP370")
  )
  plot_df$series <- factor(plot_df$series, levels = series_levels)
  
  color_map <- setNames(
    c("#1f77b4", "#2ca02c", "#d62728"),
    series_levels
  )
  
  p_with_legend <- ggplot(plot_df, aes(x = year, y = mean_aw, color = series, group = series)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = color_map, breaks = series_levels, drop = FALSE) +
    scale_x_continuous(
      limits = c(min(plot_df$year, na.rm = TRUE), max(plot_df$year, na.rm = TRUE) + 1),
      expand = c(0, 0)
    ) +
    labs(
      title = "",
      x = "Year",
      y = y_lab,
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.9),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(5, "pt"),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10, color = "black"),
      legend.position = "bottom",
      legend.box = "outside"
    )
  
  series_levels_leg <- c("Historical", "Historical (scaled)", "SSP370")
  
  legend_df <- data.frame(
    series_leg = factor(series_levels_leg, levels = series_levels_leg),
    x = 1, y = 1
  )
  
  color_map_leg <- setNames(
    as.character(color_map[series_levels]),
    series_levels_leg
  )
  
  p_legend_source <- ggplot(legend_df, aes(x = x, y = y, color = series_leg)) +
    geom_point(size = 3, show.legend = TRUE) +
    scale_color_manual(values = color_map_leg, breaks = series_levels_leg, drop = FALSE) +
    guides(color = guide_legend(
      ncol = 1, byrow = TRUE,
      keyheight = grid::unit(15, "pt")
    )) +
    theme_void(base_size = 12) +
    theme(
      legend.position   = "right",
      legend.direction  = "vertical",
      legend.title      = element_blank(),
      legend.key.height = grid::unit(12, "pt"),
      legend.spacing.y  = grid::unit(10,  "pt"),
      legend.text       = element_text(size = 10)
    )
  
  gt <- ggplotGrob(p_legend_source)
  legend_grob  <- gtable::gtable_filter(gt, "guide-box", trim = TRUE)
  legend_plot  <- cowplot::ggdraw(legend_grob)
  
  p <- p_with_legend + theme(legend.position = "none")
  
  y_hist_prev <- hist_unscaled %>% filter(year == (break_year - 1)) %>% pull(mean_aw)
  y_scal_curr <- hist_scaled   %>% filter(year ==  break_year     ) %>% pull(mean_aw)
  if (length(y_hist_prev) == 1 && length(y_scal_curr) == 1) {
    bridge_df <- data.frame(x0 = break_year - 1, y0 = y_hist_prev,
                            x1 = break_year,     y1 = y_scal_curr)
    p <- p + geom_segment(
      data = bridge_df,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      inherit.aes = FALSE,
      linewidth = 1,
      color = color_map[paste(var_pretty, "Historical (scaled)")]
    )
  }
  
  if (!is.null(save_filename)) {
    save_plot_device(
      p, save_filename,
      width_in = width_in, height_in = height_in, dpi = dpi,
      compression = compression, bg = bg
    )
    cat("Saved plot:", save_filename, "\n")
  }
  
  if (!is.null(legend_filename)) {
    save_plot_device(
      legend_plot, legend_filename,
      width_in = legend_width_in, height_in = legend_height_in, dpi = legend_dpi,
      compression = legend_compression, bg = legend_bg
    )
    cat("Saved legend:", legend_filename, "\n")
  }
  
  return(p)
}



# plot adjusted and delta maps ( final climate covariates map plot)
plot_adjusted_and_delta_maps <- function(
    adjusted_df,
    pi_df,
    map_inla,
    year,
    variable,
    
    filename_stub,
    ext = "tiff",
    width_in = 3.5, height_in = 4, dpi = 600, compression = "lzw", bg = "transparent",
    
    hist_legend_ext = NULL, hist_legend_dpi = NULL,
    hist_legend_filename  = NULL, hist_legend_width_in  = 1.4, hist_legend_height_in  = 4.0,
    hist_legend_compression = "lzw", hist_legend_bg = "transparent",
    
    filename_stub_diff,
    
    diff_legend_ext = NULL, diff_legend_dpi = NULL,
    diff_legend_filename  = NULL, diff_legend_width_in  = 1.4, diff_legend_height_in  = 4.0,
    diff_legend_compression = "lzw", diff_legend_bg = "transparent",
    
    na_fill = "grey95",
    border_col = "black", border_size = 0.2,
    
    outer_boundary_col  = "black",
    outer_boundary_size = 0.3,
    
    palette_name = "YlOrRd",
    palette_n = 100,
    palette_reverse = FALSE,
    
    palette_name_diff  = "PuBuGn",
    palette_n_diff = 100,
    palette_reverse_diff = FALSE
    
) {
  
  stopifnot("geo" %in% names(adjusted_df), "year" %in% names(adjusted_df))
  stopifnot("geo" %in% names(pi_df),       "year" %in% names(pi_df))
  stopifnot("geometry" %in% names(map_inla))
  
  if (!palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
    warning(sprintf("Unknown palette_name = '%s'. Falling back to 'YlOrRd'.", palette_name))
    palette_name <- "YlOrRd"
  }
  
  max_n <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  pal_base <- RColorBrewer::brewer.pal(max_n, palette_name)
  pal_seq <- grDevices::colorRampPalette(pal_base)(palette_n)
  
  if (isTRUE(palette_reverse)) pal_seq <- rev(pal_seq)
  
  if (!exists("palette_name_diff") || is.na(palette_name_diff) ||
      !palette_name_diff %in% rownames(RColorBrewer::brewer.pal.info)) {
    warning(sprintf(
      "Unknown palette_name_diff = '%s'. Falling back to 'PuBuGn'.",
      if (!exists("palette_name_diff")) "<missing>" else as.character(palette_name_diff)
    ))
    palette_name_diff <- "PuBuGn"
  }
  
  max_n_diff <- RColorBrewer::brewer.pal.info[palette_name_diff, "maxcolors"]
  pal_base_diff <- RColorBrewer::brewer.pal(max_n_diff, palette_name_diff)
  pal_seq_diff  <- grDevices::colorRampPalette(pal_base_diff)(palette_n)
  
  if (missing(variable) || is.null(variable)) {
    stop("Please set `variable` to a column name or a prefix (e.g., 'temperature' or 'temperature_mean').")
  }
  
  cols_use <- names(adjusted_df)[grepl(variable, names(adjusted_df), ignore.case = TRUE)]
  if (length(cols_use) == 0) stop(sprintf("No columns containing '%s' found in adjusted_df.", variable))
  
  cols_use <- intersect(cols_use, names(pi_df))
  if (length(cols_use) == 0) stop(sprintf("No overlapping '%s' columns between adjusted_df and pi_df.", variable))
  
  if (!is.null(year)) {
    adj_y <- adjusted_df %>%
      filter(.data$year == !!year) %>%
      select(geo, year, dplyr::all_of(cols_use))
    
    pi_y <- pi_df %>%
      filter(.data$year == !!year) %>%
      select(geo, year, dplyr::all_of(cols_use))
    
    if (nrow(adj_y) == 0 || nrow(pi_y) == 0)
      stop("No rows for the requested year in adjusted_df or pi_df.")
    
    df_out <- adj_y %>%
      inner_join(pi_y, by = c("geo","year"), suffix = c("", "__pi"))
    
  } else {
    adj_agg <- adjusted_df %>%
      group_by(geo) %>%
      summarise(across(all_of(cols_use), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    
    pi_agg <- pi_df %>%
      group_by(geo) %>%
      summarise(across(all_of(cols_use), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    
    df_out <- adj_agg %>%
      inner_join(pi_agg, by = "geo", suffix = c("", "__pi"))
  }
  
  diff_cols <- paste0(cols_use, "__diff")
  for (i in seq_along(cols_use)) {
    col <- cols_use[i]
    df_out[[diff_cols[i]]] <- df_out[[col]] - df_out[[paste0(col, "__pi")]]
  }
  
  year_tag <- if (is.null(year)) "aggregated" else as.character(year)
  
  adj_vals_matrix <- as.matrix(df_out %>% select(dplyr::all_of(cols_use)))
  lims_adj <- range(adj_vals_matrix, na.rm = TRUE, finite = TRUE)
  breaks_adj <- pretty(lims_adj, n = 3)
  limits_adj  <- c(min(breaks_adj), max(breaks_adj))
  
  diff_vals_matrix <- as.matrix(df_out %>% select(dplyr::all_of(diff_cols)))
  lims_diff <- range(diff_vals_matrix, na.rm = TRUE, finite = TRUE)
  breaks_diff <- pretty(lims_diff, n = 3)
  limits_diff  <- c(min(breaks_diff), max(breaks_diff))
  
  sf_y <- map_inla %>%
    select(geo, geometry) %>%
    left_join(df_out, by = "geo") %>%
    st_as_sf()
  
  map_inla <- sf::st_as_sf(map_inla, sf_column_name = "geometry")
  
  geom <- sf::st_geometry(map_inla)
  outline <- sf::st_as_sf(sf::st_union(geom))
  
  # prettify.
  prettify <- function(x) {
    gsub("_", " ", x)
  }
  
  for (col in cols_use) {
    
    col_label <- prettify(col)
    
    p_adj <- ggplot() +
      geom_sf(data = sf_y, aes(fill = .data[[col]]), color = border_col, linewidth = border_size) +
      scale_fill_gradientn(colors = pal_seq,
                           limits = limits_adj,
                           breaks = breaks_adj,
                           labels = scales::label_number(accuracy = 1),
                           oob = scales::squish,
                           na.value = na_fill, name = NULL) +
      geom_sf(data = outline,
              fill = NA,
              color = outer_boundary_col,
              linewidth = outer_boundary_size,
              inherit.aes = FALSE) +
      
      coord_sf(datum = NA) +
      labs(title = "") +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none"
      )
    
    diff_col <- paste0(col, "__diff")
    p_diff <- ggplot() +
      geom_sf(data = sf_y, aes(fill = .data[[diff_col]]),
              color = border_col,
              linewidth = border_size) +
      scale_fill_gradientn(colors = pal_seq_diff,
                           limits = limits_diff,
                           breaks = breaks_diff,
                           labels = scales::label_number(accuracy = 1),
                           oob = scales::squish,
                           na.value = na_fill,
                           name = NULL) +
      geom_sf(data = outline,
              fill = NA,
              color = outer_boundary_col,
              linewidth = outer_boundary_size,
              inherit.aes = FALSE) +
      coord_sf(datum = NA) +
      labs(title = "") +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none"
      )
    
    adj_file  <- sprintf("%s__%s__%s.%s", filename_stub,     year_tag, col_label, ext)
    diff_file <- sprintf("%s__%s__%s.%s", filename_stub_diff, year_tag, col_label, ext)
    
    save_plot_device(p_adj,  adj_file,  width_in, height_in, dpi, compression, bg)
    cat("saved:", adj_file, "\n")
    save_plot_device(p_diff, diff_file, width_in, height_in, dpi, compression, bg)
    cat("saved:", diff_file, "\n")
  }
  
  p_adj_leg <- ggplot(sf_y) +
    geom_sf(aes(fill = .data[[cols_use[2]]]), color = NA) +
    scale_fill_gradientn(colors = pal_seq, limits = limits_adj, breaks = breaks_adj,
                         labels = scales::label_number(accuracy = 1),
                         oob = scales::squish, na.value = na_fill, name = NULL) +
    theme_void(base_size = 10) + theme(legend.position = "right")
  hist_grob <- cowplot::get_legend(p_adj_leg)
  
  p_diff_leg <- ggplot(sf_y) +
    geom_sf(aes(fill = .data[[paste0(cols_use[1], "__diff")]]), color = NA) +
    scale_fill_gradientn(colors = pal_seq_diff, limits = limits_diff, breaks = breaks_diff,
                         labels = scales::label_number(accuracy = 1),
                         oob = scales::squish, na.value = na_fill, name = NULL) +
    theme_void(base_size = 10) + theme(legend.position = "right")
  diff_grob <- cowplot::get_legend(p_diff_leg)
  
  if (is.null(hist_legend_filename))
    hist_legend_filename <- sprintf("%s__%s__LEGEND.%s",
                                    filename_stub,
                                    year_tag,
                                    if (is.null(hist_legend_ext)) ext else hist_legend_ext)
  
  if (is.null(diff_legend_filename))
    diff_legend_filename <- sprintf("%s__%s__LEGEND.%s",
                                    filename_stub_diff,
                                    year_tag,
                                    if (is.null(diff_legend_ext)) ext else diff_legend_ext)
  
  if (is.null(hist_legend_dpi)) hist_legend_dpi <- dpi
  if (is.null(diff_legend_dpi)) diff_legend_dpi <- dpi
  
  save_plot_device(cowplot::ggdraw(hist_grob), hist_legend_filename,
                   hist_legend_width_in, hist_legend_height_in,
                   hist_legend_dpi, hist_legend_compression, hist_legend_bg)
  cat("saved legend:", hist_legend_filename, "\n")
  
  save_plot_device(cowplot::ggdraw(diff_grob), diff_legend_filename,
                   diff_legend_width_in, diff_legend_height_in,
                   diff_legend_dpi, diff_legend_compression, diff_legend_bg)
  cat("saved legend:", diff_legend_filename, "\n")
  
  invisible(TRUE)
}



# plotting population maps
plot_population_maps <- function(
    adjusted_df,
    map_inla,
    year,
    variable,
    
    filename_stub,
    ext = "tiff",
    
    width_in = 3.5, height_in = 4,
    dpi = 600,
    compression = "lzw",
    bg = "transparent",
    
    hist_legend_ext         = NULL,
    hist_legend_dpi         = NULL,
    hist_legend_filename    = NULL,
    hist_legend_width_in    = 1.4,
    hist_legend_height_in   = 4.0,
    hist_legend_compression = "lzw",
    hist_legend_bg          = "transparent",
    
    na_fill            = "grey95",
    border_col         = NA,
    border_size        = NA,
    
    outer_boundary_col  = "black",
    outer_boundary_size = 0.1,
    
    palette_name    = "YlGnBu",
    palette_n       = 100,
    palette_reverse = FALSE
) {
  
  stopifnot("geo"  %in% names(adjusted_df),
            "year" %in% names(adjusted_df))
  stopifnot("geometry" %in% names(map_inla))
  
  needed_cols <- c("population", "population_density")
  if (!all(needed_cols %in% names(adjusted_df))) {
    stop("adjusted_df must contain columns: 'population' and 'population_density'.")
  }
  
  allowed_vars <- c("population", "population_density", "log_population_density")
  if (!variable %in% allowed_vars) {
    stop(sprintf("variable must be one of: %s", paste(allowed_vars, collapse = ", ")))
  }
  
  if (!palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
    warning(sprintf("Unknown palette_name = '%s'. Falling back to 'YlGnBu'.",
                    palette_name))
    palette_name <- "YlGnBu"
  }
  
  max_n    <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  pal_base <- RColorBrewer::brewer.pal(max_n, palette_name)
  pal_seq  <- grDevices::colorRampPalette(pal_base)(palette_n)
  
  if (isTRUE(palette_reverse)) {
    pal_seq <- rev(pal_seq)
  }
  
  if (!is.null(year)) {
    df_out <- adjusted_df %>%
      dplyr::filter(.data$year == !!year) %>%
      dplyr::select(geo, year, population, population_density)
    
    if (nrow(df_out) == 0) {
      stop("No rows for the requested year in adjusted_df.")
    }
  } else {
    df_out <- adjusted_df %>%
      dplyr::group_by(geo) %>%
      dplyr::summarise(
        population         = mean(population,         na.rm = TRUE),
        population_density = mean(population_density, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(year = NA_real_)
  }
  
  if (variable == "population") {
    col_plot  <- "population"
    col_label <- "population"
  } else if (variable == "population_density") {
    col_plot  <- "population_density"
    col_label <- "population_density"
  } else {
    df_out <- df_out %>%
      dplyr::mutate(
        log_population_density = dplyr::if_else(
          population_density > 0,
          log(population_density),
          NA_real_
        )
      )
    col_plot  <- "log_population_density"
    col_label <- "log_population_density"
  }
  
  year_tag <- if (is.null(year)) "aggregated" else as.character(year)
  
  vals      <- df_out[[col_plot]]
  lims_adj  <- range(vals, na.rm = TRUE, finite = TRUE)
  breaks_adj <- pretty(lims_adj, n = 5)
  limits_adj <- c(min(breaks_adj), max(breaks_adj))
  
  sf_y <- map_inla %>%
    dplyr::select(geo, geometry) %>%
    dplyr::left_join(df_out, by = "geo") %>%
    sf::st_as_sf()
  
  map_inla <- sf::st_as_sf(map_inla, sf_column_name = "geometry")
  
  geom    <- sf::st_geometry(map_inla)
  outline <- sf::st_as_sf(sf::st_union(geom))
  
  p_adj <- ggplot() +
    geom_sf(
      data = sf_y,
      aes(fill = .data[[col_plot]]),
      color     = border_col,
      linewidth = border_size
    ) +
    scale_fill_gradientn(
      colors = pal_seq,
      limits = limits_adj,
      breaks = breaks_adj,
      labels = scales::label_number(accuracy = 1),
      oob    = scales::squish,
      na.value = na_fill,
      name   = NULL
    ) +
    geom_sf(
      data = outline,
      fill  = NA,
      color = outer_boundary_col,
      linewidth = outer_boundary_size,
      inherit.aes = FALSE
    ) +
    coord_sf(datum = NA) +
    labs(title = "") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "none"
    )
  
  adj_file <- sprintf(
    "%s__%s__%s.%s",
    filename_stub,
    year_tag,
    col_label,
    ext
  )
  
  save_plot_device(p_adj, adj_file, width_in, height_in, dpi, compression, bg)
  cat("saved:", adj_file, "\n")
  
  p_adj_leg <- ggplot(sf_y) +
    geom_sf(aes(fill = .data[[col_plot]]), color = NA) +
    scale_fill_gradientn(
      colors = pal_seq,
      limits = limits_adj,
      breaks = breaks_adj,
      labels = scales::label_number(accuracy = 1),
      oob    = scales::squish,
      na.value = na_fill,
      name   = NULL
    ) +
    theme_void(base_size = 10) +
    theme(legend.position = "right")
  
  hist_grob <- cowplot::get_legend(p_adj_leg)
  
  if (is.null(hist_legend_filename)) {
    hist_legend_filename <- sprintf(
      "%s__%s__LEGEND.%s",
      filename_stub,
      year_tag,
      if (is.null(hist_legend_ext)) ext else hist_legend_ext
    )
  }
  
  if (is.null(hist_legend_dpi)) hist_legend_dpi <- dpi
  
  save_plot_device(
    cowplot::ggdraw(hist_grob),
    hist_legend_filename,
    hist_legend_width_in,
    hist_legend_height_in,
    hist_legend_dpi,
    hist_legend_compression,
    hist_legend_bg
  )
  cat("saved legend:", hist_legend_filename, "\n")
  
  invisible(TRUE)
}




# plot land use maps
plot_land_use_maps <- function(
    adjusted_df,
    map_inla,
    year,
    variable,
    
    filename_stub,
    ext = "tiff",
    
    width_in = 3.5, height_in = 4,
    dpi = 600,
    compression = "lzw",
    bg = "transparent",
    
    hist_legend_ext         = NULL,
    hist_legend_dpi         = NULL,
    hist_legend_filename    = NULL,
    hist_legend_width_in    = 1.4,
    hist_legend_height_in   = 4.0,
    hist_legend_compression = "lzw",
    hist_legend_bg          = "transparent",
    
    na_fill            = "grey95",
    border_col         = NA,
    border_size        = NA,
    
    outer_boundary_col  = "black",
    outer_boundary_size = 0.1,
    
    palette_name    = "YlGnBu",
    palette_n       = 100,
    palette_reverse = FALSE,
    
    pal_seq = NULL
) {
  
  stopifnot("geo"  %in% names(adjusted_df),
            "year" %in% names(adjusted_df))
  stopifnot("geometry" %in% names(map_inla))
  
  allowed_vars <- c(
    "primary_forested_area",
    "primary_non_forested_area",
    "secondary_forested_area",
    "pasture_land",
    "rangelands",
    "urban_land"
  )
  
  if (!variable %in% allowed_vars) {
    stop(sprintf(
      "variable must be one of: %s (croplands is deliberately excluded)",
      paste(allowed_vars, collapse = ", ")
    ))
  }
  
  # Palette setup
  if (!is.null(pal_seq)) {
    if (isTRUE(palette_reverse)) {
      pal_seq <- rev(pal_seq)
    }
  } else {
    if (!palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
      warning(sprintf("Unknown palette_name = '%s'. Falling back to 'YlGnBu'.",
                      palette_name))
      palette_name <- "YlGnBu"
    }
    
    max_n    <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
    pal_base <- RColorBrewer::brewer.pal(max_n, palette_name)
    pal_seq  <- grDevices::colorRampPalette(pal_base)(palette_n)
    
    if (isTRUE(palette_reverse)) {
      pal_seq <- rev(pal_seq)
    }
  }
  
  if (!is.null(year)) {
    df_out <- adjusted_df %>%
      dplyr::filter(.data$year == !!year) %>%
      dplyr::select(geo, year, dplyr::all_of(variable))
    
    if (nrow(df_out) == 0) {
      stop("No rows for the requested year in adjusted_df.")
    }
    
  } else {
    df_out <- adjusted_df %>%
      dplyr::group_by(geo) %>%
      dplyr::summarise(
        !!variable := mean(.data[[variable]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(year = NA_real_)
  }
  
  col_plot  <- variable
  col_label <- variable
  
  year_tag <- if (is.null(year)) "aggregated" else as.character(year)
  
  vals <- df_out[[col_plot]]
  
  max_val <- max(vals, na.rm = TRUE)
  
  if (!is.finite(max_val)) {
    stop("All values for the chosen land-use variable are NA or non-finite.")
  }
  
  max_break <- ceiling(max_val / 0.1) * 0.1
  
  if (max_break == 0) max_break <- 0.1
  
  breaks_adj <- seq(0.0, max_break, by = 0.1)
  limits_adj <- c(0.0, max_break)
  
  sf_y <- map_inla %>%
    dplyr::select(geo, geometry) %>%
    dplyr::left_join(df_out, by = "geo") %>%
    sf::st_as_sf()
  
  map_inla <- sf::st_as_sf(map_inla, sf_column_name = "geometry")
  
  geom    <- sf::st_geometry(map_inla)
  outline <- sf::st_as_sf(sf::st_union(geom))
  
  p_adj <- ggplot() +
    geom_sf(
      data = sf_y,
      aes(fill = .data[[col_plot]]),
      color     = border_col,
      linewidth = border_size
    ) +
    scale_fill_gradientn(
      colors  = pal_seq,
      limits  = limits_adj,
      breaks  = breaks_adj,
      labels  = scales::label_number(accuracy = 0.1),
      oob     = scales::squish,
      na.value = na_fill,
      name    = NULL
    ) +
    geom_sf(
      data = outline,
      fill  = NA,
      color = outer_boundary_col,
      linewidth = outer_boundary_size,
      inherit.aes = FALSE
    ) +
    coord_sf(datum = NA) +
    labs(title = "") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "none"
    )
  
  adj_file <- sprintf(
    "%s__%s__%s.%s",
    filename_stub,
    year_tag,
    col_label,
    ext
  )
  
  save_plot_device(p_adj, adj_file, width_in, height_in, dpi, compression, bg)
  cat("saved:", adj_file, "\n")
  
  p_adj_leg <- ggplot(sf_y) +
    geom_sf(aes(fill = .data[[col_plot]]), color = NA) +
    scale_fill_gradientn(
      colors  = pal_seq,
      limits  = limits_adj,
      breaks  = breaks_adj,
      labels  = scales::label_number(accuracy = 0.1),
      oob     = scales::squish,
      na.value = na_fill,
      name    = NULL
    ) +
    theme_void(base_size = 10) +
    theme(legend.position = "right")
  
  hist_grob <- cowplot::get_legend(p_adj_leg)
  
  if (is.null(hist_legend_filename)) {
    hist_legend_filename <- sprintf(
      "%s__%s__LEGEND.%s",
      filename_stub,
      year_tag,
      if (is.null(hist_legend_ext)) ext else hist_legend_ext
    )
  }
  
  if (is.null(hist_legend_dpi)) hist_legend_dpi <- dpi
  
  save_plot_device(
    cowplot::ggdraw(hist_grob),
    hist_legend_filename,
    hist_legend_width_in,
    hist_legend_height_in,
    hist_legend_dpi,
    hist_legend_compression,
    hist_legend_bg
  )
  cat("saved legend:", hist_legend_filename, "\n")
  
  invisible(TRUE)
}
