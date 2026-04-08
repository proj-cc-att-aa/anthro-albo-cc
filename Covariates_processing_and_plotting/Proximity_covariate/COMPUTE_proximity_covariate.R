# ================== Spatial proximity covariate construction ==================

setwd(".../VectorNet_ 2_original")

# Source project functions
source("functions_general.R")
source("functions_plotting.R")
source("functions_risk_from_proximity.R")
source("functions_spatial_networks.R")
source("functions_matrix_operations.R")
source("functions_plot_article.R")

# Load input data
out <- attach("Data/data_INLA.RData")
map_inla <- map_inla
detach("file:Data/data_INLA.RData")

out <- load("Data/mosq_dist_until_2023.RData")
mosq_dist <- mosq_dist %>% filter(date_start>=2010) %>%
  mutate(year_idx = match(date_start, sort(unique(date_start))))

map <- map_inla

# Build the neighbourhood structure
geo_neighs <- geo_to_neighs(map)
plot_geo_neighs(geo_neighs)

plot_geo_neighs_article(geo_neighs,
                        filename = ".../Plots/A_final_plots/Networks/graph_network.tiff",
                        width_in = 3.5,
                        height_in = 4,
                        dpi = 600,
                        compression = "lzw")




# Construct precision matrices (besag neighbourhood structure)
Q_besag <-  compute_Q_besag(geo_neighs$nb_map, geo_neighs$n_geo, proper_version = TRUE, lambda = 0.5)

plot_Q_network(Q_besag, geo_neighs$coords)

plot_Q_network_article(
  Q = Q_besag,
  coords = geo_neighs$coords,
  filename = ".../Plots/A_final_plots/Networks/graph_ICAR_neigh.tiff",
  width_in = 3.5, height_in = 4, dpi = 600, compression = "lzw"
)


# Construct precision matrices (centroid distance based neighbourhood structure)
Q_besag1 <- compute_Q_besag_centroid_dist_weights(geo_neighs$nb_map, geo_neighs$coords, geo_neighs$n_geo, proper_version = FALSE, lambda = NULL)

n <- nrow(Q_besag1)

off <- Q_besag1
diag(off) <- NA

any_zero_offdiag <- any(off == 0, na.rm = TRUE)

where_zero_offdiag <- which(off == 0, arr.ind = TRUE)
head(where_zero_offdiag)

out<- plot_Q_network1(Q_besag1, geo_neighs$coords, edge_colour_by_link_weights = TRUE)

plot_Q_network_centroid_article(
  Q      = Q_besag1,
  coords = geo_neighs$coords,
  filename = ".../Plots/A_final_plots/Networks/graph_centroid_based_percentile.tiff",
  width_in = 4.5, height_in = 5, dpi = 600, compression = "lzw",
  sc_link_width = 2,
  sc_node_size  = 2,
  edge_colour_by_link_weights = TRUE,
  legend_title = "Link strength (percentile)"
)


# for proximity calculation we can't leave any data as NA as it may be neighbouring region, so we assign all the NA region as 0 for calculation of proximity
tb_data <- mosq_dist %>%
  mutate(y_to_convolve = case_when(presence==1 ~ 1,
                                   TRUE ~ 0))

y_column <- "y_to_convolve"


### raising a matrix Q to a power k (e.g.,Q^𝑘) affects the strength (or magnitude) of its eigenvalues.
exponent_vec <- c(0.5, 1, 2, 3)

## For centroid distance based proximity
Q_with_eig <- collect_Q_and_eig(Q_besag1)

# Compute proximity covariates
df_proximity <- proximity_from_Q(Q_with_eig, tb_data, y_column, map_inla, exponent_vec)

df_proximity <- df_proximity %>% rename(presence_previous_year=y)

df_surrounded <- get_surrounded_regions()     ## surrounded regions in list (inner and outer) are for those areas when a NUTS3 region completely lies inside another region (like big cities Berlin city NUTS3 is inside Berlin overall NUTS3)

proximity_col_vec <- names(df_proximity) %>% str_subset("y_proximity")


# We implemented a hierarchical spatial rule: inner regions inherit the proximity of their corresponding outer region until either region experiences its first invasion (i.e., presence = 1 in the previous year).
#This approach assumes that, prior to invasion, inner and outer regions are indistinguishable in their exposure to diffusion pressure and should share the same proximity value. Once the outer region becomes invaded, the inner region is allowed to evolve independently to reflect potential local delay, resistance, or stochastic variation in invasion dynamics. 
## Conversely, if the inner region is invaded first—a less common but plausible scenario—the inner proximity is preserved to reflect true spatial information and override prior-based assumptions.

for (j in 1:nrow(df_proximity)) {
  df_curr <- df_proximity[j,]
  if (df_curr$geo %in% df_surrounded$inner) {
    
    geo_outer <- df_surrounded %>% filter(inner == df_curr$geo) %>% pull(outer)
    df_curr_match <- df_proximity %>% filter(geo == geo_outer, year_idx == df_curr$year_idx)
    
    if (
      df_curr$presence_previous_year != 1 &&      ## The inner region mirrors the outer region only as long as it has not been invaded and the outer region has not been invaded
      df_curr_match$presence_previous_year != 1
    ) {
      for (proximity_col in proximity_col_vec) {
        df_proximity[[proximity_col]][j] <- df_curr_match[[proximity_col]]
      }
    }
  }
}# This refined rule preserves spatial diffusion continuity, avoids premature sink behavior in terminal nodes, and maintains coherence between nested regions in the absence of direct centroid-based connectivity.


## saving the data
proximity_goal = "neighbourhood"
#proximity_goal = "distance"     ## define the covariance structure (depending on Q_besag or Q_besag1 in above code)

if (proximity_goal == "distance") {
  save(df_proximity, file = ".../Data/climate/Final_covariates_data/Current_data/Final_hist_proximity_data_yearly.RData")
} else if (proximity_goal == "neighbourhood") {
  save(df_proximity, file = ".../Data/climate/Final_covariates_data/Current_data/Final_hist_ICAR_proximity_data_yearly.RData")
} else {
  warning("Unknown proximity_goal. File not saved.")
}




########## Plotting for final manuscript ####################
out_proximity <- load(".../Data/climate/Final_covariates_data/Current_data/Final_hist_proximity_data_yearly.RData")
df_proximity <- df_proximity



## This function will plot the proximity for any year but for all 4 proximity strength in 4 seperate plots (or any number as specified in value_cols), but will use the common legend and scaling (across all plots)
plot_save_proximities_for_year(
  df_proximity = df_proximity,
  map_inla     = map_inla,
  year_idx     = 5,
  
  value_cols   = c("y_proximity_exp050","y_proximity_exp100","y_proximity_exp200","y_proximity_exp300"),
  
  filename_stub = ".../Plots/A_final_plots/one_year_proximity_strength_all_values/Distance/proximity_year_2014",
  
  ext           = "tiff",
  width_in      = 3.5, height_in = 4, dpi = 600, compression = "lzw",
  
  legend_grad_filename     = ".../Plots/A_final_plots/one_year_proximity_strength_all_values/Distance/proximity_year_2014__LEGEND_GRAD.tiff",
  
  legend_presence_filename = ".../Plots/A_final_plots/one_year_proximity_strength_all_values/Distance/proximity_year_2014__LEGEND_Previous_presence.tiff",
  
  legend_grad_width_in     = 1.4, legend_grad_height_in = 5,
  legend_presence_width_in = 1.4, legend_presence_height_in = 0.9,
  legend_dpi               = 600, legend_ext = "tiff",
  
  presence_fill = "grey50",
  
  border_col = NA, border_size = NA,
  
  outer_boundary_col  = "black",
  outer_boundary_size = 0.1,
  na_fill = "grey95",
  palette_name = "Zissou1", palette_n = 100, palette_reverse = FALSE
)



## This plot code is same as the above code but now we will plot only one proximity strength, used to compare the centroid distance based proximity with neighbourhood basesd 
## Care is to be taken that since both (besag or centroid) are plotted separately, so there might be legend mismatch while comparing
plot_save_proximities_for_year(
  df_proximity = df_proximity,
  map_inla     = map_inla,
  year_idx     = 5,
  
  value_cols   = c("y_proximity_exp050"),
  
  filename_stub = ".../Plots/A_final_plots/Compare_both_proximity_method_for_oneyear/Neighbourhood_based/proximity_year_2014_neigh",
  
  ext           = "tiff",
  width_in      = 3.5, height_in = 4, dpi = 600, compression = "lzw",
  
  legend_grad_filename = ".../Plots/A_final_plots/Compare_both_proximity_method_for_oneyear/Neighbourhood_based/proximity_year_2014_neigh_LEGEND_GRAD.tiff",
  
  legend_presence_filename = ".../Plots/A_final_plots/Compare_both_proximity_method_for_oneyear/Neighbourhood_based/proximity_year_2014_neigh_LEGEND_Previous_presence.tiff",
  
  legend_grad_width_in     = 1.4, legend_grad_height_in = 5,
  legend_presence_width_in = 1.4, legend_presence_height_in = 0.9,
  legend_dpi               = 600, legend_ext = "tiff",
  
  presence_fill = "grey50",
  
  border_col = NA, border_size = NA,
  
  outer_boundary_col  = "black",
  outer_boundary_size = 0.1,
  na_fill = "grey95",
  
  palette_name = "Zissou1", palette_n = 100, palette_reverse = FALSE
)




# function plot one proximity strength for all years (1 to 15 in our case) with common legend
plot_save_proximity_all_years(
  df_proximity = df_proximity,
  map_inla     = map_inla,
  value_col    = "y_proximity_exp300",
  years        = 1:15,
  
  filename_stub = ".../Plots/A_final_plots/one_strength_prox_for_all_years/Distance/proximity_exp300",
  
  ext           = "tiff",
  width_cm = 5.7,
  height_cm = 8.5,
  dpi = 600, compression = "lzw",
  
  legend_grad_filename     = ".../Plots/A_final_plots/one_strength_prox_for_all_years/Distance/proximity_exp300__LEGEND_GRAD.tiff",
  
  legend_presence_filename = ".../Plots/A_final_plots/one_strength_prox_for_all_years/Distance/proximity_exp300__LEGEND_Previous_presence.tiff",
  
  legend_grad_width_cm     = 2.8, legend_grad_height_cm    = 8.5,
  legend_presence_width_cm  = 2.8, legend_presence_height_cm = 2.3,
  legend_dpi               = 600, legend_ext = "tiff",
  
  presence_fill = "grey50",
  
  border_col = NA, border_size = NA,
  na_fill = "grey95",
  
  outer_boundary_col  = "black",
  outer_boundary_size = 0.1,
  
  palette_name = "Zissou1", palette_n = 100, palette_reverse = FALSE
)




### NOW combined panel plot (Supp plot Fig S.13 ) for comparison in manuscript
base_dir <- ".../Plots/A_final_plots/one_strength_prox_for_all_years/Distance"

panel_files <- sprintf(
  file.path(base_dir, "proximity_exp300_y_proximity_exp300_year%02d.tiff"),
  1:15
)
legend_file <- file.path(base_dir, "proximity_exp300__LEGEND_GRAD.tiff")

panel_plots <- lapply(panel_files, function(f) {
  img <- image_read(f) |> image_trim()
  ggdraw() + draw_image(img)
})

grid_5x3 <- plot_grid(
  plotlist       = panel_plots,
  ncol           = 3,
  labels         = LETTERS[1:15],
  label_size     = 10,
  label_fontface = "bold",
  label_x        = 0.02,
  label_y        = 0.98,
  hjust          = 0,
  vjust          = 1
)

legend_plot <- {
  leg_img <- image_read(legend_file) |> image_trim()
  ggdraw() + draw_image(leg_img)
}

combined_with_legend <- plot_grid(
  grid_5x3,
  legend_plot,
  ncol       = 2,
  rel_widths = c(1, 0.08)
)

out_dir <- file.path(base_dir, "Science_journal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename    = file.path(out_dir, "proximity_15years_5x3_with_legend_3col.tiff"),
  plot        = combined_with_legend,
  width       = 18.4,
  height      = 25,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  device      = "tiff",
  bg          = "transparent"
)

ggsave(
  filename = file.path(out_dir, "proximity_15years_5x3_with_legend_3col.svg"),
  plot     = combined_with_legend,
  width       = 18.4,
  height      = 25,
  units    = "cm",
  dpi         = 600,
  device   = "svg",
  bg       = "transparent"
)
