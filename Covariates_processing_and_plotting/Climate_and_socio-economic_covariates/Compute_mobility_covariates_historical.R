# ================== Mobility covariate construction ==================

# Load functions and supporting objects
source("functions_general.R")
source("functions_plotting.R")
source("functions_risk_from_proximity.R")
source("functions_spatial_networks.R")
source("functions_matrix_operations.R")
source("functions_plot_article.R")

folder_data <- "Data"

# Load population and spatial data
out <- load(".../Current_data/Final_hist_population_data_yearly.RData")
df_population <- population

out <- attach(".../Data/data_INLA.RData")
map_inla <- map_inla
detach(".../file:Data/data_INLA.RData")

year_vec_pop <- sort(unique(df_population$year_idx))

# Compute yearly mobility summaries
df_mobility <- NULL
for (yr in year_vec_pop) {

  df_pop_year <- df_population %>%
    filter(year_idx == yr) %>%
    select(geo, population)

  map_inla_year <- map_inla %>%
    mutate(year_idx = yr)

  map_joined <- map_inla_year %>%
    left_join(df_pop_year, by = "geo")

  mobilty_flux_yearly <- mobility_from_population(map_joined)

  mobility_plot_yearly <- plot_mobility_graph(mobilty_flux_yearly, nodes = "net_flux_abs")

  #### Saving the plot for for year 2014 only
  if (yr == 5) {
    plot_mobility_graph_article(
      mobility    = mobilty_flux_yearly,
      nodes       = "net_flux_abs",
      filtering_p = 0.99,
      filename    = ".../Networks/mobility_network.tiff",
      width_in    = 3.5,
      height_in   = 4,
      dpi         = 600,
      compression = "lzw"
    )
  }

  mobilty_flux_yearly$tb_regions <- mobilty_flux_yearly$tb_regions %>%
    mutate(net_flux = flux_out - flux_in)

  df_mobility <-  df_mobility %>%
    bind_rows(tibble(geo = mobilty_flux_yearly$tb_regions$geo,
                     mobility_idx = mobilty_flux_yearly$tb_regions$mobility_idx,
                     year_idx = yr,
                     net_flux_out =  mobilty_flux_yearly$tb_regions$net_flux,
                     total_flux_in = mobilty_flux_yearly$tb_regions$flux_in))

}

# Remove residual vector names
df_mobility <- df_mobility |>
  mutate(total_flux_in = unname(total_flux_in))

df_mobility <- df_mobility |>
  mutate(net_flux_out = unname(net_flux_out))

# Save the mobility covariate
save(df_mobility, file=".../Current_data/Final_hist_mobility_data_yearly.RData")







# Plot mobility summaries
##  Line plot for total influx mobility () of some region for all years 
out_mobility <- load(".../Current_data/Final_hist_mobility_data_yearly.RData")
df_mobility <- df_mobility

out <- attach(".../Data/data_INLA.RData")
map_inla <- map_inla
detach(".../file:Data/data_INLA.RData")

top10_in_y1 <- df_mobility %>%
  filter(year_idx == 1) %>%
  arrange(desc(total_flux_in)) %>%
  slice_head(n = 10) %>%
  select(geo, mobility_idx, year_idx, total_flux_in, net_flux_out)

top10_in_y1

top10_geos <- top10_in_y1 %>%
  distinct(geo) %>%
  pull(geo)

top10_all_years_wide <- df_mobility %>%
  filter(geo %in% top10_geos) %>%
  select(geo, year_idx, total_flux_in) %>%
  left_join(map_inla %>% select(geo, LocationNa), by = "geo") %>%
  pivot_wider(names_from = year_idx, values_from = total_flux_in, names_prefix = "year_") %>%
  arrange(desc(year_1))

top10_all_years_wide


# drop Paris + Hauts-de-Seine
mob_top8_long <- top10_all_years_wide %>%
  filter(!LocationNa %in% c("Paris", "Hauts-de-Seine")) %>%
  pivot_longer(
    cols = starts_with("year_"),
    names_to = "year_idx",
    values_to = "total_flux_in"
  ) %>%
  mutate(
    year_idx = as.integer(sub("year_", "", year_idx)),
    year = 2009 + year_idx,
    total_flux_k = total_flux_in / 1e6
  ) %>%
  filter(year >= 2011, year <= 2023)



## theme
theme_science_ts <- function(legend_pos = c(0.88, 0.30)) {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial"),

      axis.title = ggplot2::element_text(size = 5),
      axis.text  = ggplot2::element_text(size = 5, colour = "black"),

      axis.line  = ggplot2::element_line(linewidth = 0.2, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.2, colour = "black"),

      panel.grid = ggplot2::element_blank(),

      legend.background     = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.key            = ggplot2::element_blank(),

      legend.title = ggplot2::element_text(
        size = 3, face = "bold",
        margin = ggplot2::margin(b = 1.5, unit = "pt")
      ),
      legend.text = ggplot2::element_text(size = 3, margin = margin(l = 1.2, unit = "pt"), colour = "black"),

      legend.key.size   = grid::unit(0.15, "cm"),
      legend.key.height = grid::unit(0.10, "cm"),

      legend.position = legend_pos
    )
}

p_mob_top8 <- ggplot2::ggplot(
  mob_top8_long,
  ggplot2::aes(
    x = year,
    y = total_flux_k,
    group = geo,
    color = LocationNa
  )
) +
  ggplot2::geom_line(linewidth = 0.3) +

  ggplot2::scale_x_continuous(
    limits = c(2011, 2023),
    breaks = seq(2011, 2023, by = 3),
    expand = c(0, 0)
  ) +

  ggplot2::scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    limits = c(1.6, 2.4),
    breaks = seq(1.6, 2.4, 0.2),
  ) +

  ggplot2::labs(
    x = "Year",
    y = "Total incoming mobility flux (millions)",
    color = "Region"
  ) +

  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1)) +

  theme_science_ts(legend_pos = c(0.15, 0.88))

p_mob_top8



ggsave(
  filename    = ".../Mobility_plot/mobility_top8_timeseries.tiff",
  plot        = p_mob_top8,
  width       = 5.7,
  height      = 5.7,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "transparent",
  device      = "tiff"
)

ggsave(
  filename = ".../Mobility_plot/mobility_top8_timeseries.svg",
  plot     = p_mob_top8,
  width    = 5.7,
  height   = 5.7,
  units    = "cm",
  device   = "svg",
  bg       = "transparent"
)
