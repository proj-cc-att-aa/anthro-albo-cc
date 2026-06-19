# Setting our working directory
setwd(".../VectorNet_ 2_original")

# Loading our packages and plotting functions

packages <- c("tidyverse", "stringr", "caret")
lapply(packages, library, character.only = TRUE)

library(dplyr)
library(raster)
library(tidyverse)
library(ggplot2)
library(mgcv)
library(dplyr)
library(rgdal)
library(sf)
library(sp)
library(tmap)
library(terra)
library(RColorBrewer)
library(scales)
library(ggnewscale)
library(data.table)

source("functions_plotting.R")

# Loading our factual model output

m_train = load(".../VectorNet_ 2_original/Output/A-Final_output/Main_model_output/main_modelvIparam110000101011001110111100_SIMP_idvI.Rdata")
m_train

train_post_samples <- p_samples_long


# Loading our counterfactual model output
m = load(".../VectorNet_ 2_original/Output/A-Final_output/Main_model_output/counterfactual_prediction_2011_2021_output.RData")
m

pred_post_samples <- p_post_pred_all

fit_tbl <- collected_results$map_full_model
pred_tbl <- map_full_model_pred_all

# Restricting our posterior samples to 2011–2021
train_post_2011_2021 <- train_post_samples %>%
  filter(year_idx <= 11)

pred_post_2011_2021 <- pred_post_samples %>%
  filter(year_idx <= 11)

# Joining our factual and counterfactual posterior draws
df_post <- train_post_2011_2021 %>%
  rename(p_fact = p) %>%
  inner_join(
    pred_post_2011_2021 %>% rename(p_cf = p),
    by = c("geo", "region_id", "year_idx", "draw")
  ) %>%
  mutate(

    p_fact = pmin(pmax(p_fact, 0), 1),
    p_cf   = pmin(pmax(p_cf,   0), 1),

  )

df_post <- df_post %>%
  dplyr::relocate(year, .after = year_idx) %>%
  dplyr::select(-data_idx.y) %>%
  dplyr::rename(data_idx = data_idx.x)

df_post


# Summarising our spatial attribution across years and posterior draws
eps <- 1e-9

spatial_draws <- df_post %>%
  group_by(draw, region_id, geo) %>%
  summarise(
    T_obs     = n(),
    
    E_fact    = mean(p_fact, na.rm = TRUE),
    E_cf      = mean(p_cf,   na.rm = TRUE),
    AF        = ifelse(E_fact > 0, (E_fact - E_cf) , NA_real_),
    .groups   = "drop"
  )


#summarising across draws for each region to get posterior mean + 95% CrI for absolute spatial burden.
spatial_abs_post <- spatial_draws %>%
  group_by(region_id, geo) %>%
  summarise(
    T_obs           = first(T_obs),
    AF_abs_mean  = mean( AF , na.rm = TRUE),
    AF_abs_med   = median( AF , na.rm = TRUE),
    AF_abs_lwr   = quantile( AF , 0.025, na.rm = TRUE),
    AF_abs_upr   = quantile( AF , 0.975, na.rm = TRUE),

    AF_abs_pct_mean = 100 * AF_abs_mean,
    AF_abs_pct_lwr  = 100 * AF_abs_lwr,
    AF_abs_pct_upr  = 100 * AF_abs_upr,

    .groups         = "drop"
  )

spatial_abs_post



###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## FIGURE 1 for article (Panel A and B)
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################

orig_data = load("Data/data_INLA.RData")
nuts3_map = map_inla

World_adm0 <- st_read('good_visual_plots_code/map_and_data_for_plots/WB_Land_10m')

presence_2021 <- fit_tbl[fit_tbl$year_idx==11, ]

# Mapping our 2021 invasion status
nuts3_map <- st_as_sf(nuts3_map, crs=4326)

map_pres_2021 <- nuts3_map %>%
  left_join(presence_2021 %>% dplyr::select(geo, y), by = "geo") %>%
  mutate(
    pres_cat = dplyr::case_when(
      is.na(y)        ~ "No data",
      y == 0          ~ "Absent",
      y == 1          ~ "Present",
      TRUE            ~ "No data"
    ),
    pres_cat = factor(pres_cat, levels = c("Present","Absent", "No data"))
  )

bb <- sf::st_bbox(nuts3_map)

dx  <- as.numeric(bb["xmax"] - bb["xmin"])
dy  <- as.numeric(bb["ymax"] - bb["ymin"])
pad <- 0.2

xlim_pad <- c(bb["xmin"] - pad*dx, bb["xmax"] + pad*dx)
ylim_pad <- c(bb["ymin"] - pad*dy, bb["ymax"] + pad*dy)

gg_presence_2021 <- ggplot() +

  geom_sf(
    data  = World_adm0,
    fill  = rgb(215/255, 215/255, 215/255),
    color = "black",
    linewidth = 0.3
  ) +

  geom_sf(
    data  = map_pres_2021,
    aes(fill = pres_cat),
    color = NA
  ) +
  scale_fill_manual(
    values = c(
      "Present"  = "darkred",
      "Absent"   = "darkblue",
      "No data"  = rgb(150/255, 150/255, 150/255)
    ),
    name  = "Invasion status\n (2021)",
    drop  = FALSE
  ) +
  coord_sf(
    xlim   = xlim_pad,
    ylim   = ylim_pad,
    expand = FALSE
  ) +
  theme_void() +
  theme(
    plot.margin     = margin(8, 8, 8, 8),
    legend.position = c(0.93, 0.22),
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 9)
  )

gg_presence_2021

# Saving our full-size invasion-status map
ext         <- "tiff"
width_in    <- 8
height_in   <- 7
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_1_plots"

ggsave(
  filename    = file.path(out_dir_save, paste0("ECDC_presence_absence_2021.", ext)),
  plot        = gg_presence_2021,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  bg          = bg
)

# Applying the journal formatting to the invasion-status map
gg_presence_2021_sci <- ggplot() +

  geom_sf(
    data  = World_adm0,
    fill  = rgb(215/255, 215/255, 215/255),
    color = "black",
    linewidth = 0.07
  ) +

  geom_sf(
    data  = map_pres_2021,
    aes(fill = pres_cat),
    color = NA,
    linewidth = 0
  ) +
  scale_fill_manual(
    values = c(
      "Present"  = "darkred",
      "Absent"   = "darkblue",
      "No data"  = rgb(150/255, 150/255, 150/255)
    ),
    name  = "Invasion status\n (2021)",
    drop  = FALSE
  ) +
  coord_sf(
    xlim   = xlim_pad,
    ylim   = ylim_pad,
    expand = FALSE
  ) +
  theme_void() +
  theme(
    text = element_text(family = "Arial"),

    legend.title    = element_text(size = 4.5, face = "bold", margin = margin(b = 1.5), hjust = 0.2),
    legend.text     = element_text(size = 3.5, margin = margin(l = 1.2, unit = "pt")),

    legend.position = c(0.88, 0.78),

    legend.justification = "center",
    legend.box          = "vertical",
    legend.box.just     = "center",

    legend.key.size = unit(0.15, "cm"),
    legend.key.width = unit(0.15, "cm"),

    legend.box.spacing = unit(0, "cm"),
    legend.margin     = margin(0, 0, 0, 0),
    plot.margin       = margin(0, 0, 0, 0)
  )

gg_presence_2021_sci

out_dir_base <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_1_plots"
out_dir_sci  <- file.path(out_dir_base, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 5.0

ggplot2::ggsave(
  filename    = file.path(out_dir_sci, "ECDC_presence_absence_2021.tiff"),
  plot        = gg_presence_2021_sci,
  device      = ragg::agg_tiff,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

ggplot2::ggsave(
  filename = file.path(out_dir_sci, "ECDC_presence_absence_2021.svg"),
  plot     = gg_presence_2021_sci,
  device   = "svg",
  width    = width_cm,
  height   = height_cm,
  units    = "cm",
  bg       = "white"
)



# Mapping our climate-attributable probability among invaded regions
map_pres_AF <- map_pres_2021 %>%
  left_join(spatial_abs_post %>% dplyr::select(geo, AF_abs_pct_mean, AF_abs_pct_lwr, AF_abs_pct_upr), by = "geo")

bb  <- st_bbox(nuts3_map)
dx  <- as.numeric(bb["xmax"] - bb["xmin"])
dy  <- as.numeric(bb["ymax"] - bb["ymin"])
pad <- 0.2

xlim_pad <- c(bb["xmin"] - pad*dx, bb["xmax"] + pad*dx)
ylim_pad <- c(bb["ymin"] - pad*dy, bb["ymax"] + pad*dy)

brks_present <- pretty(map_pres_AF$AF_abs_pct_mean[map_pres_AF$pres_cat == "Present"], 5)
brks_present <- sort(unique(c(brks_present, 0)))

gg_AF_present <- ggplot() +

  geom_sf(
    data  = World_adm0,
    fill  = rgb(215/255, 215/255, 215/255),
    color = "black",
    linewidth = 0.3
  ) +

geom_sf(
  data  = map_pres_AF %>% dplyr::filter(pres_cat == "Present"),
  aes(fill = AF_abs_pct_mean),
  color = NA
) +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    limits    = NULL,
    breaks    = brks_present,
    labels    = scales::label_number(accuracy = 1),
    name      = "Attributable\npercent (%)",
    na.value  = "grey85"
  ) +
  new_scale_fill() +

geom_sf(
  data  = map_pres_AF %>% dplyr::filter(pres_cat != "Present"),
  aes(fill = pres_cat),
  color = NA
) +
  scale_fill_manual(
    values = c(
      "Absent"  = "darkblue",
      "No data" = rgb(150/255, 150/255, 150/255)
    ),
    name = "Invasion status",
    drop = FALSE
  ) +
  coord_sf(
    xlim   = xlim_pad,
    ylim   = ylim_pad,
    expand = FALSE
  ) +
  theme_void() +
  theme(
    plot.margin     = margin(8, 8, 8, 8),
    legend.position = c(0.93, 0.22),
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 9)
  )

gg_AF_present

# Saving our full-size attribution map

ext         <- "tiff"
width_in    <- 8
height_in   <- 7
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_1_plots"

ggsave(
  filename    = file.path(out_dir_save, paste0("Attributable_percent_2021.", ext)),
  plot        = gg_AF_present,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  bg          = bg
)


# Applying journal formatting to the attribution map
gg_AF_present_sci <- ggplot() +
  geom_sf(
    data  = World_adm0,
    fill  = rgb(215/255, 215/255, 215/255),
    color = "black",
    linewidth = 0.07
  ) +

  geom_sf(
    data  = map_pres_AF %>% dplyr::filter(pres_cat == "Present"),
    aes(fill = AF_abs_pct_mean),
    color = NA,
    linewidth = 0
  ) +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    breaks    = brks_present,
    labels    = scales::label_number(accuracy = 1),
    name      = "Attributable\npercent (%)",
    na.value  = "grey85",
  ) +
  new_scale_fill() +

  geom_sf(
    data  = map_pres_AF %>% dplyr::filter(pres_cat != "Present"),
    aes(fill = pres_cat),
    color = NA,
    linewidth = 0
  ) +
  scale_fill_manual(
    values = c("Absent" = "darkblue", "No data" = rgb(150/255, 150/255, 150/255)),
    name = "Invasion\nstatus",
    drop = FALSE
  ) +
  coord_sf(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) +
  theme_void() +
  theme(
    text = element_text(family = "Arial"),

    legend.title    = element_text(size = 4.5, face = "bold", margin = margin(b = 1.5), hjust = 0.2),
    legend.text     = element_text(size = 3.5, margin = margin(l = 1.2, unit = "pt")),

    legend.position = c(0.90, 0.78),

    legend.justification = "center",
    legend.box          = "vertical",
    legend.box.just     = "center",

    legend.key.size = unit(0.15, "cm"),
    legend.key.width = unit(0.15, "cm"),
    legend.spacing.y = unit(0.2, "cm"),

    legend.box.spacing = unit(0, "cm"),
    legend.margin     = margin(0, 0, 0, 0),
    plot.margin       = margin(0, 0, 0, 0)
  ) +

  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.2),
         fill_new = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)
  )

out_dir_base <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_1_plots"
out_dir_sci  <- file.path(out_dir_base, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 5.0

ggplot2::ggsave(
  filename    = file.path(out_dir_sci, "Attributable_percent_2021.tiff"),
  plot        = gg_AF_present_sci,
  device      = ragg::agg_tiff,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

ggplot2::ggsave(
  filename = file.path(out_dir_sci, "Attributable_percent_2021.svg"),
  plot     = gg_AF_present_sci,
  device   = "svg",
  width    = width_cm,
  height   = height_cm,
  units    = "cm",
  bg       = "white"
)



###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure 1 (Panel C) for article 
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
orig_data = load("Data/data_INLA.RData")
nuts3_map = map_inla

nuts3_map <- st_as_sf(nuts3_map, crs=4326)
nuts3_map

fit_tbl

df_post

out <- load("Data/climate/Final_covariates_data/Current_data/Final_hist_population_data_yearly.RData")
df_population <- population

year_map <- unique(df_post[, c("year_idx","year")])

# Matching our observed invasion records to the historical prediction period

ecdc_tbl <- fit_tbl %>%
  left_join(year_map, by = "year_idx") %>%
  transmute(geo, year, y) %>%
  filter(
    !is.na(y),
    !is.na(year)
  )

# Restricting our population data to 2011–2021

pop_tbl <- df_population %>%
  transmute(
    geo  = as.character(geo),
    year = as.integer(year),
    population
  ) %>%
  filter(year %in% 2011:2021)

DT_post <- as.data.table(df_post[, c("geo","year","draw","p_fact","p_cf")])
DT_post[, geo  := as.character(geo)]
DT_post[, year := as.integer(year)]
setkey(DT_post, geo, year)

DT_pop <- as.data.table(df_population)
DT_pop <- DT_pop[, data.table(
  geo        = as.character(geo),
  year       = as.integer(year),
  population = as.numeric(population)
)]
DT_pop <- DT_pop[year %in% 2011:2021]
setkey(DT_pop, geo, year)

DT <- DT_post[DT_pop, nomatch = 0L]

DT[, `:=`(
  fact_invaded = (p_fact >= 0.5),
  cf_invaded   = (p_cf   >= 0.5)
)]

DT[, at_risk := (fact_invaded & !cf_invaded)]

pop_risk_draw_year <- DT[, .(
  pop_at_risk = sum(population[at_risk], na.rm = TRUE)
), by = .(draw, year)]

pop_risk_post <- pop_risk_draw_year[, .(
  pop_mean = mean(pop_at_risk, na.rm = TRUE),
  pop_lwr  = quantile(pop_at_risk, 0.025, na.rm = TRUE),
  pop_upr  = quantile(pop_at_risk, 0.975, na.rm = TRUE)
), by = year][order(year)]



pop_risk_post_hist_plot <- as.data.frame(pop_risk_post) %>%
  mutate(
    pop_mean_m = pop_mean / 1e6,
    pop_lwr_m  = pop_lwr  / 1e6,
    pop_upr_m  = pop_upr  / 1e6
  )

p_pop_risk_attributed_cont <- ggplot(pop_risk_post_hist_plot, aes(x = year)) +
  geom_ribbon(aes(ymin = pop_lwr_m, ymax = pop_upr_m), alpha = 0.2) +
  geom_line(aes(y = pop_mean_m), linewidth = 1) +
  scale_x_continuous(
    limits = range(pop_risk_post_hist_plot$year, na.rm = TRUE),
    breaks = pop_risk_post_hist_plot$year,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Year",
    y = "Expected attributable population at risk (millions)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    panel.grid  = element_blank(),
    legend.position = "none",
    axis.title  = element_text(size = 9),
    axis.text   = element_text(size = 9),
    plot.margin = margin(t = 5.5, r = 14, b = 5.5, l = 5.5)
  )

p_pop_risk_attributed_cont

# Estimating our observed population exposure from ECDC invasion records (region-years where ECDC status is observed (not NA))
DT_ecdc <- as.data.table(ecdc_tbl)
setkey(DT_ecdc, geo, year)

DT_pop <- as.data.table(pop_tbl)
setkey(DT_pop, geo, year)

DT_ecdc_obs <- copy(DT_ecdc)
DT_pop_obs  <- copy(DT_pop)

DT_ecdc_pop <- DT_ecdc_obs[DT_pop_obs, nomatch = 0L]

pop_ecdc_year <- DT_ecdc_pop[, .(
  pop_ecdc = sum(population[y == 1], na.rm = TRUE)
), by = year][order(year)]

# Combining our attributed and observed population estimates
plot_df_attr <- as.data.frame(pop_risk_post_hist_plot) %>%
  mutate(
    pop_mean_m = pop_mean / 1e6,
    pop_lwr_m  = pop_lwr  / 1e6,
    pop_upr_m  = pop_upr  / 1e6
  )

plot_df_ecdc <- as.data.frame(pop_ecdc_year) %>%
  mutate(pop_ecdc_m = pop_ecdc / 1e6)

plot_df <- plot_df_attr %>%
  left_join(plot_df_ecdc, by = "year")

y_breaks_m <- seq(
  from = floor(min(c(plot_df$pop_lwr_m, plot_df$pop_ecdc_m), na.rm = TRUE) / 50) * 50,
  to   = ceiling(max(c(plot_df$pop_upr_m, plot_df$pop_ecdc_m), na.rm = TRUE) / 50) * 50,
  by   = 50
)


p_pop_risk <- ggplot(plot_df, aes(x = year)) +

  geom_ribbon(
    aes(ymin = pop_lwr_m, ymax = pop_upr_m),
    fill  = "darkred",
    alpha = 0.20,
    show.legend = FALSE
  ) +
  geom_line(
    aes(y = pop_mean_m, color = "Climate change\nattributed"),
    linewidth = 1
  ) +

  geom_line(
    aes(y = pop_ecdc_m, color = "Absolute"),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c("Absolute" = "black", "Climate change\nattributed" = "darkred"),
    breaks = c( "Absolute", "Climate change\nattributed"),
    name   = ""
  ) +
  scale_x_continuous(
    limits = range(plot_df$year, na.rm = TRUE),
    breaks = plot_df$year,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = y_breaks_m,
    labels = function(x) paste0(x, " M"),
    expand = c(0, 0)
  ) +
  labs(
    x = "Year",
    y = "Population at risk (millions)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    panel.grid      = element_blank(),
    axis.title      = element_text(size = 9),
    axis.text       = element_text(size = 9),
    legend.position = "right",
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.text       = element_text(size = 9),

  )

p_pop_risk

# Saving our full-size historical population-at-risk plot
ext         <- "tiff"
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

width_in  <- 8.5
height_in <- 6

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_1_plots"

ggsave(
  filename    = file.path(out_dir_save, paste0("pop_at_risk_attributed_and_absolute_2011_2021.", ext)),
  plot        = p_pop_risk ,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  bg          = bg
)

ggsave(
  filename = file.path(out_dir_save, "pop_at_risk_attributed_and_absolute_2011_2021.svg"),
  plot   = p_pop_risk,
  width  = width_in,
  height = height_in,
  units  = "in",
  bg     = bg,
  device = svglite::svglite
)

# Applying journal formatting to the population-at-risk plot
p_pop_risk <- ggplot(plot_df, aes(x = year)) +
  geom_ribbon(
    aes(ymin = pop_lwr_m, ymax = pop_upr_m),
    fill  = "darkred",
    alpha = 0.20,
    show.legend = FALSE
  ) +
  geom_line(
    aes(y = pop_mean_m, color = "Climate change attributed"),
    linewidth = 0.3
  ) +
  geom_line(
    aes(y = pop_ecdc_m, color = "Absolute"),
    linewidth = 0.3
  ) +
  scale_color_manual(
    values = c("Absolute" = "black", "Climate change attributed" = "darkred"),
    breaks = c("Absolute", "Climate change attributed"),
    name   = NULL
  ) +
  scale_x_continuous(
    limits = range(plot_df$year, na.rm = TRUE),
    breaks = seq(min(plot_df$year), max(plot_df$year), by = 2),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    breaks = y_breaks_m,
    labels = function(x) paste0(x, " M"),
    expand = c(0.05, 0)
  ) +
  labs(
    x = "Year",
    y = "Population at risk (millions)"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),

    axis.title      = element_text(size = 6),
    axis.text       = element_text(size = 5, color = "black"),
    axis.line       = element_line(linewidth = 0.2),
    axis.ticks      = element_line(linewidth = 0.2, color = "black"),

    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),

    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.key.size   = unit(0.2, "cm"),
    legend.key.height = unit(0.1, "cm"),
    legend.text       = element_text(size = 5, margin = margin(l = 0.9, unit = "pt")),
    legend.margin     = margin(0, 0, 0, 0),
    legend.spacing.y  = unit(0, "pt"),

    plot.margin       = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )

out_dir_sci <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_1_plots/journal"
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 4.0
dpi       <- 600

ggplot2::ggsave(
  filename    = file.path(out_dir_sci, "pop_at_risk_journal.tiff"),
  plot        = p_pop_risk,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  dpi         = dpi,
  compression = "lzw",
  bg          = "white"
)

ggplot2::ggsave(
  filename    = file.path(out_dir_sci, "pop_at_risk_journal.svg"),
  plot        = p_pop_risk,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  bg          = "white"
)




###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Table 1 for article 
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################

table_AF_present <- map_pres_AF %>%
  sf::st_drop_geometry() %>%
  filter(y == 1) %>%
  transmute(
    Country        = CountryNam,
    ISO3_code      = CountryISO,
    NUTS3_region   = geo,
    Region_name    = LocationNa,
    `Attributable percent (AF × 100)` = AF_abs_pct_mean,
    `Lower 95% CrI`                     = AF_abs_pct_lwr,
    `Upper 95% CrI`                     = AF_abs_pct_upr
  ) %>%
  arrange(Country, Region_name)

table_AF_present

table_AF_present <- table_AF_present %>%
  mutate(
    across(
      c(`Attributable percent (AF × 100)`,
        `Lower 95% CrI`,
        `Upper 95% CrI`),
      ~ round(.x, 2)
    )
  )

table_AF_present

full_model_dir <- ".../VectorNet_ 2_original/Output/A-Final_output/Main_model_output/Attribution_table"

dir.create(full_model_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  table_AF_present,
  file = file.path(full_model_dir, "Table_AF_for_all_invaded_regions_till_2021.csv"),
  row.names = FALSE
)


###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure 2 (Panel A and C, and Panel B and D) for article
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################

orig_data = load("Data/data_INLA.RData")
nuts3_map = map_inla

nuts3_map <- st_as_sf(nuts3_map, crs=4326)
nuts3_map

# Defining our helpers for loading and identifying posterior prediction data
load_rdata_to_env <- function(path) {
  e <- new.env(parent = emptyenv())
  nms <- load(path, envir = e)
  list(env = e, names = nms)
}

pick_pred_df <- function(e) {

  objs <- ls(e)
  dfs  <- objs[vapply(objs, function(nm) inherits(e[[nm]], c("data.frame","tbl_df","tbl")), logical(1))]
  if (length(dfs) == 0) stop("No data.frame/tibble found in .RData")

  score_df <- function(df) {
    nm <- names(df)
    score <- 0L
    score <- score + ifelse("geo"  %in% nm, 10L, 0L)
    score <- score + ifelse("draw" %in% nm, 10L, 0L)
    score <- score + ifelse(any(c("year","year_idx") %in% nm), 5L, 0L)
    score <- score + ifelse(any(c("p","p_pred","prob","p_mean","p_fact","p_cf") %in% nm), 5L, 0L)
    score
  }

  scored <- tibble(obj = dfs, score = vapply(dfs, \(nm) score_df(e[[nm]]), integer(1))) %>%
    arrange(desc(score))

  df <- e[[scored$obj[1]]]

  if (!("geo" %in% names(df)) || !("draw" %in% names(df)) || !any(c("year","year_idx") %in% names(df))) {
    stop("Could not detect a prediction-like df. Candidates were: ",
         paste(scored$obj, collapse = ", "))
  }
  df
}

pick_prob_col <- function(df) {
  nm <- names(df)

  candidates <- c("p", "p_pred", "prob", "p_mean")
  col <- candidates[candidates %in% nm][1]
  if (!is.na(col)) return(col)

  stop("Could not find a probability column among: ", paste(candidates, collapse=", "),
       ". Available cols: ", paste(nm, collapse=", "))
}

# Joining our SSP and piControl posterior predictions
make_future_scenario <- function(base_dir, ssp_code) {

  file_pi  <- file.path(base_dir, sprintf(
    "future_prediction_picontrol_climate_ssp%s_socioeconomic_2022_2080_output.Rdata", ssp_code
  ))

  file_ssp <- file.path(base_dir, sprintf(
    "future_prediction_ssp%s_climate_ssp%s_socioeconomic_2022_2080_output.Rdata", ssp_code, ssp_code
  ))

  pi_env  <- load_rdata_to_env(file_pi)
  ssp_env <- load_rdata_to_env(file_ssp)

  df_pi  <- pick_pred_df(pi_env$env)
  df_ssp <- pick_pred_df(ssp_env$env)

  key_cols <- intersect(names(df_pi), names(df_ssp))

  if (!("year" %in% key_cols) && "year_idx" %in% key_cols) {

    df_pi  <- df_pi  %>% mutate(year = ifelse("year" %in% names(df_pi), year, 2021 + year_idx))
    df_ssp <- df_ssp %>% mutate(year = ifelse("year" %in% names(df_ssp), year, 2021 + year_idx))
    key_cols <- union(setdiff(key_cols, "year_idx"), "year")
  }

  join_keys <- intersect(c("geo","draw","year"), key_cols)
  if (length(join_keys) < 3) stop("Join keys missing. Need geo, draw, year in both future dfs.")

  pcol_pi  <- pick_prob_col(df_pi)
  pcol_ssp <- pick_prob_col(df_ssp)

  df_post_fact_cf <- df_ssp %>%
    transmute(geo, draw, year, p_fact = .data[[pcol_ssp]]) %>%
    left_join(
      df_pi %>% transmute(geo, draw, year, p_cf = .data[[pcol_pi]]),
      by = c("geo","draw","year")
    )

  return(df_post_fact_cf)

}


base_dir <- ".../VectorNet_ 2_original/Output/A-Final_output/Main_model_output"

future_126 <- make_future_scenario(base_dir, "126")
future_370 <- make_future_scenario(base_dir, "370")
future_585 <- make_future_scenario(base_dir, "585")

# Summarising our future factual probabilities by period and scenario

add_period <- function(df) {
  df %>%
    mutate(
      period = case_when(

        year >= 2051 & year <= 2060 ~ "2051–2060",
        year >= 2091 & year <= 2100 ~ "2091–2100",
        TRUE ~ NA_character_
      ),
      period = factor(period, levels = c("2051–2060", "2091–2100"))

    ) %>%
    filter(!is.na(period))
}

compute_spatial_pfact_by_period <- function(df_post, scenario_label = NULL) {

  df_post2 <- df_post %>%
    add_period()

  spatial_draws_period <- df_post2 %>%
    group_by(draw, geo, period) %>%
    summarise(
      T_obs  = n(),
      E_fact = mean(p_fact, na.rm = TRUE),
      .groups = "drop"
    )

  spatial_pfact_post_period <- spatial_draws_period %>%
    group_by(geo, period) %>%
    summarise(

      T_obs = first(T_obs),

      p_fact_mean = mean(E_fact, na.rm = TRUE),
      p_fact_med  = median(E_fact, na.rm = TRUE),
      p_fact_lwr  = quantile(E_fact, 0.025, na.rm = TRUE),
      p_fact_upr  = quantile(E_fact, 0.975, na.rm = TRUE),

      .groups = "drop"
    )

  if (!is.null(scenario_label)) {
    spatial_pfact_post_period <- spatial_pfact_post_period %>%
      mutate(scenario = scenario_label)
  }

  spatial_pfact_post_period
}

spatial_126_period_p <- compute_spatial_pfact_by_period(future_126, "SSP126")
spatial_370_period_p <- compute_spatial_pfact_by_period(future_370, "SSP370")
spatial_585_period_p <- compute_spatial_pfact_by_period(future_585, "SSP585")

spatial_future_period_all_p <- bind_rows(
  spatial_126_period_p,
  spatial_370_period_p,
  spatial_585_period_p
)

p_limits <- c(0, 1)
p_breaks <- seq(0, 1, by = 0.2)

pal_bwr <- "RdBu"


#Plot function: factual probability  for one SSP × one period, across all regions
plot_future_pfact_period_allregions <- function(spatial_future_period_all_p = NULL,
                                                period_label,
                                                scenario_label = "SSP126",
                                                legend_pos = c(0.93, 0.22),
                                                nuts3_map = NULL,
                                                p_limits = c(0, 1),
                                                p_breaks = seq(0, 1, 0.2)) {

  p_one <- spatial_future_period_all_p %>%
    filter(scenario == scenario_label, period == period_label) %>%
    dplyr::select(geo, p_fact_mean, p_fact_lwr, p_fact_upr)

  map_P <- nuts3_map %>%
    left_join(p_one, by = "geo")

  bb  <- sf::st_bbox(nuts3_map)
  dx  <- as.numeric(bb["xmax"] - bb["xmin"])
  dy  <- as.numeric(bb["ymax"] - bb["ymin"])
  pad <- 0.2

  xlim_pad <- c(bb["xmin"] - pad*dx, bb["xmax"] + pad*dx)
  ylim_pad <- c(bb["ymin"] - pad*dy, bb["ymax"] + pad*dy)

  ggplot() +
    geom_sf(
      data  = World_adm0,
      fill  = rgb(215/255, 215/255, 215/255),
      color = "black",
      linewidth = 0.3
    ) +
    geom_sf(
      data  = map_P,
      aes(fill = p_fact_mean),
      color = NA
    ) +
    scale_fill_distiller(

      palette   = "RdBu",
      direction = -1,
      limits    = p_limits,
      breaks    = p_breaks,
      labels    = scales::label_number(accuracy = 0.1),
      name      = "Invasion\nprobability",
      na.value  = "grey85"
    ) +
    coord_sf(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) +
    theme_void() +
    theme(
      plot.margin     = margin(8, 8, 8, 8),
      legend.position = legend_pos,
      legend.title    = element_text(size = 9, face = "bold"),
      legend.text     = element_text(size = 9)
    )
}

# SSP370 maps for  periods 2041 - 2050 
gg_P_370_2041_2050 <- plot_future_pfact_period_allregions(
  spatial_future_period_all_p = spatial_future_period_all_p,
  period_label   = "2051–2060",
  scenario_label = "SSP370",
  nuts3_map      = nuts3_map,
  p_limits       = p_limits,
  p_breaks       = p_breaks
)

gg_P_370_2041_2050

gg_P_370_2071_2080 <- plot_future_pfact_period_allregions(
  spatial_future_period_all_p = spatial_future_period_all_p,
  period_label   = "2091–2100",
  scenario_label = "SSP370",
  nuts3_map      = nuts3_map,
  p_limits       = p_limits,
  p_breaks       = p_breaks
)

gg_P_370_2071_2080


## combined maps
#scenarios <- c("SSP370")
scenarios <- c("SSP370", "SSP585")
periods   <- c("2051–2060", "2091–2100")

future_plots <- list()

for (sc in scenarios) {
  for (pr in periods) {

    pr_file <- gsub("–", "-", pr)

    key <- paste0(sc, "_", pr_file)

    future_plots[[key]] <- plot_future_pfact_period_allregions(
      spatial_future_period_all = spatial_future_period_all_p,
      period_label   = pr,
      scenario_label = sc,
      nuts3_map      = nuts3_map,
      p_limits      = p_limits,
      p_breaks      = p_breaks
    )
  }
}

ext         <- "tiff"
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

width_in  <- 8
height_in <- 7

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_2_plots"

for (nm in names(future_plots)) {

  out_file <- file.path(out_dir_save, paste0("future_invasion_probability_", nm, ".", ext))

  ggsave(
    filename    = out_file,
    plot        = future_plots[[nm]],
    width       = width_in,
    height      = height_in,
    units       = "in",
    dpi         = dpi,
    compression = compression,
    bg          = bg
  )
}


# Applying journal formatting to the future probability maps
p_limits <- c(0, 1)
p_breaks <- seq(0, 1, by = 0.2)

plot_future_pfact_period_allregions_sci <- function(spatial_future_period_all_p,
                                                    period_label,
                                                    scenario_label = "SSP126",
                                                    legend_pos = c(0.90, 0.78),
                                                    nuts3_map,
                                                    p_limits = c(0, 1),
                                                    p_breaks = seq(0, 1, 0.2),
                                                    pal_bwr = "RdBu") {

  p_one <- spatial_future_period_all_p %>%
    filter(scenario == scenario_label, period == period_label) %>%
    select(geo, p_fact_mean)

  map_P <- nuts3_map %>% left_join(p_one, by = "geo")

  bb  <- sf::st_bbox(nuts3_map)
  dx  <- as.numeric(bb["xmax"] - bb["xmin"])
  dy  <- as.numeric(bb["ymax"] - bb["ymin"])
  pad <- 0.2

  xlim_pad <- c(bb["xmin"] - pad * dx, bb["xmax"] + pad * dx)
  ylim_pad <- c(bb["ymin"] - pad * dy, bb["ymax"] + pad * dy)

  ggplot() +

    geom_sf(
      data  = World_adm0,
      fill  = rgb(215/255, 215/255, 215/255),
      color = "black",
      linewidth = 0.07
    ) +

    geom_sf(
      data  = map_P,
      aes(fill = p_fact_mean),
      color = NA,
      linewidth = 0
    ) +

    scale_fill_distiller(
      palette   = pal_bwr,
      direction = -1,
      limits    = p_limits,
      breaks    = p_breaks,
      labels    = scales::label_number(accuracy = 0.1),
      name      = "Invasion\nprobability",
      na.value  = "grey85",
      guide = guide_colorbar(

        direction     = "vertical",
        barheight     = unit(1.0, "cm"),
        barwidth      = unit(0.2, "cm"),

        ticks          = TRUE,
        ticks.colour   = "white",
        ticks.linewidth= 0.08,
        ticklength     = unit(-0.5, "mm"),

        frame.colour  = NA,

        title.position = "top",
        title.hjust    = 0.2,
        label.position = "right",
        reverse        = FALSE
      )
    ) +
    coord_sf(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) +
    theme_void() +
    theme(

      text = element_text(family = "Arial"),

      legend.title = element_text(
        size = 4.5, face = "bold",
        margin = margin(b = 2.4),
        hjust = 0.2
      ),
      legend.text  = element_text(
        size = 3.5,
        margin = margin(l = 1.2, unit = "pt")
      ),

      legend.position = legend_pos,
      legend.justification = "center",

      legend.box          = "vertical",
      legend.box.just     = "center",

      plot.margin = margin(0, 0, 0, 0)
    )
}


#scenarios <- c("SSP370")
scenarios <- c("SSP370", "SSP585")
periods   <- c("2051–2060", "2091–2100")

future_plots_sci <- list()

for (sc in scenarios) {
  for (pr in periods) {
    pr_file <- gsub("–", "-", pr)
    key <- paste0(sc, "_", pr_file)

    future_plots_sci[[key]] <- plot_future_pfact_period_allregions_sci(
      spatial_future_period_all_p = spatial_future_period_all_p,
      period_label   = pr,
      scenario_label = sc,
      nuts3_map      = nuts3_map,
      p_limits       = p_limits,
      p_breaks       = p_breaks
    )
  }
}

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_2_plots"
out_dir_sci  <- file.path(out_dir_save, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 5.0

for (nm in names(future_plots_sci)) {

  ggplot2::ggsave(
    filename    = file.path(out_dir_sci, paste0("Pfact_", nm, ".tiff")),
    plot        = future_plots_sci[[nm]],
    device      = ragg::agg_tiff,
    width       = width_cm,
    height      = height_cm,
    units       = "cm",
    dpi         = 600,
    compression = "lzw",
    bg          = "white"
  )

  ggplot2::ggsave(
    filename = file.path(out_dir_sci, paste0("Pfact_", nm, ".svg")),
    plot     = future_plots_sci[[nm]],
    device   = "svg",
    width    = width_cm,
    height   = height_cm,
    units    = "cm",
    bg       = "white"
  )
}




###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure 2 (Panel E) for article
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
load_rdata_to_env <- function(path) {
  e <- new.env(parent = emptyenv())
  nms <- load(path, envir = e)
  list(env = e, names = nms)
}


pick_pred_df <- function(e) {

  objs <- ls(e)
  dfs  <- objs[vapply(objs, function(nm) inherits(e[[nm]], c("data.frame","tbl_df","tbl")), logical(1))]
  if (length(dfs) == 0) stop("No data.frame/tibble found in .RData")

  score_df <- function(df) {
    nm <- names(df)
    score <- 0L
    score <- score + ifelse("geo"  %in% nm, 10L, 0L)
    score <- score + ifelse("draw" %in% nm, 10L, 0L)      ## this draw column will give 10 points extra to  p_post_pred_all (and it will be selected as we need all posterior samples)
    score <- score + ifelse(any(c("year","year_idx") %in% nm), 5L, 0L)
    score <- score + ifelse(any(c("p","p_pred","prob","p_mean","p_fact","p_cf") %in% nm), 5L, 0L)
    score
  }

  scored <- tibble(obj = dfs, score = vapply(dfs, \(nm) score_df(e[[nm]]), integer(1))) %>%
    arrange(desc(score))

  df <- e[[scored$obj[1]]]

  if (!("geo" %in% names(df)) || !("draw" %in% names(df)) || !any(c("year","year_idx") %in% names(df))) {
    stop("Could not detect a prediction-like df. Candidates were: ",
         paste(scored$obj, collapse = ", "))
  }
  df
}


pick_prob_col <- function(df) {

  nm <- names(df)

  candidates <- c("p", "p_pred", "prob", "p_mean")
  col <- candidates[candidates %in% nm][1]
  if (!is.na(col)) return(col)

  stop("Could not find a probability column among: ", paste(candidates, collapse=", "),
       ". Available cols: ", paste(nm, collapse=", "))
}


make_future_scenario <- function(base_dir, ssp_code) {

  file_pi  <- file.path(base_dir, sprintf(
    "future_prediction_picontrol_climate_ssp%s_socioeconomic_2022_2080_output.Rdata", ssp_code   ## note that name only has 2022_2080 (due to previous runs), but files contains predicted data from 2022 to  2100
  ))

  file_ssp <- file.path(base_dir, sprintf(
    "future_prediction_ssp%s_climate_ssp%s_socioeconomic_2022_2080_output.Rdata", ssp_code, ssp_code
  ))

  pi_env  <- load_rdata_to_env(file_pi)
  ssp_env <- load_rdata_to_env(file_ssp)

  df_pi  <- pick_pred_df(pi_env$env)
  df_ssp <- pick_pred_df(ssp_env$env)

  key_cols <- intersect(names(df_pi), names(df_ssp))

  if (!("year" %in% key_cols) && "year_idx" %in% key_cols) {

    df_pi  <- df_pi  %>% mutate(year = ifelse("year" %in% names(df_pi), year, 2021 + year_idx))
    df_ssp <- df_ssp %>% mutate(year = ifelse("year" %in% names(df_ssp), year, 2021 + year_idx))
    key_cols <- union(setdiff(key_cols, "year_idx"), "year")
  }

  join_keys <- intersect(c("geo","draw","year"), key_cols)
  if (length(join_keys) < 3) stop("Join keys missing. Need geo, draw, year in both future dfs.")

  pcol_pi  <- pick_prob_col(df_pi)
  pcol_ssp <- pick_prob_col(df_ssp)

  df_post_fact_cf <- df_ssp %>%
    transmute(geo, draw, year, p_fact = .data[[pcol_ssp]]) %>%
    left_join(
      df_pi %>% transmute(geo, draw, year, p_cf = .data[[pcol_pi]]),
      by = c("geo","draw","year")
    )

  return(df_post_fact_cf)

}


## location of .Rdata files of prediction dataset of each of pi control and ssp scenarios
base_dir <- ".../VectorNet_ 2_original/Output/A-Final_output/Main_model_output"


future_126 <- make_future_scenario(base_dir, "126")
future_370 <- make_future_scenario(base_dir, "370")
future_585 <- make_future_scenario(base_dir, "585")


# Loading our scenario-specific population projections
load_population_object <- function(rdata_path) {
  obj_names <- load(rdata_path)
  if ("population" %in% obj_names) return(get("population"))
  get(obj_names[1])
}

pop_126 <- load_population_object("Data/climate/Final_covariates_data/Future_data/Final_future_ssp126_population_data_yearly.RData")
pop_370 <- load_population_object("Data/climate/Final_covariates_data/Future_data/Final_future_ssp370_population_data_yearly.RData")
pop_585 <- load_population_object("Data/climate/Final_covariates_data/Future_data/Final_future_ssp585_population_data_yearly.RData")


# Computing our expected factual population at risk (continuous, no threshold) ----
compute_future_pop_risk_factual <- function(future_df,
                                            pop_df,
                                            scenario_label,
                                            clamp_01 = TRUE) {

  DT_fut <- as.data.table(future_df[, c("geo", "draw", "year", "p_fact")])
  DT_fut[, `:=`(
    geo  = as.character(geo),
    draw = as.integer(draw),
    year = as.integer(year)
  )]
  setkey(DT_fut, geo, year)

  DT_pop <- as.data.table(
    pop_df %>%
      transmute(
        geo        = as.character(geo),
        year       = as.integer(year),
        population = as.numeric(population)
      )
  )
  setkey(DT_pop, geo, year)

  DT <- DT_fut[DT_pop, nomatch = 0L]
  
  # double check step: # ensuring p in [0,1]
  if (clamp_01) {
    DT[, p_fact := pmin(pmax(p_fact, 0), 1)]
  }

  DT[, pop_contrib := population * p_fact]

  pop_risk_draw_year <- DT[, .(
    pop_at_risk = sum(pop_contrib, na.rm = TRUE)
  ), by = .(draw, year)]

  pop_risk_post <- pop_risk_draw_year[, .(
    pop_mean = mean(pop_at_risk, na.rm = TRUE),
    pop_lwr  = quantile(pop_at_risk, 0.025, na.rm = TRUE),
    pop_upr  = quantile(pop_at_risk, 0.975, na.rm = TRUE)
  ), by = year][order(year)]

  pop_risk_post[, scenario := scenario_label]
  pop_risk_post[]
}


pop_risk_fact_126 <- compute_future_pop_risk_factual(future_126, pop_126, "SSP126", clamp_01 = TRUE)
pop_risk_fact_370 <- compute_future_pop_risk_factual(future_370, pop_370, "SSP370", clamp_01 = TRUE)
pop_risk_fact_585 <- compute_future_pop_risk_factual(future_585, pop_585, "SSP585", clamp_01 = TRUE)

pop_risk_fact_all <- rbindlist(list(pop_risk_fact_126, pop_risk_fact_370, pop_risk_fact_585), use.names = TRUE)


# only years we want to show from 2023 (explained in mauscript)
pop_risk_fact_all_2023 <- pop_risk_fact_all %>%
  dplyr::filter(year >= 2023)


x_breaks_5y <- function(x) {
  rng <- range(x, na.rm = TRUE)
  seq(floor(rng[1] / 5) * 5, ceiling(rng[2] / 5) * 5, by = 5)
}

p_poprisk_fact <- ggplot(pop_risk_fact_all_2023, aes(x = year, y = pop_mean, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = pop_lwr, ymax = pop_upr), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = x_breaks_5y,
                     minor_breaks = NULL,
                     expand = c(0, 0)) +
  labs(x = NULL, y = "Absolute population at risk (expected)", color = "Scenario", fill = "Scenario") +
  theme_bw()

p_poprisk_fact


# Applying journal formatting to the future population plot
theme_science_small <- function(legend_pos = c(0.90, 0.38)) {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text = element_text(family = "Arial"),

      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6, color = "black"),
      axis.line  = element_line(linewidth = 0.2),
      axis.ticks = element_line(linewidth = 0.2,color = "black"),

      panel.grid = element_blank(),

      legend.title = element_text(
        size = 6, face = "bold",
        margin = margin(b = 2.5, unit = "pt")
      ),
      legend.text = element_text(size = 5),
      legend.key.size   = unit(0.15, "cm"),
      legend.key.height = unit(0.15, "cm"),
      legend.position   = legend_pos,
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
}

ssp_cols <- c(
  "SSP126" = "blue3",
  "SSP370" = "springgreen4",
  "SSP585" = "red3"
)

x_limits <- range(pop_risk_fact_all_2023$year, na.rm = TRUE)


p_pop_risk_future_sci_fact <- ggplot(

  pop_risk_fact_all_2023,
  aes(x = year)
) +

  geom_ribbon(
    aes(ymin = pop_lwr, ymax = pop_upr, fill = scenario),
    alpha = 0.12,
    linewidth = 0
  ) +

  geom_line(
    aes(y = pop_mean, color = scenario),
    linewidth = 0.4
  ) +
  scale_color_manual(values = ssp_cols, name = "Future climate\nscenario's") +
  scale_fill_manual(values = ssp_cols,  name = "Future climate\nscenario's") +

  scale_x_continuous(

    limits =  x_limits,

    breaks = seq(floor(x_limits[1]/10)*10, ceiling(x_limits[2]/10)*10, by = 10),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(

    breaks = c(180, 220, 260, 300, 340) * 1e6,
    labels = scales::label_number(scale = 1e-6, suffix = "M", accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Year",
    y = "Absolute population at risk (millions)"
  ) +
  theme_science_small(legend_pos = c(0.12, 0.90)) +
  guides(
    color = guide_legend(ncol = 1),
    fill  = guide_legend(ncol = 1)
  )

p_pop_risk_future_sci_fact


out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/supp_smoot_pop_at_risk"
out_dir_sci  <- file.path(out_dir_save, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 10
height_cm <- 7.0

ggplot2::ggsave(
  filename    = file.path(out_dir_sci, "Population_at_risk_future.tiff"),

  plot        = p_pop_risk_future_sci_fact,
  device      = ragg::agg_tiff,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

ggplot2::ggsave(
  filename = file.path(out_dir_sci, "Population_at_risk_future.svg"),

  plot     = p_pop_risk_future_sci_fact,
  device   = "svg",
  width    = width_cm,
  height   = height_cm,
  units    = "cm",
  bg       = "white"
)




###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure 3  for article - Maps and plots for attribution
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################


## first the figure 3D (attributable populate at risk)

load_rdata_to_env <- function(path) {
  e <- new.env(parent = emptyenv())
  nms <- load(path, envir = e)
  list(env = e, names = nms)
}

pick_pred_df <- function(e) {

  objs <- ls(e)
  dfs  <- objs[vapply(objs, function(nm) inherits(e[[nm]], c("data.frame","tbl_df","tbl")), logical(1))]
  if (length(dfs) == 0) stop("No data.frame/tibble found in .RData")

  score_df <- function(df) {
    nm <- names(df)
    score <- 0L
    score <- score + ifelse("geo"  %in% nm, 10L, 0L)
    score <- score + ifelse("draw" %in% nm, 10L, 0L)
    score <- score + ifelse(any(c("year","year_idx") %in% nm), 5L, 0L)
    score <- score + ifelse(any(c("p","p_pred","prob","p_mean","p_fact","p_cf") %in% nm), 5L, 0L)
    score
  }

  scored <- tibble(obj = dfs, score = vapply(dfs, \(nm) score_df(e[[nm]]), integer(1))) %>%
    arrange(desc(score))

  df <- e[[scored$obj[1]]]

  if (!("geo" %in% names(df)) || !("draw" %in% names(df)) || !any(c("year","year_idx") %in% names(df))) {
    stop("Could not detect a prediction-like df. Candidates were: ",
         paste(scored$obj, collapse = ", "))
  }
  df
}


pick_prob_col <- function(df) {

  nm <- names(df)

  candidates <- c("p", "p_pred", "prob", "p_mean")
  col <- candidates[candidates %in% nm][1]
  if (!is.na(col)) return(col)

  stop("Could not find a probability column among: ", paste(candidates, collapse=", "),
       ". Available cols: ", paste(nm, collapse=", "))
}



make_future_scenario <- function(base_dir, ssp_code) {

  file_pi  <- file.path(base_dir, sprintf(
    "future_prediction_picontrol_climate_ssp%s_socioeconomic_2022_2080_output.Rdata", ssp_code   ## note the even though names are bit different, files contains the data from 2022 to 2100
  ))

  file_ssp <- file.path(base_dir, sprintf(
    "future_prediction_ssp%s_climate_ssp%s_socioeconomic_2022_2080_output.Rdata", ssp_code, ssp_code
  ))

  pi_env  <- load_rdata_to_env(file_pi)
  ssp_env <- load_rdata_to_env(file_ssp)

  df_pi  <- pick_pred_df(pi_env$env)
  df_ssp <- pick_pred_df(ssp_env$env)

  key_cols <- intersect(names(df_pi), names(df_ssp))

  if (!("year" %in% key_cols) && "year_idx" %in% key_cols) {

    df_pi  <- df_pi  %>% mutate(year = ifelse("year" %in% names(df_pi), year, 2021 + year_idx))
    df_ssp <- df_ssp %>% mutate(year = ifelse("year" %in% names(df_ssp), year, 2021 + year_idx))
    key_cols <- union(setdiff(key_cols, "year_idx"), "year")
  }

  join_keys <- intersect(c("geo","draw","year"), key_cols)
  if (length(join_keys) < 3) stop("Join keys missing. Need geo, draw, year in both future dfs.")

  pcol_pi  <- pick_prob_col(df_pi)
  pcol_ssp <- pick_prob_col(df_ssp)

  df_post_fact_cf <- df_ssp %>%
    transmute(geo, draw, year, p_fact = .data[[pcol_ssp]]) %>%
    left_join(
      df_pi %>% transmute(geo, draw, year, p_cf = .data[[pcol_pi]]),
      by = c("geo","draw","year")
    )

  return(df_post_fact_cf)

}



base_dir <- ".../VectorNet_ 2_original/Output/A-Final_output/Main_model_output"

future_126 <- make_future_scenario(base_dir, "126")
future_370 <- make_future_scenario(base_dir, "370")
future_585 <- make_future_scenario(base_dir, "585")


load_population_object <- function(rdata_path) {
  obj_names <- load(rdata_path)
  if ("population" %in% obj_names) return(get("population"))
  get(obj_names[1])
}

# Loading our population projections for each SSP

pop_126 <- load_population_object("Data/climate/Final_covariates_data/Future_data/Final_future_ssp126_population_data_yearly.RData")
pop_370 <- load_population_object("Data/climate/Final_covariates_data/Future_data/Final_future_ssp370_population_data_yearly.RData")
pop_585 <- load_population_object("Data/climate/Final_covariates_data/Future_data/Final_future_ssp585_population_data_yearly.RData")


# Computing our continuous expected attributable population (no threshold) at risk
compute_future_pop_risk_prob <- function(future_df,
                                         pop_df,
                                         scenario_label,
                                         positive_only = TRUE,
                                         clamp_01 = TRUE) {
  
  # Keep only needed columns (critical for memory)
  DT_fut <- as.data.table(future_df[, c("geo", "draw", "year", "p_fact", "p_cf")])
  DT_fut[, `:=`(
    geo  = as.character(geo),
    draw = as.integer(draw),
    year = as.integer(year)
  )]
  setkey(DT_fut, geo, year)

  DT_pop <- as.data.table(
    pop_df %>%
      transmute(
        geo        = as.character(geo),
        year       = as.integer(year),
        population = as.numeric(population)
      )
  )
  setkey(DT_pop, geo, year)

  DT <- DT_fut[DT_pop, nomatch = 0L]
  
  # Attribution probability per row
  # Risk difference: Δp = p_fact - p_cf
  DT[, delta_p := p_fact - p_cf]
  
  # only added risk (once a region got invaded relative to picontrol, the status is positive or 0)
  if (positive_only) {
    DT[, delta_p := pmax(delta_p, 0)]
  }

  #Safe step: keep attribution prob in [0,1]
  if (clamp_01) {
    DT[, delta_p := pmin(pmax(delta_p, 0), 1)]
  }
  
  # Expected attributable population contribution
  DT[, pop_contrib := population * delta_p]

  pop_risk_draw_year <- DT[, .(
    pop_at_risk = sum(pop_contrib, na.rm = TRUE)
  ), by = .(draw, year)]

  pop_risk_post <- pop_risk_draw_year[, .(
    pop_mean = mean(pop_at_risk, na.rm = TRUE),
    pop_lwr  = quantile(pop_at_risk, 0.025, na.rm = TRUE),
    pop_upr  = quantile(pop_at_risk, 0.975, na.rm = TRUE)
  ), by = year][order(year)]

  pop_risk_post[, scenario := scenario_label]
  pop_risk_post[]
}


pop_risk_126 <- compute_future_pop_risk_prob(future_126, pop_126, "SSP126",
                                             positive_only = TRUE, clamp_01 =TRUE)
pop_risk_370 <- compute_future_pop_risk_prob(future_370, pop_370, "SSP370",
                                             positive_only = TRUE, clamp_01 = TRUE)
pop_risk_585 <- compute_future_pop_risk_prob(future_585, pop_585, "SSP585",
                                             positive_only =TRUE, clamp_01 = TRUE)

pop_risk_future_all <- rbindlist(list(pop_risk_126, pop_risk_370, pop_risk_585), use.names = TRUE, fill = TRUE)

pop_risk_future_all[, scenario := factor(scenario, levels = c("SSP126","SSP370","SSP585"))]


# Plotting our scenario-specific attributable population at risk
ssp_cols <- c(
  "SSP126" = "darkblue",
  "SSP370" = "darkgreen",
  "SSP585" = "darkred"
)

theme_clean9 <- theme_classic(base_size = 9) +
  theme(
    panel.grid        = element_blank(),
    panel.border      = element_blank(),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 9),
    axis.title        = element_text(size = 9),
    axis.text         = element_text(size = 9)
  )

p_pop_risk_future <- ggplot(pop_risk_future_all, aes(x = year)) +
  geom_ribbon(aes(ymin = pop_lwr, ymax = pop_upr, fill = scenario), alpha = 0.12) +
  geom_line(aes(y = pop_mean, color = scenario), linewidth = 1) +
  scale_color_manual(values = ssp_cols, name = "Future climate\n    scenario's") +
  scale_fill_manual(values = ssp_cols,  name = "Future climate\n    scenario's") +
  scale_x_continuous(
    limits = range(pop_risk_future_all$year, na.rm = TRUE),
    breaks = seq(min(pop_risk_future_all$year), max(pop_risk_126$year), by = 10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(40, 80,120,160,200,240) * 1e6,

    labels = scales::label_number(scale = 1e-6, suffix = "M", accuracy = 1),

    expand = c(0, 0)
  ) +
  labs(
    x = "Year",
    y = "Population at risk (in millions)"
  ) +
  theme_clean9

p_pop_risk_future


# only years we want to show from 2023 (explained in mauscript)
pop_risk_future_all_2023 <- pop_risk_future_all %>%
  dplyr::filter(year >= 2023)

x_limits <- range(pop_risk_future_all_2023$year, na.rm = TRUE)

p_pop_risk_future <- ggplot(pop_risk_future_all_2023, aes(x = year)) +
  geom_ribbon(aes(ymin = pop_lwr, ymax = pop_upr, fill = scenario), alpha = 0.12) +
  geom_line(aes(y = pop_mean, color = scenario), linewidth = 1) +
  scale_color_manual(values = ssp_cols, name = "Future climate\n    scenario's") +
  scale_fill_manual(values = ssp_cols,  name = "Future climate\n    scenario's") +
  scale_x_continuous(

    limits =  x_limits,

    breaks = seq(floor(x_limits[1]/10)*10, ceiling(x_limits[2]/10)*10, by = 10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(40, 80, 120, 160, 200, 240) * 1e6,
    labels = scales::label_number(scale = 1e-6, suffix = "M", accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(x = "Year", y = "Attributable population at risk (in millions)") +
  theme_clean9

p_pop_risk_future

ext         <- "tiff"
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

width_in  <- 8.5
height_in <- 4.8

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_3_plots"

ggsave(
  filename    = file.path(out_dir_save, paste0("pop_at_risk_future_SSP126_SSP370_SSP585_2022_2080.", ext)),
  plot        = p_pop_risk_future,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  bg          = bg
)

# Applying journal formatting to the attributable population plot
theme_science_small <- function(legend_pos = c(0.90, 0.38)) {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text = element_text(family = "Arial"),

      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6, color = "black"),
      axis.line  = element_line(linewidth = 0.2),
      axis.ticks = element_line(linewidth = 0.2,color = "black"),

      panel.grid = element_blank(),

      legend.title = element_text(
        size = 6, face = "bold",
        margin = margin(b = 2.5, unit = "pt")
      ),
      legend.text = element_text(size = 5),
      legend.key.size   = unit(0.15, "cm"),
      legend.key.height = unit(0.15, "cm"),
      legend.position   = legend_pos,
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
}

ssp_cols <- c(
  "SSP126" = "darkblue",
  "SSP370" = "darkgreen",
  "SSP585" = "darkred"
)

x_limits <- range(pop_risk_future_all_2023$year, na.rm = TRUE)

p_pop_risk_future_sci <- ggplot(

  pop_risk_future_all_2023,
  aes(x = year)
) +

  geom_ribbon(
    aes(ymin = pop_lwr, ymax = pop_upr, fill = scenario),
    alpha = 0.12,
    linewidth = 0
  ) +

  geom_line(
    aes(y = pop_mean, color = scenario),
    linewidth = 0.4
  ) +
  scale_color_manual(values = ssp_cols, name = "Future climate\nscenario's") +
  scale_fill_manual(values = ssp_cols,  name = "Future climate\nscenario's") +

  scale_x_continuous(

    limits =  x_limits,

    breaks = seq(floor(x_limits[1]/10)*10, ceiling(x_limits[2]/10)*10, by = 10),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    breaks = c(40, 80, 120, 160, 200, 240) * 1e6,
    labels = scales::label_number(scale = 1e-6, suffix = "M", accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Year",
    y = "Attributable population at risk (millions)"
  ) +
  theme_science_small(legend_pos = c(0.12, 0.90)) +
  guides(
    color = guide_legend(ncol = 1),
    fill  = guide_legend(ncol = 1)
  )

p_pop_risk_future_sci

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_3_plots"
out_dir_sci  <- file.path(out_dir_save, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 10
height_cm <- 7.0

ggplot2::ggsave(
  filename    = file.path(out_dir_sci, "Population_at_risk_future.tiff"),
  plot        = p_pop_risk_future_sci,
  device      = ragg::agg_tiff,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

ggplot2::ggsave(
  filename = file.path(out_dir_sci, "Population_at_risk_future.svg"),
  plot     = p_pop_risk_future_sci,
  device   = "svg",
  width    = width_cm,
  height   = height_cm,
  units    = "cm",
  bg       = "white"
)



# Next Figure 3A, 3B and 3C: mapping future climate attribution by period and scenario
add_period <- function(df) {
  df %>%
    mutate(
      period = case_when(
        year >= 2023 & year <= 2040 ~ "2023–2040",
        year >= 2051 & year <= 2070 ~ "2051–2070",
        year >= 2081 & year <= 2100 ~ "2081–2100",
        TRUE ~ NA_character_
      ),
      period = factor(period, levels = c("2023–2040", "2051–2070", "2081–2100"))
    ) %>%
    filter(!is.na(period))
}

# Summarising our posterior attribution within each future period
compute_spatial_by_period <- function(df_post, scenario_label = NULL) {

  df_post2 <- df_post %>%
    add_period()

  spatial_draws_period <- df_post2 %>%
    group_by(draw, geo, period) %>%
    summarise(

      T_obs = n(),

      E_fact = mean(p_fact, na.rm = TRUE),
      E_cf   = mean(p_cf,   na.rm = TRUE),

      AF = ifelse(E_fact > 0, (E_fact - E_cf), NA_real_),

      .groups = "drop"
    )

  spatial_abs_post_period <- spatial_draws_period %>%
    group_by(geo, period) %>%
    summarise(

      T_obs = first(T_obs),

      AF_abs_mean = mean(AF, na.rm = TRUE),
      AF_abs_med  = median(AF, na.rm = TRUE),
      AF_abs_lwr  = quantile(AF, 0.025, na.rm = TRUE),
      AF_abs_upr  = quantile(AF, 0.975, na.rm = TRUE),

      AF_abs_pct_mean = 100 * AF_abs_mean,
      AF_abs_pct_lwr  = 100 * AF_abs_lwr,
      AF_abs_pct_upr  = 100 * AF_abs_upr,

      .groups = "drop"
    )

  if (!is.null(scenario_label)) {
    spatial_abs_post_period <- spatial_abs_post_period %>%
      mutate(scenario = scenario_label)
  }

  spatial_abs_post_period
}

# Combining our period summaries across SSPs
spatial_126_period <- compute_spatial_by_period(future_126, "SSP126")
spatial_370_period <- compute_spatial_by_period(future_370, "SSP370")
spatial_585_period <- compute_spatial_by_period(future_585, "SSP585")


spatial_future_period_all <- bind_rows(
  spatial_126_period,
  spatial_370_period,
  spatial_585_period
)

# for common legend scale across ALL future maps (all SSPs + all periods)
af_limits <- range(spatial_future_period_all$AF_abs_pct_mean, na.rm = TRUE)

af_breaks <- pretty(spatial_future_period_all$AF_abs_pct_mean, 6)
af_breaks <- sort(unique(c(af_breaks, 0)))

# Designing common future attribution mapping function
plot_future_AF_period_allregions <- function(spatial_future_period_all = NULL,
                                             period_label,
                                             scenario_label = "SSP126",
                                             legend_pos = c(0.93, 0.22),
                                             nuts3_map = NULL,
                                             af_limits = NULL,
                                             af_breaks = NULL) {

  af_one <- spatial_future_period_all %>%
    dplyr::filter(scenario == scenario_label, period == period_label) %>%
    dplyr::select(geo, AF_abs_pct_mean, AF_abs_pct_lwr, AF_abs_pct_upr)

  map_AF <- nuts3_map %>%
    dplyr::left_join(af_one, by = "geo")

  bb  <- sf::st_bbox(nuts3_map)
  dx  <- as.numeric(bb["xmax"] - bb["xmin"])
  dy  <- as.numeric(bb["ymax"] - bb["ymin"])
  pad <- 0.2

  xlim_pad <- c(bb["xmin"] - pad*dx, bb["xmax"] + pad*dx)
  ylim_pad <- c(bb["ymin"] - pad*dy, bb["ymax"] + pad*dy)

  ggplot() +
    geom_sf(
      data  = World_adm0,
      fill  = rgb(215/255, 215/255, 215/255),
      color = "black",
      linewidth = 0.3
    ) +
    geom_sf(
      data  = map_AF,
      aes(fill = AF_abs_pct_mean),
      color = NA
    ) +
    scale_fill_distiller(
      palette   = "YlOrRd",
      direction = 1,
      limits    = af_limits,
      breaks    = af_breaks,
      labels    = scales::label_number(accuracy = 1),
      name      = "Attributable\npercent (%)",
      na.value  = "grey85"
    ) +
    coord_sf(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) +
    theme_void() +
    theme(
      plot.margin     = margin(8, 8, 8, 8),
      legend.position = legend_pos,
      legend.title    = element_text(size = 9, face = "bold"),
      legend.text     = element_text(size = 9)
    )
}

gg_AF_126_2022_2040 <- plot_future_AF_period_allregions(
  spatial_future_period_all = spatial_future_period_all,
  period_label   = "2081–2100",
  scenario_label = "SSP585",
  nuts3_map      = nuts3_map,
  af_limits      = af_limits,
  af_breaks      = af_breaks
)
gg_AF_126_2022_2040

gg_AF_126_2041_2060 <- plot_future_AF_period_allregions(
  spatial_future_period_all = spatial_future_period_all,
  period_label   = "2023–2040",
  scenario_label = "SSP585",
  nuts3_map      = nuts3_map,
  af_limits      = af_limits,
  af_breaks      = af_breaks
)
gg_AF_126_2041_2060

gg_AF_126_2061_2080 <- plot_future_AF_period_allregions(
  spatial_future_period_all = spatial_future_period_all,
  period_label   = "2081–2100",
  scenario_label = "SSP126",
  nuts3_map      = nuts3_map,
  af_limits      = af_limits,
  af_breaks      = af_breaks
)
gg_AF_126_2061_2080


# Saving the future attribution maps
scenarios <- c("SSP126", "SSP370", "SSP585")
periods   <- c("2023–2040", "2051–2070", "2081–2100")

future_plots <- list()

for (sc in scenarios) {
  for (pr in periods) {

    pr_file <- gsub("–", "-", pr)

    key <- paste0(sc, "_", pr_file)

    future_plots[[key]] <- plot_future_AF_period_allregions(
      spatial_future_period_all = spatial_future_period_all,
      period_label   = pr,
      scenario_label = sc,
      nuts3_map      = nuts3_map,
      af_limits      = af_limits,
      af_breaks      = af_breaks
    )
  }
}

ext         <- "tiff"
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

width_in  <- 8
height_in <- 7

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_4_plots"

for (nm in names(future_plots)) {

  out_file <- file.path(out_dir_save, paste0("future_attribution_", nm, ".", ext))

  ggsave(
    filename    = out_file,
    plot        = future_plots[[nm]],
    width       = width_in,
    height      = height_in,
    units       = "in",
    dpi         = dpi,
    compression = compression,
    bg          = bg
  )
}


# Applying journal formatting to the future attribution maps
af_limits <- range(spatial_future_period_all$AF_abs_pct_mean, na.rm = TRUE)

af_breaks <- pretty(spatial_future_period_all$AF_abs_pct_mean, 6)
af_breaks <- sort(unique(c(af_breaks, 0)))

plot_future_AF_period_allregions_sci <- function(spatial_future_period_all,
                                                 period_label,
                                                 scenario_label = "SSP126",
                                                 legend_pos = c(0.90, 0.78),
                                                 nuts3_map,
                                                 af_limits,
                                                 af_breaks,
                                                 pal_af = "YlOrRd") {

  af_one <- spatial_future_period_all %>%
    filter(scenario == scenario_label, period == period_label) %>%
    select(geo, AF_abs_pct_mean)

  map_AF <- nuts3_map %>% left_join(af_one, by = "geo")

  bb  <- sf::st_bbox(nuts3_map)
  dx  <- as.numeric(bb["xmax"] - bb["xmin"])
  dy  <- as.numeric(bb["ymax"] - bb["ymin"])
  pad <- 0.2

  xlim_pad <- c(bb["xmin"] - pad * dx, bb["xmax"] + pad * dx)
  ylim_pad <- c(bb["ymin"] - pad * dy, bb["ymax"] + pad * dy)

  ggplot() +

    geom_sf(
      data  = World_adm0,
      fill  = rgb(215/255, 215/255, 215/255),
      color = "black",
      linewidth = 0.07
    ) +

    geom_sf(
      data  = map_AF,
      aes(fill = AF_abs_pct_mean),
      color = NA,
      linewidth = 0
    ) +

    scale_fill_distiller(
      palette   = pal_af,
      direction = 1,
      limits    = af_limits,
      breaks    = af_breaks,
      labels    = scales::label_number(accuracy = 1),
      name      = "Attributable\npercent (%)",
      na.value  = "grey85",
      guide = guide_colorbar(
        direction      = "vertical",
        barheight      = unit(1.0, "cm"),
        barwidth       = unit(0.2, "cm"),

        ticks          = TRUE,
        ticks.colour   = "black",
        ticks.linewidth= 0.08,
        ticklength     = unit(-0.5, "mm"),

        frame.colour   = NA,

        title.position = "top",
        title.hjust    = 0.2,
        label.position = "right",
        reverse        = FALSE
      )
    ) +
    coord_sf(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) +
    theme_void() +
    theme(
      text = element_text(family = "Arial"),

      legend.title = element_text(
        size = 4.5, face = "bold",
        margin = margin(b = 2.4),
        hjust = 0.2
      ),
      legend.text  = element_text(
        size = 3.5,
        margin = margin(l = 1.2, unit = "pt")
      ),

      legend.position = legend_pos,
      legend.justification = "center",
      legend.box = "vertical",
      legend.box.just = "center",

      plot.margin = margin(0, 0, 0, 0)
    )
}

scenarios <- c("SSP126", "SSP370", "SSP585")
periods   <- c("2023–2040", "2051–2070", "2081–2100")

future_plots_AF_sci <- list()

for (sc in scenarios) {
  for (pr in periods) {
    pr_file <- gsub("–", "-", pr)
    key <- paste0(sc, "_", pr_file)

    future_plots_AF_sci[[key]] <- plot_future_AF_period_allregions_sci(
      spatial_future_period_all = spatial_future_period_all,
      period_label   = pr,
      scenario_label = sc,
      legend_pos     = c(0.90, 0.78),
      nuts3_map      = nuts3_map,
      af_limits      = af_limits,
      af_breaks      = af_breaks,
      pal_af         = "YlOrRd"
    )
  }
}

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_4_plots"
out_dir_sci  <- file.path(out_dir_save, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 5.0

for (nm in names(future_plots_AF_sci)) {

  ggplot2::ggsave(
    filename    = file.path(out_dir_sci, paste0("AF_future_", nm, ".tiff")),
    plot        = future_plots_AF_sci[[nm]],
    device      = ragg::agg_tiff,
    width       = width_cm,
    height      = height_cm,
    units       = "cm",
    dpi         = 600,
    compression = "lzw",
    bg          = "white"
  )

  ggplot2::ggsave(
    filename = file.path(out_dir_sci, paste0("AF_future_", nm, ".svg")),
    plot     = future_plots_AF_sci[[nm]],
    device   = "svg",
    width    = width_cm,
    height   = height_cm,
    units    = "cm",
    bg       = "white"
  )
}




###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure 4 (Panel A) for article 
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
future_126
future_370
future_585

# Aggregating our attribution over space and posterior draws
compute_temporal_from_fact_cf <- function(df_post_fact_cf, area_tbl, scenario_label) {

  if (!("year" %in% names(df_post_fact_cf))) {
    if ("year_idx" %in% names(df_post_fact_cf)) {
      stop("df_post has year_idx but no year. Add a year column before calling (or join a lookup).")
    } else {
      stop("df_post must contain either year or year_idx.")
    }
  }

  temporal_draws <- df_post_fact_cf %>%
    left_join(area_tbl, by = "geo") %>%
    filter(!is.na(area_m2), is.finite(area_m2), area_m2 > 0) %>%
    group_by(draw, year) %>%
    summarise(

      N_regions = n(),

      E_fact    = weighted.mean(p_fact, w = area_m2, na.rm = TRUE),
      E_cf      = weighted.mean(p_cf,   w = area_m2, na.rm = TRUE),
      AF_temp   = E_fact - E_cf,
      .groups   = "drop"
    )

  temporal_post <- temporal_draws %>%
    group_by(year) %>%
    summarise(
      N_regions = first(N_regions),
      AF_abs_mean = mean(AF_temp, na.rm = TRUE),
      AF_abs_med  = median(AF_temp, na.rm = TRUE),
      AF_abs_lwr  = quantile(AF_temp, 0.025, na.rm = TRUE),
      AF_abs_upr  = quantile(AF_temp, 0.975, na.rm = TRUE),

      AF_abs_pct_mean = 100 * AF_abs_mean,
      AF_abs_pct_lwr  = 100 * AF_abs_lwr,
      AF_abs_pct_upr  = 100 * AF_abs_upr,
      scenario = scenario_label,

      .groups = "drop"
    )
  return(temporal_post)
}


orig_data = load("Data/data_INLA.RData")
nuts3_map = map_inla

nuts3_map <- st_as_sf(nuts3_map, crs=4326)

area_tbl <- nuts3_map %>%
  st_drop_geometry() %>%
  transmute(geo, area_m2 = area)

ssp_code <- "126"
AF_SSP126 = compute_temporal_from_fact_cf(
  df_post_fact_cf = future_126,
  area_tbl = area_tbl,
  scenario_label = paste0("SSP", ssp_code)
)

ssp_code <- "370"
AF_SSP370 = compute_temporal_from_fact_cf(
  df_post_fact_cf = future_370,
  area_tbl = area_tbl,
  scenario_label = paste0("SSP", ssp_code)
)

ssp_code <- "585"
AF_SSP585 = compute_temporal_from_fact_cf(
  df_post_fact_cf = future_585,
  area_tbl = area_tbl,
  scenario_label = paste0("SSP", ssp_code)
)


AF_all <- bind_rows(AF_SSP126, AF_SSP370, AF_SSP585) %>%
  transmute(year, scenario, AF_abs_pct_mean)

# Loading and preparing our annual climate data
area_tbl <- nuts3_map %>%
  sf::st_drop_geometry() %>%
  transmute(geo, area_m2 = area)

k <- load(".../VectorNet_ 2_original/Data/climate/Final_covariates_data/Future_data/Final_ensemble_SSP1_delta_corrected_climate_data_yearly_2022_2100.Rdata")
clim_126 <- climate_summaries_region_stitched

k <- load(".../VectorNet_ 2_original/Data/climate/Final_covariates_data/Future_data/Final_ensemble_SSP3_delta_corrected_climate_data_yearly_2022_2100.Rdata")
clim_370 <- climate_summaries_region_stitched

k <- load(".../VectorNet_ 2_original/Data/climate/Final_covariates_data/Future_data/Final_ensemble_SSP5_delta_corrected_climate_data_yearly_2022_2100.Rdata")
clim_585 <- climate_summaries_region_stitched

k <- load(".../VectorNet_ 2_original/Data/climate/Final_covariates_data/Future_data/Final_ensemble_Pi_control_delta_corrected_climate_data_yearly_2022_2100.Rdata")
clim_pi <- climate_summaries_region_pi_control

# Computing our area-weighted annual temperatures
compute_area_weighted_temp <- function(clim_tbl, area_tbl) {

  clim_tbl %>%
    left_join(area_tbl, by = "geo") %>%
    filter(!is.na(area_m2), is.finite(area_m2), area_m2 > 0) %>%
    group_by(year_idx) %>%
    summarise(
      temp_aw = weighted.mean(temperature_mean, w = area_m2, na.rm = TRUE),
      .groups = "drop"
    )
}

# Calculating our SSP-minus-piControl temperature differences
compute_delta_temp <- function(clim_ssp, clim_pi, area_tbl, year_map_future, scenario_label) {

  temp_ssp <- compute_area_weighted_temp(clim_ssp, area_tbl)
  temp_pi  <- compute_area_weighted_temp(clim_pi,  area_tbl) %>%
    rename(temp_aw_pi = temp_aw)

  temp_ssp %>%
    left_join(temp_pi, by = "year_idx") %>%
    left_join(year_map_future, by = "year_idx") %>%
    transmute(
      year,
      scenario = scenario_label,
      delta_temp = temp_aw - temp_aw_pi
    )
}

# We will give a common year column we map it to calendar year (2022 is year_idx == 1)
year_map_future <- tibble(
  year_idx = sort(unique(clim_pi$year_idx)),
  year     = 2021 + year_idx
)

dT_126 <- compute_delta_temp(clim_126, clim_pi, area_tbl, year_map_future, "SSP126")
dT_370 <- compute_delta_temp(clim_370, clim_pi, area_tbl, year_map_future, "SSP370")
dT_585 <- compute_delta_temp(clim_585, clim_pi, area_tbl, year_map_future, "SSP585")

dT_all <- bind_rows(dT_126, dT_370, dT_585)

# Joining our temperature and attribution summaries

df_reg <- AF_all %>%
  left_join(dT_all, by = c("year", "scenario")) %>%
  filter(is.finite(delta_temp), is.finite(AF_abs_pct_mean))


# only years we want to show from 2023 (explained in mauscript
df_reg <- df_reg %>%
  dplyr::filter(year != 2022)


# Fitting our scenario-specific regression models
fit_tbl <- df_reg %>%
  group_by(scenario) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(AF_abs_pct_mean ~ delta_temp, data = .x)),

    grid = map(data, ~ tibble(
      delta_temp = seq(0, 8, length.out = 200)
    )),
    pred = map2(fit, grid, ~ {
      pr <- predict(.x, newdata = .y, interval = "confidence", level = 0.95)
      bind_cols(.y, as.data.frame(pr))
    })
  )

fit_tbl$fit


pred_tbl <- fit_tbl %>%
  select(scenario, pred) %>%
  unnest(pred) %>%
  rename(fit_mean = fit, lwr = lwr, upr = upr)

# Plotting our attribution response to temperature change
ssp_cols <- c(
  "SSP126" = "darkblue",
  "SSP370" = "darkgreen",
  "SSP585" = "darkred"
)

theme_clean9 <- theme_classic(base_size = 9) +
  theme(
    panel.grid        = element_blank(),
    panel.border      = element_blank(),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 9),
    axis.title        = element_text(size = 9),
    axis.text         = element_text(size = 9)
  )

p_reg <- ggplot(df_reg, aes(x = delta_temp, y = AF_abs_pct_mean,
                            color = scenario, fill = scenario)) +

  geom_ribbon(
    data = pred_tbl,
    aes(x = delta_temp, ymin = lwr, ymax = upr, fill = scenario),
    alpha = 0.12,
    color = NA,
    inherit.aes = FALSE
  ) +

  geom_line(
    data = pred_tbl,
    aes(x = delta_temp, y = fit_mean, color = scenario),
    linewidth = 1,
    inherit.aes = FALSE
  ) +

  geom_point(alpha = 0.35, size = 1) +

  scale_color_manual(values = ssp_cols, name = "Future Climate\n    scenario's") +
  scale_fill_manual(values = ssp_cols,  name = "Future Climate\n    scenario's") +

  labs(
    x = expression(Delta*" mean annual temperature ("*degree*C*")"),
    y = "Attributable percentage (AF*100)"
  ) +
  theme_clean9

p_reg

ext         <- "tiff"
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

width_in  <- 8.5
height_in <- 6

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_5_plots"

ggsave(
  filename    = file.path(out_dir_save, paste0("attribution_vs_change_in_mean_temp_future_scenario.", ext)),
  plot        = p_reg,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  bg          = bg
)


theme_science_small <- function(legend_pos = c(0.90, 0.38)) {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text       = ggplot2::element_text(family = "Arial"),
      axis.title = ggplot2::element_text(size = 7),
      axis.text  = ggplot2::element_text(size = 6, color = "black"),
      axis.line  = ggplot2::element_line(linewidth = 0.2),
      axis.ticks = ggplot2::element_line(linewidth = 0.2, color = "black"),
      panel.grid = ggplot2::element_blank(),

      legend.title = ggplot2::element_text(
        size = 6, face = "bold",
        margin = ggplot2::margin(b = 2.5, unit = "pt")
      ),
      legend.text      = ggplot2::element_text(size = 5),
      legend.key.size  = grid::unit(0.15, "cm"),
      legend.key.height= grid::unit(0.10, "cm"),
      legend.position  = legend_pos,
      legend.background    = ggplot2::element_blank(),
      legend.box.background= ggplot2::element_blank()
    )
}

x_limits <- range(c(pred_tbl$delta_temp, df_reg$delta_temp), na.rm = TRUE)

y_min <- floor(min(pred_tbl$lwr, na.rm = TRUE) / 5) * 5
y_max <- ceiling(max(pred_tbl$upr, na.rm = TRUE) / 5) * 5
y_breaks <- seq(y_min, y_max, by = 10)

lab_no_zero <- function(x) {
  lbl <- scales::label_number(accuracy = 1)(x)
  lbl[x == 0] <- ""
  lbl
}

# Applying journal formatting to the temperature-attribution plot
p_reg_sci <- ggplot2::ggplot(
  df_reg,
  ggplot2::aes(x = delta_temp, y = AF_abs_pct_mean, colour = scenario, fill = scenario)
) +
  # 95% CI band for the regression (per scenario)
  ggplot2::geom_ribbon(
    data = pred_tbl,
    ggplot2::aes(x = delta_temp, ymin = lwr, ymax = upr, fill = scenario),
    alpha = 0.12,
    linewidth = 0,
    colour = NA,
    inherit.aes = FALSE
  ) +

  ggplot2::geom_line(
    data = pred_tbl,
    ggplot2::aes(x = delta_temp, y = fit_mean, colour = scenario),
    linewidth = 0.4,
    lineend = "round",
    inherit.aes = FALSE
  ) +

  ggplot2::geom_point(alpha = 0.4, size = 0.3, shape = 16) +
  ggplot2::scale_color_manual(values = ssp_cols, name = "Future climate\nscenario's") +
  ggplot2::scale_fill_manual(values = ssp_cols,  name = "Future climate\nscenario's") +
  ggplot2::scale_x_continuous(
    limits = x_limits,
    breaks = seq(floor(x_limits[1]/2)*2,
                 ceiling(x_limits[2]/2)*2,
                 by = 2),
    labels = lab_no_zero,
    expand = expansion(mult = c(0, 0.01))) +
  ggplot2::scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = y_breaks,
    labels = scales::label_number(accuracy = 1),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = expression(Delta * "Mean annual temperature (" * degree * "C)"),
    y = "Attributable percentage (AF × 100)"
  ) +
  theme_science_small(legend_pos = c(0.2, 0.84)) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(ncol = 1),
    fill   = ggplot2::guide_legend(ncol = 1)
  )

p_reg_sci

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_5_plots"
out_dir_sci  <- file.path(out_dir_save, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 5.7

ggplot2::ggsave(
  filename    = file.path(out_dir_sci, "attribution_vs_change_in_mean_temp_future_scenario.tiff"),
  plot        = p_reg_sci,
  device      = ragg::agg_tiff,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

ggplot2::ggsave(
  filename = file.path(out_dir_sci, "attribution_vs_change_in_mean_temp_future_scenario.svg"),
  plot     = p_reg_sci,
  device   = "svg",
  width    = width_cm,
  height   = height_cm,
  units    = "cm",
  bg       = "white"
)





###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure 4 (Panel B)  for article Start
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################

pop_risk_future_all  ## loading the dependent variable of our regression model (total population at risk)

dT_all              ## contains all the aggregated (area weighted) annual temperature for all scenario (computed in Figure 5 A)


pop_risk_future_mean <- pop_risk_future_all %>%
  dplyr::select(year, scenario, pop_mean)

# Joining our annual population and temperature summaries
df_reg_pop <- pop_risk_future_mean %>%
  left_join(dT_all, by = c("year", "scenario")) %>%
  filter(is.finite(delta_temp), is.finite(pop_mean))

df_reg_pop <- df_reg_pop %>%
  dplyr::filter(year != 2022)

# Fitting our scenario-specific population regressions
fit_tbl <- df_reg_pop %>%
  group_by(scenario) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(pop_mean ~ delta_temp, data = .x)),

    grid = map(data, ~ tibble(
      delta_temp = seq(0, 8, length.out = 200)
    )),
    pred = map2(fit, grid, ~ {
      pr <- predict(.x, newdata = .y, interval = "confidence", level = 0.95)
      bind_cols(.y, as.data.frame(pr))
    })
  )

fit_tbl$fit

pred_tbl_pop <- fit_tbl %>%
  select(scenario, pred) %>%
  unnest(pred) %>%
  rename(fit_mean = fit, lwr = lwr, upr = upr)

# Plotting our population response to temperature change
ssp_cols <- c(
  "SSP126" = "darkblue",
  "SSP370" = "darkgreen",
  "SSP585" = "darkred"
)

theme_clean9 <- theme_classic(base_size = 9) +
  theme(
    panel.grid        = element_blank(),
    panel.border      = element_blank(),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 9),
    axis.title        = element_text(size = 9),
    axis.text         = element_text(size = 9)
  )

df_reg_pop_m <- df_reg_pop %>%
  dplyr::mutate(pop_mean_m = pop_mean / 1e6)

pred_tbl_pop_m <- pred_tbl_pop %>%
  dplyr::mutate(
    fit_mean_m = fit_mean / 1e6,
    lwr_m      = lwr      / 1e6,
    upr_m      = upr      / 1e6
  )

p_reg <- ggplot(df_reg_pop_m, aes(x = delta_temp, y = pop_mean_m,
                                  color = scenario, fill = scenario)) +

  geom_ribbon(
    data = pred_tbl_pop_m,
    aes(x = delta_temp, ymin = lwr_m, ymax = upr_m, fill = scenario),
    alpha = 0.12,
    color = NA,
    inherit.aes = FALSE
  ) +

  geom_line(
    data = pred_tbl_pop_m,
    aes(x = delta_temp, y = fit_mean_m, color = scenario),
    linewidth = 1,
    inherit.aes = FALSE
  ) +

  geom_point(alpha = 0.35, size = 1) +

  scale_color_manual(values = ssp_cols, name = "Future Climate\n    scenario's") +
  scale_fill_manual(values = ssp_cols,  name = "Future Climate\n    scenario's") +

  scale_y_continuous(
    labels = scales::label_number(suffix = " M", accuracy = 1),
    expand = c(0, 0)
  ) +

  labs(
    x = expression(Delta*" Annual temperature ("*degree*C*")"),
    y = "Attributable population at risk (in Millions)"
  ) +
  theme_clean9

p_reg

ext         <- "tiff"
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

width_in  <- 8.5
height_in <- 6

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_5_plots"

ggsave(
  filename    = file.path(out_dir_save, paste0("pop_attribution_vs_change_in_mean_temp_future_scenario.", ext)),
  plot        = p_reg,
  width       = width_in,
  height      = height_in,
  units       = "in",
  dpi         = dpi,
  compression = compression,
  bg          = bg
)


theme_science_small <- function(legend_pos = c(0.90, 0.38)) {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text       = ggplot2::element_text(family = "Arial"),
      axis.title = ggplot2::element_text(size = 6),
      axis.text  = ggplot2::element_text(size = 6, color = "black"),
      axis.line  = ggplot2::element_line(linewidth = 0.2),
      axis.ticks = ggplot2::element_line(linewidth = 0.2, color = "black"),
      panel.grid = ggplot2::element_blank(),

      legend.title = ggplot2::element_text(
        size = 6, face = "bold",
        margin = ggplot2::margin(b = 2.5, unit = "pt")
      ),
      legend.text      = ggplot2::element_text(size = 5),
      legend.key.size  = grid::unit(0.15, "cm"),
      legend.key.height= grid::unit(0.10, "cm"),
      legend.position  = legend_pos,
      legend.background    = ggplot2::element_blank(),
      legend.box.background= ggplot2::element_blank()
    )
}

x_limits <- range(c(pred_tbl_pop_m$delta_temp, df_reg_pop_m$delta_temp), na.rm = TRUE)

lab_no_zero <- function(x) {
  lbl <- scales::label_number(accuracy = 1)(x)
  lbl[x == 0] <- ""
  lbl
}

# Applying journal formatting to the population regression plot
p_reg_sci_pop <- ggplot2::ggplot(
  df_reg_pop_m,
  ggplot2::aes(x = delta_temp, y = pop_mean_m, colour = scenario, fill = scenario)
) +

  ggplot2::geom_ribbon(
    data = pred_tbl_pop_m,
    ggplot2::aes(x = delta_temp, ymin = lwr_m, ymax = upr_m, fill = scenario),
    alpha = 0.12,
    linewidth = 0,
    colour = NA,
    inherit.aes = FALSE
  ) +

  ggplot2::geom_line(
    data =  pred_tbl_pop_m,
    ggplot2::aes(x = delta_temp, y = fit_mean_m, colour = scenario),
    linewidth = 0.4,
    lineend = "round",
    inherit.aes = FALSE
  ) +

  ggplot2::geom_point(alpha = 0.4, size = 0.3, shape = 16) +
  ggplot2::scale_color_manual(values = ssp_cols, name = "Future climate\nscenario's") +
  ggplot2::scale_fill_manual(values = ssp_cols,  name = "Future climate\nscenario's") +
  ggplot2::scale_x_continuous(
    limits = x_limits,
    breaks = seq(floor(x_limits[1]/2)*2,
                 ceiling(x_limits[2]/2)*2,
                 by = 2),
    labels = lab_no_zero,
    expand = expansion(mult = c(0, 0.01))) +
  ggplot2::scale_y_continuous(

    breaks = c(40, 80, 120, 160, 200, 240, 280),
    labels = scales::label_number(suffix = "M", accuracy = 1),
    expand = c(0, 0)

  ) +
  ggplot2::labs(
    x = expression(Delta * "Mean annual temperature (" * degree * "C)"),
    y = "Attributable population at risk (in Millions)"
  ) +
  theme_science_small(legend_pos = c(0.2, 0.84)) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(ncol = 1),
    fill   = ggplot2::guide_legend(ncol = 1)
  )

p_reg_sci_pop

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_5_plots"
out_dir_sci  <- file.path(out_dir_save, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 5.7

ggplot2::ggsave(
  filename    = file.path(out_dir_sci, "pop_at_risk_vs_change_in_mean_temp_future_scenario.tiff"),
  plot        = p_reg_sci_pop,
  device      = ragg::agg_tiff,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

ggplot2::ggsave(
  filename = file.path(out_dir_sci, "pop_at_risk_vs_change_in_mean_temp_future_scenario.svg"),
  plot     = p_reg_sci_pop,
  device   = "svg",
  width    = width_cm,
  height   = height_cm,
  units    = "cm",
  bg       = "white"
)





###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure 5  for article
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################

# assuming these objects below are already in memory exactly --
# nuts3_map
# World_adm0
# future_126, future_370, future_585
# clim_126, clim_370, clim_585, clim_pi
# pop_126, pop_370, pop_585

# Defining our common future periods
make_period <- function(year) {
  dplyr::case_when(
    year >= 2023 & year <= 2040 ~ "2023–2040",
    year >= 2051 & year <= 2070 ~ "2051–2070",
    year >= 2081 & year <= 2100 ~ "2081–2100",
    TRUE ~ NA_character_
  )
}

# Standardising our climate tables
prep_clim_tbl <- function(clim_tbl, scenario_name) {
  clim_tbl %>%
    transmute(
      geo      = geo,
      year     = as.integer(year),
      temp_summer = temperature_summer_mean,
      scenario = scenario_name
    )
}

clim_126_s <- prep_clim_tbl(clim_126, "SSP126")
clim_370_s <- prep_clim_tbl(clim_370, "SSP370")
clim_585_s <- prep_clim_tbl(clim_585, "SSP585")
clim_pi_s  <- prep_clim_tbl(clim_pi,  "piControl")

# Standardising the population tables
prep_pop_tbl <- function(pop_tbl, scenario_name) {
  pop_tbl %>%
    transmute(
      geo        = geo,
      year       = as.integer(year),
      population = population,
      scenario   = scenario_name
    )
}

pop_126_s <- prep_pop_tbl(pop_126, "SSP126")
pop_370_s <- prep_pop_tbl(pop_370, "SSP370")
pop_585_s <- prep_pop_tbl(pop_585, "SSP585")

HR_T <- 1.55     ## published estimate by Zia et al. hazard ratio per +1°C summer warming

# Calculating our factual-to-counterfactual outbreak hazard ratios
make_hr_future_tbl <- function(clim_fact_tbl, clim_cf_tbl, scenario_name) {
  
  clim_fact_tbl %>%
    select(geo, year, temp_summer) %>%
    rename(temp_fact = temp_summer) %>%
    left_join(
      clim_cf_tbl %>%
        select(geo, year, temp_summer) %>%
        rename(temp_cf = temp_summer),
      by = c("geo", "year")
    ) %>%
    filter(year >= 2023) %>%
    mutate(
      scenario = scenario_name,
      hr_f_cf  = HR_T^(temp_fact - temp_cf),
      log_hr   = log(hr_f_cf),
      period   = make_period(year)
    )
}

hr_126 <- make_hr_future_tbl(clim_126_s, clim_pi_s, "SSP126")
hr_370 <- make_hr_future_tbl(clim_370_s, clim_pi_s, "SSP370")
hr_585 <- make_hr_future_tbl(clim_585_s, clim_pi_s, "SSP585")

hr_future_all <- bind_rows(hr_126, hr_370, hr_585)


# Extracting one area value for each region
area_tbl <- pop_126 %>%
  group_by(geo) %>%
  summarise(
    area_km_2 = first(area_km_2),
    .groups = "drop"
  )

make_decade <- function(year) {
  dplyr::case_when(
    year >= 2031 & year <= 2040 ~ "2031–2040",
    year >= 2041 & year <= 2050 ~ "2041–2050",
    year >= 2051 & year <= 2060 ~ "2051–2060",
    year >= 2061 & year <= 2070 ~ "2061–2070",
    year >= 2071 & year <= 2080 ~ "2071–2080",
    year >= 2081 & year <= 2090 ~ "2081–2090",
    year >= 2091 & year <= 2100 ~ "2091–2100",
    TRUE ~ NA_character_
  )
}


# Aggregating our area-weighted hazard estimates by decade
hr_future_decade <- hr_future_all %>%
  mutate(decade = make_decade(year)) %>%
  left_join(area_tbl, by = "geo") %>%
  filter(
    !is.na(decade),
    !is.na(area_km_2),
    area_km_2 > 0
  ) %>%
  group_by(scenario, year, decade) %>%
  summarise(
    HR_aw_year = exp(weighted.mean(log_hr, w = area_km_2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  group_by(scenario, decade) %>%
  summarise(
    HR_aw_mean = mean(HR_aw_year, na.rm = TRUE),
    
    HR_aw_lwr  = quantile(HR_aw_year, 0.025, na.rm = TRUE),
    HR_aw_upr  = quantile(HR_aw_year, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    decade = factor(
      decade,
      levels = c(
        "2031–2040",
        "2041–2050",
        "2051–2060",
        "2061–2070",
        "2071–2080",
        "2081–2090",
        "2091–2100"
      )
    )
  )

hr_future_decade



# Designing  theme
theme_science_small <- function(legend_pos = c(0.90, 0.38)) {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text = element_text(family = "Arial"),
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6, color = "black"),
      axis.line  = element_line(linewidth = 0.2),
      axis.ticks = element_line(linewidth = 0.2, color = "black"),
      panel.grid = element_blank(),
      legend.title = element_text(
        size = 6, face = "bold",
        margin = margin(b = 2.5, unit = "pt")
      ),
      legend.text = element_text(size = 5),
      legend.key.size   = unit(0.15, "cm"),
      legend.key.height = unit(0.15, "cm"),
      legend.position   = legend_pos,
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
}

ssp_cols <- c(
  "SSP126" = "darkblue",
  "SSP370" = "darkgreen",
  "SSP585" = "darkred"
)

# Offset scenario-specific points and error bars within each decade
pd <- position_dodge(width = 0.45)

y_min <- floor(min(hr_future_decade$HR_aw_lwr, na.rm = TRUE) * 10) / 10
y_max <- ceiling(max(hr_future_decade$HR_aw_upr, na.rm = TRUE) * 10) / 10
y_breaks <- pretty(c(y_min, y_max), n = 5)

decade_labels <- c(
  "2031–2040" = "2031–2040",
  "2041–2050" = "2041–2050",
  "2051–2060" = "2051–2060",
  "2061–2070" = "2061–2070",
  "2071–2080" = "2071–2080",
  "2081–2090" = "2081–2090",
  "2091–2100" = "2091–2100"
)

# Plotting the decade-level outbreak hazard estimates
p_hr_future_decade_sci <- ggplot(
  hr_future_decade,
  aes(x = decade, y = HR_aw_mean, color = scenario)
) +
  
  geom_errorbar(
    aes(ymin = HR_aw_lwr, ymax = HR_aw_upr),
    width = 0.12,
    linewidth = 0.4,
    position = pd
  ) +
  
  geom_point(
    position = pd,
    size = 1.2,
    stroke = 0.2
  ) +
  
  scale_color_manual(values = ssp_cols, name = "Future climate\nscenarios") +
  
  scale_x_discrete(
    labels = decade_labels,
    expand = expansion(add = c(0.4, 0.0))
  ) +
  
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = y_breaks,
    labels = scales::label_number(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Decade",
    y = "Outbreak hazard estimate"
  ) +
  theme_science_small(legend_pos = c(0.015, 0.98)) +
  theme(
    legend.justification = c(0, 1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  guides(
    color = guide_legend(ncol = 1)
  )

p_hr_future_decade_sci

# Saving our outbreak hazard estimates
out_dir_hr <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_6_plots/journal"
dir.create(out_dir_hr, recursive = TRUE, showWarnings = FALSE)

ggsave(
  
  filename    = file.path(out_dir_hr, "Outbreak_HR_decade_points.tiff"),
  
  plot        = p_hr_future_decade_sci,
  device      = ragg::agg_tiff,
  width       = 10,
  height      = 7,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

ggsave(
  
  filename = file.path(out_dir_hr, "Outbreak_HR_decade_points.svg"),
  
  plot        = p_hr_future_decade_sci,
  device   = "svg",
  width       = 10,
  height   =   7,
  dpi         = 600,
  units    = "cm",
  bg       = "white"
)





###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure 6 (Panel A, B, C) for article
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################

## assuming these objects are already in memory
hr_future_all
future_126
clim_126_s
clim_pi_s
pop_126_s
pop_126


# Defining the common future periods for spatial summaries
make_period <- function(year) {
  dplyr::case_when(
    year >= 2023 & year <= 2040 ~ "2023–2040",
    year >= 2051 & year <= 2070 ~ "2051–2070",
    year >= 2081 & year <= 2100 ~ "2081–2100",
    TRUE ~ NA_character_
  )
}

# Summarising our spatial outbreak opportunity index
compute_deltaR_map <- function(df_post, clim_fact_tbl, clim_cf_tbl, scenario_name) {
  
  inp <- make_deltaR_inputs(clim_fact_tbl, clim_cf_tbl, scenario_name)
  
  df_post %>%
    transmute(
      geo   = geo,
      draw  = draw,
      year  = as.integer(year),
      p_fact = p_fact,
      p_cf   = p_cf
    ) %>%
    filter(year >= 2023) %>%
    
    left_join(inp$hr_tbl %>% select(geo, year, hr_f_cf), by = c("geo", "year")) %>%
    
    mutate(
      deltaR = p_fact * hr_f_cf - p_cf,
      period = make_period(year)
    ) %>%
    filter(!is.na(period)) %>%
    
    group_by(geo, year, period) %>%
    summarise(
      deltaR_mean = mean(deltaR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    
    group_by(geo, period) %>%
    summarise(
      deltaR_period_mean = mean(deltaR_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      scenario = scenario_name,
      deltaR_period_mean_pct = 100 * deltaR_period_mean
    )
}

spatial_deltaR_126 <- compute_deltaR_map(future_126, clim_126_s, clim_pi_s, "SSP126")
spatial_deltaR_370 <- compute_deltaR_map(future_370, clim_370_s, clim_pi_s, "SSP370")
spatial_deltaR_585 <- compute_deltaR_map(future_585, clim_585_s, clim_pi_s, "SSP585")

spatial_deltaR_period_all <- bind_rows(
  spatial_deltaR_126,
  spatial_deltaR_370,
  spatial_deltaR_585
)

spatial_deltaR_period_all

# Designing the spatial outbreak opportunity maps
plot_future_deltaR_period_allregions_sci <- function(
    spatial_deltaR_period_all,
    period_label,
    scenario_label = "SSP126",
    legend_pos = c(0.90, 0.78),
    nuts3_map,
    delta_limits,
    delta_breaks,
    pal_delta = "YlOrRd"
) {
  
  delta_one <- spatial_deltaR_period_all %>%
    filter(scenario == scenario_label, period == period_label) %>%
    select(geo, deltaR_period_mean)
  
  map_deltaR <- nuts3_map %>% left_join(delta_one, by = "geo")
  
  bb  <- sf::st_bbox(nuts3_map)
  dx  <- as.numeric(bb["xmax"] - bb["xmin"])
  dy  <- as.numeric(bb["ymax"] - bb["ymin"])
  pad <- 0.2
  
  xlim_pad <- c(bb["xmin"] - pad * dx, bb["xmax"] + pad * dx)
  ylim_pad <- c(bb["ymin"] - pad * dy, bb["ymax"] + pad * dy)
  
  ggplot() +
    
    geom_sf(
      data = World_adm0,
      fill = rgb(215/255, 215/255, 215/255),
      color = "black",
      linewidth = 0.07
    ) +
    
    geom_sf(
      data  = map_deltaR,
      aes(fill = deltaR_period_mean),
      color = NA,
      linewidth = 0
    ) +
    
    scale_fill_distiller(
      palette   = pal_delta,
      direction = 1,
      limits    = delta_limits,
      breaks    = delta_breaks,
      labels    = scales::label_number(accuracy = 1),
      name      = "Outbreak\nRisk Index",
      na.value  = "grey85",
      guide = guide_colorbar(
        direction      = "vertical",
        barheight      = unit(1.0, "cm"),
        barwidth       = unit(0.2, "cm"),
        
        ticks          = TRUE,
        ticks.colour   = "black",
        ticks.linewidth= 0.08,
        ticklength     = unit(-0.5, "mm"),
        frame.colour   = NA,
        title.position = "top",
        title.hjust    = 0.2,
        label.position = "right",
        reverse        = FALSE
      )
    ) +
    coord_sf(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) +
    theme_void() +
    theme(
      text = element_text(family = "Arial"),
      
      legend.title = element_text(
        size = 4.5, face = "bold",
        margin = margin(b = 2.4),
        hjust = 0.2
      ),
      legend.text = element_text(
        size = 3.5,
        margin = margin(l = 1.2, unit = "pt")
      ),
      legend.position = legend_pos,
      legend.justification = "center",
      legend.box = "vertical",
      legend.box.just = "center",
      plot.margin = margin(0, 0, 0, 0)
    )
}

delta_limits <- range(spatial_deltaR_period_all$deltaR_period_mean, na.rm = TRUE)
delta_breaks <- pretty(spatial_deltaR_period_all$deltaR_period_mean, n = 6)

scenarios <- c("SSP126", "SSP370", "SSP585")
periods   <- c("2023–2040", "2051–2070", "2081–2100")

# Building our outbreak opportunity maps with a common scale
future_plots_deltaR_sci <- list()

for (sc in scenarios) {
  for (pr in periods) {
    pr_file <- gsub("–", "-", pr)
    key <- paste0(sc, "_", pr_file)
    
    future_plots_deltaR_sci[[key]] <- plot_future_deltaR_period_allregions_sci(
      spatial_deltaR_period_all = spatial_deltaR_period_all,
      period_label   = pr,
      scenario_label = sc,
      legend_pos     = c(0.90, 0.78),
      nuts3_map      = nuts3_map,
      delta_limits   = delta_limits,
      delta_breaks   = delta_breaks,
      pal_delta = "RdPu"
    )
  }
}

future_plots_deltaR_sci[["SSP126_2081-2100"]]

future_plots_deltaR_sci[["SSP126_2051-2070"]]

# Saving our common-scale outbreak opportunity maps
out_dir_deltaR_map_common <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figure_7_plots/journal"
dir.create(out_dir_deltaR_map_common, recursive = TRUE, showWarnings = FALSE)

for (nm in names(future_plots_deltaR_sci)) {
  ggsave(
    filename    = file.path(out_dir_deltaR_map_common, paste0("Outbreak_deltaR_map_common_legend", nm, ".tiff")),
    plot        = future_plots_deltaR_sci[[nm]],
    device      = ragg::agg_tiff,
    width       = 5.7,
    height      = 5.0,
    units       = "cm",
    dpi         = 600,
    compression = "lzw",
    bg          = "white"
  )
  
  ggsave(
    filename = file.path(out_dir_deltaR_map_common, paste0("Outbreak_deltaR_map_common_legend", nm, ".svg")),
    plot     = future_plots_deltaR_sci[[nm]],
    device   = "svg",
    width    = 5.7,
    height   = 5.0,
    units    = "cm",
    bg       = "white"
  )
}

# Generating our outbreak opportunity maps with period-specific scales
make_deltaR_maps_for_one_period<- function(
    period_label,
    spatial_deltaR_period_all,
    scenarios = c("SSP126", "SSP370", "SSP585"),
    nuts3_map,
    legend_pos = c(0.90, 0.78),
    pal_delta = "YlOrBr"
) {
  
  delta_period_values <- spatial_deltaR_period_all %>%
    filter(period == period_label, scenario %in% scenarios) %>%
    pull(deltaR_period_mean)
  
  delta_limits <- range(delta_period_values, na.rm = TRUE)
  
  if (period_label == "2023–2040") {
    
    delta_breaks <- seq(
      from = ceiling(delta_limits[1]),
      to   = floor(delta_limits[2]),
      by   = 1
    )
    
    if (length(delta_breaks) < 2) {
      delta_breaks <- seq(
        from = floor(delta_limits[1]),
        to   = ceiling(delta_limits[2]),
        by   = 1
      )
    }
    
    delta_breaks <- sort(unique(delta_breaks))
    
  } else {
    
    delta_breaks <- pretty(delta_period_values, n = 6)
    delta_breaks <- sort(unique(delta_breaks))
    
  }
  
  plots <- list()
  
  for (sc in scenarios) {
    plots[[sc]] <- plot_future_deltaR_period_allregions_sci(
      spatial_deltaR_period_all = spatial_deltaR_period_all,
      period_label   = period_label,
      scenario_label = sc,
      legend_pos     = legend_pos,
      nuts3_map      = nuts3_map,
      delta_limits   = delta_limits,
      delta_breaks   = delta_breaks,
      pal_delta         = pal_delta
    )
  }
  
  return(plots)
}

# Building the period-specific outbreak opportunity maps
future_plots_deltaR_period_specific <- list()

for (pr in periods) {
  
  pr_file <- gsub("–", "-", pr)
  
  period_maps <- make_deltaR_maps_for_one_period(
    period_label = pr,
    spatial_deltaR_period_all = spatial_deltaR_period_all,
    scenarios = scenarios,
    nuts3_map = nuts3_map,
    legend_pos = c(0.90, 0.78),
    pal_delta = "RdPu"
  )
  
  for (sc in scenarios) {
    key <- paste0(sc, "_", pr_file)
    future_plots_deltaR_period_specific[[key]] <- period_maps[[sc]]
  }
}

future_plots_deltaR_period_specific[["SSP126_2023-2040"]]
future_plots_deltaR_period_specific[["SSP370_2023-2040"]]
future_plots_deltaR_period_specific[["SSP585_2023-2040"]]

# Saving our period-specific outbreak opportunity maps
out_dir_deltaR_map <- ".../VectorNet_ 2_original/Plots/A_final_plots/Main_article_plot/Figures_extended_1_plots/journal"
dir.create(out_dir_deltaR_map, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 5.0

for (nm in names(future_plots_HR_period_specific)) {
  
  ggplot2::ggsave(
    filename    = file.path(out_dir_deltaR_map, paste0("Outbreak_deltaR_map_", nm, ".tiff")),
    plot        = future_plots_deltaR_period_specific[[nm]],
    device      = ragg::agg_tiff,
    width       = width_cm,
    height      = height_cm,
    units       = "cm",
    dpi         = 600,
    compression = "lzw",
    bg          = "white"
  )
  
  ggplot2::ggsave(
    filename = file.path(out_dir_deltaR_map, paste0("Outbreak_deltaR_map_", nm, ".svg")),
    plot     = future_plots_deltaR_period_specific[[nm]],
    device   = "svg",
    width    = width_cm,
    height   = height_cm,
    units    = "cm",
    bg       = "white"
  )
}





###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure S1 (A)  for supplementary section 
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
orig_data = load("Data/data_INLA.RData")
invasion_data  = tb_inla
nuts3_map = map_inla

World_adm0 <- st_read('good_visual_plots_code/map_and_data_for_plots/WB_Land_10m')

latest_year <- max(invasion_data$date_start)
latest_year

invasion_year_region <- invasion_data %>%
  group_by(geo) %>%
  summarize(

    latest_year = max(date_start, na.rm = TRUE),

    presence_valid = all(presence == 0 | lag(presence) == 1, na.rm = TRUE),

    invasion_year = if (all(is.na(presence))) {
      NA
    } else if (!latest_year %in% date_start[presence == 1]) {
      2024
    } else if (!presence_valid) {
      max(date_start[presence == 1 & lag(presence) == 0], na.rm = TRUE)
    } else {
      min(date_start[presence == 1], na.rm = TRUE)
    }
  ) %>%
  ungroup()

examp_year <- invasion_year_region %>% filter(geo == "HU110")
print(examp_year)

invasion_year_region <- invasion_year_region %>%
  dplyr::select(-presence_valid)

print(invasion_year_region)

# Joining our invasion years to the NUTS3 geometries
nuts3_map$YearObserved <- invasion_year_region$invasion_year[match(nuts3_map$geo, invasion_year_region$geo)]

nuts3_map <- st_as_sf(nuts3_map)

palette_info <- brewer.pal.info["BrBG", ]
palette_info


# Designing our final invasion-year map
nuts3_map_till_2023 <- nuts3_map %>%
  dplyr::filter(is.na(YearObserved) | YearObserved <= 2023)

nuts3_map_cont <- nuts3_map_till_2023 %>%
  dplyr::filter(!is.na(YearObserved))

nuts3_map_na <- nuts3_map %>%
  dplyr::filter(is.na(YearObserved))   # "No data" layer

nuts3_map_2024 <- nuts3_map %>%
  dplyr::filter(YearObserved == 2024) # "Absent" layer


bb  <- sf::st_bbox(nuts3_map)
dx  <- as.numeric(bb["xmax"] - bb["xmin"])
dy  <- as.numeric(bb["ymax"] - bb["ymin"])
pad <- 0.2

xlim_pad <- c(bb["xmin"] - pad*dx, bb["xmax"] + pad*dx)
ylim_pad <- c(bb["ymin"] - pad*dy, bb["ymax"] + pad*dy)

breaks_year <- c(2010, 2012, 2014, 2016, 2018, 2020, 2022)
year_labels <- function(x) ifelse(x == 2010, "\u2264 2010", as.character(x))

legend_pos <- c(0.94, 0.25)

c <- ggplot() +
  geom_sf(
    data  = World_adm0,
    fill  = rgb(215/255, 215/255, 215/255),
    color = "black",
    linewidth = 0.3
  ) +

  geom_sf(
    data  = nuts3_map_cont,
    aes(fill = YearObserved),
    color = NA
  ) +
  scale_fill_distiller(
    palette   = "BrBG",
    direction = 1,
    limits    = c(2010, 2023),
    breaks    = breaks_year,
    labels    = year_labels,
    name      = "Invasion\nyear",
    guide     = guide_colorbar(
      order = 1,
      barheight = unit(35, "mm"),
      ticks = TRUE
    )
  ) +

  ggnewscale::new_scale_fill() +

  geom_sf(
    data = nuts3_map_na,
    aes(fill = "No data"),
    color = NA,
    show.legend = TRUE
  ) +
  geom_sf(
    data = nuts3_map_2024,
    aes(fill = "Absent"),
    color = NA,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = c("No data" = rgb(150/255, 150/255, 150/255),
               "Absent"  = "darkblue"),
    breaks = c("No data", "Absent"),
    name   = NULL,
    guide  = guide_legend(
      order = 2,
      override.aes = list(color = NA)
    )
  ) +

  coord_sf(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) +
  theme_void() +
  theme(
    plot.margin     = margin(8, 8, 8, 8),
    legend.position = legend_pos,
    legend.box      = "vertical",
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 9),
    legend.background = element_rect(fill = NA, color = NA)
  )

ext         <- "tiff"
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

width_in  <- 8
height_in <- 7

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Inkscape_final_arranged_plots/supp_tiff_600_dpi/Plots_S1"
dir.create(out_dir_save, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(out_dir_save, paste0("map_obs_EU_invasion_year.", ext))

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    filename    = outfile,
    plot        = p_map_obs_EU,
    device      = ragg::agg_tiff,
    width       = width_in,
    height      = height_in,
    units       = "in",
    dpi         = dpi,
    compression = compression,
    background  = bg
  )
} else {
  ggsave(
    filename    = outfile,
    plot        = p_map_obs_EU,
    device      = "tiff",
    width       = width_in,
    height      = height_in,
    units       = "in",
    dpi         = dpi,
    compression = compression,
    bg          = bg,
    type        = "cairo"
  )
}


# Applying journal formatting to the invasion-year map
bb  <- sf::st_bbox(nuts3_map)
dx  <- as.numeric(bb["xmax"] - bb["xmin"])
dy  <- as.numeric(bb["ymax"] - bb["ymin"])
pad <- 0.2
xlim_pad <- c(bb["xmin"] - pad * dx, bb["xmax"] + pad * dx)
ylim_pad <- c(bb["ymin"] - pad * dy, bb["ymax"] + pad * dy)

breaks_year  <- c(2010, 2012, 2014, 2016, 2018, 2020, 2022)
year_labels  <- function(x) ifelse(x == 2010, "\u2264 2010", as.character(x))

gg_invasion_year_sci <- ggplot() +
  geom_sf(
    data      = World_adm0,
    fill      = rgb(215/255, 215/255, 215/255),
    color     = "black",
    linewidth = 0.07
  ) +

  geom_sf(
    data      = nuts3_map_cont,
    aes(fill  = YearObserved),
    color     = NA,
    linewidth = 0
  ) +
  scale_fill_distiller(
    palette   = "BrBG",
    direction = 1,
    limits    = c(2010, 2023),
    breaks    = breaks_year,
    labels    = year_labels,
    name      = "Invasion\nyear",
    guide     = guide_colorbar(
      order          = 1,
      title.position = "top",
      title.hjust    = 0.2,
      barheight      = unit(0.8, "cm"),
      barwidth       = unit(0.15, "cm"),
      ticks          = TRUE
    )
  ) +

  ggnewscale::new_scale_fill() +

  geom_sf(
    data        = nuts3_map_na,
    aes(fill    = "No data"),
    color       = NA,
    linewidth   = 0,
    show.legend = TRUE
  ) +
  geom_sf(
    data        = nuts3_map_2024,
    aes(fill    = "Absent"),
    color       = NA,
    linewidth   = 0,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "Absent"  = "darkblue",
      "No data" = rgb(150/255, 150/255, 150/255)
    ),
    breaks = c("Absent", "No data"),
    name   = "Invasion\nstatus",
    drop   = FALSE,
    guide  = guide_legend(
      order          = 2,
      title.position = "top",
      title.hjust    = 0.2,
      ncol           = 1,
      override.aes   = list(color = NA)
    )
  ) +

  coord_sf(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) +
  theme_void() +
  theme(
    text = element_text(family = "Arial"),

    legend.title = element_text(
      size   = 4.5,
      face   = "bold",
      margin = margin(b = 1.5),
      hjust  = 0.2
    ),
    legend.text = element_text(
      size   = 3.0,
      margin = margin(l = 1.2, unit = "pt")
    ),

    legend.position       = c(0.90, 0.78),
    legend.justification  = "center",
    legend.box            = "vertical",
    legend.box.just       = "center",

    legend.key.size       = unit(0.15, "cm"),
    legend.key.width      = unit(0.15, "cm"),
    legend.spacing.y      = unit(0.25, "cm"),
    legend.box.spacing    = unit(0, "cm"),
    legend.margin         = margin(0, 0, 0, 0),

    plot.margin           = margin(0, 0, 0, 0),

    legend.background     = element_rect(fill = NA, color = NA)
  )

out_dir_base <- ".../VectorNet_ 2_original/Plots/A_final_plots/Inkscape_final_arranged_plots/supp_tiff_600_dpi/Plots_S1"
out_dir_sci  <- file.path(out_dir_base, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 5.0

ggsave(
  filename    = file.path(out_dir_sci, "Invasion_year_map.tiff"),
  plot        = gg_invasion_year_sci,
  device      = ragg::agg_tiff,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

ggsave(
  filename = file.path(out_dir_sci, "Invasion_year_map.svg"),
  plot     = gg_invasion_year_sci,
  device   = "svg",
  width    = width_cm,
  height   = height_cm,
  units    = "cm",
  bg       = "white"
)




###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
## Figure S1(B)  for supplementary section 
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################

cum_inv_EU_nuts = table(invasion_year_region$invasion_year)

cum_inv_EU_nuts = data.frame(Year = as.numeric(names(cum_inv_EU_nuts)),
                             Count = cumsum(cum_inv_EU_nuts))

cum_inv_EU_nuts

cum_inv_EU_nuts <- rbind(
  cum_inv_EU_nuts[1:9, ],
  data.frame(Year = 2019, Count = 427),
  cum_inv_EU_nuts[10:nrow(cum_inv_EU_nuts), ]
)

cum_inv_EU_nuts2 <- cum_inv_EU_nuts %>%
  dplyr::filter(Year <= 2023)

# Computing our cumulative NUTS3 invasion counts

p2 <- ggplot(cum_inv_EU_nuts2, aes(x = Year, y = Count)) +
  geom_col() +
  scale_y_continuous(
    name   = "Cumulative NUTS3 regions",
    breaks = c(0, 100, 200, 300, 400, 500, 600),
    labels = c(0, 100, 200, 300, 400, 500, 600),
    limits = c(0, 600),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    name   = "Year",
    breaks = seq(min(cum_inv_EU_nuts2$Year, na.rm = TRUE),
                 max(cum_inv_EU_nuts2$Year, na.rm = TRUE),
                 by = 2),
    expand = c(0.01, 0.01)
  ) +
  theme_classic(base_size = 9) +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_text(size = 9, face = "bold"),
    axis.text    = element_text(size = 9),
    plot.margin  = margin(8, 8, 8, 8)
  )

print(p2)

ext         <- "tiff"
dpi         <- 600
compression <- "lzw"
bg          <- "transparent"

width_in  <- 8.5
height_in <- 6

out_dir_save <- ".../VectorNet_ 2_original/Plots/A_final_plots/Inkscape_final_arranged_plots/supp_tiff_600_dpi/Plots_S1"
dir.create(out_dir_save, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(out_dir_save, paste0("cum_invaded_nuts3_bar.", ext))

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    filename    = outfile,
    plot        = p2,
    device      = ragg::agg_tiff,
    width       = width_in,
    height      = height_in,
    units       = "in",
    dpi         = dpi,
    compression = compression,
    background  = bg
  )
} else {
  ggsave(
    filename    = outfile,
    plot        = p2,
    device      = "tiff",
    width       = width_in,
    height      = height_in,
    units       = "in",
    dpi         = dpi,
    compression = compression,
    bg          = bg,
    type        = "cairo"
  )
}


# Applying journal formatting to the cumulative invasion plot
p2_sci <- ggplot(cum_inv_EU_nuts2, aes(x = Year, y = Count)) +
  geom_col(
    width     = 0.95,
    linewidth = 0.07
  ) +
  scale_y_continuous(
    name   = "Cumulative NUTS3 regions",
    breaks = c(0, 100, 200, 300, 400, 500, 600),
    labels = c(0, 100, 200, 300, 400, 500, 600),
    limits = c(0, 600),
    expand = c(0.0, 0.0)
  ) +
  scale_x_continuous(
    name   = "Year",
    breaks = seq(min(cum_inv_EU_nuts2$Year, na.rm = TRUE),
                 max(cum_inv_EU_nuts2$Year, na.rm = TRUE),
                 by = 2),
    expand = c(0.01, 0.01)
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),

    axis.title = element_text(size =  7, color = "black"),
    axis.text  = element_text(size = 6, color = "black"),

    axis.line  = element_line(linewidth = 0.20),
    axis.ticks = element_line(linewidth = 0.20),
    axis.ticks.length = unit(1.5, "pt"),

    panel.grid  = element_blank(),

    panel.border = element_blank()
  )

print(p2_sci)

out_dir_base <- ".../VectorNet_ 2_original/Plots/A_final_plots/Inkscape_final_arranged_plots/supp_tiff_600_dpi/Plots_S1"
out_dir_sci  <- file.path(out_dir_base, "journal")
dir.create(out_dir_sci, recursive = TRUE, showWarnings = FALSE)

width_cm  <- 5.7
height_cm <- 5.0

ggsave(
  filename    = file.path(out_dir_sci, "Cumulative_invaded_NUTS3.tiff"),
  plot        = p2_sci,
  device      = ragg::agg_tiff,
  width       = width_cm,
  height      = height_cm,
  units       = "cm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

ggsave(
  filename = file.path(out_dir_sci, "Cumulative_invaded_NUTS3.svg"),
  plot     = p2_sci,
  device   = "svg",
  width    = width_cm,
  height   = height_cm,
  units    = "cm",
  bg       = "white"
)




################################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################
# Figure S15   for supplementary section 
###############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################

m_train = load(".../VectorNet_ 2_original/Output/A-Final_output/Main_model_output/main_modelvIparam110000101011001110111100_SIMP_idvI.Rdata")
m_train

spde_results<- collected_results_spde$tb_region_results

orig_data = load("Data/data_INLA.RData")
nuts3_map = map_inla
nuts3_sf <- sf::st_as_sf(nuts3_map, sf_column_name = "geometry")

outline <- sf::st_as_sf(sf::st_union(sf::st_geometry(nuts3_sf)))

spde_df <- spde_results %>%
  select(geo, year_idx, spde_integrated)

years <- sort(unique(spde_df$year_idx))

vals_for_scale <- spde_df %>%
  filter(year_idx %in% years) %>%
  pull(spde_integrated)

lims_raw <- range(vals_for_scale, na.rm = TRUE, finite = TRUE)

step  <- 2
min_b <- floor(lims_raw[1] / step) * step
max_b <- ceiling(lims_raw[2] / step) * step
breaks <- seq(min_b, max_b, by = step)
lims   <- c(min_b, max_b)

palette_name <- "Rushmore1"; palette_n <- 100; palette_reverse <- FALSE
pal <- wesanderson::wes_palette(palette_name, palette_n, type = "continuous")
if (isTRUE(palette_reverse)) pal <- rev(pal)

na_fill = "grey95"
border_col = NA; border_size = NA
outer_boundary_col  = "black"; outer_boundary_size = 0.1


# Designing our yearly spatial-effect maps
make_spde_map <- function(yy) {

  sf_year <- nuts3_sf %>%
    dplyr::select(geo, geometry) %>%
    dplyr::left_join(spde_df %>% dplyr::filter(year_idx == yy), by = "geo")

  ggplot() +
    geom_sf(
      data = sf_year,
      aes(fill = spde_integrated),
      color = border_col, linewidth = border_size
    ) +
    scale_fill_gradientn(
      colors = pal,
      limits = lims,
      breaks = breaks,
      labels = scales::label_number(accuracy = 0.1),
      oob = scales::squish,
      na.value = na_fill,
      name = "SPDE effect"
    ) +
    geom_sf(
      data = outline,
      fill = NA,
      color = outer_boundary_col,
      linewidth = outer_boundary_size,
      inherit.aes = FALSE
    ) +
    coord_sf(datum = NA, expand = FALSE, clip = "on") +
    theme_void(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
}


# Designing our shared spatial-effect legend
make_spde_legend_plot <- function() {

  p_leg_src <- ggplot() +
    geom_sf(
      data = nuts3_sf %>% dplyr::slice(1),
      aes(fill = 0),
      color = NA
    ) +
    scale_fill_gradientn(
      colors = pal,
      limits = lims,
      breaks = breaks,
      labels = scales::label_number(accuracy = 0.1),
      oob = scales::squish,
      na.value = na_fill,
      name = NULL
    ) +
    theme_void(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 8)
    ) +
    guides(fill = guide_colorbar(
      barheight = unit(40, "mm"),
      barwidth  = unit(4.5, "mm"),
      ticks.colour = "white",
      frame.colour = NA
    ))

  leg_grob <- cowplot::get_legend(p_leg_src)
  cowplot::ggdraw(leg_grob)
}

legend_plot <- make_spde_legend_plot()

base_dir <- ".../VectorNet_ 2_original/Plots/A_final_plots/Spatial_effects_plot"

out_dir <- file.path(base_dir, "journal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

map_w_cm <- 5.7
map_h_cm <- 5.8

legend_w_cm <- 3.4
legend_h_cm <- 7.0

for (yy in years) {
  p_map <- make_spde_map(yy)

  ggsave(
    filename    = file.path(out_dir, sprintf("spde_integrated_year%02d_map.tiff", yy)),
    plot        = p_map,
    width       = map_w_cm, height = map_h_cm, units = "cm",
    dpi         = 600,
    compression = "lzw",
    device      = "tiff",
    bg          = "transparent"
  )
}

ggsave(file.path(out_dir, "spde_integrated__LEGEND.tiff"),
       legend_plot, width = legend_w_cm, height = legend_h_cm, units = "cm",
       device = "tiff", dpi = 600, compression = "lzw", bg = "transparent")

for (yy in years) {
  p_map <- make_spde_map(yy)

  ggsave(
    filename    = file.path(out_dir, sprintf("svg/spde_integrated_year%02d_map.svg", yy)),
    plot        = p_map,
    width       = map_w_cm,
    height = map_h_cm,
    units = "cm",
    dpi         = 600,

    device = "svg",
    bg          = "transparent"
  )
}

ggsave(file.path(out_dir, "svg/spde_integrated__LEGEND.svg"),
       legend_plot, width = legend_w_cm, height = legend_h_cm, units = "cm",
       device = "svg", dpi = 600,bg = "transparent")

