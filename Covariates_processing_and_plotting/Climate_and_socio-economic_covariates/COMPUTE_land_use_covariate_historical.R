# ================== Land-use covariate harmonisation ==================

setwd(".../VectorNet_ 2_original")

# Load functions and supporting objects
source("functions_general.R")
source("functions_plotting.R")
source("functions_plot_article.R")

# Load the spatial framework
out <- attach(".../Data/data_INLA.RData")
map_inla <- map_inla
detach(".../file:Data/data_INLA.RData")

# Assemble annual land-use tables
tb_list <- NULL
count <- 0

for (year in 2006:2023) {

  print(sprintf("Processing year: %d", year))

  land_use <-
    read.csv(paste0(".../land_use_", year, ".csv"))%>%
    as_tibble()

  if ("X" %in% colnames(land_use)) {
    land_use <- land_use %>% select(-X)
  }

  print(land_use)

  vals <- tibble(Location=land_use$Location,
                 LocationNa=map_inla$LocationNa) %>%
    mutate(id=1:n())

  vals$is_equal <- NA

  for (j in 1:nrow(vals)) {
    vals$is_equal[j] <- identical(vals$Location[j],vals$LocationNa[j])
  }

  vals %>% filter(!is_equal)

  land_use <- land_use %>% mutate(
    Location = case_when(
      Location == "Cote-d?Or" ~ "Cote-d’Or",
      Location == "Cotes-d?Armor" ~ "Cotes-d’Armor",
      TRUE ~ Location
    )
  )

  tb_land_use <- tibble(
    geo=map_inla$geo,
    LocationNa=land_use$Location,
    primary_forested_area = land_use$primary.foreseted.area,
    primary_non_forested_area = land_use$primary.non.foreseted.area,
    secondary_forested_area = land_use$secondary.foreseted.area,
    croplands = land_use$croplands,
    pasture_land = land_use$pasture.land,
    rangelands = land_use$rangelands,
    urban_land = land_use$urban.land
  )

  count <- count + 1
  tb_list[[count]] <- tb_land_use

}


## we need yearly data as covariates 
tb_list <- lapply(seq_along(tb_list), function(i) {
  tb_list[[i]] %>% mutate(year_idx = i)
})

combined_df <- bind_rows(tb_list)

combined_df <- combined_df %>%
  relocate(year_idx, .after = geo)

print(combined_df)

unique(combined_df$year_idx)


#checking the dataset created 
combined_df[combined_df$year_idx == 18,]
tb_list[[18]]

combined_df <- combined_df %>%
  mutate(year = 2006 + year_idx - 1) %>%
  relocate(year, .after = year_idx)

unique(combined_df$year)

# Add a column with labels based on year and scenario (Since ISIMIP3b and LUH (Land use Harmonization) has only historical data till 2014 and after that ssp370 dataset is added))
combined_df <- combined_df %>%
  mutate(scenario = ifelse(year <= 2014, "historical", "SSP370")) %>%
  relocate(scenario, .after = year)

print(combined_df)




# Estimate scenario offsets and rescale post-2014 series
process_column <- function(data, column_name) {

  region_models <- data %>%
    group_by(geo) %>%
    group_map(~ glm(as.formula(paste(column_name, "~ year + scenario")), data = .)) %>%
    setNames(sort(unique(data$geo)))

  region_models$UKI34      ## note that glm models contain coefficient (along with glm description) in alphabetical order of geo 

  glm_data <- data %>% filter(geo == "UKI34")

  glm_model <- glm(as.formula(paste(column_name, "~ year + scenario")), data = glm_data)

  glm_model

  region_coeff <- NULL
  individual_coeff <- NULL

  region_coeff <- coef(region_models$UKI34)
  individual_coeff <- coef(glm_model)

  # Compare coefficients
  if (all.equal(region_coeff, individual_coeff, tolerance = 1e-6)) {
    message(paste0("The glm coefficient in region-wise approach matches with each individual region fitting for region for covariate ", column_name))
  } else {
    message(paste0("The glm coefficient does NOT match for region for covariate", column_name))
  }

  region_coefficients <- sapply(region_models, function(model) coef(model)["scenarioSSP370"])

  names(region_coefficients) <- gsub("\\.scenarioSSP370", "", names(region_coefficients))

  region_coefficients_df <- data.frame(
    geo = names(region_coefficients),
    region_coefficient = unlist(region_coefficients)
  )

  rownames(region_coefficients_df) <- NULL

  return(region_coefficients_df)
}



columns_to_process <- names(combined_df)[
  grepl("forested|non|secondary|croplands|pasture|rangelands|urban", names(combined_df))
]


# Initialize adjusted_data with the original combined_df
adjusted_data <- combined_df


# Loop through each column to process
for (column_name in columns_to_process) {
  region_coefficients_df <- process_column(combined_df, column_name)

  # Join region coefficients to the main data
  adjusted_data <- adjusted_data %>%
    left_join(region_coefficients_df, by = "geo")

  J = adjusted_data[adjusted_data$geo == "UKI34", ]
  J$region_coefficient

  # Adjust SSP data for each region
  adjusted_data <- adjusted_data %>%
    group_by(geo) %>%
    mutate(
      !!paste0("scaled_",column_name) := if_else(
        scenario == "SSP370",
        !!sym(column_name) - as.numeric(region_coefficient),
        !!sym(column_name)
      )
    )%>%
    ungroup()        ## ungrouping it because we don't want to carry the grouping information to further lines of code

  adjusted_data <- adjusted_data %>%
    select(-region_coefficient)

  # Debug print (optional)
  print(paste("Adjustment of ", column_name, "completed."))
}

print(adjusted_data)




############ Plotting European aggregates
create_aggregated_plot <- function(adjusted_data, map_inla, variable) {

  combined_mean_with_area <- adjusted_data %>%
    left_join(map_inla %>% select(geo, area), by = "geo")

  aggregated_data <- combined_mean_with_area %>%
    group_by(year) %>%
    summarise(
      weighted_mean = sum(get(variable) * area, na.rm = TRUE) /
        sum(area, na.rm = TRUE)
    )

  # Determine the scaled variable name dynamically  (add scaled_ to the main varible name)
  scaled_variable <- paste0("scaled_", variable)

  combined_scaled_mean_with_area <- adjusted_data %>%
    left_join(map_inla %>% select(geo, area), by = "geo")

  # Aggregate data for all regions using area for scaled weighted mean
  aggregated_data_scaled <- combined_scaled_mean_with_area %>%
    group_by(year) %>%
    summarise(
      weighted_scaled_mean = sum(get(scaled_variable) * area, na.rm = TRUE) /
        sum(area, na.rm = TRUE)
      )

  # Remove '_mean' from variable name and capitalize the first letter
  clean_variable_name <- function(var_name) {
    gsub(" Mean$", "", tools::toTitleCase(gsub("_", " ", var_name)))
  }

  aggregated_data_combined <- aggregated_data %>%
    rename(mean = weighted_mean) %>%
    mutate(variable = ifelse(year >= 2015,
                             paste(clean_variable_name(variable), "SSP370", sep = " "),
                             paste(clean_variable_name(variable), "Historical", sep = " "))) %>%
    bind_rows(
      aggregated_data_scaled %>%
        rename(mean = weighted_scaled_mean) %>%
        mutate(variable = paste(clean_variable_name(variable), "Historical (scaled)", sep = " "))
    )

  historical_name <- paste(clean_variable_name(variable), "Historical", sep = " ")
  ssp370_name <- paste(clean_variable_name(variable), "SSP370", sep = " ")
  combined_name <- paste(clean_variable_name(variable), "Historical & SSP370", sep = " ")

  # Use precomputed strings in mutate()
  aggregated_data_combined <- aggregated_data_combined %>%
    mutate(group = ifelse(
      variable %in% c(historical_name, ssp370_name),
      combined_name,
      variable
    ))

  # Define color mapping explicitly
  variable_name <- clean_variable_name(variable)
  historical_name <- paste(variable_name, "Historical", sep = " ")
  ssp370_name <- paste(variable_name, "SSP370", sep = " ")
  scaled_name <- paste(variable_name, "Historical (scaled)", sep = " ")

  color_mapping <- setNames(
    c("blue", "red", "green"),
    c(historical_name, ssp370_name, scaled_name)
  )

  # Plot the combined data
  ggplot(aggregated_data_combined, aes(x = year, y = mean, color = variable, group = group)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = paste("Total", clean_variable_name(variable)),
      x = "Year",
      y = paste(clean_variable_name(variable), "fraction (of EU)"),
      color = "Climate scenario's"
    ) +
    scale_color_manual(values = color_mapping) +
    scale_x_continuous(
      limits = c(min(aggregated_data_combined$year), max(aggregated_data_combined$year) + 2),
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(5, "pt"),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 12)
    )
}

# Plotting the function
create_aggregated_plot(adjusted_data, map_inla, "primary_forested_area")
create_aggregated_plot(adjusted_data, map_inla, "primary_non_forested_area")
create_aggregated_plot(adjusted_data, map_inla, "urban_land")
create_aggregated_plot(adjusted_data, map_inla, "croplands")
create_aggregated_plot(adjusted_data, map_inla, "pasture_land")
create_aggregated_plot(adjusted_data, map_inla, "rangelands")




##!!!!!! Plotting for the article - 
## Run below code only if we want to save it for the article

legend_name <- "legend"       ## name of legend

land_use_covariate = "croplands"
y_axis_name <- "Croplands"

y_axis_name <- paste0(y_axis_name, "\n(all categories)")

# # for urban_land areas
# land_use_covariate = "urban_land"
# y_axis_name <- "Urban lands"
# 
# land_use_covariate = "primary_forested_area"
# y_axis_name <- "Primary forested areas"
# 
# land_use_covariate = "pasture_land"
# y_axis_name <- "Pasture lands"
# 
# land_use_covariate = "rangelands"
# y_axis_name <- "Range lands"
# 
# land_use_covariate = "secondary_forested_area"
# y_axis_name <- "Secondary forested areas"

create_aggregated_plot_article(adjusted_data,
                               map_inla,
                               variable = land_use_covariate,
                               y_axis_label = y_axis_name,
                               break_year = 2015,
                               save_filename = paste0(".../land_use_covariates_scaled", land_use_covariate, ".tiff"),
                               width_in = 6, height_in = 4, dpi = 600,
                               compression = "lzw", bg = "transparent",
                               legend_width_in = 2.8,
                               legend_height_in = 0.8,
                               legend_dpi = 600,

                               legend_filename <- NULL,

                               legend_compression = "lzw",
                               legend_bg = "transparent" )


# Retain the scaled covariates used in the analysis
adjusted_data_final <- adjusted_data %>%
  select(geo, year_idx, year, LocationNa, starts_with("scaled"), )

adjusted_data_final <- adjusted_data_final %>%
  rename_with(~ gsub("^scaled_", "", .), starts_with("scaled_"))

adjusted_data_final <- adjusted_data_final %>%
  filter(!(year >= 2006 & year <= 2009)) %>%
  mutate(year_idx = year_idx - min(year_idx) + 1)

adjusted_data_final

land_use <- adjusted_data_final

# Save the harmonised covariate
save("land_use",
     file=paste0(".../Current_data/Final_hist_land_use_scaled_data_yearly.RData"))





##### Plotting the covariate maps
source("functions_plot_article.R")

out <- load(".../Current_data/Final_hist_land_use_scaled_data_yearly.RData")

## The loaded object is `land_use`; rename to `land_use_factual` for clarity
land_use_factual <- land_use

land_use_factual <- land_use_factual %>%
  dplyr::filter(year_idx > 1) %>%
  dplyr::mutate(year_idx = year_idx - 1)

unique(land_use_factual$year_idx)

land_use_factual <- land_use_factual %>%
  dplyr::filter(year_idx <= 11)          # keep only years 2011–2021 (consistent with counterfactual dataset)

unique(land_use_factual$year_idx)

## 4) Add a calendar year column from year_idx (1 → 2011, ..., 11 → 2021)
land_use_factual <- land_use_factual %>%
  dplyr::mutate(year = 2011 + year_idx - 1) %>%       # map 1→2011, 2→2012, ..., 11→2021
  dplyr::relocate(year, .after = geo)

unique(land_use_factual$year)


# ---- Year selection ----
# <- numeric year to plot
# year      <- 2017        ## for any particular year
# year <- NULL             ## for aggregation across all years
year <- NULL

#variable <- "primary_forested_area"   # change as needed
#variable <- "urban_land"
#variable <- "primary_non_forested_area"
#variable <- "secondary_forested_area"
#variable <- "pasture_land"
variable <- "rangelands"


base_dir <- ".../map_plot_covariates"

filename_stub_hist <- paste0(
  base_dir, "/",
  variable, "/",        # folder "primary_forested_area", "pasture_land", etc.
  variable
)

dir.create(dirname(filename_stub_hist), recursive = TRUE, showWarnings = FALSE)

## Palette length & direction
palette_n            <- 100
palette_reverse_case <- FALSE


## Palette choices for different land-use types
if (variable == "primary_forested_area") {
  pal_seq <- c(
    "#E5E5E5",
    grDevices::colorRampPalette(
      c("gray97", "chartreuse4"),
      bias = 1
    )(palette_n)
  )

} else if (variable == "secondary_forested_area") {
  pal_seq <- c(
    "#E5E5E5",
    grDevices::colorRampPalette(
      c("gray97", "olivedrab3"),
      bias = 1
    )(palette_n)
  )

} else if (variable == "primary_non_forested_area") {
  pal_seq <- c(
    "#E5E5E5",
    grDevices::colorRampPalette(
      c("gray97", "darkseagreen4"),
      bias = 1
    )(palette_n)
  )

} else if (variable %in% c("pasture_land")) {
  pal_seq <- c(
    "#E5E5E5",
    grDevices::colorRampPalette(
      c("gray97", "burlywood3"),
      bias = 1
    )(palette_n)
  )

} else if (variable %in% c("rangelands")) {
  pal_seq <- c(
    "#E5E5E5",
    grDevices::colorRampPalette(
      c("gray97", "navajowhite4"),
      bias = 1
    )(palette_n)
  )

}else if (variable == "urban_land") {
  pal_seq <- c(
    "#E5E5E5",
    grDevices::colorRampPalette(
      c("gray97", "saddlebrown"),
      bias = 1
    )(palette_n)
  )

} else {
  # Fallback for any other land-use-like variable
  pal_seq <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(9, "YlGnBu")
  )(palette_n)
}


## plotting
plot_land_use_maps(
  adjusted_df = land_use_factual,
  map_inla    = map_inla,
  year        = year,
  variable    = variable,

  filename_stub = filename_stub_hist,
  ext           = "tiff",
  width_in      = 3.5,
  height_in     = 4,
  dpi           = 600,
  compression   = "lzw",
  bg            = "transparent",

  hist_legend_ext         = "tiff",
  hist_legend_dpi         = 600,
  hist_legend_filename    = NULL,
  hist_legend_width_in    = 1.4,
  hist_legend_height_in   = 5.0,
  hist_legend_compression = "lzw",
  hist_legend_bg          = "transparent",

  na_fill            = "grey95",
  border_col         = NA,
  border_size        = NA,
  outer_boundary_col = "black",
  outer_boundary_size = 0.1,

  palette_name    = "custom",
  palette_n       = palette_n,
  palette_reverse = palette_reverse_case,
  pal_seq         = pal_seq
)
## If a warning appears that means we are getting different legends for each seasonal aggregation of climate variable.