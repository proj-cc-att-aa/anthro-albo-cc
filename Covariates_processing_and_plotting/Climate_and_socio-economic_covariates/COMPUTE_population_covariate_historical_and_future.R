# ================== Population covariate harmonisation ==================

setwd(".../VectorNet_ 2_original")

# Load functions and supporting objects
source("functions_general.R")
source("functions_plotting.R")
source("functions_plot_article.R")

# Load the spatial framework
out <- attach(".../Data/data_INLA.RData")
map_inla <- map_inla
detach(".../file:Data/data_INLA.RData")

# Assemble annual population tables
tb_list <- NULL
count <- 0

for (year in 2010:2022) {

  print(sprintf("Processing year: %d", year))

  population <-
    read.csv(paste0(".../population_", year, ".csv"))%>%
    as_tibble()

  if ("X" %in% colnames(population)) {
    population <- population %>% select(-X)
  }

  print(population)

  population <- population %>%
    select(Location, `total.population`)

  vals <- tibble(Location=population$Location,
                 LocationNa=map_inla$LocationNa) %>%
    mutate(id=1:n())

  vals$is_equal <- NA

  for (j in 1:nrow(vals)) {
    vals$is_equal[j] <- identical(vals$Location[j],vals$LocationNa[j])
  }

  vals %>% filter(!is_equal)

  population <- population %>% mutate(
    Location = case_when(
      Location == "Cote-d?Or" ~ "Cote-d’Or",
      Location == "Cotes-d?Armor" ~ "Cotes-d’Armor",
      TRUE ~ Location
    )
  )

  ## Since we have directly the annual mean data unlike the daily data
  tb_population <- tibble(
    geo=map_inla$geo,
    LocationNa=population$Location,
    population = population$total.population
  )

  count <- count + 1
  tb_list[[count]] <- tb_population

}

tb_list <- lapply(seq_along(tb_list), function(i) {
  tb_list[[i]] %>% mutate(year_idx = i)
})




combined_df <- bind_rows(tb_list)

combined_df <- combined_df %>%
  relocate(year_idx, .after = geo)

print(combined_df)

unique(combined_df$year_idx)

#checking the dataset created 
combined_df[combined_df$year_idx == 13,]
tb_list[[13]]

combined_df <- combined_df %>%
  mutate(year = 2010 + year_idx - 1) %>%
  relocate(year, .after = year_idx)

unique(combined_df$year)



# Extrapolate one additional year by region (for year 2023)
if (!("year" %in% names(combined_df))) {
  combined_df <- combined_df %>% mutate(year = 2010 - 1 + year_idx)
}
combined_df <- combined_df %>% mutate(year = as.integer(year))

last_year     <- max(combined_df$year, na.rm = TRUE)
new_year      <- last_year + 1
new_year_idx  <- max(combined_df$year_idx, na.rm = TRUE) + 1

static_geo <- combined_df %>%
  group_by(geo) %>%
  summarise(
    LocationNa = first(LocationNa),
    .groups = "drop"
  )

fit_pop_glm_all_years <- function(d) {
  d <- d %>%
    arrange(year) %>%
    filter(is.finite(population), population > 0)

  if (nrow(d) < 3) return(NULL)

  glm(population ~ year, data = d, family = gaussian(link = "log"))
}

predict_pop_1yr <- function(mod, year_val) {
  if (is.null(mod)) {
    return(tibble(population = NA_real_, pop_lwr = NA_real_, pop_upr = NA_real_))
  }
  pr <- predict(mod, newdata = data.frame(year = year_val), type = "link", se.fit = TRUE)

  tibble(
    population = exp(pr$fit),
    pop_lwr    = exp(pr$fit - 1.96 * pr$se.fit),
    pop_upr    = exp(pr$fit + 1.96 * pr$se.fit)
  )
}

df_pop_2023 <- combined_df %>%
  group_by(geo) %>%
  nest() %>%
  mutate(
    model = map(data, fit_pop_glm_all_years),
    pred  = map(model, ~ predict_pop_1yr(.x, new_year))
  ) %>%
  select(geo, pred) %>%
  unnest(pred) %>%
  left_join(static_geo, by = "geo") %>%
  mutate(
    year     = new_year,
    year_idx = new_year_idx,
    is_extrapolated = TRUE
  ) %>%
  select(geo, year_idx, year, LocationNa, population, pop_lwr, pop_upr, is_extrapolated)

# Optional safety cap: limit 2022→2023 change 
CAP_LO <- 0.85
CAP_HI <- 1.15

last_pop <- combined_df %>%
  group_by(geo) %>%
  filter(year == last_year) %>%
  summarise(pop_last = first(population), .groups = "drop")

df_pop_2023 <- df_pop_2023 %>%
  left_join(last_pop, by = "geo") %>%
  mutate(
    population = if_else(!is.finite(population), pop_last, population),

    population = pmin(pmax(population, pop_last * CAP_LO), pop_last * CAP_HI)
  ) %>%
  select(-pop_last)

geo_levels <- combined_df %>% distinct(geo) %>% pull(geo)

df_pop_2023_clean <- df_pop_2023 %>%
  select(all_of(names(combined_df)))

combined_df_ext <- bind_rows(combined_df, df_pop_2023_clean) %>%
  mutate(
    geo = factor(geo, levels = geo_levels)
    ) %>%
  mutate(geo = as.character(geo))

unique(combined_df_ext$year)
combined_df_ext %>% filter(year == 2023) %>% count()

population <- combined_df_ext

# Derive population density and area
map_inla <- map_inla %>%
  mutate(area_km_2 = area / 1e6)

nuts_area <- map_inla %>%
  st_drop_geometry() %>%
  select(geo, area_km_2)

population <- population %>%
  left_join(nuts_area, by = "geo") %>%
  mutate(population_density = population / area_km_2)

population <- population %>%
  dplyr::rename(area = area_km_2)

population

# Save the harmonised covariate
save("population",
     file=paste0(".../Current_data/Final_hist_population_data_yearly.RData"))





#### PLOTTING POPULATION DATA MAP ##################
pop_globpop_out <- load(".../Current_data/Final_hist_population_data_yearly.RData")
pop_globpop <- population
pop_globpop

variable <- "population_density"

base_dir <- ".../map_plot_covariates"

filename_stub_hist <- paste0(
  base_dir, "/",
  variable, "/",
  variable
)

dir.create(dirname(filename_stub_hist), recursive = TRUE, showWarnings = FALSE)


# Palette name, length & direction
palette_name <- "Blues"
palette_n <- 100
palette_reverse_case <- FALSE


### double check step
stopifnot(all(c("geo", "year", "population", "population_density") %in% names(pop_globpop)))



## plotting
plot_population_maps(
  adjusted_df = pop_globpop,
  map_inla    = map_inla,

  year        = NULL,            # aggregated mean across years  (Output names will use "__aggregated__")
  variable    = variable,
  filename_stub = filename_stub_hist,

  ext         = "tiff",
  width_in    = 3.5,
  height_in   = 4,
  dpi         = 600,
  compression = "lzw",
  bg          = "transparent",

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

  palette_name    = palette_name,
  palette_n       = palette_n,
  palette_reverse = palette_reverse_case
)



##switch variable cleanly 
# If we  want all 3 products (population, density, log-density)
vars_to_plot <- c("population", "population_density", "log_population_density")


for (v in vars_to_plot) {

  variable <- v
  filename_stub_hist <- file.path(base_dir, variable, variable)
  dir.create(dirname(filename_stub_hist), recursive = TRUE, showWarnings = FALSE)

  plot_population_maps(
    adjusted_df = pop_globpop,
    map_inla    = map_inla,

    year        = NULL,
    variable    = variable,
    filename_stub = filename_stub_hist,
    ext         = "tiff",
    width_in    = 3.5,
    height_in   = 4,
    dpi         = 600,
    compression = "lzw",
    bg          = "transparent",
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
    palette_name    = palette_name,
    palette_n       = palette_n,
    palette_reverse = palette_reverse_case
  )
}


## loop over multiple years
# This will save one map + one legend per year.
 years_to_plot <- sort(unique(pop_globpop$year))

for (yy in years_to_plot) {
  plot_population_maps(
    adjusted_df = pop_globpop,
    map_inla    = map_inla,

    year        = yy,
    variable    = variable,
    filename_stub = filename_stub_hist,
    ext         = "tiff",
    width_in    = 3.5,
    height_in   = 4,
    dpi         = 600,
    compression = "lzw",
    bg          = "transparent",
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
    palette_name    = palette_name,
    palette_n       = palette_n,
    palette_reverse = palette_reverse_case
  )
}
