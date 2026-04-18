# Delta-method alignment of future climate covariates

# Setup -------------------------------------------------------------------
setwd(".../VectorNet_ 2_original")    

packages <- c("tidyverse", "stringr", "caret")
lapply(packages, library, character.only = TRUE)


## loading the required functions
library(dplyr)
source("functions_general.R")
source("functions_plotting.R")



#Delta method 

##Logic ---- ###
#Instead of extracting 2021 obs and 2021 model as anchors, we compute simple means over overlap years 2015–2021:
#obsBase = mean(obs, 2015–2021) per geo
##sspBase = mean(ensemble SSP of all models, 2015–2021) per geo

#Then for each future year (≥2022) we do:
#temp: obsBase + (ssp - sspBase)
#precip: obsBase * (ssp / sspBase)
#RH: logit(obsBase) + (logit(ssp) - logit(sspBase))



## Code alignment
# we load all 5 model folders for a given SSP (e.g., SSP1),
# we compute the ensemble mean (per geo × year_idx) across models
# we compute overlap baselines (for observational overlap as well as ssp's overlap (used for subtraction from each year ssp's data )) using simple means over 2015–2021 (not a single year to avoid artificial jumps),
# we stitch future years (≥2022) onto the observational scale:
# temperature: additive
# precipitation: multiplicative ratio (as per the paper of delta method)
# humidity: logit-delta (bounded)

setwd(".../VectorNet_ 2_original/Data")

covariates_folder_data <- "climate/Final_covariates_data/Current_data"
factual_file <- file.path(covariates_folder_data,
                          "Final_20CRV3-ERA5_factual_climate_data_yearly.RData")

future_root <- "...VectorNet_ 2_original/Data/climate/Final_covariates_data/Future_data"

out_root <- file.path(future_root, "ENSEMBLE_STITCHED")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Load factual obs (2010..2023)
load(factual_file)  #climate_summaries_region
df_climate <- climate_summaries_region %>%
  filter(year_idx > 9) %>%           # we drop 1..9 so that we start at 2010
  mutate(year_idx = year_idx - 9)    # we reindex so 2010 -> 1, ..., 2023 -> 14

df_obs <- df_climate %>%
  mutate(year = 2010 + year_idx - 1) %>%
  filter(year <= 2021)               # we keep only the training-era observational window


# function to load all model files for an SSP and build ensemble mean
load_ssp_ensemble <- function(ssp_folder, ssp_tag) {
  
  model_dirs <- list.dirs(ssp_folder, recursive = FALSE, full.names = TRUE)
  
  model_files <- map_chr(model_dirs, function(d) {
    ff <- list.files(d, pattern = "\\.RData$", full.names = TRUE)
    if (length(ff) == 0) stop("No .RData file found in: ", d)
    if (length(ff) > 1) {
      ff2 <- ff[str_detect(basename(ff), "Final_.*_climate_data_yearly_2015_2100\\.RData")]
      if (length(ff2) >= 1) return(ff2[1])
    }
    ff[1]
  })
  
  lst <- map2(model_files, model_dirs, function(f, d) {
    load(f)  
    climate_summaries_region %>%
      mutate(model = basename(d),    ## create a column which give model tags to all files
             ssp   = ssp_tag)        ## give corresponding ssp tag
  })
  
  df_all <- bind_rows(lst)  ## bind all the 5 model list (each model has data from 2015 to 2100) created above the only thing separating them is model tag (for any ssp's)
  
  key_order <- df_all %>%
    distinct(geo, year_idx) %>%     ## of all distinct geo and year_idx (the exact order is recorded)
    mutate(.ord = row_number())
  
  # we compute our ensemble mean per geo × year_idx across models
  num_cols <- names(df_all)[sapply(df_all, is.numeric)]       ##  all numeric columns  (except location)
  num_cols <- setdiff(num_cols, c("year_idx")) # we do not average year_idx itself (seperate out year_idx column)
  
  # computing ensemble across all models and then re-apply the same ordering
  df_ens <- df_all %>%
    group_by(geo, year_idx) %>%     ## group by each region and each year (now we have five models for each region year combination, so we take ensemble mean across 5 models)
    summarise(
      LocationNa = first(LocationNa),
      across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(key_order, by = c("geo", "year_idx")) %>%
    arrange(.ord) %>%                         # we restore df_all ordering (unique geo and year_idx column in actual data )
    select(-.ord) %>%
    mutate(
      year = 2015 + year_idx - 1,    ## add one year column too....year_idx = 1 == year 2015 in our ssp files
      ssp  = ssp_tag
    )
  
  df_ens
}

# Overlap-mean stitcher (simple means over overlap window)
#    - temperature: additive
#    - precipitation: multiplicative ratio
#    - humidity: logit-delta
stitch_by_overlap_means <- function(df_obs, df_ssp_ens,
                                    overlap_years = 2015:2021,
                                    future_year_min = 2022,
                                    eps_pr = 1e-6,
                                    eps_rh = 1e-4) {
  
  # df_obs         = df_obs
  # df_ssp_ens     = df_ssp_ens
  # overlap_years  = overlap_years
  # future_year_min = 2022
  # eps_pr = 1e-6
  # eps_rh = 1e-4
  
  covar_cols <- names(df_obs)[sapply(df_obs, is.numeric)]
  covar_cols <- setdiff(covar_cols, c("year_idx", "year"))
  
  
  # helper transforms for RH
  logit   <- function(p) log(p / (1 - p))
  expit   <- function(z) 1 / (1 + exp(-z))
  clamp01 <- function(x, eps) pmin(pmax(x, eps), 1 - eps)
  
  rh_cols <- covar_cols[str_detect(covar_cols, regex("humidity", ignore_case = TRUE))]
  rh_scale <- setNames(rep(100, length(rh_cols)), rh_cols)
  if (length(rh_cols) > 0) {
    for (cc in rh_cols) {
      mx <- max(df_obs[[cc]], df_ssp_ens[[cc]], na.rm = TRUE)
      if (is.finite(mx) && mx <= 1.5) rh_scale[[cc]] <- 1      # if max value of rh in tha column is less than 1.5 then put rh_scale list as 1 otherwise 100 as default percentage
    }
  }
  
  
  ## observational baseline means over overlap years
  obs_base <- df_obs %>%
    filter(year %in% overlap_years) %>%
    group_by(geo) %>%
    summarise(
      across(all_of(covar_cols), ~ mean(.x, na.rm = TRUE)),    ## mean of observational basleine for overlap years
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(.x, "_obsBase"), all_of(covar_cols))   ## rename each column with observation in that extracted baseline year data (mean over 2015 - 2021)
  
  # We extract the model reference (2015-2021) from the future file
  ssp_base <- df_ssp_ens %>%
    filter(year %in% overlap_years) %>% 
    group_by(geo) %>%
    summarise(
      across(all_of(covar_cols), ~ mean(.x, na.rm = TRUE)),   ## mean of model baseline for overlap years
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(.x, "_sspBase"), all_of(covar_cols))
  
  # stitch future years onto observational baseline
  df_fut <- df_ssp_ens %>%
    filter(year >= future_year_min) %>%
    
    left_join(obs_base, by = "geo") %>%   # we attach obs overlap means to every future row (per geo)
    left_join(ssp_base, by = "geo")       # we attach SSP overlap means to every future row (per geo)
  ## The order does not matter here, because left_join(..., by = "geo") matches rows by the key geo, not by row position.
  
  for (cc in covar_cols) {
    
    obs_cc <- paste0(cc, "_obsBase")    #baseline observed data
    ssp_cc <- paste0(cc, "_sspBase")    #model reference data as per scenario
    
    # temperature: additive stitch
    if (str_detect(cc, regex("temperature", ignore_case = TRUE))) {
      
      df_fut[[cc]] <- df_fut[[obs_cc]] + (df_fut[[cc]] - df_fut[[ssp_cc]])
      
      # precipitation: multiplicative ratio stitch
    } else if (str_detect(cc, regex("precipitation", ignore_case = TRUE))) {
      
      denom <- pmax(df_fut[[ssp_cc]], eps_pr)
      df_fut[[cc]] <- df_fut[[obs_cc]] * (df_fut[[cc]] / denom)
      
      # humidity: logit-delta (bounded)
    } else if (str_detect(cc, regex("humidity", ignore_case = TRUE))) {
      
      sc <- if (!is.null(rh_scale[[cc]])) rh_scale[[cc]] else 100
      
      h_obs <- clamp01(df_fut[[obs_cc]] / sc, eps_rh)
      h_ref <- clamp01(df_fut[[ssp_cc]] / sc, eps_rh)
      h_mod <- clamp01(df_fut[[cc]]     / sc, eps_rh)
      
      z_obs <- logit(h_obs)
      z_ref <- logit(h_ref)
      z_mod <- logit(h_mod)
      
      z_corr <- z_obs + (z_mod - z_ref)
      df_fut[[cc]] <- sc * expit(z_corr)
      
    } else {
      
      # !caution step
      df_fut[[cc]] <- df_fut[[cc]]
    }
  }
  
  # dropping helper baseline columns and reindexing year_idx so 2022 -> 1, ...
  df_out <- df_fut %>%
    dplyr::select(-ends_with("_obsBase"), -ends_with("_sspBase")) %>%
    mutate(year_idx = year - future_year_min + 1L) %>%      ## reindexing year_idx again from 1
    dplyr::select(geo, year_idx, year, LocationNa, all_of(covar_cols), ssp)
  
  df_out
}


#  Run for SSPs 
#    - we build ensemble mean first
#    - then we stitch using overlap years 2015..2021

overlap_years <- 2015:2021

ssp_folders <- c("SSP1", "SSP3", "SSP5")

future_corrected_ens <- list()

for (ssp_tag in ssp_folders) {
  
  #ssp_tag <-  "SSP1"
  
  ssp_folder <- file.path(future_root, ssp_tag)
  
  # -ensemble mean across models (2015..2100) ----
  df_ssp_ens <- load_ssp_ensemble(ssp_folder = ssp_folder, ssp_tag = ssp_tag)
  
  # saving ensemble file (optional, but useful for debugging later )
  ens_outfile <- file.path(out_root, paste0("EnsembleMean_", ssp_tag, "_climate_data_yearly_2015_2100.RData"))
  climate_summaries_region_ensemble <- df_ssp_ens
  save(climate_summaries_region_ensemble, file = ens_outfile)
  
  # stitch future years (2022..2100) using overlap means 2015..2021 ----
  df_ssp_stitched <- stitch_by_overlap_means(
    df_obs         = df_obs,
    df_ssp_ens     = df_ssp_ens,
    overlap_years  = overlap_years,
    future_year_min = 2022
  )
  
  future_corrected_ens[[ssp_tag]] <- df_ssp_stitched
  stitched_outfile <- file.path(out_root, paste0("Final_ensemble_", ssp_tag, "_delta_corrected_climate_data_yearly_2022_2100.RData"))
  climate_summaries_region_stitched <- df_ssp_stitched
  save(climate_summaries_region_stitched, file = stitched_outfile)
}




## Plotting

# factual aggregated series (<= 2021)
df_factual_summer <- df_obs %>%
  group_by(year) %>%
  summarise(
    temp_summer = mean(temperature_summer_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(scenario = "Factual (20CRv3/ERA5)")

agg_future_one <- function(df_future_stitched, sc_label) {
  df_future_stitched %>%
    group_by(year) %>%
    summarise(
      temp_summer = mean(temperature_summer_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(scenario = sc_label)
}

df_ssp1 <- agg_future_one(future_corrected_ens[["SSP1"]], "SSP1 ensemble (overlap-mean (2015-2021))")
df_ssp3 <- agg_future_one(future_corrected_ens[["SSP3"]], "SSP3 ensemble (overlap-mean (2015-2021))")
df_ssp5 <- agg_future_one(future_corrected_ens[["SSP5"]], "SSP5 ensemble (overlap-mean (2015-2021))")

df_plot <- bind_rows(df_factual_summer, df_ssp1, df_ssp3, df_ssp5) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "Factual (20CRv3/ERA5)",
        
        "SSP1 ensemble (overlap-mean (2015-2021))",
        "SSP3 ensemble (overlap-mean (2015-2021))",
        "SSP5 ensemble (overlap-mean (2015-2021))"
      )
    )
  )

p <- ggplot(df_plot, aes(x = year, y = temp_summer, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 2021, linetype = "dashed", linewidth = 0.4) +
  labs(
    x = "Year",
    y = "Mean temperature (averaged over NUTS3)",
    colour = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.height = grid::unit(0.55, "cm"),
    legend.key.width  = grid::unit(0.90, "cm"),
    axis.line  = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3)
  )

print(p)





# Parallel workflow for Pi_control, fully consistent with what we did for SSPs:
## we load counterfactual observational (20CRV3–ERA5 counterfactual) and keep 2010–2021
## we load Pi_control future files from Pi_control/ (5 model folders)
## we compute ensemble mean across the 5 models (per geo × year_idx)
## we compute overlap means over 2015–2021
## we stitch Pi_control 2022–2100 onto the counterfactual observational baseline
##  we save outputs in the same ENSEMBLE_STITCHED folder with Pi_control names


setwd("...VectorNet_ 2_original/Data")

covariates_folder_data <- "climate/Final_covariates_data/Current_data"
cf_file <- file.path(covariates_folder_data,
                     "Final_20CRV3-ERA5_counterfactual_climate_data_yearly.RData")

future_root <- ".../Future_data"
pi_folder   <- file.path(future_root, "Pi_control")

out_root <- file.path(future_root, "ENSEMBLE_STITCHED")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Load counterfactual observational data (2010..2021)
load(cf_file)  # 

df_climate_cf <- climate_summaries_region %>%
  filter(year_idx > 9) %>%            # we drop 1..9 so that we start at 2010
  mutate(year_idx = year_idx - 9)     # we reindex so 2010 -> 1, ..., 2021 -> 12

df_cf_obs <- df_climate_cf %>%
  mutate(year = 2010 + year_idx - 1) %>%
  filter(year <= 2021)               # we keep only the observational counter factual window


# function to load all model files for Pi_control and build ensemble mean
load_pi_control_ensemble <- function(pi_folder, scen_tag = "Pi_control") {
  
  model_dirs <- list.dirs(pi_folder, recursive = FALSE, full.names = TRUE)
  
  model_files <- map_chr(model_dirs, function(d) {
    ff <- list.files(d, pattern = "\\.RData$", full.names = TRUE)
    if (length(ff) == 0) stop("No .RData file found in: ", d)
    
    if (length(ff) > 1) {
      ff2 <- ff[str_detect(basename(ff), "Final_.*_climate_data_yearly_2015_2100\\.RData")]
      if (length(ff2) >= 1) return(ff2[1])
    }
    ff[1]
  })
  
  lst <- map2(model_files, model_dirs, function(f, d) {
    load(f)  
    climate_summaries_region %>%
      mutate(
        model = basename(d),
        scenario = scen_tag
      )
  })
  
  df_all <- bind_rows(lst)
  
  key_order <- df_all %>%
    distinct(geo, year_idx) %>%
    mutate(.ord = row_number())
  
  num_cols <- names(df_all)[sapply(df_all, is.numeric)]
  num_cols <- setdiff(num_cols, "year_idx")
  
  df_ens <- df_all %>%
    group_by(geo, year_idx) %>%
    summarise(
      LocationNa = first(LocationNa),
      across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(key_order, by = c("geo", "year_idx")) %>%
    arrange(.ord) %>%                 # restore df_all ordering
    select(-.ord) %>%
    mutate(
      year = 2015 + year_idx - 1,
      scenario = scen_tag
    )
  
  df_ens
}


# Overlap-mean stitcher (same logic as SSPs)
#    - temperature: additive
#    - precipitation: multiplicative ratio
#    - humidity: logit-delta
stitch_by_overlap_means <- function(df_obs, df_future_ens,
                                    overlap_years = 2015:2021,
                                    future_year_min = 2022,
                                    eps_pr = 1e-6,
                                    eps_rh = 1e-4) {
  
  covar_cols <- names(df_obs)[sapply(df_obs, is.numeric)]
  covar_cols <- setdiff(covar_cols, c("year_idx", "year"))
  
  logit   <- function(p) log(p / (1 - p))
  expit   <- function(z) 1 / (1 + exp(-z))
  clamp01 <- function(x, eps) pmin(pmax(x, eps), 1 - eps)
  
  rh_cols <- covar_cols[str_detect(covar_cols, regex("humidity", ignore_case = TRUE))]
  rh_scale <- setNames(rep(100, length(rh_cols)), rh_cols)
  if (length(rh_cols) > 0) {
    for (cc in rh_cols) {
      mx <- max(df_obs[[cc]], df_future_ens[[cc]], na.rm = TRUE)
      if (is.finite(mx) && mx <= 1.5) rh_scale[[cc]] <- 1
    }
  }
  
  # observational baseline means over overlap years ----
  obs_base <- df_obs %>%
    filter(year %in% overlap_years) %>%
    group_by(geo) %>%
    summarise(
      across(all_of(covar_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(.x, "_obsBase"), all_of(covar_cols))
  
  # model baseline means over overlap years ----
  fut_base <- df_future_ens %>%
    filter(year %in% overlap_years) %>%
    group_by(geo) %>%
    summarise(
      across(all_of(covar_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(.x, "_futBase"), all_of(covar_cols))
  
  # stitch future years onto observational baseline ----
  df_fut <- df_future_ens %>%
    filter(year >= future_year_min) %>%
    left_join(obs_base, by = "geo") %>%
    left_join(fut_base, by = "geo")
  
  for (cc in covar_cols) {
    
    obs_cc <- paste0(cc, "_obsBase")
    ref_cc <- paste0(cc, "_futBase")
    
    if (str_detect(cc, regex("temperature", ignore_case = TRUE))) {
      
      df_fut[[cc]] <- df_fut[[obs_cc]] + (df_fut[[cc]] - df_fut[[ref_cc]])
      
    } else if (str_detect(cc, regex("precipitation", ignore_case = TRUE))) {
      
      denom <- pmax(df_fut[[ref_cc]], eps_pr)
      df_fut[[cc]] <- df_fut[[obs_cc]] * (df_fut[[cc]] / denom)
      
    } else if (str_detect(cc, regex("humidity", ignore_case = TRUE))) {
      
      sc <- if (!is.null(rh_scale[[cc]])) rh_scale[[cc]] else 100
      
      h_obs <- clamp01(df_fut[[obs_cc]] / sc, eps_rh)
      h_ref <- clamp01(df_fut[[ref_cc]] / sc, eps_rh)
      h_mod <- clamp01(df_fut[[cc]]     / sc, eps_rh)
      
      z_obs <- logit(h_obs)
      z_ref <- logit(h_ref)
      z_mod <- logit(h_mod)
      
      z_corr <- z_obs + (z_mod - z_ref)
      df_fut[[cc]] <- sc * expit(z_corr)
      
    } else {
      
      df_fut[[cc]] <- df_fut[[cc]]
    }
  }
  
  df_out <- df_fut %>%
    select(-ends_with("_obsBase"), -ends_with("_futBase")) %>%
    mutate(year_idx = year - future_year_min + 1L) %>%
    select(geo, year_idx, year, LocationNa, all_of(covar_cols), scenario)
  
  df_out
}


overlap_years <- 2015:2021
scenario_tag  <- "Pi_control"      ##  pi control

# ensemble mean across models (2015..2100)
df_pi_ens <- load_pi_control_ensemble(pi_folder = pi_folder, scen_tag = scenario_tag)

ens_outfile <- file.path(out_root, "EnsembleMean_Pi_control_climate_data_yearly_2015_2100.RData")
climate_summaries_region_ensemble <- df_pi_ens
save(climate_summaries_region_ensemble, file = ens_outfile)

# stitch future years (2022..2100) using overlap means 2015..2021 ----
df_pi_stitched <- stitch_by_overlap_means(
  df_obs          = df_cf_obs,        # counterfactual observational baseline
  df_future_ens   = df_pi_ens,        # ensemble mean Pi_control
  overlap_years   = overlap_years,
  future_year_min = 2022
)

# save stitched output (2022..2100)
stitched_outfile <- file.path(out_root,
                              "Final_ensemble_Pi_control_delta_corrected_climate_data_yearly_2022_2100.RData")
climate_summaries_region_pi_control <- df_pi_stitched
save(climate_summaries_region_pi_control, file = stitched_outfile)




# Quick plot check 
df_cf_summer <- df_cf_obs %>%
  group_by(year) %>%
  summarise(
    temp_summer = mean(temperature_summer_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(scenario = "Counterfactual obs (20CRv3/ERA5)")

df_pi_summer <- df_pi_stitched %>%
  group_by(year) %>%
  summarise(
    temp_summer = mean(temperature_summer_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(scenario = "Pi_control ensemble (overlap-mean stitched)")


df_plot <- bind_rows(df_cf_summer, df_pi_summer)


p <- ggplot(df_plot, aes(x = year, y = temp_summer, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 2021, linetype = "dashed", linewidth = 0.4) +
  labs(
    x = "Year",
    y = "Mean summer temperature (averaged over NUTS3)",
    colour = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.height = grid::unit(0.55, "cm"),
    legend.key.width  = grid::unit(0.90, "cm"),
    axis.line  = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3)
  )

print(p)
