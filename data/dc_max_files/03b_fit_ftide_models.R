# -------------------------------------------------------------------------
# 03_fit_ftide_models.R
#
# Fits TideHarmonics::ftide models to the combined tide-gauge dataset.
#
# Inputs:
#   processed/02_combined_dataset/da.rda
#
# Saves internal analysis objects to:
#   rds/03_ftide/
#
# Saves large residual datasets as partitioned parquet to:
#   processed/03_ftide/residuals_all_time/
##   processed/03_ftide/residuals_annual/      (annual branch not used here)
#
## Saves final human-facing tables to:      (not used here)
##   outputs/tables/03_ftide/
#
## Saves figures to:                        (not used here)
##   figs/03_ftide/
# -------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(tibble)
library(TideHarmonics)
library(arrow)
library(here)

source(here::here("data/dc_max_files/fit_ftide_helpers.R"))
## Not sent, and only used by the blocks commented out below
## source(here::here("data/dc_max_files/calculate_surge_thresholds.R"))
## source(here::here("data/dc_max_files/calculate_annual_surge_maxima.R"))
## source(here::here("data/dc_max_files/calculate_site_issue_flags.R"))

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

step_id <- "03_ftide"

inst_paper <- "data/dc_max_files"

input_rda <- here::here(inst_paper,
    "processed", "02_combined_dataset", "da.rda"
)

processed_ftide_dir <- here::here(inst_paper,
  "processed", step_id
)

residuals_all_time_dir <- file.path(
  processed_ftide_dir, "residuals_all_time"
)

## residuals_annual_dir <- file.path(
##   processed_ftide_dir, "residuals_annual"
## )

rds_dir <- here::here(inst_paper, "rds", step_id)

## outputs_tables_dir <- here::here(inst_paper,
##   "outputs", "tables", step_id
## )

## figs_dir <- here::here(inst_paper, "figs", step_id)

dir.create(processed_ftide_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(residuals_all_time_dir, recursive = TRUE, showWarnings = FALSE)
## dir.create(residuals_annual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
## dir.create(outputs_tables_dir, recursive = TRUE, showWarnings = FALSE)
## dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Settings
# -------------------------------------------------------------------------

min_obs_all_time <- 24 * 30
## min_obs_year     <- 24 * 30

use_smsl  <- TRUE
use_nodal <- TRUE

hcn_use <- TideHarmonics::hc37

## Only used by the blocks commented out below
## min_prediction_coverage_warn <- 0.80
## max_rmse_warn <- 0.30
## min_obs_year_for_annual_max_flag <- 24 * 180

# -------------------------------------------------------------------------
# Load combined dataset
# -------------------------------------------------------------------------

if (!file.exists(input_rda)) {
  stop("Could not find combined dataset: ", input_rda)
}

load(input_rda)

required_cols <- c("site", "time", "latitude", "longitude", "sea_level", "qc")
missing_cols <- setdiff(required_cols, names(da))

if (length(missing_cols) > 0) {
  stop(
    "The dataset is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

da <- da |>
  mutate(
    site = as.character(site),
    time = as.POSIXct(time, tz = "UTC")
  ) |>
  arrange(site, time) |>
  distinct(site, time, .keep_all = TRUE)

raw_data_by_year <- da |>
  mutate(year = lubridate::year(time)) |>
  group_by(site, year) |>
  summarise(
    n_rows = n(),
    n_nonmissing = sum(!is.na(sea_level)),
    n_missing = sum(is.na(sea_level)),
    prop_missing = n_missing / n_rows,
    .groups = "drop"
  )

saveRDS(
  raw_data_by_year,
  file.path(rds_dir, "raw_data_by_year.rds")
)

# -------------------------------------------------------------------------
# Fit one model per site using the full available time series
# -------------------------------------------------------------------------

site_names <- sort(unique(da$site))

all_time_summary_list <- vector("list", length(site_names))
all_time_skipped_list <- vector("list", length(site_names))
all_time_residuals_list <- vector("list", length(site_names))

for (i in seq_along(site_names)) {
  site_nm <- site_names[i]
  message("Fitting all-time model: ", site_nm)
  df_site <- da |>
    filter(site == site_nm)
  res <- fit_ftide_one_group(
    df = df_site,
    min_obs_required = min_obs_all_time,
    hcn_use = hcn_use,
    use_nodal = use_nodal,
    use_smsl = use_smsl
  )
  if (isTRUE(res$ok)) {
    all_time_residuals_list[[i]] <- res$residuals |>
      mutate(
        site = site_nm,
        model_type = "site_all_time"
      ) |>
      select(
        site, model_type, time, latitude, longitude, sea_level, qc,
        tide_pred, residual,
        skew_surge_date, daily_max_sea_level, daily_max_tide_pred,
        surge_estimate
      )
    all_time_summary_list[[i]] <- res$summary |>
      mutate(
        site = site_nm,
        model_type = "site_all_time"
      )
    all_time_skipped_list[[i]] <- NULL
  } else {
    all_time_summary_list[[i]] <- NULL
    all_time_skipped_list[[i]] <- tibble(
      site = site_nm,
      model_type = "site_all_time",
      reason = res$reason
    )
  }
  rm(df_site, res)
  gc()
}

all_time_summary <- bind_rows(all_time_summary_list)
all_time_skipped <- bind_rows(all_time_skipped_list)
all_time_residuals_combined <- bind_rows(all_time_residuals_list)

saveRDS(
  all_time_summary,
  file.path(rds_dir, "all_time_summary.rds")
)

saveRDS(
  all_time_skipped,
  file.path(rds_dir, "all_time_skipped.rds")
)


arrow::write_dataset(
  all_time_residuals_combined,
  path = residuals_all_time_dir,
  format = "parquet",
  partitioning = c("site"),
  existing_data_behavior = "overwrite"
)

##==============================================================================
## # -------------------------------------------------------------------------
## # Fit one model per site per calendar year
## # -------------------------------------------------------------------------

## da_year <- da |>
##   mutate(year = lubridate::year(time))

## site_year_keys <- da_year |>
##   distinct(site, year) |>
##   arrange(site, year)

## year_summary_list <- vector("list", nrow(site_year_keys))
## year_skipped_list <- vector("list", nrow(site_year_keys))
## annual_residuals_list <- vector("list", nrow(site_year_keys))

## for (i in seq_len(nrow(site_year_keys))) {
##   site_nm <- site_year_keys$site[i]
##   yr      <- site_year_keys$year[i]
##   message("Fitting annual model: ", site_nm, " / ", yr)
##   df_group <- da_year |>
##     filter(site == site_nm, year == yr)
##   res <- fit_ftide_one_group(
##     df = df_group,
##     min_obs_required = min_obs_year,
##     hcn_use = hcn_use,
##     use_nodal = use_nodal,
##     use_smsl = use_smsl
##   )
##   if (isTRUE(res$ok)) {
##     annual_residuals_list[[i]] <- res$residuals |>
##       mutate(
##         site = site_nm,
##         year = yr,
##         model_type = "site_year"
##       ) |>
##       select(
##         site, year, model_type, time, latitude, longitude, sea_level, qc,
##         tide_pred, residual,
##         skew_surge_date, daily_max_sea_level, daily_max_tide_pred,
##         surge_estimate
##       )
##     year_summary_list[[i]] <- res$summary |>
##       mutate(
##         site = site_nm,
##         year = yr,
##         model_type = "site_year"
##       )
##     year_skipped_list[[i]] <- NULL
##   } else {
##     year_summary_list[[i]] <- NULL
##     year_skipped_list[[i]] <- tibble(
##       site = site_nm,
##       year = yr,
##       model_type = "site_year",
##       reason = res$reason
##     )
##   }
##   rm(df_group, res)
##   gc()
## }

## annual_summary <- bind_rows(year_summary_list)
## annual_skipped <- bind_rows(year_skipped_list)
## annual_residuals_combined <- bind_rows(annual_residuals_list)

## saveRDS(
##   annual_summary,
##   file.path(rds_dir, "annual_summary.rds")
## )

## saveRDS(
##   annual_skipped,
##   file.path(rds_dir, "annual_skipped.rds")
## )


## arrow::write_dataset(
##   annual_residuals_combined,
##   path = residuals_annual_dir,
##   format = "parquet",
##   partitioning = c("site", "year"),
##   existing_data_behavior = "overwrite"
## )

# -------------------------------------------------------------------------
# Calculate surge thresholds and annual maxima
# -------------------------------------------------------------------------

## annual_site_surge_thresholds <- calculate_surge_thresholds(
##   annual_residuals_combined
## )

## all_time_site_surge_thresholds <- calculate_surge_thresholds(
##   all_time_residuals_combined
## )

## saveRDS(
##   annual_site_surge_thresholds,
##   file.path(rds_dir, "site_surge_thresholds_from_annual_residuals.rds")
## )

## saveRDS(
##   all_time_site_surge_thresholds,
##   file.path(rds_dir, "site_surge_thresholds_from_all_time_residuals.rds")
## )

## annual_surge_max_from_annual_residuals <- calculate_annual_surge_maxima(
##   residuals_df = annual_residuals_combined,
##   raw_data_by_year = raw_data_by_year,
##   site_surge_thresholds = annual_site_surge_thresholds,
##   min_obs_year_for_annual_max_flag = min_obs_year_for_annual_max_flag,
##   min_prediction_coverage_warn = min_prediction_coverage_warn
## )

##==============================================================================
## THIS WAS PROVIDED
## annual_surge_max_from_all_time_residuals <- calculate_annual_surge_maxima(
##   residuals_df = all_time_residuals_combined,
##   raw_data_by_year = raw_data_by_year,
##   site_surge_thresholds = all_time_site_surge_thresholds,
##   min_obs_year_for_annual_max_flag = min_obs_year_for_annual_max_flag,
##   min_prediction_coverage_warn = min_prediction_coverage_warn
## )

## saveRDS(
##   annual_surge_max_from_all_time_residuals,
##   file.path(rds_dir, "annual_surge_max_from_all_time_residuals.rds")
## )
##==============================================================================

## saveRDS(
##   annual_surge_max_from_annual_residuals,
##   file.path(rds_dir, "annual_surge_max_from_annual_residuals.rds")
## )


# -------------------------------------------------------------------------
# Calculate final site issue flags
# -------------------------------------------------------------------------

## site_issue_flags <- calculate_site_issue_flags(
##   all_time_summary = all_time_summary,
##   annual_surge_max = annual_surge_max_from_all_time_residuals,
##   min_prediction_coverage_warn = min_prediction_coverage_warn,
##   max_rmse_warn = max_rmse_warn
## )

## saveRDS(
##   site_issue_flags,
##   file.path(rds_dir, "site_issue_flags.rds")
## )

# -------------------------------------------------------------------------
# Save final human-facing output tables only
# -------------------------------------------------------------------------

## readr::write_csv(
##   site_issue_flags,
##   file.path(outputs_tables_dir, "site_issue_flags.csv")
## )

## readr::write_csv(
##   annual_surge_max_from_all_time_residuals,
##   file.path(outputs_tables_dir, "annual_surge_max_from_all_time_residuals.csv")
## )

## readr::write_csv(
##   annual_surge_max_from_annual_residuals,
##   file.path(outputs_tables_dir, "annual_surge_max_from_annual_residuals.csv")
## )

## readr::write_csv(
##   all_time_site_surge_thresholds,
##   file.path(outputs_tables_dir, "site_surge_thresholds_from_all_time_residuals.csv")
## )

## readr::write_csv(
##   annual_site_surge_thresholds,
##   file.path(outputs_tables_dir, "site_surge_thresholds_from_annual_residuals.csv")
## )

message("03_fit_ftide_models.R complete.")
