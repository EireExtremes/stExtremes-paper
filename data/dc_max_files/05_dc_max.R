# -------------------------------------------------------------------------
# 05_dc_max.R
#
# Converts annual surge maxima summaries into dc_max format for downstream
# spatial extreme-value and dependence modelling.
#
# The resulting dc_max tables contain one annual maximum value per
# site-year combination, formatted for later stExtremes / INLA workflows.
#
# Inputs:
#   rds/03_ftide/annual_surge_max_from_all_time_residuals.rds
##   rds/03_ftide/annual_surge_max_from_annual_residuals.rds  (not used here)
#
# Saves internal analysis objects to:
#   rds/05_dc_max/
#
## Saves final human-facing tables to:      (not used here)
##   outputs/tables/05_dc_max/
# -------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(here)

inst_paper <- "data/dc_max_files"

source(here::here(inst_paper, "make_dc_max.R"))

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

step_id <- "05_dc_max"

rds_ftide_dir <- here::here(inst_paper, "rds", "03_ftide")

rds_dc_max_dir <- here::here(inst_paper,
  "rds", step_id
)

## outputs_dc_max_dir <- here::here(inst_paper,
##   "outputs", "tables", step_id
## )

dir.create(rds_dc_max_dir, recursive = TRUE, showWarnings = FALSE)
## dir.create(outputs_dc_max_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Import annual maxima objects
# -------------------------------------------------------------------------

annual_surge_max_from_all_time_residuals <- readRDS(
  file.path(
    rds_ftide_dir,
    "annual_surge_max_from_all_time_residuals.rds"
  )
)

## annual_surge_max_from_annual_residuals <- readRDS(
##   file.path(
##     rds_ftide_dir,
##     "annual_surge_max_from_annual_residuals.rds"
##   )
## )

# -------------------------------------------------------------------------
# Create dc_max tables
# -------------------------------------------------------------------------

dc_max_all_time_residual <- make_dc_max(
  annual_results = annual_surge_max_from_all_time_residuals,
  max_col = max_residual_surge
)

dc_max_all_time_skew <- make_dc_max(
  annual_results = annual_surge_max_from_all_time_residuals,
  max_col = max_skew_surge
)

## dc_max_annual_residual <- make_dc_max(
##   annual_results = annual_surge_max_from_annual_residuals,
##   max_col = max_residual_surge
## )

## dc_max_annual_skew <- make_dc_max(
##   annual_results = annual_surge_max_from_annual_residuals,
##   max_col = max_skew_surge
## )

# -------------------------------------------------------------------------
# Save internal R objects
# -------------------------------------------------------------------------

saveRDS(
    dc_max_all_time_residual,
    file.path(
        rds_dc_max_dir,
        "dc_max_all_time_residual.rds"
    )
)

saveRDS(
    dc_max_all_time_skew,
    file.path(
        rds_dc_max_dir,
        "dc_max_all_time_skew.rds"
    )
)

## saveRDS(
##   dc_max_annual_residual,
##   file.path(
##     rds_dc_max_dir,
##     "dc_max_annual_residual.rds"
##   )
## )

## saveRDS(
##   dc_max_annual_skew,
##   file.path(
##     rds_dc_max_dir,
##     "dc_max_annual_skew.rds"
##   )
## )

# -------------------------------------------------------------------------
# Save human-readable tables
# -------------------------------------------------------------------------

## write_csv(
##   dc_max_all_time_residual,
##   file.path(
##     outputs_dc_max_dir,
##     "dc_max_all_time_residual.csv"
##   )
## )

## write_csv(
##   dc_max_all_time_skew,
##   file.path(
##     outputs_dc_max_dir,
##     "dc_max_all_time_skew.csv"
##   )
## )

## write_csv(
##   dc_max_annual_residual,
##   file.path(
##     outputs_dc_max_dir,
##     "dc_max_annual_residual.csv"
##   )
## )

## write_csv(
##   dc_max_annual_skew,
##   file.path(
##     outputs_dc_max_dir,
##     "dc_max_annual_skew.csv"
##   )
## )

message("05_dc_max.R complete.")
