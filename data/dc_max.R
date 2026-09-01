##==============================================================================
## Checks for dc_max
here::i_am("data/dc_max.R")
library(tidyverse)
library(stExtremes)

load(here::here("data/dc_max_all_time_residual.rda"))
## Avoid confusion
dc_max_new <- dc_max
str(dc_max_new)
dc_max_new

data(area, package = "stExtremes") # sea domain (study area polygon)
data(barrier, package = "stExtremes") # land (the physical barrier)

## Re-attach lon/lat columns from the geometry (kept after transforms)
## Common lat/lon projection
proj <- st_crs("+proj=longlat +datum=WGS84")
dc_max_new <- st_as_sf(dc_max_new, coords = c("lon", "lat"), crs = proj)
add_coords <- function(x, names = c("lon", "lat")) {
    st_coordinates(x) |>
        as_tibble() |>
        setNames(names) |>
        bind_cols(x) |>
        st_as_sf()
}
dc_max_new <- add_coords(dc_max_new)

## Map
ggplot() +
    geom_sf(data = area, fill = "lightblue") +
    geom_sf(data = barrier, fill = "grey90") +
    geom_sf(data = dc_max_new, col = 2, size = 4) +
    labs(x = "Longitude", y = "Latitude")

dc_max_new |>
    filter(lon > -2)
## Potsmouth should be removed

## Sample size
dc_max_new |>
    filter(!is.na(max)) |>
    group_by(country, site_name) |>
    summarise(n_years = n_distinct(year), .groups = "drop") |>
    ggplot(aes(x = forcats::fct_reorder(site_name, n_years, .fun = max),
        y = n_years)) +
    geom_col() +
    labs(x = "", y = "Number of years") +
    coord_flip() +
    facet_wrap(~ country, scales = "free_y")

## Time series
dc_max_new |>
    filter(!is.na(max)) |>
    ggplot(aes(x = year, y = max, color = site_name)) +
    geom_line() +
    geom_point() +
    labs(x = "Year", y = "Surge Maxima") +
    theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 6)) +
    guides(color = guide_legend(nrow = 3))

## Check number of observations per year
with(subset(dc_max_new, !is.na(max)), table(year))
barplot(with(subset(dc_max_new, !is.na(max)), table(year)))
abline(h = 5, col = 2, lty = 2, lwd = 2)
## NOTE: the numer of observations is too small before 1961 (wich the previous
## dataset used). As the number of years forward have increased, I'll adopt "at
## least 5 data points per year" as a criterion. This way, the new data should
## start in 1968 onwards. Also, 2026 is incomplete data (it's a block year), so
## I'll stop in 2025. So
2025 - 1968

##==============================================================================
## Create this new dataset

## Load again to be sure
load(here::here("data/dc_max_all_time_residual.rda"))
## Avoid confusion
dc_max_new <- dc_max
str(dc_max_new)
dc_max_new

## Filter Portsmouth
dc_max_new <- dc_max_new |>
    filter(site_name != "Portsmouth")
str(dc_max_new)

## Filter by years
dc_max_new <- dc_max_new |>
    filter(year >= 1968 & year <= 2025)
str(dc_max_new)

## Fix year_idx
dc_max_new <- dc_max_new |>
    arrange(country, site_name, year) |>
    mutate(year_idx = year - min(year) + 1)
str(dc_max_new)

with(dc_max_new, table(year))
with(dc_max_new, table(year_idx))

##------------------------------------------------------------------------------
## Check again
dc_max_new_geo <- st_as_sf(dc_max_new, coords = c("lon", "lat"), crs = proj)
add_coords <- function(x, names = c("lon", "lat")) {
    st_coordinates(x) |>
        as_tibble() |>
        setNames(names) |>
        bind_cols(x) |>
        st_as_sf()
}
dc_max_new_geo <- add_coords(dc_max_new_geo)

## Map
ggplot() +
    geom_sf(data = area, fill = "lightblue") +
    geom_sf(data = barrier, fill = "grey90") +
    geom_sf(data = dc_max_new_geo, col = 2, size = 4) +
    labs(x = "Longitude", y = "Latitude")

## Sample size
dc_max_new |>
    filter(!is.na(max)) |>
    group_by(country, site_name) |>
    summarise(n_years = n_distinct(year), .groups = "drop") |>
    ggplot(aes(x = forcats::fct_reorder(site_name, n_years, .fun = max),
        y = n_years)) +
    geom_col() +
    labs(x = "", y = "Number of years") +
    coord_flip() +
    facet_wrap(~ country, scales = "free_y")

## Time series
dc_max_new |>
    filter(!is.na(max)) |>
    ggplot(aes(x = year, y = max, color = site_name)) +
    geom_line() +
    geom_point() +
    labs(x = "Year", y = "Surge Maxima") +
    theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 6)) +
    guides(color = guide_legend(nrow = 3))

## Check number of observations per year
with(subset(dc_max_new, !is.na(max)), table(year))
barplot(with(subset(dc_max_new, !is.na(max)), table(year)))
abline(h = 5, col = 2, lty = 2, lwd = 2)

## Export
write.table(dc_max_new, file = here::here("data/dc_max.csv"),
    sep = ",", dec = ".", row.names = FALSE)
