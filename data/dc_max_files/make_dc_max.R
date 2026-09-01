# -------------------------------------------------------------------------
# make_dc_max.R
#
# Helper for converting annual surge maxima tables into dc_max format.
# -------------------------------------------------------------------------

ireland_tide_gauge_sites <- function() {
  c(
    "AranmoreIsland-Leabgarrow",
    "ArklowHarbour",
    "BallycottonHarbour",
    "BallyglassHarbour",
    "CastletownberePort",
    "DingleHarbour",
    "DublinPort",
    "Dundalk",
    "DunmoreEastHarbour",
    "Fenit",
    "GalwayPort",
    "Howth",
    "Inishmore",
    "KillybegsPort",
    "KilrushLough",
    "MalinHead-PortmorePier",
    "PortOriel",
    "RingaskiddyNMCI",
    "RoonaghPier",
    "Rosslare",
    "SkerriesHarbour",
    "Sligo",
    "UnionHallHarbor",
    "WexfordHarbour"
  )
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  x[1]
}

make_dc_max <- function(annual_results,
                        max_col,
                        irl_sites = ireland_tide_gauge_sites()) {
  
  max_col <- rlang::ensym(max_col)
  
  year_min <- min(annual_results$year, na.rm = TRUE)
  year_max <- max(annual_results$year, na.rm = TRUE)
  
  site_lookup <- annual_results |>
    dplyr::group_by(site) |>
    dplyr::summarise(
      lat = first_non_missing(latitude),
      lon = first_non_missing(longitude),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      country = dplyr::if_else(site %in% irl_sites, "IRL", "GBR"),
      site_name = site
    ) |>
    dplyr::select(country, site_name, lat, lon)
  
  annual_results |>
    dplyr::transmute(
      site_name = site,
      year = as.integer(year),
      max = !!max_col
    ) |>
    dplyr::right_join(
      tidyr::expand_grid(
        site_name = sort(unique(annual_results$site)),
        year = year_min:year_max
      ),
      by = c("site_name", "year")
    ) |>
    dplyr::left_join(site_lookup, by = "site_name") |>
    dplyr::mutate(
      year_idx = year - year_min + 1
    ) |>
    dplyr::select(
      country,
      site_name,
      year,
      lat,
      lon,
      max,
      year_idx
    ) |>
    dplyr::arrange(country, site_name, year)
}
