##==============================================================================
## Exploratory data analysis for the revised dataset (dc_max.csv)
##
## Produces every data table and figure for the Data section of
## the manuscript and for the referee responses. Each block is tagged:
##
##   [PAPER]    -- candidate for the manuscript Data section
##   [APPENDIX] -- candidate for an appendix / supplement
##   [LETTER]   -- numbers quoted in the response letter only
##
## Referee anchors (see chats/TODO.md, item 2):
##  - Referee 2, minor ("p. 5, top; Figure 2"): "How many daily
##    observations have to be present for the annual max to be considered
##    not missing?" -> blocks 6-7 (within-year coverage from the raw
##    record).
##  - Referee 2, exploratory figures: the old EDA figure was judged
##    uninformative; replaced by the study-area map (block 3), the
##    record-length figure (block 4) and the glyph map (block 5), which
##    give the maxima geographical context.
##  - Referee 2, busy figures: everything here is one message per figure.
##
## Inputs:
##  - data/dc_max.csv (the modelled dataset)
##  - data/dc_max_files/processed/02_combined_dataset/da.rda
##    (raw record; gitignored, local-only -- coverage blocks are skipped
##    with a message when it is absent)
##  - rds/pinned_geometry.rds (mesh figure only)
##
## Outputs: figures in figs/ (png + pdf), tables in rds/ with the
## explore_data_ prefix. The name map_study_area is pinned: the
## fig-study-area chunk in the manuscript includes it.
##==============================================================================

##==============================================================================
## Setup

here::i_am("00_explore_data.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(sf)
})

## The stExtremes package provides the area and barrier datasets
library(stExtremes)

figs_dir <- here::here("figs")
if(!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)
rds_dir <- here::here("rds")
if(!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)

theme_set(theme_bw(base_size = 14))

save_fig <- function(plot, name, w, h) {
    ggsave(sprintf("%s/%s.png", figs_dir, name), plot,
        width = w, height = h, dpi = 300, device = png)
    ggsave(sprintf("%s/%s.pdf", figs_dir, name), plot,
        width = w, height = h, dpi = 300, device = cairo_pdf)
}

## Neutral, print-friendly palette. Land and sea are greys, so the black
## series in the glyph map stay legible against the background. Country
## is never encoded by colour: it is already given by position on the map
## and by the facet strips everywhere else.
col_sea <- "white"
col_land <- "grey85"
col_coast <- "grey50"
col_ink <- "grey15"
col_bar <- "grey45"

##==============================================================================
## Data

proj <- st_crs("+proj=longlat +datum=WGS84")
kmproj <- st_crs("+proj=utm +zone=30 +datum=WGS84 +units=km")

dc_max <- read_csv(here::here("data/dc_max.csv"),
    show_col_types = FALSE)

## One row per site, for maps
sites <- dc_max |>
    distinct(country, site_name, lat, lon)
sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = proj)

data(area, package = "stExtremes")
data(barrier, package = "stExtremes")

year_min <- min(dc_max$year)
year_max <- max(dc_max$year)
n_years <- year_max - year_min + 1

##==============================================================================
## 1. Dataset summary [PAPER, Data section]
##
## One-row overview quoted in the text (and complementing tbl-data in
## the manuscript). Maxima are in METRES; the text must say so (TODO
## item 1).

tbl_overview <- dc_max |>
    summarise(
        sites = n_distinct(site_name),
        sites_irl = n_distinct(site_name[country == "IRL"]),
        sites_gbr = n_distinct(site_name[country == "GBR"]),
        years = n_years,
        cells = n(),
        observed = sum(!is.na(max)),
        observed_pct = 100 * observed / cells,
        max_min_m = min(max, na.rm = TRUE),
        max_max_m = max(max, na.rm = TRUE))
tbl_overview
saveRDS(tbl_overview, file.path(rds_dir, "explore_data_overview_tbl.rds"))

##==============================================================================
## 2. Per-site record summary [APPENDIX]
##
## 44 rows: observed years, first/last year, range of the maxima. Too
## long for the text; the paper quotes the extremes (58 years at Milford
## down to 3 at RoonaghPier).

tbl_sites <- dc_max |>
    filter(!is.na(max)) |>
    group_by(country, site_name) |>
    summarise(
        n_years = n(),
        first = min(year),
        last = max(year),
        max_min_m = min(max),
        max_max_m = max(max),
        .groups = "drop") |>
    arrange(country, desc(n_years))
tbl_sites
saveRDS(tbl_sites, file.path(rds_dir, "explore_data_sites_tbl.rds"))

##==============================================================================
## 3. Study-area map with the gauges [PAPER]
##
## The fig-study-area chunk in the manuscript waits for exactly this file
## name. Sea white, land barrier in grey, gauges as dark dots.
##
## This map is the reader's gazetteer, so it names things the other
## figures only imply: the two countries (as the glyph map does), all 44
## gauges, and the three exposure regions used as the spread covariate in
## the results.
##
## The exposure rule is the one in 01_explore_models.R and
## 05_explore_spread_shape_inla.R, and it is a pair of half-planes in
## lon/lat rather than a hydrographic boundary:
##
##   lon < -8.0                         -> Atlantic
##   lat < 52.0 & lon > -5.5 (else)     -> Bristol Channel
##   otherwise                          -> Irish Sea
##
## Because it is that simple, it draws exactly: one meridian at -8.0
## spanning the map, and an L joining the parallel at 52.0 to the
## meridian at -5.5. Dashed, so they cannot be mistaken for coastline.
## Keep the definition here in step with those two scripts; it is
## duplicated rather than shared because they build model matrices and
## this builds a picture.
##
## Gauge names are the raw site_name strings, NOT prettified, so they
## match @tbl-sites in the appendix exactly and a reader can move between
## the two. ggrepel places them: 44 labels on one panel cannot be
## hand-set, and the seed pins the result so the figure is stable across
## runs.

exposure_lines <- tribble(
    ~x,    ~xend, ~y,    ~yend,
    -8.0,  -8.0,  51.0,  55.98,   # Atlantic boundary, full height
    -5.5,  -2.38, 52.0,  52.0,    # Bristol Channel, northern edge
    -5.5,  -5.5,  51.0,  52.0)    # Bristol Channel, western edge

## Region names sit in open water, clear of every gauge and of the
## country labels. Checked against the rendered panel, not eyeballed from
## the coordinates.
exposure_lab <- tibble(
    lon = c(-10.4, -4.75, -4.50),
    lat = c(55.0, 54.6, 51.36),
    lab = c("Atlantic", "Irish Sea", "Bristol Channel"))

## Same two country labels as the glyph map, in the same open interior.
## IRL is nudged off the -8.0 meridian, which is an exposure boundary
## here and would otherwise run straight through the word.
country_lab <- tibble(
    lon = c(-7.5, -3.5),
    lat = c(53.2, 52.3),
    lab = c("IRL", "GBR"))

p_map <- ggplot() +
    geom_sf(data = area, fill = col_sea, colour = NA) +
    geom_sf(data = barrier, fill = col_land, colour = col_coast,
        linewidth = 0.2) +
    geom_segment(data = exposure_lines,
        aes(x = x, xend = xend, y = y, yend = yend),
        linetype = "dashed", colour = col_bar, linewidth = 0.4) +
    geom_text(data = exposure_lab, aes(x = lon, y = lat, label = lab),
        colour = col_bar, size = 3.6, fontface = "italic") +
    geom_text(data = country_lab, aes(x = lon, y = lat, label = lab),
        colour = col_coast, size = 5, fontface = "bold") +
    geom_sf(data = sites_sf, size = 1.4, colour = col_ink) +
    ggrepel::geom_text_repel(data = sites,
        aes(x = lon, y = lat, label = site_name),
        size = 2.3, colour = col_ink, segment.colour = col_coast,
        segment.size = 0.2, min.segment.length = 0,
        box.padding = 0.3, point.padding = 0.1, max.overlaps = Inf,
        seed = 20260730) +
    coord_sf(xlim = c(-11.2, -2.38), ylim = c(51.0, 55.98),
        expand = FALSE) +
    labs(x = "Longitude", y = "Latitude")
p_map
save_fig(p_map, "map_study_area", w = 7.5, h = 7.5)

##==============================================================================
## 4. Record length per site [PAPER]
##
## Replacement for the old Figure 2 (the barplot the referee's coverage
## question is anchored to). Number of OBSERVED years per site; the
## within-year coverage behind each bar is quantified in blocks 6-7.

p_length <- dc_max |>
    filter(!is.na(max)) |>
    count(country, site_name, name = "n_years") |>
    ggplot(aes(x = fct_reorder(site_name, n_years), y = n_years)) +
    geom_col(fill = col_bar) +
    coord_flip() +
    facet_wrap(~ country, scales = "free_y") +
    labs(x = NULL, y = "Observed years (annual maxima)")
p_length
save_fig(p_length, "data_record_length", w = 9, h = 6)

##==============================================================================
## 5. Annual maxima on the map, one glyph per gauge [PAPER]
##
## Based on chats/graph-overlap.R (time series overlapped on a map), but
## simplified: instead of one patchwork inset per site (unmanageable for
## 44 gauges), each series is drawn directly in map coordinates as a
## small glyph anchored at its gauge. All glyphs share BOTH scales, so
## levels and timing are comparable across sites:
##  - horizontally, every glyph spans 1968-2025 -> gaps line up;
##  - vertically, every glyph rises from 0 m at the gauge latitude, so
##    the HEIGHT above the baseline reads directly as surge magnitude.
## Anchoring at the baseline (not the glyph centre) is what makes each
## series visually belong to its dot. Glyph width is kept near the median
## gauge spacing (about 40 km, i.e. 0.6 deg of longitude here) so glyphs
## stay next to their own gauge instead of sprawling inland.
## Lines bridge missing years; the record-length figure and the coverage
## tile carry the missingness message, this figure carries level and
## timing.
##
## Overlap: 27 of the 44 gauges fall into 10 clusters close enough for
## their glyph boxes to collide (the largest is five gauges on the Irish
## east coast between Dundalk and Dublin Port). Drawn at their true
## positions the crowded ones hide each other, so the glyph anchors are
## separated by the short force-directed pass below and a leader line
## joins each glyph back to its gauge. The gauge dots stay at the true
## coordinates; only the glyphs move.

glyph_w <- 0.55   # glyph width, degrees longitude
glyph_h <- 0.35   # glyph height, degrees latitude
y_lo <- 0         # baseline is 0 m, so height = surge magnitude
y_hi <- max(dc_max$max, na.rm = TRUE)

## Push apart anchors whose boxes overlap. Each box spans [x - w/2,
## x + w/2] horizontally and [y, y + h] vertically. Overlapping pairs are
## separated along whichever axis needs the smaller move, and a weak
## spring pulls every anchor back towards its gauge so glyphs settle near
## home rather than drifting across the map.
##
## The spring must stay weak: it balances the repulsion at equilibrium,
## so a strong one (0.02 was tried) locks in residual overlaps of about
## 10% of a box. At 0.001 the crowded pairs settle tangent instead, the
## worst residual overlap being 0.3% of a box, i.e. invisible, and no
## glyph moves more than about 0.45 degrees from its gauge.
separate_anchors <- function(x, y, w, h, n_iter = 3000, push = 0.5,
    spring = 0.001) {
    x0 <- x
    y0 <- y
    n <- length(x)
    for(it in seq_len(n_iter)) {
        hit <- FALSE
        for(i in seq_len(n - 1)) {
            for(j in seq((i + 1), n)) {
                dx <- x[j] - x[i]
                dy <- y[j] - y[i]
                ox <- w - abs(dx)
                oy <- h - abs(dy)
                if(ox <= 0 || oy <= 0) next
                hit <- TRUE
                if(ox / w < oy / h) {
                    s <- push * ox * if(dx == 0) 1 else sign(dx)
                    x[i] <- x[i] - s / 2
                    x[j] <- x[j] + s / 2
                } else {
                    s <- push * oy * if(dy == 0) 1 else sign(dy)
                    y[i] <- y[i] - s / 2
                    y[j] <- y[j] + s / 2
                }
            }
        }
        x <- x + spring * (x0 - x)
        y <- y + spring * (y0 - y)
        if(!hit) break
    }
    tibble(ax = x, ay = y)
}

## Row order is irrelevant to the separation, so bind straight onto
## sites. Do NOT reorder one side only: dplyr::arrange() sorts in the C
## locale and base order() in the system locale, and the two disagree on
## names such as PortEllen and Portbury, which silently pairs each site
## with another site's anchor.
anchors <- bind_cols(sites,
    separate_anchors(x = sites$lon, y = sites$lat,
        w = glyph_w, h = glyph_h))

## How far the separation had to move things, in km-equivalent degrees
anchors |>
    mutate(shift = sqrt((ax - lon)^2 + (ay - lat)^2)) |>
    summarise(moved = sum(shift > 1e-6), max_shift_deg = max(shift))

glyphs <- dc_max |>
    filter(!is.na(max)) |>
    left_join(select(anchors, site_name, ax, ay), by = "site_name") |>
    mutate(
        gx = ax - glyph_w / 2 +
            glyph_w * (year - year_min) / (year_max - year_min),
        gy = ay + glyph_h * (max - y_lo) / (y_hi - y_lo))

## Scale key, drawn in the empty Atlantic corner
key_x <- -10.8
key_y <- 51.15

## Reference labels, so the reader can orient without a separate map.
##
## All positions below are HAND-SET. The glyphs occupy most of the sea and
## reach inland from the coast-facing gauges, so these were chosen to sit
## in open interior: the Irish midlands and mid-Wales. If glyph_w or
## glyph_h is ever changed, a series may move onto one, in which case
## nudge the coordinates here rather than shrinking the glyphs.
region_lab <- tibble(
    lon = c(-8.0, -3.5),
    lat = c(53.0, 52.3),
    lab = c("IRL", "GBR"))

## Four ports give a geographical anchor: one on each side of the Irish
## Sea, one on the Atlantic coast, and the Bristol Channel, which is
## where the largest maxima are. The label sits away from its gauge in
## open ground and is joined to it by a dotted leader. Dotted
## deliberately: the solid grey leaders in this figure already mean "this
## glyph was displaced from its gauge", and reusing them here would
## conflate two different things.
## The label positions are not eyeballed. They were found by searching a
## 0.05-degree grid for the point nearest each gauge whose label box
## clears every separated glyph box, the two region labels and the scale
## key, using the label extents implied by size 3.6 on this panel. Redo
## that search if glyph_w, glyph_h or the label size changes; two of the
## first hand-picked positions turned out to sit on the very series they
## were naming.
port_lab <- tribble(
    ~site_name,   ~lab,        ~lx,   ~ly,
    "DublinPort", "Dublin",   -6.75, 53.70,
    "GalwayPort", "Galway",   -8.80, 53.00,
    "Liverpool",  "Liverpool", -3.00, 53.10,
    "Avonmouth",  "Bristol",   -2.40, 51.35) |>
    left_join(select(sites, site_name, lon, lat), by = "site_name")
stopifnot(!anyNA(port_lab$lon))

p_glyph <- ggplot() +
    geom_sf(data = area, fill = col_sea, colour = NA) +
    geom_sf(data = barrier, fill = col_land, colour = col_coast,
        linewidth = 0.2) +
    ## Leader line from the true gauge position to its displaced glyph
    geom_segment(data = anchors,
        aes(x = lon, y = lat, xend = ax, yend = ay),
        colour = "grey55", linewidth = 0.2) +
    ## Baseline (0 m) under each glyph, so levels are read off a common
    ## reference rather than guessed
    geom_segment(data = anchors,
        aes(x = ax - glyph_w / 2, xend = ax + glyph_w / 2,
            y = ay, yend = ay),
        colour = "grey55", linewidth = 0.25) +
    geom_line(data = glyphs,
        aes(x = gx, y = gy, group = site_name),
        colour = "black", linewidth = 0.45) +
    geom_point(data = sites, aes(x = lon, y = lat), size = 0.7,
        colour = col_ink) +
    ## Country labels in the open interior
    geom_text(data = region_lab, aes(x = lon, y = lat, label = lab),
        size = 5, fontface = "bold", colour = "grey35") +
    ## Port leaders and labels. The label is drawn on a translucent white
    ## patch so it stays readable if a series happens to pass behind it.
    ## linewidth = 0 removes the patch border: label.size did the same
    ## job but is deprecated from ggplot2 3.5.0.
    geom_segment(data = port_lab,
        aes(x = lon, y = lat, xend = lx, yend = ly),
        colour = "grey55", linewidth = 0.2, linetype = "dotted") +
    geom_label(data = port_lab, aes(x = lx, y = ly, label = lab),
        size = 3.6, colour = col_ink, fill = alpha("white", 0.75),
        linewidth = 0, label.padding = unit(0.12, "lines")) +
    ## Key: same geometry as a glyph, with the two scales labelled
    annotate("segment", x = key_x, xend = key_x + glyph_w,
        y = key_y, yend = key_y, colour = col_ink, linewidth = 0.3) +
    annotate("segment", x = key_x, xend = key_x,
        y = key_y, yend = key_y + glyph_h, colour = col_ink,
        linewidth = 0.3) +
    annotate("text", x = key_x, y = key_y - 0.12, label = year_min,
        size = 2.6, hjust = 0.5) +
    annotate("text", x = key_x + glyph_w, y = key_y - 0.12,
        label = year_max, size = 2.6, hjust = 0.5) +
    annotate("text", x = key_x - 0.08, y = key_y, label = "0",
        size = 2.6, hjust = 1) +
    annotate("text", x = key_x - 0.08, y = key_y + glyph_h,
        label = sprintf("%.1f m", y_hi), size = 2.6, hjust = 1) +
    labs(x = "Longitude", y = "Latitude")
p_glyph
save_fig(p_glyph, "data_glyph_map", w = 8, h = 8)

##==============================================================================
## 6. Site-by-year data presence [APPENDIX]
##
## The missingness pattern behind Figure 2, cell by cell. Sites ordered
## by latitude within country, so spatial neighbours are adjacent rows.
## Block 7e is the same figure recoloured by within-year coverage and is
## strictly more informative, so only one of the two belongs in the
## paper. This one is kept because it needs only dc_max.csv, while
## 7e needs the local-only raw record.

p_presence <- dc_max |>
    mutate(site_name = fct_reorder(site_name, lat)) |>
    ggplot(aes(x = year, y = site_name, fill = !is.na(max))) +
    geom_tile(colour = "grey85", linewidth = 0.1) +
    scale_fill_manual(name = NULL, values = c(`TRUE` = col_ink,
        `FALSE` = "white"), labels = c(`TRUE` = "observed",
        `FALSE` = "missing")) +
    facet_grid(country ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Year", y = NULL) +
    theme(axis.text.y = element_text(size = 7))
p_presence
save_fig(p_presence, "data_presence_tile", w = 9, h = 8)

##==============================================================================
## 7. Within-year coverage from the raw record [APPENDIX + LETTER]
##
## The direct answer to Referee 2's "daily observations" question. No
## minimum-coverage rule was applied upstream (see
## chats/chat_dc_max_provenance_and_coverage.md), so the honest answer
## is the coverage DISTRIBUTION: distinct calendar days with at least
## one record, per observed site-year. Day counts are comparable across
## the mixed sampling intervals (5, 6, 15, 60 min); raw record counts
## are not. Needs the local-only raw record; skipped when absent.

da_file <- here::here("data/dc_max_files/processed",
    "02_combined_dataset/da.rda")
cov_file <- file.path(rds_dir, "explore_data_cov_year.rds")

if(!file.exists(cov_file) && file.exists(da_file)) {
    message("## Computing within-year coverage from da.rda ...")
    load(da_file)  # provides da (~35M rows; takes a couple of minutes)
    cov_year <- da |>
        distinct(site, time) |>  # exact duplicate rows exist in 2023
        mutate(year = year(time), day = as_date(time)) |>
        group_by(site, year) |>
        summarise(
            n_rec = n(),
            n_days = n_distinct(day),
            n_months = n_distinct(month(time)),
            .groups = "drop")
    saveRDS(cov_year, cov_file)
    rm(da)
}

if(file.exists(cov_file)) {
    cov_year <- readRDS(cov_file)
    ## Coverage of the OBSERVED site-years in the modelled dataset
    obs_cov <- dc_max |>
        filter(!is.na(max)) |>
        left_join(cov_year, by = c("site_name" = "site", "year"))
    stopifnot(!anyNA(obs_cov$n_days))
    ## Objects inside this block need an explicit print(): R only
    ## auto-prints at the top level of a script, not inside braces.
    ## 7a. Quantiles of days-with-data [LETTER, one sentence in PAPER]
    tbl_cov_quant <- tibble(
        prob = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1),
        n_days = quantile(obs_cov$n_days, prob))
    print(tbl_cov_quant)
    ## 7b. Site-years below candidate thresholds [LETTER]
    tbl_cov_thr <- tibble(min_days = c(30, 60, 90, 120, 180, 274)) |>
        mutate(
            n_below = map_int(min_days, \(t) sum(obs_cov$n_days < t)),
            pct_below = 100 * n_below / nrow(obs_cov))
    print(tbl_cov_thr)
    ## 7c. Do sparse years bias the maxima low? [APPENDIX + LETTER]
    ## Yes: mean maximum rises monotonically with coverage.
    tbl_cov_bias <- obs_cov |>
        mutate(cov_bin = cut(n_days, breaks = c(0, 30, 90, 180, 300, 366),
            include.lowest = TRUE)) |>
        group_by(cov_bin) |>
        summarise(
            n = n(),
            mean_max_m = mean(max),
            median_max_m = median(max),
            .groups = "drop")
    print(tbl_cov_bias)
    saveRDS(list(quantiles = tbl_cov_quant, thresholds = tbl_cov_thr,
        bias = tbl_cov_bias),
        file.path(rds_dir, "explore_data_cov_tables.rds"))
    ## 7d. Coverage histogram [APPENDIX]
    p_cov <- ggplot(obs_cov, aes(x = n_days)) +
        geom_histogram(binwidth = 15, boundary = 0, fill = col_bar,
            colour = "white", linewidth = 0.2) +
        labs(x = "Days with data in the year", y = "Site-years")
    p_cov
    save_fig(p_cov, "data_coverage_hist", w = 7, h = 4.5)
    ## 7e. Presence tile recoloured by coverage [PAPER candidate]
    ## Stronger than block 6: shows WHERE the thin years are, not just
    ## that they exist. If this goes in the paper, block 6 is redundant
    ## and moves out.
    p_cov_tile <- dc_max |>
        left_join(cov_year, by = c("site_name" = "site", "year")) |>
        mutate(
            n_days = if_else(is.na(max), NA, n_days),
            site_name = fct_reorder(site_name, lat)) |>
        ggplot(aes(x = year, y = site_name, fill = n_days)) +
        geom_tile(colour = "grey92", linewidth = 0.1) +
        ## Sequential grey: light = thin year, dark = full year. Missing
        ## cells are white, separated from the lightest grey by the tile
        ## borders.
        scale_fill_gradient(name = "Days\nwith data", low = "grey80",
            high = "grey10", na.value = "white") +
        facet_grid(country ~ ., scales = "free_y", space = "free_y") +
        labs(x = "Year", y = NULL) +
        theme(axis.text.y = element_text(size = 7))
    p_cov_tile
    save_fig(p_cov_tile, "data_coverage_tile", w = 9, h = 8)
} else {
    message("## da.rda not found: coverage blocks 7a-7e skipped. ",
        "The raw record is local-only (gitignored); see ",
        "chats/chat_dc_max_provenance_and_coverage.md.")
}

##==============================================================================
## 8. Mesh over the study area [PAPER or APPENDIX]
##
## Speaks to Referee 2's mesh-quality point. Uses the PINNED geometry
## (never rebuild it inline; see 01_explore_models.R). Combined with
## the study-area map this closes the "data and mesh figures" half of
## TODO item 2.

geom_file <- file.path(rds_dir, "pinned_geometry.rds")

if(file.exists(geom_file) && requireNamespace("inlabru", quietly = TRUE)) {
    mesh <- readRDS(geom_file)$mesh
    area_km <- st_transform(area, kmproj)
    sites_km <- st_transform(sites_sf, kmproj)
    p_mesh <- ggplot() +
        inlabru::gg(mesh) +
        geom_sf(data = area_km, fill = NA, colour = "grey45",
            linewidth = 0.4) +
        geom_sf(data = sites_km, colour = col_ink, size = 1) +
        coord_sf(datum = kmproj) +
        labs(x = "Easting (km)", y = "Northing (km)")
    p_mesh
    save_fig(p_mesh, "map_mesh", w = 7, h = 7)
} else {
    message("## pinned_geometry.rds (or inlabru) not available: ",
        "mesh figure skipped.")
}

message("00_explore_data.R complete.")
