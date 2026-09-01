##==============================================================================
## Barrier vs non-barrier comparison: spatial fields and return levels.
##
## Consumes the cached spatio-temporal fits from 01_explore_models.R
## (the barrier model m0st and its non-barrier twin m0st_nb) and the pinned
## geometry, then EXTRACTS and PLOTS:
##   1. the posterior-mean SPDE field u(s, t) on a grid, for selected years;
##   2. bGEV return levels at the gauge sites (curves) and on the grid (maps);
##   3. barrier vs non-barrier difference maps (stationary minus barrier).
##
## The fits are read straight from the cached objects; only posterior summaries
## (summary.fixed, summary.random, summary.fitted.values and the hyperparameter
## marginals) are used, so the pruned fits from 01 work here unchanged.
##
## Extraction and plotting are kept separate: every quantity is first written
## to a tidy object in rds, so figures can be redrawn without refitting or
## re-projecting. Nothing here refits a model or rebuilds the mesh; it reads
## the same pinned geometry used for fitting, so the field lines up with the
## fitted values (a sanity check below confirms this).
##==============================================================================

here::i_am("02_compare_barrier.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(sf)
    library(INLA)
})

## The stExtremes package provides create_indices(), marginals_st_gev() and
## return_level_bgev2() (the new-parametrisation return level: it takes the
## INLA bGEV outputs q_alpha and s_beta and converts internally).
library(stExtremes)

theme_set(theme_bw(base_size = 14))

rds_dir <- here::here("rds")
figs_dir <- here::here("figs")
if(!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

## Common projections (match the fitting script)
proj <- st_crs("+proj=longlat +datum=WGS84")
kmproj <- st_crs("+proj=utm +zone=30 +datum=WGS84 +units=km")

## Cache helper: reuse a saved result if present, otherwise compute and save.
## R evaluates `expr` lazily, so a cached run skips the projection and
## return-level work entirely. The cache is keyed by file name only: delete
## the relevant file to force a fresh extraction if you change n_pix,
## sel_years or map_periods below.
cache_rds <- function(file, expr) {
    if(file.exists(file)) return(readRDS(file))
    obj <- expr
    saveRDS(obj, file)
    obj
}

##==============================================================================
## Geometry, data and polygons
##==============================================================================
geom <- readRDS(file.path(rds_dir, "pinned_geometry.rds"))
mesh <- geom$mesh
n_data <- length(geom$y)
n_years <- length(unique(geom$year_idx))
stopifnot(length(geom$year_idx) == n_data)

## Data, in the SAME row order used to build the geometry, so fitted values
## and the stacked field line up with the sites/years below.
dc_max <- read_csv(here::here("data", "dc_max.csv"),
    show_col_types = FALSE)
stopifnot(isTRUE(all.equal(dc_max$max, geom$y)))

dc_sf <- st_as_sf(dc_max, coords = c("lon", "lat"), crs = proj) |>
    st_transform(kmproj)

## Domain (for masking the grid) and land (drawn grey, masks the field).
data(area, package = "stExtremes")
data(barrier, package = "stExtremes")
area_km <- st_union(st_transform(area, kmproj))
barrier_sf <- st_sf(geometry = st_union(st_transform(barrier, kmproj)))

## Gauge sites (unique), in km, for overlay on the maps.
sites_km <- dc_sf |>
    mutate(x = st_coordinates(geometry)[, 1],
        y = st_coordinates(geometry)[, 2]) |>
    st_drop_geometry() |>
    distinct(country, site_name, x, y)

##==============================================================================
## Load the two spatio-temporal fits
##==============================================================================
mv <- sprintf("mesh_%d", mesh$n)
st_models <- c(m0st = "barrier", m0st_nb = "stationary")
fit_files <- file.path(rds_dir,
    sprintf("explore_dc_max_bGEV_%s_%s.rds", names(st_models), mv))
stopifnot(all(file.exists(fit_files)))
fits <- map(fit_files, readRDS)
names(fits) <- names(st_models)

## Readable facet labels (ASCII only).
mod_labs <- c(m0st = "barrier, st", m0st_nb = "stationary, st")

##==============================================================================
## Per-model bGEV parameters and fixed effects
##==============================================================================
## bGEV mixing probabilities: must match control.bgev$q.mix in the fit.
p_a <- 0.05
p_b <- 0.2

## Posterior-mean spread and tail (natural scale, via marginals_st_gev) and the
## fixed-effect intercept of the bGEV location q_alpha(s, t).
model_pars <- imap(fits, \(fit, nm) {
    mg <- marginals_st_gev(fit)
    spread <- inla.zmarginal(as.matrix(mg$spread), silent = TRUE)$mean
    tail <- inla.zmarginal(as.matrix(mg$tail), silent = TRUE)$mean
    fx <- fit$summary.fixed
    tibble(model = nm, type = unname(st_models[[nm]]),
        beta_m = fx["m", "mean"], spread = spread, tail = tail)
}) |>
    list_rbind()
model_pars

##==============================================================================
## Prediction grid over the domain and the mesh -> grid projector
##==============================================================================
## The grid is restricted to the data-informed coastal zone: sea cells (inside
## area_km) that also lie within mask_km of a gauge. Cells farther offshore --
## the open NW Atlantic in particular -- carry no information, so the latent
## field reverts to its zero-mean prior there and the maps show spurious
## low/negative return levels. Masking to the gauge neighbourhood removes those
## artefacts from every grid map below. mask_km is the only knob: lower it for
## a tighter coastal band, raise it to show more. (Distance to the coastline
## via barrier_sf is an alternative criterion, but distance to a gauge better
## reflects where the data actually constrain the field.)
mask_km <- 75
bb <- st_bbox(area_km)
n_pix <- 140
gx <- seq(bb[["xmin"]], bb[["xmax"]], length.out = n_pix)
gy <- seq(bb[["ymin"]], bb[["ymax"]], length.out = n_pix)
grid_sf <- expand_grid(x = gx, y = gy) |>
    st_as_sf(coords = c("x", "y"), crs = kmproj, remove = FALSE)
sites_sf <- st_as_sf(sites_km, coords = c("x", "y"), crs = kmproj)
gauge_zone <- st_union(st_buffer(sites_sf, mask_km))
inside <- lengths(st_intersects(grid_sf, area_km)) > 0 &
    lengths(st_intersects(grid_sf, gauge_zone)) > 0
grid_df <- st_drop_geometry(grid_sf)[inside, , drop = FALSE]
A_grid <- inla.spde.make.A(mesh, loc = as.matrix(grid_df[, c("x", "y")]))
stopifnot(ncol(A_grid) == mesh$n)
## Grid-cache tag: the resolution and the mask define the grid, so embedding
## them in the cached-extract filenames rebuilds those extracts automatically
## when either changes, instead of silently reloading a stale grid.
gv <- sprintf("%s_n%d_msk%d", mv, n_pix, mask_km)

##==============================================================================
## Field extraction helper (robust to INLA's internal field ordering)
##==============================================================================
## The grouped ST field is one long vector; create_indices() gives, for each
## entry, its spatial node (w) and its year group (w.group). field_year()
## returns the field for a given year group in node order (length mesh$n),
## without assuming how the entries are laid out.
w.index <- create_indices(geom$year_idx, mesh)
stopifnot(length(w.index$w) == mesh$n * n_years)

field_year <- function(w, g) {
    sel <- w.index$w.group == g
    out <- numeric(mesh$n)
    out[w.index$w[sel]] <- w[sel]
    out
}

## Year-group (year_idx) lookup and the snapshots to map.
year_lookup <- dc_max |>
    distinct(year, year_idx) |>
    arrange(year)
sel_years <- c(1970, 1985, 2000, 2015, 2025)
sel_idx <- year_lookup$year_idx[match(sel_years, year_lookup$year)]
stopifnot(!anyNA(sel_idx))
last_idx <- max(geom$year_idx)
last_year <- year_lookup$year[match(last_idx, year_lookup$year_idx)]

##==============================================================================
## Sanity check: rebuild q_alpha(s, t) at the gauges and match the fitted
## location. This validates the fixed-effect + field reconstruction used for
## the grid maps below. The estimation stack is est-only, so the fitted rows
## are 1:n_data in data order.
##==============================================================================
idat <- seq_len(n_data)
fit_check <- imap(fits, \(fit, nm) {
    w <- fit$summary.random$w$mean
    u_obs <- as.vector(geom$Ast %*% w)
    pr <- filter(model_pars, model == nm)
    qa <- pr$beta_m + u_obs
    fitted_loc <- fit$summary.fitted.values$mean[idat]
    tibble(model = nm,
        max_abs_diff = max(abs(qa - fitted_loc), na.rm = TRUE))
}) |>
    list_rbind()
fit_check

##==============================================================================
## 1. Spatial field u(s, t) on the grid, selected years, both models
##==============================================================================
field_grid <- cache_rds(
    file.path(rds_dir, sprintf("compare_field_grid_%s.rds", gv)),
    {
        imap(fits, \(fit, nm) {
            w <- fit$summary.random$w$mean
            stopifnot(length(w) == mesh$n * n_years)
            map(seq_along(sel_idx), \(k) {
                u <- as.vector(A_grid %*% field_year(w, sel_idx[k]))
                tibble(model = nm, year = sel_years[k],
                    x = grid_df$x, y = grid_df$y, u = u)
            }) |>
                list_rbind()
        }) |>
            list_rbind() |>
            mutate(model = factor(model, levels = names(st_models)))
    })

##==============================================================================
## 1b. Spatial-field posterior SD u(s, t) on the grid, selected years
##==============================================================================
## Companion to the posterior-mean field above: the marginal posterior SD of
## the latent field at each mesh node, projected with the same A_grid. This is
## a linear interpolation of the node-wise marginal SDs (the standard
## inla.mesh.project approach), not a full variance propagation, so it shows
## WHERE the field is well or poorly constrained -- tight near the gauges,
## widening offshore -- rather than an exact predictive SD. Referee 2 asked for
## the uncertainty to accompany the mean field.
field_sd_grid <- cache_rds(
    file.path(rds_dir, sprintf("compare_field_sd_grid_%s.rds", gv)),
    {
        imap(fits, \(fit, nm) {
            w_sd <- fit$summary.random$w$sd
            stopifnot(length(w_sd) == mesh$n * n_years)
            map(seq_along(sel_idx), \(k) {
                u_sd <- as.vector(A_grid %*% field_year(w_sd, sel_idx[k]))
                tibble(model = nm, year = sel_years[k],
                    x = grid_df$x, y = grid_df$y, u_sd = u_sd)
            }) |>
                list_rbind()
        }) |>
            list_rbind() |>
            mutate(model = factor(model, levels = names(st_models)))
    })

##==============================================================================
## 2. Return levels at the gauge sites, fitted location at last_year
##==============================================================================
periods <- 2:100
rl_sites <- cache_rds(
    file.path(rds_dir, sprintf("compare_rlevels_sites_%s.rds", mv)),
    {
        imap(fits, \(fit, nm) {
            pr <- filter(model_pars, model == nm)
            loc_all <- fit$summary.fitted.values$mean[idat]
            dc_max |>
                mutate(loc = loc_all) |>
                filter(year == last_year) |>
                select(country, site_name, loc) |>
                mutate(rl = map(loc, \(l)
                    tibble(period = periods,
                        rl = map_dbl(periods, \(tt)
                            return_level_bgev2(tt, l, pr$spread, pr$tail,
                                alpha = 0.5, beta = 0.5,
                                p_a = p_a, p_b = p_b))))) |>
                select(-loc) |>
                unnest(rl) |>
                mutate(model = nm, type = unname(st_models[[nm]]),
                    .before = 1)
        }) |>
            list_rbind() |>
            mutate(model = factor(model, levels = names(st_models)))
    })

##==============================================================================
## 3. Return-level maps on the grid, at last_year, selected return periods
##==============================================================================
map_periods <- c(10, 50, 100)
rl_grid <- cache_rds(
    file.path(rds_dir, sprintf("compare_rlevels_grid_%s.rds", gv)),
    {
        imap(fits, \(fit, nm) {
            pr <- filter(model_pars, model == nm)
            w <- fit$summary.random$w$mean
            u <- as.vector(A_grid %*% field_year(w, last_idx))
            qa <- pr$beta_m + u
            map(map_periods, \(tt) {
                tibble(model = nm, period = tt,
                    x = grid_df$x, y = grid_df$y,
                    rl = map_dbl(qa, \(l)
                        return_level_bgev2(tt, l, pr$spread, pr$tail,
                            alpha = 0.5, beta = 0.5,
                            p_a = p_a, p_b = p_b)))
            }) |>
                list_rbind()
        }) |>
            list_rbind() |>
            mutate(model = factor(model, levels = names(st_models)))
    })

##==============================================================================
## Differences (stationary minus barrier) at last_year
##==============================================================================
## Compare the models on the IDENTIFIED location surface q_alpha = beta_m + u,
## not the raw field u. beta_m and u are only jointly identified, and the two
## models split the overall level differently between them (barrier
## beta_m = 0.76, stationary 0.27), so the stationary field sits ~0.5 m higher
## purely to compensate. A raw-field difference is dominated by that bookkeeping
## offset; q_alpha (and the return levels) are pinned to the data, so their
## differences are the real barrier-vs-stationary effect. With q.location = 0.5
## the location is the bGEV median, i.e. the 2-yr return level.
bm <- setNames(model_pars$beta_m, model_pars$model)
qalpha_grid <- field_grid |>
    mutate(qalpha = bm[as.character(model)] + u)

qalpha_diff <- qalpha_grid |>
    filter(year == last_year) |>
    select(model, x, y, qalpha) |>
    pivot_wider(names_from = model, values_from = qalpha) |>
    transmute(x, y, diff = m0st_nb - m0st)

rl_diff <- rl_grid |>
    filter(period == 50) |>
    select(model, x, y, rl) |>
    pivot_wider(names_from = model, values_from = rl) |>
    transmute(x, y, diff = m0st_nb - m0st)

##==============================================================================
## Plots
##==============================================================================
## Shared map layers: land in grey, gauges as small points, square aspect, and
## the panel cropped to the domain. The grid rasters carry km coordinates, so
## coord_sf(crs = kmproj) keeps them aligned with the sf layers.
add_map_layers <- function(p) {
    p +
        geom_sf(data = barrier_sf, fill = "grey85", colour = "grey55",
            linewidth = 0.2, inherit.aes = FALSE) +
        geom_point(data = sites_km, aes(x, y), inherit.aes = FALSE,
            size = 1, colour = "black") +
        coord_sf(crs = kmproj, datum = kmproj,
            xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
            expand = FALSE) +
        labs(x = NULL, y = NULL)
}

save_fig <- function(plot, name, w, h) {
    ggsave(sprintf("%s/%s.png", figs_dir, name), plot,
        width = w, height = h, dpi = 300, device = png)
    ggsave(sprintf("%s/%s.pdf", figs_dir, name), plot,
        width = w, height = h, dpi = 300, device = cairo_pdf)
}

##------------------------------------------------------------------------------
## Figure 1: posterior-mean field, models x years
p_field <- ggplot(field_grid, aes(x, y, fill = u)) +
    geom_raster() +
    facet_grid(model ~ year, labeller = labeller(model = mod_labs)) +
    scale_fill_gradient2(name = "u(s,t)", low = "dodgerblue3", mid = "white",
        high = "firebrick", midpoint = 0)
    ## scale_fill_viridis_c(name = "u(s,t)", option = "D")
p_field <- add_map_layers(p_field)
save_fig(p_field, "compare_field_mean", 13, 11)

## Figure 1 (SD): posterior-SD field, models x years
## Companion uncertainty map for the posterior-mean field above. The SD is
## non-negative, so a sequential (viridis) scale is used instead of the
## diverging scale of the mean-field figure.
p_field_sd <- ggplot(field_sd_grid, aes(x, y, fill = u_sd)) +
    geom_raster() +
    facet_grid(model ~ year, labeller = labeller(model = mod_labs)) +
    scale_fill_viridis_c(name = "SD u(s,t)", option = "D")
p_field_sd <- add_map_layers(p_field_sd)
save_fig(p_field_sd, "compare_field_sd", 13, 11)

## Figure 1b: location field q_alpha = m + u
## qalpha_grid is built in the Differences section above. Putting the intercept
## back in places both models on a common level, so the panels are directly
## comparable; q_alpha is the bGEV median (the 2-yr return level).
p_qalpha <- ggplot(qalpha_grid, aes(x, y, fill = qalpha)) +
    geom_raster() +
    facet_grid(model ~ year, labeller = labeller(model = mod_labs)) +
    scale_fill_viridis_c(name = "q_alpha (m)")
p_qalpha <- add_map_layers(p_qalpha)
save_fig(p_qalpha, "compare_qalpha_mean", 13, 11)

##------------------------------------------------------------------------------
## Figure 2: location (q_alpha) difference (stationary minus barrier), last year
p_qalpha_diff <- ggplot(qalpha_diff, aes(x, y, fill = diff)) +
    geom_raster() +
    scale_fill_gradient2(name = "stationary - barrier", low = "dodgerblue3",
        mid = "white", high = "firebrick", midpoint = 0) +
    ggtitle(sprintf("Location (q_alpha) difference, %d", last_year))
p_qalpha_diff <- add_map_layers(p_qalpha_diff)
save_fig(p_qalpha_diff, "compare_qalpha_diff", 11, 6)

##------------------------------------------------------------------------------
## Figure 3: return-level maps, models x return periods, last year
p_rl_map <- ggplot(rl_grid, aes(x, y, fill = rl)) +
    geom_raster() +
    facet_grid(model ~ period,
        labeller = labeller(model = mod_labs,
            period = \(p) paste0(p, "-yr"))) +
    scale_fill_viridis_c(name = "Return\nlevel (m)", option = "C")
p_rl_map <- add_map_layers(p_rl_map)
save_fig(p_rl_map, "compare_rlevel_map", 11, 11)

##------------------------------------------------------------------------------
## Figure 4: return-level difference (stationary minus barrier), 50-yr period
p_rl_diff <- ggplot(rl_diff, aes(x, y, fill = diff)) +
    geom_raster() +
    scale_fill_gradient2(name = "stationary - barrier", low = "dodgerblue3",
        mid = "white", high = "firebrick", midpoint = 0) +
    ggtitle(sprintf("50-yr return-level difference, %d", last_year))
p_rl_diff <- add_map_layers(p_rl_diff)
save_fig(p_rl_diff, "compare_rlevel_diff", 11, 6)

##------------------------------------------------------------------------------
## Figure 5: per-site return-level curves, barrier vs stationary
## Sites ordered by country then name so the facets group geographically.
rl_sites <- rl_sites |>
    mutate(site_name = fct_reorder(site_name, as.integer(factor(country))))
p_rl_sites <- ggplot(rl_sites,
        aes(period, rl, colour = type)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ site_name, scales = "free_y") +
    ## scale_x_log10() +
    scale_colour_manual(name = "Spatial model",
        values = c(barrier = "firebrick", stationary = "dodgerblue3")) +
    labs(x = "Return period (years)", y = "Return level (m)") +
    theme(legend.position = "top",
        axis.text = element_text(size = 7))
save_fig(p_rl_sites, "compare_rlevel_sites", 14, 10)
