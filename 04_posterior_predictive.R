##==============================================================================
## Posterior-predictive return levels for the spatio-temporal models.
##
## Extends the point-estimate maps and curves of 02_compare_barrier.R with
## FULL posterior-predictive uncertainty. It propagates the joint posterior of
## the latent field u(s, t), the intercept m, and the bGEV spread and tail
## through return_level_bgev2() via inla.posterior.sample(), for the two
## spatio-temporal models only (m0st = barrier, m0st_nb = stationary), and
## produces:
##   * return-level MAPS on the grid at the last year: posterior median, the
##     95% credible-interval width, and a 2.5/50/97.5% quantile panel for the
##     50-yr level -- the uncertainty companions to 02's compare_rlevel_map;
##   * per-site return-level CURVES with 95% credible ribbons, in two figures:
##     one with all 44 sites and one with a few key coastal cities.
##
## WHY A NEW SCRIPT, AND WHY REFIT RATHER THAN UN-PRUNE 01
## ------------------------------------------------------
## inla.posterior.sample() needs misc$configs (config = TRUE), which
## 01_explore_models.R strips before caching (prune_fits <- TRUE). Two ways
## to get it back: (a) set prune_fits <- FALSE in 01, or (b) refit the two ST
## models with config = TRUE where the samples are drawn. We take (b):
##   - un-pruning bloats each cached ST fit and slows every downstream reader
##     (02 loads those files but needs only summaries);
##   - refitting confines the heavy config object to RAM during sampling and
##     writes only small tidy posterior summaries to rds.
##
## The refit is CHEAP because it does not re-estimate anything. We read the
## converged internal hyperparameter mode from the cached fit (fit$mode$theta,
## which survives pruning) and pass it back with
## control.mode = list(theta = ..., restart = FALSE): INLA reuses that mode and
## skips the mode-FINDING iterations, but still builds the ccd integration grid
## around it. We deliberately do NOT use fixed = TRUE: the bGEV spread and tail
## are hyperparameters (theta[1], theta[2] here), so fixing theta would freeze
## them across all draws and collapse exactly the return-level uncertainty this
## script exists to quantify. restart = FALSE keeps the full hyperparameter
## integration (spread, tail, range, sigma, rho all vary across draws) while
## reproducing the cached fit's mode, which also removes the small drift a fresh
## re-optimisation would introduce.
##
## To guarantee the refit is byte-identical to 01's model spec, we SOURCE 01 (it
## is written to be source()-able) and reuse its control lists, fms, stacks,
## priors and pinned geometry verbatim; a sanity check below confirms the refit
## reproduces the cached point estimates. 01 and 02 are left untouched.
##
## This is still HEAVIER than 02 (two refits with config = TRUE plus n_samples
## posterior draws), so it lives outside the fast, refit-free 02 and is not read
## by the manuscript yet. On a re-run it reloads the cached summaries and only
## redraws, without sourcing 01 or resampling.
##==============================================================================

here::i_am("04_posterior_predictive.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(sf)
    library(INLA)
})

## The stExtremes package provides create_indices(), marginals_st_gev(), return_level_bgev2().
library(stExtremes)

theme_set(theme_bw(base_size = 14))

rds_dir <- here::here("rds")
figs_dir <- here::here("figs")
if(!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

## Common projections (match the fitting script)
proj <- st_crs("+proj=longlat +datum=WGS84")
kmproj <- st_crs("+proj=utm +zone=30 +datum=WGS84 +units=km")

## Number of posterior draws. 1000 gives stable 2.5/97.5% quantiles; lower it
## (e.g. 500) if RAM is tight -- the field sample matrix is the largest object.
n_samples <- 1000L
pp_seed <- 20240601L

## The two ST models (same keys and labels as 02).
st_models <- c(m0st = "barrier", m0st_nb = "stationary")
mod_labs <- c(m0st = "barrier, st", m0st_nb = "stationary, st")
type_of <- c(m0st = "barrier", m0st_nb = "stationary")

## Key Irish coastal cities for the reduced per-site figure. site_name is the
## gauge label in dc_max.csv; city is the display label. Edit freely.
key_tbl <- tibble(
    site_name = c("DublinPort", "GalwayPort", "RingaskiddyNMCI",
        "Rosslare", "Sligo", "Dundalk"),
    city = c("Dublin", "Galway", "Cork (Ringaskiddy)",
        "Rosslare", "Sligo", "Dundalk"))

##==============================================================================
## Geometry, data and polygons (loaded directly, like 02; no 01 source needed
## on the redraw path)
##==============================================================================
geom <- readRDS(file.path(rds_dir, "pinned_geometry.rds"))
mesh <- geom$mesh
n_data <- length(geom$y)
n_years <- length(unique(geom$year_idx))

dc_max <- read_csv(here::here("data", "dc_max.csv"),
    show_col_types = FALSE)
stopifnot(isTRUE(all.equal(dc_max$max, geom$y)))

dc_sf <- st_as_sf(dc_max, coords = c("lon", "lat"), crs = proj) |>
    st_transform(kmproj)

data(area, package = "stExtremes")
data(barrier, package = "stExtremes")
area_km <- st_union(st_transform(area, kmproj))
barrier_sf <- st_sf(geometry = st_union(st_transform(barrier, kmproj)))
bb <- st_bbox(area_km)

## Unique gauge sites (km), for the site projector and the map overlay.
sites_km <- dc_sf |>
    mutate(x = st_coordinates(geometry)[, 1],
        y = st_coordinates(geometry)[, 2]) |>
    st_drop_geometry() |>
    distinct(country, site_name, x, y)
n_site <- nrow(sites_km)
stopifnot(all(key_tbl$site_name %in% sites_km$site_name))

## Year lookup and the last (mapped) year.
year_lookup <- dc_max |>
    distinct(year, year_idx) |>
    arrange(year)
last_idx <- max(geom$year_idx)
last_year <- year_lookup$year[match(last_idx, year_lookup$year_idx)]

## Field index and a matrix version of 02's field_year(): place the group-g
## entries of each posterior draw (columns) into mesh-node order.
w.index <- create_indices(geom$year_idx, mesh)
stopifnot(length(w.index$w) == mesh$n * n_years)

field_year_mat <- function(w_mat, g) {
    sel <- w.index$w.group == g
    node_of <- w.index$w[sel]
    out <- matrix(0, mesh$n, ncol(w_mat))
    out[node_of, ] <- w_mat[sel, ]
    out
}

##==============================================================================
## Prediction grid (identical construction and mask to 02) and cache tags
##==============================================================================
mask_km <- 75
n_pix <- 140
mv <- sprintf("mesh_%d", mesh$n)
gv <- sprintf("%s_n%d_msk%d", mv, n_pix, mask_km)

## Return periods: maps at a few key periods (match 02), curves on a full range.
map_periods <- c(10, 50, 100)
site_periods <- 2:100

##==============================================================================
## Plot helpers (shared with 02)
##==============================================================================
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

## Row-wise posterior summary across draws (columns): mean and 2.5/50/97.5%.
row_summ <- function(mat) {
    q <- t(apply(mat, 1, quantile, probs = c(0.025, 0.5, 0.975),
        names = FALSE))
    tibble(mean = rowMeans(mat), lo = q[, 1], med = q[, 2], hi = q[, 3])
}

##==============================================================================
## Sampling (only when the summaries are not already cached)
##==============================================================================
## Two tidy caches: the grid return-level summaries and the per-site curve
## summaries. Present -> skip the refit and the sampling entirely and just
## redraw. Delete them (or change n_pix / mask_km / n_samples) to resample.
grid_cache <- file.path(rds_dir, sprintf("pp_rlevel_grid_%s.rds", gv))
site_cache <- file.path(rds_dir, sprintf("pp_rlevel_sites_%s.rds", mv))

if(file.exists(grid_cache) && file.exists(site_cache)) {
    message("## Loading cached posterior-predictive summaries.")
    pp_grid <- readRDS(grid_cache)
    pp_sites <- readRDS(site_cache)
} else {
    ## Building blocks for the config = TRUE refit. Sourcing 01 reuses its
    ## deterministic INLA setup (threads, BLAS, ccd), pinned geometry, priors,
    ## stacks and control lists, and loads the cached (pruned) fits into lmn for
    ## the mode reuse and the sanity check. If the ST fit caches are absent this
    ## triggers 01's full fit, so run 01 first on a fresh machine.
    message("## Sourcing 01_explore_models.R for the fitting setup ...")
    source(here::here("01_explore_models.R"))
    ## bGEV constants, taken from the sourced control so they cannot drift from
    ## the fit: q.location / q.spread are the location/spread quantile levels
    ## and q.mix the bGEV mixing probabilities.
    alpha <- inla_control_bgev$q.location
    beta <- inla_control_bgev$q.spread
    p_a <- inla_control_bgev$q.mix[1]
    p_b <- inla_control_bgev$q.mix[2]
    ## Mesh -> grid and mesh -> site projectors. The grid is masked to the
    ## data-informed coastal band exactly as in 02 (sea cells within mask_km of
    ## a gauge), so the maps line up with 02's point-estimate versions.
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
    A_sites <- inla.spde.make.A(mesh,
        loc = as.matrix(sites_km[, c("x", "y")]))
    stopifnot(ncol(A_grid) == mesh$n, ncol(A_sites) == mesh$n)
    ## Location/period index vectors for the site curves: every (period, site)
    ## pair, so one vectorised return_level_bgev2() call covers a whole draw.
    per_vec <- rep(site_periods, each = n_site)
    site_vec <- rep(seq_len(n_site), times = length(site_periods))
    ## Refit one ST model at the CACHED converged hyperparameter mode, with
    ## config = TRUE so the fit can be sampled. control.mode restart = FALSE
    ## reuses fit$mode$theta and skips the mode search but keeps the ccd
    ## integration, so spread, tail, range, sigma and rho still vary across
    ## draws. Every other control is 01's, so the spec is identical.
    refit_at_mode <- function(key) {
        inla(
            formula = fms[[key]]$form,
            data = st_inla_data,
            family = inla_family,
            control.family = list(hyper = hyper_bgev,
                control.bgev = inla_control_bgev),
            control.predictor = inla_pred_st,
            control.fixed = inla_control_fixed,
            control.compute = inla_control_compute,
            control.inla = inla_control_inla,
            control.mode = list(theta = lmn[[key]]$mode$theta,
                restart = FALSE),
            verbose = inla_verbose,
            safe = inla_safe)
    }
    ## Pull the field draws (nodes x groups), intercept and bGEV spread/tail out
    ## of an inla.posterior.sample() list into plain matrices/vectors.
    extract_draws <- function(smp) {
        ln <- rownames(smp[[1]]$latent)
        w_rows <- grep("^w:", ln)
        m_row <- grep("^m:", ln)
        stopifnot(length(w_rows) == mesh$n * n_years, length(m_row) == 1L)
        hn <- names(smp[[1]]$hyperpar)
        sp_nm <- grep("spread for BGEV", hn, value = TRUE)
        tl_nm <- grep("^tail for BGEV", hn, value = TRUE)
        stopifnot(length(sp_nm) == 1L, length(tl_nm) == 1L)
        list(
            W = vapply(smp, \(s) s$latent[w_rows, 1],
                numeric(length(w_rows))),
            m = vapply(smp, \(s) s$latent[m_row, 1], numeric(1L)),
            spread = vapply(smp, \(s) s$hyperpar[[sp_nm]], numeric(1L)),
            tail = vapply(smp, \(s) s$hyperpar[[tl_nm]], numeric(1L)))
    }
    ## Refit at the cached mode, sample, and reduce to the grid and site
    ## return-level summaries. One model at a time keeps the large draw matrix
    ## (mesh$n * n_years x n_samples) from co-existing.
    pp_list <- imap(st_models, \(type, key) {
        message(sprintf(
            "## Refitting %s (%s) at the cached mode and sampling ...",
            key, type))
        fit <- refit_at_mode(key)
        ## Sanity: reusing the mode reproduces 01's cached intercept to within
        ## numerical noise. The ccd grid is rebuilt from a finite-difference
        ## Hessian around the reused mode, so a gap of a few percent of the
        ## posterior SD is expected and harmless; only a gap of a large fraction
        ## of an SD would mean the mode was not actually reused, so the check is
        ## scaled by the cached posterior SD rather than an absolute tolerance.
        b_refit <- fit$summary.fixed["m", "mean"]
        b_cache <- lmn[[key]]$summary.fixed["m", "mean"]
        sd_cache <- lmn[[key]]$summary.fixed["m", "sd"]
        pct_sd <- 100 * abs(b_refit - b_cache) / sd_cache
        message(sprintf(
            "##   intercept refit = %.5f, cached = %.5f (%.1f%% of post. SD)",
            b_refit, b_cache, pct_sd))
        if(pct_sd > 10) {
            warning(sprintf(
                "%s: refit intercept is %.1f%% of an SD from the cache.",
                key, pct_sd))
        }
        smp <- inla.posterior.sample(n_samples, fit, seed = pp_seed)
        dr <- extract_draws(smp)
        rm(smp)
        ## Field at the last year for every draw, then q_alpha = m + u.
        w_year <- field_year_mat(dr$W, last_idx)
        u_grid <- as.matrix(A_grid %*% w_year)
        u_site <- as.matrix(A_sites %*% w_year)
        q_grid <- sweep(u_grid, 2, dr$m, "+")
        q_site <- sweep(u_site, 2, dr$m, "+")
        ## Grid return levels: per period, one vectorised call per draw over all
        ## grid cells (all periods >= 2 land in the closed-form bGEV branch).
        grid_sum <- map(map_periods, \(tt) {
            rl <- vapply(seq_len(n_samples),
                \(i) return_level_bgev2(tt, q_grid[, i],
                    dr$spread[i], dr$tail[i], alpha = alpha, beta = beta,
                    p_a = p_a, p_b = p_b),
                numeric(nrow(q_grid)))
            bind_cols(
                tibble(model = key, period = tt,
                    x = grid_df$x, y = grid_df$y),
                row_summ(rl))
        }) |>
            list_rbind()
        ## Site curves: one vectorised call per draw over every (period, site).
        rl_site <- vapply(seq_len(n_samples),
            \(i) return_level_bgev2(per_vec, q_site[site_vec, i],
                dr$spread[i], dr$tail[i], alpha = alpha, beta = beta,
                p_a = p_a, p_b = p_b),
            numeric(length(per_vec)))
        site_sum <- bind_cols(
            tibble(model = key,
                country = sites_km$country[site_vec],
                site_name = sites_km$site_name[site_vec],
                period = per_vec),
            row_summ(rl_site))
        list(grid = grid_sum, site = site_sum)
    })
    pp_grid <- map(pp_list, "grid") |> list_rbind()
    pp_sites <- map(pp_list, "site") |> list_rbind()
    saveRDS(pp_grid, grid_cache)
    saveRDS(pp_sites, site_cache)
}

## Common factor levels / labels for the plots.
pp_grid <- pp_grid |>
    mutate(model = factor(model, levels = names(st_models)),
        ci_width = hi - lo)
pp_sites <- pp_sites |>
    mutate(model = factor(model, levels = names(st_models)),
        type = type_of[as.character(model)])

##==============================================================================
## Maps
##==============================================================================
## Map 1: posterior-median return level, models x periods (companion to 02's
## compare_rlevel_map, which shows the point estimate).
p_med <- ggplot(pp_grid, aes(x, y, fill = med)) +
    geom_raster() +
    facet_grid(model ~ period,
        labeller = labeller(model = mod_labs,
            period = \(p) paste0(p, "-yr"))) +
    scale_fill_viridis_c(name = "Median\nRL (m)", option = "C")
p_med <- add_map_layers(p_med)
save_fig(p_med, "pp_rlevel_median_map", 11, 8)

## Map 2: 95% credible-interval width (97.5% - 2.5%) -- the uncertainty map.
p_ciw <- ggplot(pp_grid, aes(x, y, fill = ci_width)) +
    geom_raster() +
    facet_grid(model ~ period,
        labeller = labeller(model = mod_labs,
            period = \(p) paste0(p, "-yr"))) +
    scale_fill_viridis_c(name = "95% CI\nwidth (m)", option = "D")
p_ciw <- add_map_layers(p_ciw)
save_fig(p_ciw, "pp_rlevel_ciwidth_map", 11, 8)

## Map 3: the 50-yr level as a 2.5 / 50 / 97.5% quantile panel, so the credible
## band is visible spatially for both models.
pp_grid_50 <- pp_grid |>
    filter(period == 50) |>
    select(model, x, y, `2.5%` = lo, `50%` = med, `97.5%` = hi) |>
    pivot_longer(c(`2.5%`, `50%`, `97.5%`),
        names_to = "quantile", values_to = "rl") |>
    mutate(quantile = factor(quantile, levels = c("2.5%", "50%", "97.5%")))
p_q50 <- ggplot(pp_grid_50, aes(x, y, fill = rl)) +
    geom_raster() +
    facet_grid(model ~ quantile, labeller = labeller(model = mod_labs)) +
    scale_fill_viridis_c(name = "50-yr\nRL (m)", option = "C")
p_q50 <- add_map_layers(p_q50)
save_fig(p_q50, "pp_rlevel_quantile_map_50yr", 11, 8)

##==============================================================================
## Per-site return-level curves with 95% credible ribbons
##==============================================================================
## Figure A: all 44 sites, ordered by country then name so the facets group
## geographically. Barrier vs stationary overlaid, each with its ribbon.
pp_sites_all <- pp_sites |>
    mutate(site_name = fct_reorder(site_name,
        as.integer(factor(country))))
p_sites_all <- ggplot(pp_sites_all,
        aes(period, med, colour = type, fill = type)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ site_name, scales = "free_y") +
    scale_colour_manual(name = "Spatial model",
        values = c(barrier = "firebrick", stationary = "dodgerblue3")) +
    scale_fill_manual(name = "Spatial model",
        values = c(barrier = "firebrick", stationary = "dodgerblue3")) +
    labs(x = "Return period (years)", y = "Return level (m)") +
    theme(legend.position = "top",
        axis.text = element_text(size = 7))
save_fig(p_sites_all, "pp_rlevel_sites_all", 14, 10)

## Figure B: a few key coastal cities, larger panels so the ribbons read.
pp_sites_key <- pp_sites |>
    inner_join(key_tbl, by = "site_name") |>
    mutate(city = factor(city, levels = key_tbl$city))
p_sites_key <- ggplot(pp_sites_key,
        aes(period, med, colour = type, fill = type)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ city, scales = "free_y") +
    scale_colour_manual(name = "Spatial model",
        values = c(barrier = "firebrick", stationary = "dodgerblue3")) +
    scale_fill_manual(name = "Spatial model",
        values = c(barrier = "firebrick", stationary = "dodgerblue3")) +
    labs(x = "Return period (years)", y = "Return level (m)",
        title = sprintf("Posterior-predictive return levels, %d", last_year)) +
    theme(legend.position = "top")
save_fig(p_sites_key, "pp_rlevel_sites_key", 11, 7)

message("## Done. Wrote 5 figures to figs and 2 summaries to rds.")
