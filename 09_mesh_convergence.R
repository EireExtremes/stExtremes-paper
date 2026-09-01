##==============================================================================
## Mesh convergence (TODO item 5; Referee 2, Sec. 4.2)
##
## Referee 2: "There is a lot of space devoted to claims about
## constructing a good mesh, but no investigation into whether the mesh
## so constructed is actually good. Or even how one might decide how good
## a mesh is."
##
## The check: refit the selected model on a coarser and a finer mesh and
## report how far the hyperparameters and the 50-year return levels move.
## If they barely move, the production mesh is fine enough for the
## quantities reported, which is the operational meaning of "good enough"
## here.
##
## Three meshes, controlled by the single divisor that sets max.edge in
## 01_explore_models.R:
##
##   divisor 10   318 nodes   coarse
##   divisor 15   577 nodes   production
##   divisor 22  1098 nodes   fine
##
## roughly half and double the production node count.
##
## TWO DELIBERATE SIMPLIFICATIONS, both to keep the runtime sane.
##
## 1. Only m0st_nb, the selected model whose results the paper reports.
##    Adding m0st would double the cost. Set `models` below to both if it
##    is wanted later; nothing else needs changing.
##
## 2. Return levels are computed at the POSTERIOR MEANS of the location,
##    spread and tail, not from joint posterior draws. A convergence
##    check asks whether the answer moves with the mesh, not what its
##    uncertainty is, so the point estimate is the right tool and it
##    avoids config = TRUE and the sampling step entirely.
##
## THE PRODUCTION MESH IS REFITTED HERE rather than read from the cache.
## That is not redundant. The cached fits are from INLA 26.5.21 and this
## script runs under 26.6.8, so comparing an alternative mesh against the
## cache would confound the mesh with the software version. All three
## fits here share one version, so the comparison is clean.
##
## The alternative geometries are BUILT ONCE AND CACHED, exactly as the
## production geometry is, so the mesh-dependent floating point work is
## not repeated and the check is reproducible.
##
## Inputs:  data/dc_max.csv, and 01 for the model settings
## Outputs: rds/mesh_geom_d<divisor>.rds        (pinned geometries)
##          rds/mesh_fit_<model>_d<divisor>.rds (fits, pruned)
##          rds/mesh_conv_hyper.rds, mesh_conv_rlevel.rds
##          figs/mesh_conv_rlevel.{png,pdf}
##
## Runtime about 4 hours: the coarse fit is quicker than production and
## the fine one is roughly twice as slow. Every geometry and every fit is
## cached, so an interrupted run resumes where it stopped.
##
## The conclusion is written at the end of this file after the run.
##==============================================================================

##==============================================================================
## Setup

here::i_am("09_mesh_convergence.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(sf)
    library(INLA)
})

library(stExtremes)

rds_dir <- here::here("rds")
figs_dir <- here::here("figs")
if(!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

theme_set(theme_bw(base_size = 14))

save_fig <- function(plot, name, w, h) {
    ggsave(sprintf("%s/%s.png", figs_dir, name), plot,
        width = w, height = h, dpi = 300, device = png)
    ggsave(sprintf("%s/%s.pdf", figs_dir, name), plot,
        width = w, height = h, dpi = 300, device = cairo_pdf)
}

col_ink <- "grey15"
col_bar <- "grey45"

## Which models, and which meshes. The production divisor must stay in
## the list: it is the baseline everything else is measured against.
models <- c("m0st_nb")
divisors <- c(coarse = 10, production = 15, fine = 22)
ref_div <- 15
ref_period <- 50

##==============================================================================
## Model settings, taken from 01 so they cannot drift
##
## Sourcing 01 also loads its cached fits, which is unnecessary here but
## harmless, and it is the only way to be sure the priors and controls
## are identical to the ones behind the reported results.

message("## Sourcing 01_explore_models.R for the settings ...")
source(here::here("01_explore_models.R"))

## Capture what is needed under names that cannot be shadowed by the
## per-mesh objects built below.
set_hyper <- hyper_bgev
set_bgev <- inla_control_bgev
set_fixed <- inla_control_fixed
set_inla <- inla_control_inla
set_rho <- prior_rho
set_range <- prior_range
set_sigma <- prior_sigma
alpha_loc <- set_bgev$q.location
beta_spr <- set_bgev$q.spread
p_a <- set_bgev$q.mix[1]
p_b <- set_bgev$q.mix[2]
s_blend <- set_bgev$beta.ab

## Lean compute: a convergence check needs the fit and WAIC, nothing
## else. CPO is the slow part and config the heavy part, so both are off.
conv_compute <- list(waic = TRUE, dic = TRUE, cpo = FALSE, config = FALSE,
    return.marginals.predictor = FALSE)

##==============================================================================
## Data and geometry

dc_csv <- read_csv(here::here("data/dc_max.csv"),
    show_col_types = FALSE)
dc_km <- st_as_sf(dc_csv, coords = c("lon", "lat"), crs = proj) |>
    st_transform(kmproj)
sites_km <- dc_km |>
    mutate(x = st_coordinates(geometry)[, 1],
        y = st_coordinates(geometry)[, 2]) |>
    st_drop_geometry() |>
    distinct(country, site_name, x, y)

data(area, package = "stExtremes")
area_km2 <- st_transform(area, kmproj)
yr_idx <- dc_csv$year_idx
n_yr <- length(unique(yr_idx))
last_yr <- max(yr_idx)

## Build one geometry for a given divisor, mirroring the build branch of
## 01 exactly and varying only max.edge. No barrier triangles: the
## stationary model does not use them, and skipping them keeps this
## simple.
build_geom <- function(dv) {
    f <- file.path(rds_dir, sprintf("mesh_geom_d%d.rds", dv))
    if(file.exists(f)) return(readRDS(f))
    message(sprintf("## Building geometry for divisor %d ...", dv))
    r2 <- st_bbox(area_km2) |>
        matrix(nrow = 2) |>
        apply(1, diff)
    rr <- mean(r2)
    max_edge <- (rr / dv) * c(1, 5)
    bound_outer <- max_edge[2] * 1.1
    buf <- (max_edge[1] / 2) / 1.5
    dom <- st_buffer(st_simplify(area_km2, TRUE, buf), buf)
    m0 <- fmesher::fm_mesh_2d_inla(boundary = dom,
        max.edge = max_edge[1] / 2, cutoff = max_edge[1] / 4)
    mm <- fmesher::fm_mesh_2d_inla(loc = m0$loc[, 1:2],
        max.edge = max_edge[2], offset = bound_outer,
        cutoff = max_edge[1] / 2)
    fmesher::fm_crs(mm) <- kmproj
    loc <- st_coordinates(dc_km)
    out <- list(divisor = dv, mesh = mm,
        As = inla.spde.make.A(mm, loc = loc),
        Ast = inla.spde.make.A(mm, loc = loc, group = yr_idx,
            n.group = n_yr))
    saveRDS(out, f)
    out
}

geoms <- map(set_names(divisors, names(divisors)), build_geom)
map_dbl(geoms, \(g) g$mesh$n)

##==============================================================================
## Fit one model on one geometry

fit_on_mesh <- function(key, g) {
    f <- file.path(rds_dir,
        sprintf("mesh_fit_%s_d%d.rds", key, g$divisor))
    if(file.exists(f)) {
        message(sprintf("## Cached: %s on %d nodes", key, g$mesh$n))
        return(readRDS(f))
    }
    message(sprintf("## Fitting %s on %d nodes (divisor %d) ...",
        key, g$mesh$n, g$divisor))
    spde <- inla.spde2.pcmatern(g$mesh, prior.range = set_range,
        prior.sigma = set_sigma, constr = TRUE)
    widx <- create_indices(yr_idx, g$mesh)
    n_d <- length(dc_csv$max)
    stk <- inla.stack(data = list(y = dc_csv$max),
        effects = list(widx, data.frame(m = rep(1, n_d))),
        A = list(g$Ast, 1), tag = "est")
    dat <- inla.stack.data(stk)
    dat$spread_x <- dat$tail_x <- matrix(nrow = n_d, ncol = 0)
    form <- inla.mdata(y, spread_x, tail_x) ~ -1 + m +
        f(w, model = spde, group = w.group, vb.correct = inla_vb_random,
            control.group = list(model = "ar1", hyper = set_rho))
    fit <- try(inla(form, data = dat, family = "bgev",
        control.family = list(hyper = set_hyper, control.bgev = set_bgev),
        control.predictor = list(A = inla.stack.A(stk), compute = TRUE,
            link = 1),
        control.fixed = set_fixed,
        control.compute = conv_compute,
        control.inla = set_inla,
        control.mode = list(theta = fms[[key]]$mode, restart = TRUE),
        verbose = FALSE, safe = TRUE), silent = FALSE)
    if(inherits(fit, "try-error")) {
        warning(sprintf("%s on divisor %d failed.", key, g$divisor))
        return(NULL)
    }
    ## Keep only what the summaries below need, so the cache stays small.
    keep <- list(divisor = g$divisor, model = key, nodes = g$mesh$n,
        waic = fit$waic$waic, dic = fit$dic$dic,
        secs = unname(fit$cpu.used[["Total"]]),
        hyper = marginals_st_gev(fit) |>
            imap(\(mg, nm) {
                z <- inla.zmarginal(as.matrix(mg), silent = TRUE)
                tibble(param = nm, mean = z$mean, lo = z$quant0.025,
                    hi = z$quant0.975)
            }) |>
            list_rbind(),
        m = fit$summary.fixed["m", "mean"],
        w = fit$summary.random$w$mean)
    saveRDS(keep, f)
    rm(fit)
    gc()
    keep
}

runs <- cross_join(tibble(key = models), tibble(nm = names(divisors))) |>
    pmap(\(key, nm) fit_on_mesh(key, geoms[[nm]])) |>
    compact()

##==============================================================================
## 1. Hyperparameters across meshes [APPENDIX]

hyper_tab <- map(runs, \(r) mutate(r$hyper, model = r$model,
    nodes = r$nodes)) |>
    list_rbind() |>
    filter(param %in% c("m", "spread", "tail", "rangeM", "sigmaM", "rho"))
hyper_tab |>
    select(model, nodes, param, mean, lo, hi) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    arrange(model, param, nodes) |>
    print(n = 30)
saveRDS(hyper_tab, file.path(rds_dir, "mesh_conv_hyper.rds"))

tibble(model = map_chr(runs, "model"), nodes = map_dbl(runs, "nodes"),
    WAIC = round(map_dbl(runs, "waic"), 1),
    secs = round(map_dbl(runs, "secs")))

##==============================================================================
## 2. Return levels across meshes [PAPER]
##
## Posterior means only, as explained in the header. The field is stored
## as (node, group) stacked by group, so the last year is the final block
## of mesh$n values.

rlevel_of <- function(r) {
    g <- geoms[[which(divisors == r$divisor)]]
    idx <- ((last_yr - 1) * g$mesh$n + 1):(last_yr * g$mesh$n)
    A <- inla.spde.make.A(g$mesh,
        loc = as.matrix(sites_km[, c("x", "y")]))
    q <- as.vector(A %*% r$w[idx]) + r$m
    hp <- \(nm) r$hyper$mean[r$hyper$param == nm]
    tibble(model = r$model, nodes = r$nodes, site = sites_km$site_name,
        rl = return_level_bgev2(ref_period, q = q, sb = hp("spread"),
            xi = hp("tail"), alpha = alpha_loc, beta = beta_spr,
            p_a = p_a, p_b = p_b, s = s_blend))
}

rlevel_tab <- map(runs, rlevel_of) |> list_rbind()
saveRDS(rlevel_tab, file.path(rds_dir, "mesh_conv_rlevel.rds"))

ref_nodes <- runs[[which(map_dbl(runs, "divisor") == ref_div)]]$nodes
shift_tab <- rlevel_tab |>
    left_join(rlevel_tab |>
        filter(nodes == ref_nodes) |>
        select(model, site, ref_rl = rl), by = c("model", "site")) |>
    filter(nodes != ref_nodes) |>
    group_by(model, nodes) |>
    summarise(max_abs_m = max(abs(rl - ref_rl)),
        med_abs_m = median(abs(rl - ref_rl)),
        max_pct = 100 * max(abs(rl - ref_rl) / ref_rl), .groups = "drop") |>
    mutate(across(where(is.numeric), \(x) round(x, 4)))
shift_tab
saveRDS(shift_tab, file.path(rds_dir, "mesh_conv_shift.rds"))

##==============================================================================
## 3. Figure [PAPER]

p_mesh <- rlevel_tab |>
    mutate(site = fct_reorder(site, rl),
        Mesh = factor(nodes, levels = sort(unique(nodes)),
            labels = paste(sort(unique(nodes)), "nodes"))) |>
    ggplot(aes(x = site, y = rl, shape = Mesh)) +
    geom_point(size = 1.8, colour = col_ink) +
    scale_shape_manual(values = c(1, 16, 3)) +
    coord_flip() +
    labs(x = NULL, y = sprintf("%g-year return level (m)", ref_period)) +
    theme(legend.position = "top", axis.text.y = element_text(size = 7))
p_mesh
save_fig(p_mesh, "mesh_conv_rlevel", w = 8, h = 9)

message("09_mesh_convergence.R complete.")

##==============================================================================
## DECISION RULE, written before the run
##
## Read `shift_tab`, column max_pct.
##
## (a) under about 2%. The production mesh has converged for the reported
##     quantities. Write it up in one short paragraph plus the figure,
##     and close Referee 2 Sec. 4.2.
## (b) 2 to 5%. Report the number honestly and say which direction the
##     bias runs. Still defensible, since it is small against the
##     credible intervals, which are 0.4 to 0.7 m wide.
## (c) above 5%, or the coarse and fine meshes disagree in opposite
##     directions from production. The mesh is NOT converged and the
##     production results would need refitting on the finer mesh. Check
##     the hyperparameters first: a moving Matern range with stable
##     return levels means the field is being reparametrised rather than
##     the answer changing.
##
## Whatever the outcome, do not claim the mesh is "good" in the abstract.
## The defensible claim is narrower: the quantities this paper reports
## are insensitive to halving or doubling the node count.
##==============================================================================
