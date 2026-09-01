##==============================================================================
## Final model run for the spatio-temporal bGEV paper.
##
## Reproducible, source()-able fitting of the retained models on any machine
## (laptop or server), using the revised data set dc_max.csv. INLA settings
## are the final, deterministic ones (full Laplace, ccd integration, single
## inner thread). Plots and exploratory data analysis are removed; only model
## fitting, caching and diagnostics remain.
##
## Changes from the exploratory version:
##  - All models carrying the global linear time trend (trendc) are removed.
##    A single coastline-wide linear trend in the location was judged
##    misleading, so only the intercept, spatial and spatio-temporal models
##    (barrier and stationary) are kept.
##  - INLA settings tightened for the final run (see control.inla below).
##  - Fitted objects are pruned before caching: the per-configuration latent
##    approximations stored by config = TRUE (the bulk of the ~1.4 GB
##    spatio-temporal files) are dropped AFTER the group CV is computed, and
##    the per-node latent marginals are not stored. This cuts the saved-file
##    and in-RAM size by an order of magnitude without changing any reported
##    quantity. Set prune_fits <- FALSE to keep the full objects (needed only
##    if a downstream step draws posterior samples from these fits).
##  - A hold-out cross-validation block is appended: leave selected sites out
##    (prediction at ungauged sites) and selected years out (temporal AR1
##    interpolation), refit the spatio-temporal models and compare observed
##    versus predicted. Run for the barrier and stationary ST models on the
##    SAME held-out sites/years so the two are directly comparable.
##
## Reproducibility strategy:
##  1. The mesh and every geometry-derived object are built ONCE and reused
##     from a saved file. This removes the GEOS / PROJ / fmesher floating
##     point differences that change the mesh node count (574 vs 577) and
##     hence every downstream estimate. Build the file once on one machine,
##     then COPY it to the others. Do NOT regenerate it per machine.
##  2. Single-threaded linear algebra (BLAS B = 1, INLA inner threads = 1)
##     removes threading nondeterminism; outer parallelism gives the speed.
##  3. int.strategy = "ccd" integrates over the hyperparameters and is far
##     more stable across machines than the "eb" used originally.
##==============================================================================

## Anchor the project root. Edit this path to match where you save the file.
here::i_am("01_explore_models.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(sf)
    library(INLA)
})

## The stExtremes package provides tri_on_mesh(), create_indices(), map_tail()
library(stExtremes)

rds_dir <- here::here("rds")
if(!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)

## Reproducible parallelism.
## INLA nested threads "A:B": A = parallel hyperparameter (theta) evaluations
## for speed, B = threads per evaluation. B MUST be 1: it keeps the per-theta
## linear algebra deterministic, which is what pins the results across
## machines. A can be set to the cores you have WITHOUT breaking cross-machine
## reproducibility, because the ccd integration combines the per-theta
## evaluations by a fixed weighted sum that does not depend on evaluation
## order. The default below is "8:1"; raising A only costs memory (each
## parallel evaluation holds a copy of the latent field), so lower A if the
## large spatio-temporal models run out of RAM, and raise it if you have more
## cores. Set "1:1" for a fully serial reference run (identical results).
n_threads <- "8:1"
inla.setOption(num.threads = n_threads)

## Single-threaded BLAS removes a major nondeterminism source and barely
## affects speed here (the heavy work is inside INLA, not R-level BLAS).
if(requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
} else {
    message("## RhpcBLASctl not installed: BLAS stays multi-threaded.")
}

## Record the full computational stack. Diff this file across machines: the
## geometry libraries and BLAS are the uncontrolled variables, not the R
## package versions.
env_info_file <- file.path(rds_dir,
    sprintf("env_info_%s.rds", Sys.info()[["nodename"]]))
if (!file.exists(env_info_file)) {
    env_info <- list(
        when = Sys.time(),
        node = Sys.info()[["nodename"]],
        r = R.version.string,
        blas = sessionInfo()$BLAS,
        lapack = sessionInfo()$LAPACK,
        geos_proj_gdal = sf::sf_extSoftVersion(),
        pkgs = vapply(
            c("INLA", "fmesher", "INLAspacetime", "sf", "Matrix"),
            \(p) as.character(packageVersion(p)), character(1)))
    saveRDS(env_info, env_info_file)
}
env_info <- readRDS(env_info_file)

## Common projections
proj <- st_crs("+proj=longlat +datum=WGS84")
kmproj <- st_crs("+proj=utm +zone=30 +datum=WGS84 +units=km")

## Pinned geometry.
## Built once and reused. The build branch is the only place that touches the
## floating-point geometry pipeline; everything afterwards is deterministic.
geom_file <- file.path(rds_dir, "pinned_geometry.rds")

if(!file.exists(geom_file)) {
    message("## Building pinned geometry (run once, then copy the file) ...")
    ## data(dc_max, package = "stExtremes")
    ## load(here::here("data/dc_max_all_time_residual.rda"))
    dc_max <- read_csv(here::here("data/dc_max.csv"))
    dc_max <- dc_max |>
        st_as_sf(coords = c("lon", "lat"), crs = proj)
    data(area, package = "stExtremes")
    area_km <- st_transform(area, kmproj)
    dc_max_km <- st_transform(dc_max, kmproj)
    ## Mesh parameters from the size of the domain
    r2 <- st_bbox(area_km) |>
        matrix(nrow = 2) |>
        apply(1, diff)
    rr <- mean(r2)
    max.edge <- (rr / 15) * c(1, 5)
    bound.outer <- max.edge[2] * 1.1
    buffer_size <- (max.edge[1] / 2) / 1.5
    domain_buffered <- st_buffer(
        st_simplify(area_km, TRUE, buffer_size), buffer_size)
    mesh0 <- fmesher::fm_mesh_2d_inla(
        boundary = domain_buffered,
        max.edge = max.edge[1] / 2,
        cutoff = max.edge[1] / 4)
    mesh <- fmesher::fm_mesh_2d_inla(
        loc = mesh0$loc[, 1:2],
        max.edge = max.edge[2],
        offset = bound.outer,
        cutoff = max.edge[1] / 2)
    fmesher::fm_crs(mesh) <- kmproj
    tri_barrier <- tri_on_mesh(mesh, domain_buffered, inverse = TRUE)
    As <- inla.spde.make.A(mesh, loc = st_coordinates(dc_max_km))
    Ast <- inla.spde.make.A(mesh,
        loc = st_coordinates(dc_max_km),
        group = dc_max_km$year_idx,
        n.group = length(unique(dc_max_km$year_idx)))
    saveRDS(
        list(mesh = mesh, tri_barrier = tri_barrier, As = As, Ast = Ast,
            y = dc_max_km$max, year_idx = dc_max_km$year_idx),
        geom_file)
}

geom <- readRDS(geom_file)
mesh <- geom$mesh
tri_barrier <- geom$tri_barrier
n_data <- length(geom$y)
n_years <- length(unique(geom$year_idx))

## Projector rows must sum to the number of observations
stopifnot(isTRUE(all.equal(n_data, sum(geom$As))))
stopifnot(isTRUE(all.equal(n_data, sum(geom$Ast))))

## Priors and SPDE barrier model
prior_sigma <- c(1, 0.01)            # P(sigma > 1) = 0.01
prior_range <- c(50, 0.01)           # P(range < 50 km) = 0.01
barrier_model <- INLAspacetime::barrierModel.define(
    mesh,
    barrier.triangles = tri_barrier,
    prior.range = prior_range,
    prior.sigma = prior_sigma,
    constr = TRUE)

## Non-barrier counterpart: a stationary Matern SPDE on the SAME mesh and the
## SAME PC priors. Used for the barrier vs non-barrier comparison so the only
## difference between the two is whether spatial correlation respects the
## coastline. Everything else (priors, stacks, indices, A, bGEV) is held fixed.
stationary_model <- inla.spde2.pcmatern(
    mesh,
    prior.range = prior_range,
    prior.sigma = prior_sigma,
    constr = TRUE)

## PC prior for the AR1 temporal correlation. pc.cor1 has its base (null)
## model at rho = 1, i.e. it shrinks the AR1 field towards "constant in time"
## (the parsimonious choice for a latent field) and penalises departures from
## it. pc.cor1 is only valid when alpha > sqrt((1 - u) / 2); at u = 0 that floor
## is sqrt(0.5) ~ 0.707, so the original P(rho > 0) = 0.9 can be weakened only
## to about that bound. We use P(rho > 0) = 0.8 (20% of the prior mass below 0,
## vs 10% before) -- a genuine relaxation that stays comfortably valid. 0.7 is
## INVALID (below the floor) and is rejected by inla; ~0.72 is the vaguest
## pc.cor1 allows at u = 0. To shrink towards independence instead of
## persistence (which has no such floor near 0), switch to "pc.cor0", e.g.
## param = c(0.5, 0.5).
prior_rho <- list(theta = list(prior = "pc.cor0", param = c(0.5, 0.5)))

## bGEV priors and controls
hyper_spread <- list(initial = 1, fixed = FALSE, prior = "loggamma",
    param = c(3, 3))
tail <- 0.1
tail_interval <- c(0, 0.5)
tail_intern <- map_tail(tail, tail_interval, inverse = TRUE)
hyper_tail <- list(initial = tail_intern, prior = "pc.gevtail",
    param = c(7, tail_interval), fixed = FALSE)
hyper_bgev <- list(spread = hyper_spread, tail = hyper_tail)

inla_control_bgev <- list(
    q.location = 0.5,                # median for the location
    q.spread = 0.5,                  # quantile level for the spread
    q.mix = c(0.05, 0.2),
    beta.ab = 5)

## Indices and stacks (data and stacks differ for spatial and spatio-temporal)
w.index <- create_indices(geom$year_idx, mesh)

stk.s <- inla.stack(
    data = list(y = geom$y),
    effects = list(
        w = 1:barrier_model$f$n,
        ## Intercept only. The second A block is the identity, so this
        ## data.frame MUST have n_data rows: rep(1, n_data), not a scalar 1.
        data.frame(m = rep(1, n_data))),
    A = list(geom$As, 1),
    tag = "est.s")

stk.st <- inla.stack(
    data = list(y = geom$y),
    effects = list(
        w.index,
        ## Intercept only. The second A block is the identity, so this
        ## data.frame MUST have n_data rows: rep(1, n_data), not a scalar 1.
        data.frame(m = rep(1, n_data))),
    A = list(geom$Ast, 1),
    tag = "est.st")

s_inla_data <- inla.stack.data(stk.s)
st_inla_data <- inla.stack.data(stk.st)
s_inla_data$spread_x <- s_inla_data$tail_x <-
    st_inla_data$spread_x <- st_inla_data$tail_x <-
        matrix(nrow = n_data, ncol = 0)

inla_pred_s <- list(A = inla.stack.A(stk.s), compute = TRUE, link = 1)
inla_pred_st <- list(A = inla.stack.A(stk.st), compute = TRUE, link = 1)

## Exposure covariate on the SPREAD.
##
## 05_explore_spread_shape_inla.R shows that a single spread for the
## whole coast is the one clearly violated assumption of these models:
## freeing the spread gauge by gauge is worth about 173 WAIC, while
## freeing the tail is worse than pooling. The spread tracks basin
## geometry rather than exposure to the open Atlantic, the Bristol
## Channel being the largest because a converging estuary amplifies
## surge. INLA's bGEV takes fixed covariates on the spread through the
## second argument of inla.mdata, so this is a data change, not a new
## likelihood: the spread becomes spread * exp(spread_x %*% beta).
##
## The grouping is a stated rule, not a fitted one, and is deliberately
## coarse. Atlantic is the reference level, so beta1 is the Bristol
## Channel contrast and beta2 the Irish Sea contrast, both against the
## Atlantic. Their priors are the bGEV defaults, normal(0, 300).
dc_csv <- read_csv(here::here("data/dc_max.csv"),
    show_col_types = FALSE)

## The pinned geometry stores y in the row order of that file. Verify it
## rather than assume it: a silent misalignment here would attach every
## gauge's covariate to the wrong rows.
stopifnot(nrow(dc_csv) == n_data,
    identical(is.na(dc_csv$max), is.na(geom$y)),
    isTRUE(all.equal(dc_csv$max[!is.na(dc_csv$max)],
        geom$y[!is.na(geom$y)])))

exposure <- with(dc_csv, case_when(
    lon < -8.0 ~ "Atlantic",
    lat < 52.0 & lon > -5.5 ~ "Bristol Channel",
    TRUE ~ "Irish Sea"))
spread_x_exp <- model.matrix(~ exposure)[, -1, drop = FALSE]

s_inla_data_exp <- s_inla_data
st_inla_data_exp <- st_inla_data
s_inla_data_exp$spread_x <- st_inla_data_exp$spread_x <- spread_x_exp

## INLA controls (final run: deterministic, full Laplace)
inla_family <- "bgev"
inla_vb_random <- FALSE
inla_verbose <- TRUE
## safe = TRUE lets INLA detect a failed optimisation and rerun it from a
## perturbed start. It is kept ON because the hold-out refits run on reduced
## data, where a warm start can occasionally stall; the rerun recovers it. It
## does NOT affect reproducibility: the result is still pinned by the shared
## geometry, single-threaded BLAS and inner threads = 1.
inla_safe <- TRUE

inla_control_fixed <- list(prec.intercept = 100, prec = 100)
## config = TRUE is required by inla.group.cv (it reads the per-configuration
## latent approximations), so it stays ON for the fit. Those configs and the
## per-node latent-field marginals are the bulk of the object size, so the fit
## is PRUNED after the CV is computed (see prune_inla() and the fit loop). We do
## NOT set return.marginals = FALSE here: in this INLA build that also drops the
## hyperparameter marginals that marginals_st_gev() needs. Instead INLA returns
## everything and prune_inla() strips only the heavy latent-field marginals,
## leaving marginals.hyperpar intact.
inla_control_compute <- list(cpo = TRUE, waic = TRUE, dic = TRUE,
    config = TRUE, return.marginals.predictor = FALSE)

## Lean compute settings for the hold-out refits below. They only need the
## fitted values (the posterior of the median annual maximum at each
## site-year), so CPO / WAIC / DIC and config sampling are switched off (CPO is
## the slow part and config is the heavy part). The latent-field marginals are
## left at their default and stripped by prune_inla() before caching, so the
## saved hold-out objects are still small.
inla_control_compute_pred <- list(cpo = FALSE, waic = FALSE, dic = FALSE,
    config = FALSE, return.marginals.predictor = FALSE)

## Final integration / approximation settings.
## int.strategy = "ccd": integrates over the 5 hyperparameters and is far more
## stable across machines than the "eb" the instability traced back to. With 5
## hyperparameters a full grid is infeasible and "eb" is too crude, so ccd is
## the practical middle ground. strategy = "laplace" is the full Laplace
## approximation for the latent field -- more accurate than simplified.laplace,
## but slower on the large spatio-temporal field. If a fit is too slow,
## "simplified.laplace" is the documented fallback and shifts the estimates
## only marginally. cmin = 0 is set explicitly (matters on flat surfaces); h,
## tolerance and restart are the tightened final-run values.
inla_control_inla <- list(
    int.strategy = "ccd",
    strategy = "laplace",
    cmin = 0,
    h = 0.005,
    tolerance = 1e-8,
    restart = 1L)

## Single entry point for every bGEV fit (main fits and hold-out refits). The
## shared controls live here; callers vary only the formula, the data/stack,
## the starting mode, the predictor control and (for the hold-out) the leaner
## compute list.
fit_bgev <- function(form, data, mode, pred,
                     compute = inla_control_compute,
                     safe = inla_safe) {
    inla(
        formula = form,
        data = data,
        family = inla_family,
        control.family = list(hyper = hyper_bgev,
            control.bgev = inla_control_bgev),
        control.predictor = pred,
        control.fixed = inla_control_fixed,
        control.compute = compute,
        control.inla = inla_control_inla,
        control.mode = list(theta = mode, restart = TRUE),
        verbose = inla_verbose,
        safe = safe)
}

## Prune the parts of a fitted object that are not needed once the group CV and
## the summaries have been extracted. misc$configs (stored by config = TRUE) is
## the per-configuration Gaussian approximation of the full latent field and is
## what inflates the spatio-temporal objects to ~1.4 GB; the per-node latent
## marginals are the next largest. Removing them leaves every quantity this
## script reports (DIC/WAIC/CPO, summary.fixed/random/hyperpar, the
## hyperparameter marginals, summary.fitted.values) untouched, so the results
## are unchanged -- only the file and the in-RAM footprint shrink.
##
## TRADE-OFF: a pruned object can no longer feed inla.posterior.sample() or a
## fresh inla.group.cv(). The group CV is computed and cached before pruning,
## so this script is unaffected; if a DOWNSTREAM step needs posterior samples
## from these fits (e.g. posterior-predictive return levels), either set
## prune_fits <- FALSE here or refit the single selected model with
## config = TRUE in that step. The pruning is purely a storage choice and does
## not affect reproducibility of the fit itself.
prune_fits <- TRUE

prune_inla <- function(fit) {
    if(is.null(fit)) return(fit)
    fit$misc$configs <- NULL
    fit$marginals.random <- NULL
    fit$marginals.linear.predictor <- NULL
    fit$marginals.fitted.values <- NULL
    fit
}

## Model definitions. Each entry carries its formula, a short label, the stack
## type ("s" or "st") and a starting mode for the hyperparameters.
##
## The trend models (m1, m1s, m1st and their _nb twins) are removed: a single
## coastline-wide linear trend in the location was judged misleading. Dropping
## trendc is a fixed-effect change only, so the starting modes of the retained
## models are unchanged (trendc carried no hyperparameter). Seven models:
## the intercept (m0), the barrier spatial (m0s) and spatio-temporal (m0st)
## models, their stationary non-barrier counterparts m0s_nb and m0st_nb,
## which replace barrier_model by a stationary Matern on the SAME mesh and the
## SAME PC priors, and the two spread-covariate versions of the
## spatio-temporal pair (m0st_exp, m0st_nb_exp). The barrier vs non-barrier
## pair isolates the single effect of whether spatial correlation respects the
## coastline; the _exp pair isolates the single effect of letting the spread
## depend on the basin.
fms <- list(
    m0 = list(
        form = inla.mdata(y, spread_x, tail_x) ~ -1 + m,
        label = "y ~ 1",
        type = "s",
        mode = c(-0.5, -3.15)),
    m0s = list(
        form = inla.mdata(y, spread_x, tail_x) ~ -1 + m +
            f(w, model = barrier_model, vb.correct = inla_vb_random),
        label = "y ~ 1 + s",
        type = "s",
        mode = c(-0.7, -4.4, 5.7, -1.3)),
    m0st = list(
        form = inla.mdata(y, spread_x, tail_x) ~ -1 + m +
            f(w, model = barrier_model, group = w.group,
                vb.correct = inla_vb_random,
                control.group = list(model = "ar1", hyper = prior_rho)),
        label = "y ~ 1 + st",
        type = "st",
        mode = c(-1.1, -1.6, 5.9, -0.4, 3.9)),
    ## Non-barrier (stationary Matern) counterparts. Formulas, stacks and
    ## starting modes match their barrier twins; only the SPDE object differs
    ## (stationary_model replaces barrier_model).
    m0s_nb = list(
        form = inla.mdata(y, spread_x, tail_x) ~ -1 + m +
            f(w, model = stationary_model, vb.correct = inla_vb_random),
        label = "y ~ 1 + s (no barrier)",
        type = "s",
        mode = c(-0.7, -4.4, 5.7, -1.3)),
    m0st_nb = list(
        form = inla.mdata(y, spread_x, tail_x) ~ -1 + m +
            f(w, model = stationary_model, group = w.group,
                vb.correct = inla_vb_random,
                control.group = list(model = "ar1", hyper = prior_rho)),
        label = "y ~ 1 + st (no barrier)",
        type = "st",
        mode = c(-1.1, -1.6, 5.9, -0.4, 3.9)),
    ## Spread covariate versions of the two spatio-temporal models. The
    ## formulas are identical to m0st and m0st_nb; the only difference is
    ## that spread_x carries the exposure dummies instead of zero
    ## columns, which is why they need the _exp data lists. Fitting them
    ## here rather than in a separate script keeps them in the same
    ## comparison table as every other model.
    ##
    ## Starting modes: the two spread coefficients enter the
    ## hyperparameter vector immediately AFTER the tail and BEFORE the
    ## f() terms, i.e. spread, tail, beta1, beta2, range, stdev, rho.
    ## This ordering was checked against a fitted model, not assumed. The
    ## twin's mode is reused with two zeros inserted, zero meaning no
    ## exposure effect.
    m0st_exp = list(
        form = inla.mdata(y, spread_x, tail_x) ~ -1 + m +
            f(w, model = barrier_model, group = w.group,
                vb.correct = inla_vb_random,
                control.group = list(model = "ar1", hyper = prior_rho)),
        label = "y ~ 1 + st, spread ~ exposure",
        type = "st",
        spread = "exposure",
        mode = c(-1.1, -1.6, 0, 0, 5.9, -0.4, 3.9)),
    m0st_nb_exp = list(
        form = inla.mdata(y, spread_x, tail_x) ~ -1 + m +
            f(w, model = stationary_model, group = w.group,
                vb.correct = inla_vb_random,
                control.group = list(model = "ar1", hyper = prior_rho)),
        label = "y ~ 1 + st, spread ~ exposure (no barrier)",
        type = "st",
        spread = "exposure",
        mode = c(-1.1, -1.6, 0, 0, 5.9, -0.4, 3.9)))

nm <- length(fms)

## Result file names tied to the (now fixed) mesh
mv <- sprintf("mesh_%d", mesh$n)
fnms <- file.path(rds_dir,
    sprintf("explore_dc_max_bGEV_%s_%s.rds", names(fms), mv))

## Fit, group-cross-validate, prune and cache each model in one pass.
## A failed fit is NOT cached, so it is retried on the next run (the original
## script cached NULL on failure, which silently froze the failure). The fit
## and its group CV are cached together: the CV is computed from the full
## (config = TRUE) object, then the object is pruned (see prune_inla) before it
## is written, so the saved file is small. On a re-run the pruned fit and the
## cached CV are loaded together and nothing is recomputed.
lmn <- vector("list", nm)
names(lmn) <- names(fms)
lmn_cv <- vector("list", nm)
names(lmn_cv) <- names(fms)
fnms_cv <- gsub("\\.rds$", "_cv.rds", fnms)

t_start <- Sys.time()
for(i in seq_len(nm)) {
    message(sprintf("## Fitting %s (%s) ==========",
        names(fms)[i], fms[[i]]$label))
    if(file.exists(fnms[i])) {
        lmn[[i]] <- readRDS(fnms[i])
        if(file.exists(fnms_cv[i])) lmn_cv[[i]] <- readRDS(fnms_cv[i])
        next
    }
    is_st <- fms[[i]]$type == "st"
    ## Models without a spread entry get the zero-column spread_x. The
    ## isTRUE() guard keeps that working for the entries that predate the
    ## exposure models and carry no such field.
    use_exp <- isTRUE(fms[[i]]$spread == "exposure")
    fit_data <- if(is_st) {
        if(use_exp) st_inla_data_exp else st_inla_data
    } else {
        if(use_exp) s_inla_data_exp else s_inla_data
    }
    fit <- try(
        fit_bgev(
            form = fms[[i]]$form,
            data = fit_data,
            mode = fms[[i]]$mode,
            pred = if(is_st) inla_pred_st else inla_pred_s),
        silent = FALSE)
    if(inherits(fit, "try-error")) {
        warning(sprintf("Model %s failed; not cached, will retry next run.",
            names(fms)[i]))
        next
    }
    ## Group cross-validation. Dependence-aware and the metric that stayed
    ## stable across runs (leave-one-out CPO / lCPO below is unreliable under
    ## strong latent dependence). It reads the config = TRUE approximations, so
    ## it MUST run before the fit is pruned.
    cv <- try(inla.group.cv(fit, num.level.sets = 5), silent = FALSE)
    if(!inherits(cv, "try-error")) {
        lmn_cv[[i]] <- cv
        saveRDS(cv, fnms_cv[i])
    }
    if(prune_fits) fit <- prune_inla(fit)
    lmn[[i]] <- fit
    saveRDS(fit, fnms[i])
}
t_fit <- difftime(Sys.time(), t_start, units = "mins")

## Diagnostics.
labels <- vapply(fms, \(x) x$label, character(1))

## Model comparison criteria (matches the table layout reported before)
diag_tab <- tibble(
    Model = names(fms),
    Formula = labels,
    DIC = map_dbl(lmn, \(x) x$dic$dic %||% NA_real_),
    WAIC = map_dbl(lmn, \(x) x$waic$waic %||% NA_real_),
    lCPO = map_dbl(lmn, \(x)
        if(is.null(x)) NA_real_ else -sum(log(x$cpo$cpo), na.rm = TRUE)),
    LGOCV = map_dbl(lmn_cv, \(x)
        if(is.null(x)) NA_real_ else -sum(log(x$cv), na.rm = TRUE)))

## Convergence and timing: a 3-hour run or a non-zero mode status means the
## optimiser struggled and the criteria above are not trustworthy.
conv_tab <- tibble(
    Model = names(fms),
    mode_status = map_dbl(lmn, \(x) x$mode$mode.status %||% NA_real_),
    cpo_fail = map_dbl(lmn, \(x)
        if(is.null(x)) NA_real_ else sum(x$cpo$failure > 0, na.rm = TRUE)),
    secs = map_dbl(lmn, \(x)
        if(is.null(x)) NA_real_ else unname(x$cpu.used[["Total"]])))

## Key hyperparameters on the natural (interpretable) scale.
## summary.hyperpar reports "Theta1/Theta2 for w" on INLA's internal (log)
## scale, so reading their means directly is wrong. marginals_st_gev() back
## transforms each marginal (range = exp(Theta1), sigma = exp(Theta2); rho and
## the bGEV spread/tail are identity) and renames them to rangeM, sigmaM, rho,
## spread, tail. We then summarise the transformed marginals.
hp_summary <- function(fit) {
    if(is.null(fit)) return(tibble())
    marginals_st_gev(fit) |>
        imap(\(marg, nm) {
            z <- inla.zmarginal(as.matrix(marg), silent = TRUE)
            tibble(param = nm, mean = z$mean, sd = z$sd)
        }) |>
        list_rbind()
}

## Long form: posterior mean and sd (natural scale) for every fixed effect
## and hyperparameter of each model.
param_summary <- imap(lmn, \(fit, m) {
    s <- hp_summary(fit)
    if(nrow(s) == 0) NULL else mutate(s, Model = m, .before = 1)
}) |>
    compact() |>
    list_rbind()

## Compact table with the key bGEV and spatio-temporal hyperparameters
param_names <- c("m", "spread", "tail", "rangeM", "sigmaM", "rho")
param_tab <- param_summary |>
    filter(param %in% param_names) |>
    select(Model, param, mean) |>
    pivot_wider(names_from = param, values_from = mean) |>
    select(Model, any_of(param_names))

## Full hyperparameter summaries kept for exact cross-machine comparison
hyperpars <- map(lmn, \(x) if(is.null(x)) NULL else x$summary.hyperpar)


fl <- file.path(rds_dir, sprintf("explore_diagnostics_%s.rds", env_info$node))
if(!file.exists(fl)) {
    diagnostics <- list(
        env = env_info,
        mesh_n = mesh$n,
        fit_minutes = as.numeric(t_fit),
        diag_tab = diag_tab,
        conv_tab = conv_tab,
        param_tab = param_tab,
        hyperpars = hyperpars)
    saveRDS(diagnostics, fl)
} else {
    diagnostics <- readRDS(fl)
}

## Summary of posterior for all (hyper)parameters. NULL-guarded so a failed
## (uncached) fit stays NULL instead of erroring the whole summary.
res_marg <- map(lmn,
    \(x) if(is.null(x)) NULL else summary(marginals_st_gev(x)))

## Print the results so they appear when the script is source()d
print(env_info)
message(sprintf("## Mesh nodes: %d   (must match across machines)", mesh$n))
## message(sprintf("## Total fitting time: %.1f minutes", as.numeric(t_fit)))
print(diag_tab)
print(conv_tab)
print(param_tab)
print(res_marg)

## The PIT values are only meaningful if cpo_fail == 0. See below:
## INLA computes CPO without refitting, by correcting the full-data posterior
## for the removal of one point. When an observation is strongly informed by
## its neighbours -- exactly the case under an AR1-in-time plus spatial field
## -- removing it shifts the posterior materially, the approximation's
## assumptions break, and INLA sets failure > 0. The spatial-only models have
## weaker dependence (no long temporal chaining per site), so few or no
## failures; the intercept model has no latent field, so none. For the
## spatio-temporal models a large fraction of the site-years are flagged, and
## that is precisely why lCPO is unreliable for those rows: log(cpo) blows up
## at the flagged points and inflates the sum, while the dependence-aware
## LGOCV stays well behaved. Read LGOCV, not lCPO, for the ST models.

## ok <- c("m0", "m0s", "m0s_nb")
## pit_df <- imap(lmn[ok], \(fit, m) tibble(Model = m, pit = fit$cpo$pit)) |>
##     list_rbind()
## ggplot(pit_df, aes(pit)) +
##     geom_histogram(breaks = seq(0, 1, 0.1), boundary = 0) +
##     facet_wrap(~ Model) +
##     labs(x = "PIT", y = "Count")

## For the ST models (m0st, m0st_nb) a PIT histogram does not work: the PIT
## values at the flagged points are as unreliable as their CPO, so a histogram
## over all points is contaminated, and restricting to failure == 0 leaves a
## non-random subset -- the least-influential points -- which tends to look
## deceptively uniform and tells you little. The dependence-aware route for
## those two is the hold-out (leave sites or years and predict) calibration
## check implemented below, not PIT.


##==============================================================================
## Hold-out cross-validation for the spatio-temporal models.
##
## Two leave-out experiments, each refitting the ST models on reduced data and
## comparing observed vs predicted at the held-out cells:
##
##  * SITE hold-out (leave whole sites out): emulates prediction at an UNGAUGED
##    location, which is the paper's stated aim and the test where the barrier
##    and stationary fields genuinely differ -- they borrow strength across the
##    coastline in different ways, so leaving a site out and predicting it from
##    its neighbours is where a barrier-vs-stationary gap should appear (most
##    sharply for sites that face their neighbours across water).
##  * YEAR hold-out (leave whole years out): tests temporal interpolation by
##    the shared AR1. Because the barrier and stationary models use the SAME
##    AR1 in time, expect only a small difference here; read it as a temporal
##    calibration check rather than a barrier discriminator.
##
## Method: INLA-idiomatic missing-response prediction. We copy the ST stack
## data, set the held-out responses to NA, refit, and read the posterior of the
## fitted values at those rows. With control.bgev$q.location = 0.5 and
## control.predictor$link = 1, summary.fitted.values$mean is the posterior mean
## of the MEDIAN annual maximum at each site-year -- the natural point
## prediction to compare against the observed annual maximum. The 0.025/0.975
## fitted quantiles are a credible interval for that median, NOT a posterior
## predictive interval for a new observation (a full predictive interval needs
## return_level_bgev2() over posterior samples; see the note at the end). The
## SAME held-out sites/years are used for the barrier and stationary refits, so
## the two are directly comparable.
##
## Only the ST models are refit (m0st = barrier, m0st_nb = stationary): the
## hold-out targets the spatial/temporal prediction the location field drives.
##==============================================================================

## Recover the site/year identity of every stacked row. The pinned geometry was
## built by reading this CSV in order and never reordering, so row k of the CSV
## matches element k of geom$y / geom$year_idx. Assert that before relying on
## it; if the data file is ever re-exported in a different order this stops the
## script rather than silently mislabelling the hold-out.
dc_id <- read_csv(
    here::here("data", "dc_max.csv"),
    show_col_types = FALSE)
stopifnot(
    nrow(dc_id) == n_data,
    identical(as.integer(dc_id$year_idx), as.integer(geom$year_idx)),
    isTRUE(all.equal(dc_id$max, geom$y)))

## Hold-out configuration. Seeds make the random site/year choice reproducible;
## the tag goes into the cache file names so a changed configuration does not
## silently reuse old fits.
holdout_seed <- 20240601L
n_sites_per_country <- 2L            # sites held out per country
n_years_holdout <- 4L                # interior years held out
min_years_site <- 10L                # only hold out well-observed sites
min_sites_year <- 20L                # only hold out well-observed years
ho_tag <- sprintf("ho_spc%d_y%d_seed%d",
    n_sites_per_country, n_years_holdout, holdout_seed)

## Observed cells only (the NA cells are predictions, never observations).
obs_id <- dc_id |>
    filter(!is.na(max))

## Count observed years per site and observed sites per year.
site_counts <- obs_id |>
    count(country, site_name, name = "n_obs")
year_counts <- obs_id |>
    count(year, name = "n_obs")

## Select held-out sites: stratified by country, restricted to sites with at
## least min_years_site observed years so there are enough points to score.
set.seed(holdout_seed)
holdout_sites <- site_counts |>
    filter(n_obs >= min_years_site) |>
    group_by(country) |>
    slice_sample(n = n_sites_per_country) |>
    ungroup() |>
    arrange(country, site_name)
stopifnot(nrow(holdout_sites) ==
    n_sites_per_country * n_distinct(site_counts$country))

## Select held-out years: interior only (never the first or last observed year,
## to avoid AR1 boundary effects), restricted to years with at least
## min_sites_year observed sites. If this stopifnot trips, lower min_sites_year.
yr_lo <- min(obs_id$year)
yr_hi <- max(obs_id$year)
set.seed(holdout_seed + 1L)
holdout_years <- year_counts |>
    filter(year > yr_lo, year < yr_hi, n_obs >= min_sites_year) |>
    slice_sample(n = n_years_holdout) |>
    arrange(year)
stopifnot(nrow(holdout_years) == n_years_holdout)

## Held-out masks over the stacked rows (observed cells in the chosen sites /
## years). Computed ONCE and shared by the barrier and stationary refits so the
## comparison is on identical held-out cells. Sites are matched on the
## country+name pair in case a station name is reused across countries.
ho_site_keys <- paste(holdout_sites$country, holdout_sites$site_name)
is_holdout_site <-
    (paste(dc_id$country, dc_id$site_name) %in% ho_site_keys) & !is.na(geom$y)
is_holdout_year <- (dc_id$year %in% holdout_years$year) & !is.na(geom$y)

## Build the two reduced data sets: copy the full ST stack data and blank the
## held-out responses. Everything else (stack, projector A, indices, priors) is
## unchanged, so the latent structure the model fits is identical -- only the
## observations differ.
d_site <- st_inla_data
d_site$y <- geom$y
d_site$y[is_holdout_site] <- NA_real_

d_year <- st_inla_data
d_year$y <- geom$y
d_year$y[is_holdout_year] <- NA_real_

ho_data <- list(site = d_site, year = d_year)
ho_held <- list(site = is_holdout_site, year = is_holdout_year)

## Row indices of the est.st block within summary.fitted.values. Reading the
## fitted values at these rows puts them back in the original (geom / dc_id)
## order, so the held-out masks above line up directly.
idx_st <- inla.stack.index(stk.st, "est.st")$data

## The ST models to refit, paired with their fms entry (formula + start mode).
ho_models <- tibble(
    model = c("barrier", "stationary"),
    spde_key = c("m0st", "m0st_nb"))

## One row per (experiment, model). Each row is a single reduced-data refit.
ho_specs <- expand_grid(
    experiment = c("site", "year"),
    ho_models)

## Run the hold-out refits, caching each one. A failed fit is skipped (warned)
## rather than aborting the whole loop.
ho_preds <- vector("list", nrow(ho_specs))
for(j in seq_len(nrow(ho_specs))) {
    ex <- ho_specs$experiment[j]
    key <- ho_specs$spde_key[j]
    mlab <- ho_specs$model[j]
    spec <- fms[[key]]
    held <- ho_held[[ex]]
    fnm_ho <- file.path(rds_dir,
        sprintf("holdout_%s_%s_%s_%s.rds", ex, key, ho_tag, mv))
    message(sprintf("## Hold-out: experiment = %s, model = %s (%s) ====",
        ex, mlab, key))
    if(file.exists(fnm_ho)) {
        fit_ho <- readRDS(fnm_ho)
    } else {
        fit_ho <- try(
            fit_bgev(
                form = spec$form,
                data = ho_data[[ex]],
                mode = spec$mode,
                pred = inla_pred_st,
                compute = inla_control_compute_pred),
            silent = FALSE)
        if(inherits(fit_ho, "try-error")) {
            warning(sprintf(
                "Hold-out fit failed: %s / %s; skipping.", ex, mlab))
            next
        }
        ## Only summary.fitted.values is read below, so prune the rest before
        ## caching (config is already off; this drops any residual marginals).
        if(prune_fits) fit_ho <- prune_inla(fit_ho)
        saveRDS(fit_ho, fnm_ho)
    }
    ff <- fit_ho$summary.fitted.values[idx_st, ]
    ho_preds[[j]] <- tibble(
        experiment = ex,
        model = mlab,
        country = dc_id$country[held],
        site_name = dc_id$site_name[held],
        year = dc_id$year[held],
        observed = dc_id$max[held],
        predicted = ff[["mean"]][held],
        pred_lo = ff[["0.025quant"]][held],
        pred_hi = ff[["0.975quant"]][held])
}
ho_pred <- ho_preds |>
    compact() |>
    list_rbind()

## If every hold-out refit failed, ho_pred is an empty 0-column tibble; give it
## the expected columns so the summaries below return empty results (with a
## message) instead of erroring on missing grouping columns.
if(nrow(ho_pred) == 0) {
    message("## All hold-out refits failed; hold-out metrics are empty.")
    ho_pred <- tibble(
        experiment = character(), model = character(),
        country = character(), site_name = character(),
        year = integer(), observed = double(), predicted = double(),
        pred_lo = double(), pred_hi = double())
}

## Point-prediction skill per experiment and model. cor() is over the held-out
## cells; bias is signed (predicted - observed).
ho_metrics <- ho_pred |>
    group_by(experiment, model) |>
    summarise(
        n = n(),
        bias = mean(predicted - observed),
        rmse = sqrt(mean((predicted - observed)^2)),
        mae = mean(abs(predicted - observed)),
        cor = cor(predicted, observed),
        .groups = "drop")

## Paired barrier-vs-stationary comparison on identical held-out cells. A
## positive diff_abserr means the barrier model was closer on that cell. The
## pivot needs both model columns, so guard against a failed refit leaving one
## (or both) absent: in that case skip the comparison rather than crash.
ho_paired <- ho_pred |>
    select(experiment, country, site_name, year, observed, model, predicted) |>
    pivot_wider(names_from = model, values_from = predicted)

if(all(c("barrier", "stationary") %in% names(ho_paired))) {
    ho_paired <- ho_paired |>
        mutate(
            abserr_barrier = abs(barrier - observed),
            abserr_stationary = abs(stationary - observed),
            diff_abserr = abserr_stationary - abserr_barrier)
    ho_paired_summary <- ho_paired |>
        group_by(experiment) |>
        summarise(
            n = n(),
            mae_barrier = mean(abserr_barrier, na.rm = TRUE),
            mae_stationary = mean(abserr_stationary, na.rm = TRUE),
            mean_diff = mean(diff_abserr, na.rm = TRUE),
            pct_barrier_better = mean(diff_abserr > 0, na.rm = TRUE),
            .groups = "drop")
} else {
    message("## Paired comparison skipped: both barrier and stationary ",
        "hold-out fits are required.")
    ho_paired_summary <- tibble()
}

## Persist the hold-out results alongside the main diagnostics.
fl <- file.path(rds_dir, sprintf("holdout_results_%s.rds", env_info$node))
if(!file.exists(fl)) {
    holdout <- list(
        config = list(
            seed = holdout_seed,
            n_sites_per_country = n_sites_per_country,
            n_years_holdout = n_years_holdout,
            min_years_site = min_years_site,
            min_sites_year = min_sites_year,
            holdout_sites = holdout_sites,
            holdout_years = holdout_years,
            tag = ho_tag),
        pred = ho_pred,
        metrics = ho_metrics,
        paired = ho_paired,
        paired_summary = ho_paired_summary)
    saveRDS(holdout, fl)
} else {
    holdout <- readRDS(fl)
}

## Print the hold-out summaries when the script is source()d.
message("## Held-out sites:")
print(holdout_sites)
message("## Held-out years:")
print(holdout_years)
message("## Hold-out skill (per experiment and model):")
print(ho_metrics)
message("## Barrier vs stationary on identical held-out cells:")
print(ho_paired_summary)

## Observed-vs-predicted plot (kept commented: this script does not draw).
ggplot(ho_pred, aes(observed, predicted)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_linerange(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.3) +
    geom_point(alpha = 0.6) +
    facet_grid(experiment ~ model) +
    labs(x = "Observed annual maximum",
        y = "Predicted median annual maximum")

## Posterior-predictive extension (optional, heavier). The fitted-value
## interval above is for the MEDIAN. For a calibrated predictive interval of a
## NEW annual maximum at a held-out cell, refit with config = TRUE, draw with
## inla.posterior.sample(), map each sample's (location, spread, tail) hypers
## through return_level_bgev2() in bgev.R, and take predictive quantiles. That
## turns this point-prediction check into a coverage check.
