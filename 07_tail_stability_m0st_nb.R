##==============================================================================
## Does the selected model's tail posterior reproduce? (TODO item 7,
## diagnostic B)
##
## This is the question that decides whether the tail issue touches the
## paper at all. Everything found so far concerns a small pooled model
## with no latent field. The selected model is m0st_nb, and its return
## levels are a headline result of the revision, built by pushing joint
## posterior draws of the tail through the bGEV return-level function.
## If its tail posterior does not reproduce, those intervals do not
## either.
##
## What is already known, and why a clean answer is expected:
##
##   model            tail CrI                width   internal mode
##   pooled (05, 06)  (~0.002, ~0.499)        0.497   -5.01
##   m0st             (0.0148, 0.0763)        0.062   -2.28
##   m0st_nb          (0.0182, 0.0643)        0.046   -2.28
##
## The pooled model sits against the xi -> 0 boundary, where the internal
## transform logit(2 * xi) degenerates and the log-posterior goes flat.
## The spatio-temporal models sit well clear of it and identify the tail
## an order of magnitude better. 06_tail_integration_check.R showed the
## pooled model's MODE is perfectly stable and only the integration over
## the flat direction wanders, so a model with real curvature there
## should be reproducible. This script is the confirmation, not a fishing
## expedition.
##
## The design. Refit m0st_nb n_refit times under two starting regimes:
##
##   warm  control.mode set to the cached fit's theta, restart = TRUE.
##         This is exactly what 01_explore_models.R does, so it tests
##         the pipeline as it actually runs.
##   cold  no control.mode at all, so INLA finds its own mode from the
##         prior initial values. This is what the pooled fits in 05 did,
##         and it is the harsher and more honest test: a warm start could
##         make the answer reproducible for trivial reasons.
##
## Everything else matches 01_explore_models.R exactly.
##
## COST. m0st_nb takes about 61 minutes of INLA CPU. The defaults below
## are 3 refits x 2 regimes = 6 fits, so roughly 6 hours. To shorten it,
## drop `starts` to "warm" only (halves it) or set n_refit to 2. Group
## cross-validation and config sampling are switched OFF: they are the
## expensive extras in 01 and nothing here needs them.
##
## This script does NOT write to the rds model cache. It writes only
## its own diagnostic files, so it cannot disturb the cached fits that
## the paper reads.
##
## Inputs:  rds/pinned_geometry.rds
##          rds/explore_dc_max_bGEV_m0st_nb_mesh_<n>.rds (warm start)
## Outputs: rds/explore_tailstab_runs.rds
##          rds/explore_tailstab_summary.rds
##          figs/tailstab_marginals.{png,pdf}
##
## Read the decision rule at the end of this file before interpreting.
##==============================================================================

##==============================================================================
## Setup

here::i_am("07_tail_stability_m0st_nb.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(sf)
    library(INLA)
})

library(stExtremes)

## Threading pinned exactly as in 01_explore_models.R. B = 1 is what
## makes the per-theta linear algebra deterministic; A only affects
## speed. If any variation shows up below it is therefore NOT threading.
inla.setOption(num.threads = "8:1")
if(requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
}

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

col_ink <- "grey15"
col_bar <- "grey45"

## Run size. See the COST note in the header before raising these.
n_refit <- 3
starts <- c("warm", "cold")

##==============================================================================
## Geometry and model, mirroring 01_explore_models.R
##
## The values below are copied from that script and MUST stay in step
## with it. If any prior or control there changes, this diagnostic stops
## testing the model the paper actually fits.

geom <- readRDS(file.path(rds_dir, "pinned_geometry.rds"))
mesh <- geom$mesh
n_data <- length(geom$y)
mesh_n <- mesh$n

prior_sigma <- c(1, 0.01)
prior_range <- c(50, 0.01)
stationary_model <- inla.spde2.pcmatern(mesh, prior.range = prior_range,
    prior.sigma = prior_sigma, constr = TRUE)

prior_rho <- list(theta = list(prior = "pc.cor0", param = c(0.5, 0.5)))

hyper_spread <- list(initial = 1, fixed = FALSE, prior = "loggamma",
    param = c(3, 3))
tail_interval <- c(0, 0.5)
hyper_tail <- list(initial = map_tail(0.1, tail_interval, inverse = TRUE),
    prior = "pc.gevtail", param = c(7, tail_interval), fixed = FALSE)
hyper_bgev <- list(spread = hyper_spread, tail = hyper_tail)

inla_control_bgev <- list(q.location = 0.5, q.spread = 0.5,
    q.mix = c(0.05, 0.2), beta.ab = 5)
inla_control_fixed <- list(prec.intercept = 100, prec = 100)
inla_control_inla <- list(int.strategy = "ccd", strategy = "laplace",
    cmin = 0, h = 0.005, tolerance = 1e-8, restart = 1L)

w.index <- create_indices(geom$year_idx, mesh)

stk.st <- inla.stack(
    data = list(y = geom$y),
    effects = list(w.index, data.frame(m = rep(1, n_data))),
    A = list(geom$Ast, 1),
    tag = "est.st")

st_inla_data <- inla.stack.data(stk.st)
st_inla_data$spread_x <- st_inla_data$tail_x <-
    matrix(nrow = n_data, ncol = 0)
inla_pred_st <- list(A = inla.stack.A(stk.st), compute = TRUE, link = 1)

form_m0st_nb <- inla.mdata(y, spread_x, tail_x) ~ -1 + m +
    f(w, model = stationary_model, group = w.group, vb.correct = FALSE,
        control.group = list(model = "ar1", hyper = prior_rho))

## Starting mode for the warm regime, taken from the cached fit rather
## than hard-coded, so it cannot drift away from what 01 produced.
cached_file <- file.path(rds_dir,
    sprintf("explore_dc_max_bGEV_m0st_nb_mesh_%d.rds", mesh_n))
stopifnot(file.exists(cached_file))
cached <- readRDS(cached_file)
theta_warm <- cached$mode$theta

## CPO, config and predictor marginals are all off: this diagnostic only
## needs the hyperparameter posterior, and those three are what make the
## spatio-temporal fits slow and heavy.
inla_control_compute <- list(waic = TRUE, dic = TRUE, cpo = FALSE,
    config = FALSE, return.marginals.predictor = FALSE)

fit_once <- function(start) {
    ctrl_mode <- if(start == "warm") {
        list(theta = theta_warm, restart = TRUE)
    } else {
        list(restart = TRUE)
    }
    inla(formula = form_m0st_nb, data = st_inla_data, family = "bgev",
        control.family = list(hyper = hyper_bgev,
            control.bgev = inla_control_bgev),
        control.predictor = inla_pred_st,
        control.fixed = inla_control_fixed,
        control.compute = inla_control_compute,
        control.inla = inla_control_inla,
        control.mode = ctrl_mode,
        verbose = FALSE,
        safe = TRUE)
}

##==============================================================================
## Run

row_tail <- "tail for BGEV observations"
row_spread <- "spread for BGEV observations"

intern_tail <- function(f) {
    th <- f$mode$theta
    i <- grep("tail", names(th))
    if(length(i) == 0) return(c(NA_real_, NA_real_))
    sdv <- NA_real_
    if(!is.null(f$misc$cov.intern)) sdv <- sqrt(diag(f$misc$cov.intern))[i]
    c(unname(th[i]), unname(sdv))
}

specs <- expand_grid(start = starts, refit = seq_len(n_refit))

out <- pmap(specs, \(start, refit) {
    message(sprintf("## m0st_nb, %s start, refit %d of %d",
        start, refit, n_refit))
    t0 <- Sys.time()
    f <- try(fit_once(start), silent = TRUE)
    if(inherits(f, "try-error")) {
        return(list(row = tibble(start = start, refit = refit,
            failed = TRUE), marg = NULL))
    }
    it <- intern_tail(f)
    row <- tibble(start = start, refit = refit, failed = FALSE,
        spread = f$summary.hyperpar[row_spread, "mean"],
        tail = f$summary.hyperpar[row_tail, "mean"],
        tail_sd = f$summary.hyperpar[row_tail, "sd"],
        tail_lo = f$summary.hyperpar[row_tail, "0.025quant"],
        tail_hi = f$summary.hyperpar[row_tail, "0.975quant"],
        waic = f$waic$waic,
        tail_intern_mode = it[1], tail_intern_sd = it[2],
        secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
    ## Keep the marginal itself: two fits can share a posterior mean and
    ## still have different distributions, which the summary would hide.
    marg <- as_tibble(f$marginals.hyperpar[[row_tail]]) |>
        mutate(start = start, refit = refit)
    list(row = row, marg = marg)
})

runs <- map(out, "row") |> list_rbind()
margs <- map(out, "marg") |> compact() |> list_rbind()

runs |>
    mutate(across(where(is.numeric), \(x) round(x, 4))) |>
    print(n = 20)
saveRDS(list(runs = runs, margs = margs,
    cached_tail = cached$summary.hyperpar[row_tail, ]),
    file.path(rds_dir, "explore_tailstab_runs.rds"))

##==============================================================================
## Summary
##
## The reference is the cached fit that the paper currently reports:
## tail mean 0.0400, CrI (0.0182, 0.0643). Refits should land on it.

stab <- runs |>
    filter(!failed) |>
    group_by(start) |>
    summarise(n = n(),
        tail_min = min(tail), tail_max = max(tail),
        tail_range = max(tail) - min(tail),
        ci_lo_range = max(tail_lo) - min(tail_lo),
        ci_hi_range = max(tail_hi) - min(tail_hi),
        ci_width_med = median(tail_hi - tail_lo),
        spread_range = max(spread) - min(spread),
        waic_range = max(waic) - min(waic),
        intern_mode_med = median(tail_intern_mode),
        secs_med = median(secs),
        .groups = "drop") |>
    mutate(across(where(is.numeric), \(x) round(x, 4)))
stab
saveRDS(stab, file.path(rds_dir, "explore_tailstab_summary.rds"))

##==============================================================================
## Figure: the tail marginals themselves, overlaid

p_marg <- margs |>
    ggplot(aes(x = x, y = y, group = interaction(start, refit))) +
    geom_line(colour = col_bar, linewidth = 0.4) +
    geom_vline(xintercept = cached$summary.hyperpar[row_tail, "mean"],
        colour = col_ink, linewidth = 0.4, linetype = "dashed") +
    facet_wrap(~ start) +
    labs(x = expression(paste("Tail ", xi)), y = "Posterior density")
p_marg
save_fig(p_marg, "tailstab_marginals", w = 9, h = 4.5)

message("07_tail_stability_m0st_nb.R complete.")

##==============================================================================
## DECISION RULE, written before the run
##
## Read `stab`, and look at the overlaid marginals rather than only the
## summaries.
##
## (a) tail_range small (say under 0.005) in BOTH regimes, marginals
##     visually indistinguishable, and all runs near the cached mean of
##     0.0400.
##     The issue is closed. The instability belongs to the pooled
##     diagnostic model and does not touch the paper. Record that in
##     05_explore_spread_shape_inla.R and in the response letter, keep
##     the return levels as they stand, and treat item 6.3 as an ordinary
##     prior-sensitivity question on an identified parameter.
##
## (b) warm stable, cold not.
##     The pipeline is reproducible because 01 always warm-starts from a
##     cached mode, but the posterior is not intrinsically well
##     determined. That is acceptable for the paper as long as it is
##     stated: the reported fit is conditional on a pinned starting mode.
##     Say so in the reproducibility appendix rather than leaving it
##     implicit.
##
## (c) both unstable, but the marginals overlap heavily and only the mean
##     moves.
##     Same reading as the pooled model: report the interval and the
##     mode, never the mean, and check whether the return-level intervals
##     move materially, since they integrate over the whole marginal
##     rather than using its mean.
##
## (d) both unstable AND the marginals differ in shape or location.
##     This is the serious branch. The return levels and their credible
##     intervals are built from joint posterior draws that include the
##     tail, so they inherit the problem. Item 6.3 would then have to
##     come first and the return-level uncertainty would need to be
##     re-examined before anything is resubmitted. Do not paper over it
##     by pinning the mode.
##
## In every branch, compare intern_mode_med against the pooled model's
## -5.01. If it stays near -2.28 the model is clear of the boundary and
## any variation has some other cause, which is worth knowing before
## reaching for a fix.
##==============================================================================
