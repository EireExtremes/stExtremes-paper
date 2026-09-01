##==============================================================================
## Tail prior sensitivity (TODO item 6.3; Referee 1, point 3)
##
## Referee 1: "I did not find any sensitivity analysis exploring
## alternative priors for the tail parameter, the robustness of return
## levels and trends to prior assumptions, or the consequences of
## assuming a common tail parameter across all sites."
##
## Three things are being asked. This script answers the first two. The
## third is already answered by 05_explore_spread_shape_inla.R, where
## letting the tail vary gauge by gauge is 21 WAIC WORSE than holding it
## common, so a shared tail is supported rather than merely assumed. Do
## not redo it here; cite it.
##
## The question is well posed only because of what came before. In the
## pooled diagnostic model the tail is not identified at all: its
## credible interval spans the whole admissible range under every
## integration setting tested (06_tail_integration_check.R), so asking
## about prior sensitivity there would just be measuring the prior. The
## spatio-temporal models are different. m0st_nb has a tail interval of
## (0.0182, 0.0643), an order of magnitude tighter, and sits clear of the
## xi -> 0 boundary. On that model the question has an answer worth
## having.
##
## WHAT IS VARIED. The tail carries a PC prior, pc.gevtail(lambda, low,
## high), whose base model is the Gumbel case xi = 0. lambda sets how
## hard the prior shrinks towards it: larger means more shrinkage. The
## production setting is lambda = 7 on (0, 0.5).
##
##   lambda = 1    very weak shrinkage, heavy tails barely penalised
##   lambda = 3    weak
##   lambda = 7    production, the reference run
##   lambda = 15   strong shrinkage towards Gumbel
##
## plus one interval variant, lambda = 7 on (0, 1.0), to show that the
## upper bound is not doing any work: the posterior sits near 0.04, far
## below either ceiling. It is wrapped in a guard because INLA may not
## accept every interval.
##
## The spread prior is deliberately held fixed. The referee asked about
## the tail, and varying two priors at once would make it impossible to
## attribute any movement. Add a loggamma variant later if wanted.
##
## WHAT IS MEASURED. Not the hyperparameter. The referee asks about the
## robustness of RETURN LEVELS, so the comparison is made on the same
## quantity the paper reports: joint posterior draws of the latent field,
## the intercept, the spread and the tail pushed through
## return_level_bgev2(), exactly as 04_posterior_predictive.R does.
##
## The headline number is a RATIO, not a difference: the shift in the
## 50-year return level across priors, divided by the width of its 95%
## credible interval under the production prior. A shift that is small
## relative to the uncertainty already reported is the definition of
## robust, and it is a fairer statement than any absolute tolerance.
##
## PLANNED PLACEMENT of the outputs (each block is tagged below):
##
##   [PAPER]    one short paragraph in the results, plus fig
##              tailprior_rlevel_shift: the 50-year return level at the
##              gauges under each prior, with the production credible
##              band drawn behind it. One figure carries the whole
##              message: the bands overlap, the medians barely move.
##   [APPENDIX] tbl tailprior_hyper (spread and tail posteriors under
##              each prior) and tbl tailprior_rlevel (return levels at
##              all gauges for 10, 50 and 100 years). The detail a
##              sceptical reader wants, out of the main text.
##   [LETTER]   the ratio above, as a single sentence answering Referee 1
##              point 3, together with the pointer to item 3 for the
##              common-tail half of the question.
##
## COST. Each fit is a full re-estimation of m0st_nb, about 61 minutes,
## because a new prior moves the mode and it cannot be reused. Five
## settings is therefore roughly 5 hours plus sampling. The runs are
## cached one at a time, so the script can be stopped and resumed.
##
## PREREQUISITE. 01_explore_models.R must have completed, since this
## sources it for the model spec and the cached mode.
##
## Read the decision rule at the end of this file before interpreting.
##==============================================================================

##==============================================================================
## Setup

here::i_am("08_tail_prior_sensitivity.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(sf)
    library(INLA)
})

library(stExtremes)

rds_dir <- here::here("rds")
figs_dir <- here::here("figs")
if(!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)
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

## Same sampling settings as 04, so the return levels here and there are
## directly comparable rather than merely similar.
n_samples <- 1000L
pp_seed <- 20240601L
periods <- c(10, 50, 100)
ref_period <- 50

## The model under test. m0st_nb is the selected model and the one whose
## return levels the paper reports. Adding m0st doubles the runtime.
key_model <- "m0st_nb"

## Prior settings. The first is the production one and is the reference
## against which everything else is measured.
prior_grid <- tibble(
    label = c("lambda = 7 (production)", "lambda = 1", "lambda = 3",
        "lambda = 15", "lambda = 7, interval (0, 1)"),
    lambda = c(7, 1, 3, 15, 7),
    low = c(0, 0, 0, 0, 0),
    high = c(0.5, 0.5, 0.5, 0.5, 1.0)) |>
    mutate(tag = sprintf("l%g_h%g", lambda, high))
ref_label <- "lambda = 7 (production)"

##==============================================================================
## Model spec, taken from 01 so it cannot drift
##
## Sourcing 01 reuses its deterministic INLA setup (threads, BLAS, ccd),
## the pinned geometry, the priors, the stacks and the control lists, and
## loads the cached fits into lmn. If the fit caches are absent this
## triggers 01's full fit, so run 01 first.

message("## Sourcing 01_explore_models.R for the fitting setup ...")
source(here::here("01_explore_models.R"))

alpha <- inla_control_bgev$q.location
beta <- inla_control_bgev$q.spread
p_a <- inla_control_bgev$q.mix[1]
p_b <- inla_control_bgev$q.mix[2]
s_blend <- inla_control_bgev$beta.ab

n_years <- length(unique(geom$year_idx))
last_idx <- max(geom$year_idx)

## Gauge locations in km, built exactly as in 04_posterior_predictive.R
## and 02_compare_barrier.R. 01 does not define these: it only needs
## the projector matrices stored in the pinned geometry, not the site
## coordinates themselves.
dc_max_csv <- read_csv(
    here::here("data", "dc_max.csv"),
    show_col_types = FALSE)
stopifnot(isTRUE(all.equal(dc_max_csv$max, geom$y)))
sites_km <- st_as_sf(dc_max_csv, coords = c("lon", "lat"), crs = proj) |>
    st_transform(kmproj) |>
    mutate(x = st_coordinates(geometry)[, 1],
        y = st_coordinates(geometry)[, 2]) |>
    st_drop_geometry() |>
    distinct(country, site_name, x, y)

A_sites <- inla.spde.make.A(mesh,
    loc = as.matrix(sites_km[, c("x", "y")]))
n_site <- nrow(sites_km)
stopifnot(ncol(A_sites) == mesh$n)

##==============================================================================
## Fit under one prior and reduce to return levels

## Rebuild hyper_bgev with a different tail prior, leaving the spread
## prior exactly as 01 set it.
make_hyper <- function(lambda, low, high) {
    list(spread = hyper_spread,
        tail = list(initial = map_tail(0.1, c(low, high), inverse = TRUE),
            prior = "pc.gevtail", param = c(lambda, low, high),
            fixed = FALSE))
}

## A new prior moves the mode, so restart = TRUE and a real optimisation.
## The cached mode is still supplied as a starting point, which is only a
## speed-up and does not pin the answer.
fit_under_prior <- function(lambda, low, high) {
    inla(formula = fms[[key_model]]$form,
        data = st_inla_data,
        family = inla_family,
        control.family = list(hyper = make_hyper(lambda, low, high),
            control.bgev = inla_control_bgev),
        control.predictor = inla_pred_st,
        control.fixed = inla_control_fixed,
        control.compute = inla_control_compute,
        control.inla = inla_control_inla,
        control.mode = list(theta = lmn[[key_model]]$mode$theta,
            restart = TRUE),
        verbose = inla_verbose,
        safe = inla_safe)
}

## Same extraction as 04: field, intercept, spread and tail per draw.
extract_draws <- function(smp) {
    ln <- rownames(smp[[1]]$latent)
    w_rows <- grep("^w:", ln)
    m_row <- grep("^m:", ln)
    stopifnot(length(w_rows) == mesh$n * n_years, length(m_row) == 1L)
    hn <- names(smp[[1]]$hyperpar)
    sp_nm <- grep("spread for BGEV", hn, value = TRUE)
    tl_nm <- grep("^tail for BGEV", hn, value = TRUE)
    list(W = vapply(smp, \(s) s$latent[w_rows, 1],
            numeric(length(w_rows))),
        m = vapply(smp, \(s) s$latent[m_row, 1], numeric(1L)),
        spread = vapply(smp, \(s) s$hyperpar[[sp_nm]], numeric(1L)),
        tail = vapply(smp, \(s) s$hyperpar[[tl_nm]], numeric(1L)))
}

## Field for the last year: rows of W are (node, group) stacked by group.
field_last_year <- function(W) {
    idx <- ((last_idx - 1) * mesh$n + 1):(last_idx * mesh$n)
    W[idx, , drop = FALSE]
}

## One prior -> hyperparameter summary + per-site return-level summary.
run_one <- function(label, lambda, low, high, tag) {
    cache <- file.path(rds_dir,
        sprintf("tailprior_%s_%s.rds", key_model, tag))
    if(file.exists(cache)) {
        message(sprintf("## Cached: %s", label))
        return(readRDS(cache))
    }
    message(sprintf("## Fitting %s under %s ...", key_model, label))
    fit <- try(fit_under_prior(lambda, low, high), silent = FALSE)
    if(inherits(fit, "try-error")) {
        warning(sprintf("%s failed; skipped.", label))
        return(NULL)
    }
    hyper <- as_tibble(fit$summary.hyperpar, rownames = "param") |>
        filter(str_detect(param, "BGEV")) |>
        transmute(label = label, param, mean, sd,
            lo = `0.025quant`, hi = `0.975quant`)
    smp <- inla.posterior.sample(n_samples, fit, seed = pp_seed)
    dr <- extract_draws(smp)
    rm(smp)
    ## q_alpha at every gauge for every draw, then return levels. One
    ## vectorised call per (draw, period) since spread and tail are
    ## scalar within a draw.
    u_sites <- as.matrix(A_sites %*% field_last_year(dr$W))
    q_draw <- sweep(u_sites, 2, dr$m, "+")
    rl <- map(periods, \(pp) {
        m_rl <- vapply(seq_len(n_samples), \(k)
            return_level_bgev2(pp, q = q_draw[, k], sb = dr$spread[k],
                xi = dr$tail[k], alpha = alpha, beta = beta,
                p_a = p_a, p_b = p_b, s = s_blend),
            numeric(n_site))
        tibble(label = label, period = pp,
            site = sites_km$site_name,
            med = apply(m_rl, 1, median),
            lo = apply(m_rl, 1, quantile, 0.025),
            hi = apply(m_rl, 1, quantile, 0.975))
    }) |>
        list_rbind()
    out <- list(hyper = hyper, rlevel = rl)
    saveRDS(out, cache)
    rm(fit, dr)
    gc()
    out
}

res <- pmap(prior_grid, run_one) |> compact()

hyper_tab <- map(res, "hyper") |> list_rbind()
rlevel_tab <- map(res, "rlevel") |> list_rbind()
saveRDS(list(hyper = hyper_tab, rlevel = rlevel_tab),
    file.path(rds_dir, "tailprior_summary.rds"))

##==============================================================================
## 1. Hyperparameter posteriors under each prior [APPENDIX]

hyper_tab |>
    mutate(across(where(is.numeric), \(x) round(x, 4))) |>
    arrange(param, label) |>
    print(n = 40)

##==============================================================================
## 2. The headline robustness number [PAPER + LETTER]
##
## For the reference return period, at every gauge: how far does the
## posterior median move across priors, relative to the width of the 95%
## credible interval under the production prior? A ratio below about 0.1
## means the prior moves the answer by less than a tenth of the
## uncertainty already reported.

ref_band <- rlevel_tab |>
    filter(label == ref_label, period == ref_period) |>
    transmute(site, ref_med = med, ref_width = hi - lo)

shift_tab <- rlevel_tab |>
    filter(period == ref_period) |>
    left_join(ref_band, by = "site") |>
    mutate(shift = abs(med - ref_med), rel = shift / ref_width) |>
    group_by(label) |>
    summarise(max_shift_m = max(shift), med_shift_m = median(shift),
        max_rel = max(rel), med_rel = median(rel), .groups = "drop") |>
    mutate(across(where(is.numeric), \(x) round(x, 4))) |>
    arrange(max_rel)
shift_tab
saveRDS(shift_tab, file.path(rds_dir, "tailprior_shift.rds"))

##==============================================================================
## 3. Return-level table, all gauges and periods [APPENDIX]

rlevel_tab |>
    mutate(across(c(med, lo, hi), \(x) round(x, 3))) |>
    arrange(period, site, label) |>
    print(n = 30)

##==============================================================================
## 4. The figure that carries the message [PAPER]
##
## Gauges ordered by their production median. The grey band is the
## production 95% credible interval; the points are the medians under
## each prior. If the points sit inside the band the answer is robust,
## and that is legible at a glance without reading a single number.

plot_dat <- rlevel_tab |>
    filter(period == ref_period) |>
    left_join(ref_band, by = "site") |>
    mutate(site = fct_reorder(site, ref_med))

band <- plot_dat |>
    filter(label == ref_label) |>
    distinct(site, lo, hi)

p_shift <- ggplot(plot_dat, aes(x = site)) +
    geom_linerange(data = band, aes(ymin = lo, ymax = hi),
        colour = "grey80", linewidth = 2.5) +
    geom_point(aes(y = med, shape = label), size = 1.8, colour = col_ink) +
    scale_shape_manual(name = NULL, values = c(16, 1, 2, 3, 4)) +
    coord_flip() +
    labs(x = NULL,
        y = sprintf("%g-year return level (m)", ref_period)) +
    theme(legend.position = "top", axis.text.y = element_text(size = 7))
p_shift
save_fig(p_shift, "tailprior_rlevel_shift", w = 9, h = 9)

message("08_tail_prior_sensitivity.R complete.")

##==============================================================================
## DECISION RULE, written before the run
##
## Read `shift_tab`, column max_rel, and look at the figure.
##
## (a) max_rel below about 0.1 for every prior.
##     The return levels are robust to the tail prior. Write the [PAPER]
##     paragraph as a positive result, put the tables in the appendix,
##     and answer Referee 1 point 3 in two sentences: the return levels
##     move by less than a tenth of their own credible interval across a
##     fifteen-fold change in the prior rate, and a common tail is
##     supported by WAIC rather than assumed.
##
## (b) max_rel between roughly 0.1 and 0.5.
##     Real but modest sensitivity. Report it honestly, name the gauges
##     where it is largest (expect the data-sparse Atlantic ones), and
##     say the ranking of sites is unchanged even where the level moves.
##     Do NOT bury this in the appendix; a referee who asked for it will
##     look.
##
## (c) max_rel above about 0.5 anywhere.
##     The prior is doing a substantial share of the work in the
##     headline result. That has to go in the main text, and the
##     conclusions about return levels need hedging accordingly. Check
##     first whether it is confined to lambda = 1, which is a
##     deliberately extreme setting, before treating it as the finding.
##
## (d) the lambda = 7 interval (0, 1) run differs materially from the
##     production run.
##     Then the upper bound IS doing work, which would be surprising
##     given the posterior sits near 0.04. Suspect a fitting problem
##     rather than a real effect, and check the mode and the internal sd
##     before drawing any conclusion.
##
## In every branch: the hyperparameter table is the diagnostic, not the
## result. A tail posterior that tracks lambda closely while the return
## levels do not move is the expected and reportable pattern, since the
## return level at these periods depends far more on the location field
## and the spread than on the tail.
##==============================================================================

##==============================================================================
## RESULT (run 2026-07-28)
##
## Outcome (a): the return levels are robust to the tail prior.
##
##   prior         max shift (m)  median shift (m)  max shift / CrI width
##   lambda = 3        0.024          0.0055               0.028
##   lambda = 15       0.027          0.0072               0.035
##   lambda = 1        0.026          0.0041               0.039
##
## Over a fifteen-fold change in the PC prior rate, the 50-year return
## level moves by at most 2.7 cm at any gauge, and by at most 3.9% of the
## width of its own 95% credible interval. The median gauge moves under a
## centimetre. That is comfortably inside branch (a) of the rule above,
## so this can be written as a positive result.
##
## The spread is untouched, 0.266 under every prior, which is the
## expected control: the tail prior should not move it and does not.
##
## FOUR PRIORS, NOT FIVE. The interval variant, lambda = 7 on (0, 1),
## failed. INLA produced NaN values in the log-likelihood and then
## aborted on the tabulate-Qfunc assertion, twice, including the safe
## rerun. The reason is substantive rather than a bug in this script:
## widening the admissible range to 1 lets the sampler visit shapes at
## which the bGEV mean does not exist (xi >= 1) and the variance has long
## since ceased to (xi >= 1/2), so the likelihood degenerates. The
## comparison over lambda at the fixed (0, 0.5) interval is complete and
## is the part the referee asked about, so nothing is lost. Do not report
## the failed variant as evidence that the upper bound "matters"; it
## shows only that the wider interval cannot be fitted.
##
## TWO THINGS THE HYPERPARAMETER TABLE SHOWS THAT THE HEADLINE DOES NOT.
##
## 1. The tail posterior is NOT monotone in lambda, and its width varies
##    fourfold:
##
##      lambda = 1   mean 0.0427  sd 0.0180
##      lambda = 3   mean 0.0392  sd 0.0209
##      lambda = 7   mean 0.0461  sd 0.0052
##      lambda = 15  mean 0.0354  sd 0.0093
##
##    A stronger PC prior should shrink the tail towards the Gumbel case
##    monotonically. It does not, and lambda = 7 returns both the largest
##    mean and much the narrowest interval.
##
## 2. The production-prior refit here does not reproduce the cached
##    m0st_nb fit, although it is the same model, the same prior and the
##    same settings: tail mean 0.0461 against 0.0400, and a credible
##    interval of width 0.021 against 0.046.
##
##    This is the diagnostic B question arriving early and unasked, and
##    the answer is that the m0st_nb tail posterior is NOT perfectly
##    reproducible either. The assumption recorded in chats/TODO.md, that
##    the instability belongs to the pooled diagnostic model alone, is
##    therefore too strong.
##
## WHY THAT DOES NOT UNDERMINE ANYTHING REPORTED. The published 50-year
## return levels from 04_posterior_predictive.R and those from the
## production refit here differ by at most 0.020 m over the 44 gauges,
## median 0.006 m, at most 3.1% of the published credible-interval width.
## So the tail posterior wobbles and the return levels do not follow it.
## Both facts should be stated together: the second is what makes the
## first tolerable, and reporting either alone would mislead.
##
## The reason is structural. At a 50-year return period the return level
## is governed mainly by the location field and the spread, both of which
## are stable here, and the tail contributes a small correction. That is
## also why the prior sensitivity above is so mild.
##==============================================================================
