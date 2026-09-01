##==============================================================================
## Exploratory analysis of site-wise spread and shape (TODO item 3)
##
## Asks whether the constant spread and constant shape assumption of the
## fitted models is supported by the data. The fitted models carry
## spatio-temporal structure on the location quantile only. The spread
## s_beta and the tail xi are single numbers shared by all 44 gauges.
##
## Each block is tagged:
##
##   [PAPER]    -- candidate for the manuscript
##   [APPENDIX] -- candidate for an appendix / supplement
##   [LETTER]   -- numbers quoted in the response letter only
##
## Referee anchors (see chats/TODO.md, item 3):
##  - Referee 2, major point 1: the referee would be "very surprised" if
##    constant s_beta and xi were not obviously violated, and notes there
##    is not even an exploratory analysis.
##  - Referee 1, model validation: the same structural assumption, raised
##    from the validation side.
##
## Method. Independent GEV fits by maximum likelihood at every gauge with
## enough years, then a likelihood-ratio decomposition that asks whether
## the between-site variation exceeds sampling noise. The naive plot of
## per-site estimates is NOT enough on its own: with 15 to 58 years per
## gauge the estimates scatter widely even when the true values are
## identical, so a formal comparison is needed before claiming the
## assumption is violated.
##
## Why GEV and not bGEV. The two differ only below the mixing quantile
## p_a = 0.05, so the spread and the shape are effectively the same
## quantity in both. The package bGEV likelihood returns 1e10 for xi < 0,
## which is a hard barrier: the optimiser then parks on xi = 0 with a
## singular Hessian and no standard error. The unconstrained GEV shows
## where the data actually put the shape, which is the whole question
## here. Estimates are reported in the (q_alpha, s_beta, xi)
## parametrisation of the paper via old_to_new(), with alpha = 0.5 and
## beta = 0.5 to match control.bgev in 01_explore_models.R.
##
## Inputs:
##  - data/dc_max.csv
##  - rds/explore_dc_max_bGEV_m0st_nb_mesh_*.rds (model comparison
##    overlay; skipped with a message when absent)
##
## Outputs: figures in figs/ (png + pdf), tables in rds/ under
## the explore_ss_ prefix.
##
## Runtime about 4 minutes, all of it in the profile likelihoods of
## block 2, which refit every gauge at each step of the outer search. No
## INLA fitting happens here; the model is only read for the overlay.
##
## The conclusion is written out at the end of this file.
##==============================================================================

##==============================================================================
## Setup

here::i_am("05_explore_spread_shape.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(sf)
})

## The stExtremes package provides nllik_gev(), old_to_new(), marginals_st_gev()
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

## Same neutral palette as 00_explore_data.R
col_ink <- "grey15"
col_bar <- "grey45"

## bGEV settings of the fitted models, from 01_explore_models.R
alpha_loc <- 0.5    # control.bgev$q.location
beta_spr <- 0.5     # control.bgev$q.spread

## A gauge needs enough years before a three-parameter fit means
## anything. The shape is the first thing to become unidentifiable, so
## the threshold is deliberately not low. 15 years keeps 32 of the 44
## gauges and, importantly, keeps both countries represented.
min_years <- 15

##==============================================================================
## Data

dc_max <- read_csv(here::here("data/dc_max.csv"),
    show_col_types = FALSE) |>
    filter(!is.na(max))

sites <- dc_max |>
    count(country, site_name, name = "n_years") |>
    left_join(distinct(dc_max, site_name, lat, lon), by = "site_name")

fit_sites <- sites |>
    filter(n_years >= min_years)

series <- dc_max |>
    filter(site_name %in% fit_sites$site_name) |>
    group_by(site_name) |>
    summarise(y = list(max), .groups = "drop") |>
    deframe()

##==============================================================================
## Likelihood helpers
##
## nllik_gev() from the package, wrapped so that an invalid scale or a
## support violation returns a large finite value instead of Inf or NaN,
## which Nelder-Mead cannot recover from.

nll_site <- function(par, x) {
    if(par[2] <= 0) return(1e10)
    v <- nllik_gev(par, x)
    if(!is.finite(v)) 1e10 else v
}

start_par <- function(x) {
    c(mean(x) - 0.45 * sd(x), max(sd(x) * 0.78, 1e-3), 0.05)
}

## Two passes: Nelder-Mead is restarted from its own optimum, which is
## the standard guard against a premature stop on a flat surface.
fit_gev_site <- function(x) {
    o <- optim(start_par(x), nll_site, x = x, hessian = TRUE,
        control = list(maxit = 10000, reltol = 1e-12))
    optim(o$par, nll_site, x = x, hessian = TRUE,
        control = list(maxit = 10000, reltol = 1e-12))
}

##==============================================================================
## 1. Per-site fits [APPENDIX + LETTER]

fits <- map(series, fit_gev_site)

se_of <- function(o) {
    s <- tryCatch(sqrt(diag(solve(o$hessian))), error = \(e) rep(NA_real_, 3))
    ifelse(is.finite(s), s, NA_real_)
}

site_tab <- imap(fits, \(o, nm) {
    np <- old_to_new(o$par, alpha = alpha_loc, beta = beta_spr)
    s <- se_of(o)
    tibble(site_name = nm, conv = o$convergence, nll = o$value,
        mu = o$par[1], sigma = o$par[2], xi = o$par[3],
        se_sigma = s[2], se_xi = s[3],
        q_alpha = np$q, s_beta = np$s)
}) |>
    list_rbind() |>
    left_join(fit_sites, by = "site_name") |>
    ## s_beta is a fixed multiple of sigma at fixed xi, so its standard
    ## error is rescaled by the same factor
    mutate(se_s_beta = se_sigma * s_beta / sigma,
        lo_xi = xi - 1.96 * se_xi, hi_xi = xi + 1.96 * se_xi,
        lo_s = s_beta - 1.96 * se_s_beta,
        hi_s = s_beta + 1.96 * se_s_beta)

stopifnot(all(site_tab$conv == 0))
site_tab |>
    select(country, site_name, n_years, s_beta, se_s_beta, xi, se_xi) |>
    arrange(xi) |>
    print(n = 40)
saveRDS(site_tab, file.path(rds_dir, "explore_ss_site_tab.rds"))

##==============================================================================
## 2. Likelihood-ratio decomposition [PAPER + LETTER]
##
## Four nested models, all with a free location at every gauge:
##
##   A0  common spread, common shape        k + 2 parameters
##   A1  free spread,   common shape       2k + 1
##   A2  common spread, free shape         2k + 1
##   B   free spread,   free shape         3k
##
## Fitted by profiling: the parameters shared across gauges are the only
## ones in the outer optimisation, and each gauge's own parameters are
## optimised inside. This keeps every optimisation low-dimensional
## instead of throwing 96 parameters at Nelder-Mead at once.

k <- length(series)

## A0: outer over (sigma, xi), inner over mu at each gauge
prof_a0 <- function(par) {
    sigma <- par[1]
    xi <- par[2]
    if(sigma <= 0) return(1e10)
    sum(map_dbl(series, \(x) optimize(\(mu) nll_site(c(mu, sigma, xi), x),
        interval = range(x) + c(-3, 3) * sd(x))$objective))
}

## A1: outer over xi, inner over (mu, sigma) at each gauge
prof_a1 <- function(xi) {
    sum(map_dbl(series, \(x) optim(start_par(x)[1:2],
        \(p) nll_site(c(p, xi), x),
        control = list(maxit = 5000, reltol = 1e-10))$value))
}

## A2: outer over sigma, inner over (mu, xi) at each gauge
prof_a2 <- function(sigma) {
    if(sigma <= 0) return(1e10)
    sum(map_dbl(series, \(x) optim(c(start_par(x)[1], 0.05),
        \(p) nll_site(c(p[1], sigma, p[2]), x),
        control = list(maxit = 5000, reltol = 1e-10))$value))
}

nll_b <- sum(map_dbl(fits, \(o) o$value))

o_a0 <- optim(c(median(site_tab$sigma), median(site_tab$xi)), prof_a0,
    control = list(maxit = 5000, reltol = 1e-12))
o_a1 <- optimize(prof_a1, interval = c(-0.9, 0.9), tol = 1e-8)
o_a2 <- optimize(prof_a2, interval = c(1e-3, 3 * max(site_tab$sigma)),
    tol = 1e-8)

lrt <- tibble(
    Comparison = c("Both constant vs both free",
        "Shape constant vs shape free (spread free)",
        "Spread constant vs spread free (shape free)"),
    Deviance = c(2 * (o_a0$value - nll_b),
        2 * (o_a1$objective - nll_b),
        2 * (o_a2$objective - nll_b)),
    df = c(2 * (k - 1), k - 1, k - 1)) |>
    mutate(p_value = pchisq(Deviance, df, lower.tail = FALSE),
        Deviance = round(Deviance, 1))
lrt
saveRDS(lrt, file.path(rds_dir, "explore_ss_lrt.rds"))

## Common values under A0, for the overlay and the text
common_sigma <- o_a0$par[1]
common_xi <- o_a0$par[2]
common_s_beta <- old_to_new(c(0, common_sigma, common_xi),
    alpha = alpha_loc, beta = beta_spr)$s
tibble(common_sigma, common_xi, common_s_beta)

##==============================================================================
## 3. What the fitted model says [PAPER + LETTER]
##
## The selected stationary spatio-temporal model reports one spread and
## one tail for the whole coast. Read them to overlay on the per-site
## estimates.

geom_file <- file.path(rds_dir, "pinned_geometry.rds")
mod_spread <- NA_real_
mod_tail <- NA_real_

if(file.exists(geom_file) && requireNamespace("INLA", quietly = TRUE)) {
    mesh_n <- readRDS(geom_file)$mesh$n
    f <- file.path(rds_dir,
        sprintf("explore_dc_max_bGEV_m0st_nb_mesh_%d.rds", mesh_n))
    if(file.exists(f)) {
        marg <- marginals_st_gev(readRDS(f))
        zz <- map(marg[c("spread", "tail")],
            \(m) INLA::inla.zmarginal(as.matrix(m), silent = TRUE))
        mod_spread <- zz$spread$mean
        mod_tail <- zz$tail$mean
    }
}
tibble(mod_spread, mod_tail)

if(is.na(mod_tail)) {
    message("## m0st_nb fit not available: model overlay omitted.")
}

##==============================================================================
## 4. Spread and shape against latitude [PAPER]
##
## The referee expects windward and leeward differences. Latitude is the
## cleanest single axis here, and the exposure grouping in block 5 is the
## more direct test. The solid line is the common value fitted under A0
## and the dashed line is the posterior mean of the fitted model.

plot_par <- function(dat, est, lo, hi, ylab, common, model) {
    p <- ggplot(dat, aes(x = lat, y = {{ est }})) +
        geom_hline(yintercept = common, colour = col_ink, linewidth = 0.4) +
        geom_linerange(aes(ymin = {{ lo }}, ymax = {{ hi }}),
            colour = col_bar, linewidth = 0.3) +
        geom_point(aes(shape = country), size = 2, colour = col_ink) +
        scale_shape_manual(name = NULL, values = c(GBR = 16, IRL = 17)) +
        labs(x = "Latitude (degrees north)", y = ylab)
    if(!is.na(model)) {
        p <- p + geom_hline(yintercept = model, colour = col_ink,
            linewidth = 0.4, linetype = "dashed")
    }
    p
}

p_spread <- plot_par(site_tab, s_beta, lo_s, hi_s,
    expression(paste("Spread ", s[beta], " (m)")),
    common_s_beta, mod_spread)
p_spread
save_fig(p_spread, "ss_spread_lat", w = 7, h = 4.5)

p_shape <- plot_par(site_tab, xi, lo_xi, hi_xi,
    expression(paste("Shape ", xi)), common_xi, mod_tail)
p_shape
save_fig(p_shape, "ss_shape_lat", w = 7, h = 4.5)

##==============================================================================
## 5. Exposure grouping [PAPER or APPENDIX]
##
## A coarse windward and leeward proxy. Prevailing storms arrive from the
## southwest, so the Atlantic-facing gauges are the exposed ones, the
## Irish Sea basin is sheltered by both islands, and the Bristol Channel
## is a funnel that amplifies surge. The rule is stated here rather than
## fitted, and is only meant to be indicative.

site_tab <- site_tab |>
    mutate(exposure = case_when(
        lon < -8.0 ~ "Atlantic",
        lat < 52.0 & lon > -5.5 ~ "Bristol Channel",
        TRUE ~ "Irish Sea"))

site_tab |>
    count(exposure, country) |>
    pivot_wider(names_from = country, values_from = n, values_fill = 0)

exposure_tab <- site_tab |>
    group_by(exposure) |>
    summarise(n = n(),
        spread_med = median(s_beta), spread_min = min(s_beta),
        spread_max = max(s_beta),
        xi_med = median(xi), xi_min = min(xi), xi_max = max(xi),
        .groups = "drop")
exposure_tab
saveRDS(exposure_tab, file.path(rds_dir, "explore_ss_exposure_tab.rds"))

p_exposure <- site_tab |>
    select(site_name, exposure, s_beta, xi) |>
    pivot_longer(c(s_beta, xi), names_to = "par", values_to = "value") |>
    mutate(par = factor(par, levels = c("s_beta", "xi"),
        labels = c("Spread (m)", "Shape"))) |>
    ggplot(aes(x = exposure, y = value)) +
    geom_boxplot(outlier.shape = NA, colour = col_bar, fill = NA) +
    geom_jitter(width = 0.15, height = 0, size = 1.6, colour = col_ink) +
    facet_wrap(~ par, scales = "free_y") +
    labs(x = NULL, y = NULL)
p_exposure
save_fig(p_exposure, "ss_exposure", w = 8, h = 4.5)

##==============================================================================
## 6. Where the shape sits relative to the model support [PAPER + LETTER]
##
## The fitted models use the bGEV with the tail confined to [0, 0.5] by
## its PC prior. Count how many per-site estimates fall outside that
## interval, and how many have an interval that excludes zero.

shape_pos <- site_tab |>
    summarise(
        n = n(),
        n_negative = sum(xi < 0),
        n_ci_below_zero = sum(hi_xi < 0),
        n_ci_above_zero = sum(lo_xi > 0),
        n_ci_covers_zero = sum(lo_xi <= 0 & hi_xi >= 0))
shape_pos
saveRDS(shape_pos, file.path(rds_dir, "explore_ss_shape_pos.rds"))

message("05_explore_spread_shape.R complete.")

##==============================================================================
## CONCLUSION
##
## 1. The constant assumption is rejected, and the referee is right.
##    Joint likelihood-ratio test of "both constant" against "both free":
##    deviance 305.0 on 62 degrees of freedom, p = 8e-34. This is not a
##    borderline result and it is not an artefact of small samples: the
##    test already accounts for the sampling noise that makes 32 separate
##    estimates scatter.
##
## 2. The violation is almost entirely in the SPREAD, not the shape.
##    Decomposing the joint test:
##      spread constant vs free (shape free): deviance 166.4 on 31 df,
##        p = 2e-20
##      shape constant vs free (spread free): deviance  53.3 on 31 df,
##        p = 0.008
##    So a single spread for the whole coast is badly wrong, while a
##    single shape is only mildly questionable. This matters for what to
##    do next: the actionable change is a covariate on the spread.
##
## 3. The spread varies by a factor of four, from 0.19 m at
##    CastletownberePort to 0.78 m at Avonmouth, and it tracks basin
##    geometry rather than exposure to the open Atlantic. Medians are
##    0.28 m for the Atlantic gauges, 0.36 m in the Irish Sea, and 0.45 m
##    in the Bristol Channel. The referee expected a windward and leeward
##    contrast. What the data show is funnel amplification: the most
##    sheltered basin has the LARGEST spread, because a converging
##    estuary amplifies surge. An exposure indicator built on that
##    geometry is therefore a better covariate than a windward and
##    leeward one.
##
## 4. CAVEAT on the magnitude, not on the heterogeneity. These are
##    independent per-site fits with a constant location, so all
##    within-site temporal variation is absorbed into the spread. The
##    fitted model instead attributes part of that variation to the
##    AR(1) field. That is why the common spread here (0.45 m) exceeds
##    the model posterior mean (0.27 m). The two numbers are not
##    comparable and should not be quoted side by side as a discrepancy.
##    The BETWEEN-SITE variation is unaffected by this and is the real
##    finding.
##
## 5. Separate and more serious: the shape sits where the model cannot
##    follow. Of the 32 gauges, 30 have a negative point estimate and 19
##    have a 95% interval lying entirely below zero. Not one is
##    significantly positive. The common shape is -0.149. The fitted
##    models use the bGEV, whose upper tail is Frechet by construction
##    and requires xi > 0, with the PC prior confining the tail to
##    (0, 0.5); the reported posterior mean is 0.040, i.e. pinned against
##    the lower edge of the admissible range.
##
##    So the data prefer a BOUNDED upper tail while the likelihood can
##    only represent an unbounded one. The practical consequence is that
##    long-period return levels are extrapolated with a heavier tail than
##    the data support, which is conservative and therefore safe for
##    coastal risk, but it must be stated rather than left implicit. This
##    is a limitation of the chosen likelihood family, not of the spatial
##    model, and it is not fixed by letting the tail vary in space.
##
## In one line: the constant-spread assumption should be relaxed with a
## covariate, the constant-shape assumption is defensible, and the
## positivity constraint on the tail deserves an explicit paragraph.
##==============================================================================
