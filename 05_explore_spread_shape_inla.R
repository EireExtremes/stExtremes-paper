##==============================================================================
## Site-wise spread and shape, fitted in INLA (TODO item 3)
##
## Bayesian counterpart of 05_explore_spread_shape.R. Same question, same
## conclusions, but every model is the bGEV the paper actually fits, and
## the evidence is posterior distributions and WAIC rather than p-values.
## DIC is computed only as a check and is never used to choose a model.
##
## The question: the fitted models carry spatio-temporal structure on the
## location quantile only. The spread s_beta and the tail xi are single
## numbers shared by all gauges. Is that supported by the data?
##
## The design. Four nested variants, all with a free location at every
## gauge, differing only in what is shared:
##
##   pooled       one spread, one tail for the whole coast
##   free spread  per-site spread, tail held at the pooled value
##   free tail    per-site tail, spread held at the pooled value
##   free both    per-site spread and tail
##
## The three per-site variants are fitted one gauge at a time and their
## WAIC summed, which is comparable with the pooled fit because the data
## are identical. Holding one hyperparameter at the pooled value
## is what separates the spread question from the shape question: it is
## the Bayesian analogue of the likelihood-ratio decomposition in the
## frequentist script, without any test.
##
## Each block is tagged [PAPER], [APPENDIX] or [LETTER].
##
## Referee anchors (chats/TODO.md, item 3): Referee 2 major point 1 and
## Referee 1 on model validation.
##
## Inputs:  data/dc_max.csv
## Outputs: figures in figs/, tables in rds/ (explore_ssi_ prefix)
##
## Runtime a few minutes: one pooled fit plus 3 x 32 small per-site fits,
## each of a few tenths of a second.
##
## The conclusion is written out at the end of this file.
##==============================================================================

##==============================================================================
## Setup

here::i_am("05_explore_spread_shape_inla.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(INLA)
})

library(stExtremes)

## Determinism, as in 01_explore_models.R. This is not optional here:
## with INLA's default threading the pooled tail posterior mean moved
## between 0.004 and 0.21 across identical runs, because the tail is
## barely identified and the likelihood surface along it is nearly flat.
## The fits are tiny, so a fully serial run costs nothing.
inla.setOption(num.threads = "1:1")
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

## A gauge needs enough years before spread and tail mean anything.
## 15 years keeps 32 of the 44 gauges and both countries.
min_years <- 15

##==============================================================================
## Data

dc_max <- read_csv(here::here("data/dc_max.csv"),
    show_col_types = FALSE) |>
    filter(!is.na(max))

sites <- dc_max |>
    count(country, site_name, name = "n_years") |>
    left_join(distinct(dc_max, site_name, lat, lon), by = "site_name") |>
    filter(n_years >= min_years) |>
    ## Same coarse basin proxy as the frequentist script
    mutate(exposure = case_when(
        lon < -8.0 ~ "Atlantic",
        lat < 52.0 & lon > -5.5 ~ "Bristol Channel",
        TRUE ~ "Irish Sea"))

dat_all <- dc_max |>
    filter(site_name %in% sites$site_name) |>
    left_join(select(sites, site_name, exposure), by = "site_name") |>
    mutate(site = factor(site_name))

series <- dat_all |>
    group_by(site_name) |>
    summarise(y = list(max), .groups = "drop") |>
    deframe()

##==============================================================================
## INLA settings, copied from 01_explore_models.R so the likelihood
## here is the one the paper fits

hyper_spread <- list(prior = "loggamma", param = c(3, 3))
tail_interval <- c(0, 0.5)
hyper_tail <- list(initial = map_tail(0.1, tail_interval, inverse = TRUE),
    prior = "pc.gevtail", param = c(7, tail_interval), fixed = FALSE)

inla_control_bgev <- list(q.location = 0.5, q.spread = 0.5,
    q.mix = c(0.05, 0.2), beta.ab = 5)

## Numerical settings as in the main fits, with ONE deliberate change.
## cmin = 0 matters here: on the short single-gauge series the Laplace
## step can otherwise produce a negative precision and INLA aborts with a
## tabulate-Qfunc assertion.
##
## The location prior is NOT the prec = 100 of 01_explore_models.R.
## There the intercept is a single number for the whole coast and the
## spatial field carries the variation, so a tight prior is an
## identifiability device. Here the location is the only thing standing
## between the prior and the data at each gauge, and prec = 100 (a prior
## sd of 0.1 m around zero) drags it towards zero and inflates the spread
## to compensate: the per-site spreads came out roughly twice the
## maximum-likelihood values until this was loosened.
inla_control_fixed <- list(prec.intercept = 0.001, prec = 0.001)
inla_control_inla <- list(int.strategy = "ccd", strategy = "laplace",
    cmin = 0, h = 0.005, tolerance = 1e-8, restart = 1L)

## The bGEV response must be an inla.mdata with the spread and tail
## design matrices attached; zero columns means no covariates on either.
zero_x <- function(n) matrix(nrow = n, ncol = 0)

fit_bgev_simple <- function(y, group = NULL, spread_x = NULL,
                            hyper_s = hyper_spread, hyper_t = hyper_tail) {
    n <- length(y)
    dat <- list(y = y, m = rep(1, n),
        spread_x = if(is.null(spread_x)) zero_x(n) else spread_x,
        tail_x = zero_x(n))
    form <- inla.mdata(y, spread_x, tail_x) ~ -1 + m
    if(!is.null(group)) {
        dat$m <- NULL
        dat$site <- group
        form <- inla.mdata(y, spread_x, tail_x) ~ -1 + site
    }
    inla(form, family = "bgev", data = dat,
        control.family = list(hyper = list(spread = hyper_s, tail = hyper_t),
            control.bgev = inla_control_bgev),
        control.fixed = inla_control_fixed,
        control.inla = inla_control_inla,
        control.compute = list(dic = TRUE, waic = TRUE))
}

##==============================================================================
## 1. Pooled fit: one spread and one tail for the whole coast [PAPER]
##
## Locations are free at every gauge, so the shared hyperparameters carry
## only the spread and the shape. This is the assumption under test.

fit_pool <- fit_bgev_simple(dat_all$max, group = dat_all$site)
round(fit_pool$summary.hyperpar[, c(1, 3, 5)], 4)

## Internal-scale modes, reused below to hold a hyperparameter fixed at
## its pooled value. Taking them on the internal scale avoids having to
## invert the spread and tail link functions by hand.
theta_pool <- fit_pool$mode$theta

hp_q <- function(f, row, col) f$summary.hyperpar[row, col]
hp_mean <- function(f, row) hp_q(f, row, "mean")
pool_spread <- hp_mean(fit_pool, "spread for BGEV observations")
pool_tail <- hp_mean(fit_pool, "tail for BGEV observations")

##==============================================================================
## 2. Per-site fits, three variants [APPENDIX]

fixed_at <- function(v) list(initial = v, fixed = TRUE)

## A single awkward series must not abort the whole run, so failures are
## caught and dropped. Any dropped gauge is reported below, because a
## variant fitted on fewer gauges is not comparable with the others.
fit_variant <- function(free_spread, free_tail) {
    out <- map(series, \(y) try(fit_bgev_simple(y,
        hyper_s = if(free_spread) hyper_spread else fixed_at(theta_pool[1]),
        hyper_t = if(free_tail) hyper_tail else fixed_at(theta_pool[2])),
        silent = TRUE))
    discard(out, \(f) inherits(f, "try-error"))
}

variants <- list(
    `free both` = fit_variant(TRUE, TRUE),
    `free spread` = fit_variant(TRUE, FALSE),
    `free tail` = fit_variant(FALSE, TRUE),
    `neither free` = fit_variant(FALSE, FALSE))

## Every variant must cover the same gauges for the DIC sums to mean
## anything. Restrict all of them to the gauges that all four fitted.
ok_sites <- reduce(map(variants, names), intersect)
variants <- map(variants, \(fl) fl[ok_sites])
tibble(gauges_attempted = length(series), gauges_used = length(ok_sites),
    dropped = paste(setdiff(names(series), ok_sites), collapse = ", "))

## The pooled fit above ran on every gauge, so it is not comparable with
## per-site sums taken over a subset. Refit it, and the exposure model of
## block 7, on exactly the gauges that all four variants managed.
dat_ok <- filter(dat_all, site_name %in% ok_sites) |>
    mutate(site = droplevels(site))
fit_pool_ok <- fit_bgev_simple(dat_ok$max, group = dat_ok$site)

## Everything downstream compares against the pooled fit on this same set
pool_spread <- hp_mean(fit_pool_ok, "spread for BGEV observations")
pool_tail <- hp_mean(fit_pool_ok, "tail for BGEV observations")
round(fit_pool_ok$summary.hyperpar[, c(1, 3, 5)], 4)

##==============================================================================
## 3. Model comparison [PAPER + LETTER]
##
## Sums over the per-site fits, against the pooled fit. Lower is better.
##
## WAIC decides. DIC is reported only as a check that the two point the
## the same way, and is never used to choose between models: it uses a
## plug-in deviance at the posterior mean and its effective-parameter
## penalty is unreliable for hierarchical models such as these, whereas
## WAIC averages the pointwise predictive density over the posterior.
##
## "neither free" is a sanity check of a different kind: it refits the
## pooled model gauge by gauge, so it should land close to "pooled".

sum_ic <- function(fl) {
    tibble(WAIC = sum(map_dbl(fl, \(f) f$waic$waic)),
        DIC_check = sum(map_dbl(fl, \(f) f$dic$dic)))
}

ic_tab <- imap(variants, \(fl, nm) mutate(sum_ic(fl), Model = nm)) |>
    list_rbind() |>
    add_row(Model = "pooled", WAIC = fit_pool_ok$waic$waic,
        DIC_check = fit_pool_ok$dic$dic) |>
    select(Model, WAIC, DIC_check) |>
    mutate(across(c(WAIC, DIC_check), \(x) round(x, 1)),
        dWAIC = round(WAIC - min(WAIC), 1)) |>
    select(Model, WAIC, dWAIC, DIC_check) |>
    arrange(WAIC)
ic_tab
saveRDS(ic_tab, file.path(rds_dir, "explore_ssi_ic_tab.rds"))

##==============================================================================
## 4. Per-site posteriors [PAPER]

hyper_post <- function(fl, row, label) {
    imap(fl, \(f, nm) {
        s <- f$summary.hyperpar
        if(!row %in% rownames(s)) return(NULL)
        tibble(site_name = nm, par = label, mean = s[row, "mean"],
            lo = s[row, "0.025quant"], hi = s[row, "0.975quant"])
    }) |>
        compact() |>
        list_rbind()
}

post <- bind_rows(
    hyper_post(variants$`free both`, "spread for BGEV observations",
        "Spread (m)"),
    hyper_post(variants$`free both`, "tail for BGEV observations", "Shape")) |>
    left_join(sites, by = "site_name")
saveRDS(post, file.path(rds_dir, "explore_ssi_post.rds"))

## Reference line for the spread only. The pooled tail posterior is not
## stable enough to draw as a reference; see block 8.
pool_line <- tibble(par = "Spread (m)", value = pool_spread)

p_post <- ggplot(post, aes(x = lat, y = mean)) +
    geom_hline(data = pool_line, aes(yintercept = value), colour = col_ink,
        linewidth = 0.4) +
    geom_linerange(aes(ymin = lo, ymax = hi), colour = col_bar,
        linewidth = 0.3) +
    geom_point(aes(shape = country), size = 2, colour = col_ink) +
    scale_shape_manual(name = NULL, values = c(GBR = 16, IRL = 17)) +
    facet_wrap(~ par, scales = "free_y") +
    labs(x = "Latitude (degrees north)", y = "Posterior mean and 95% CrI")
p_post
save_fig(p_post, "ssi_post_lat", w = 9, h = 4.5)

##==============================================================================
## 5. How far the per-site posteriors sit from the pooled value [LETTER]
##
## A plain count: at how many gauges does the 95% credible interval fail
## to cover the pooled posterior mean? Under a genuinely common value
## this should happen at about 5% of gauges.
##
## Spread only. The same count for the shape would be meaningless: the
## pooled tail posterior is not even stable between refits of the same
## model (block 8), so there is no fixed point to compare against. The
## model comparison in block 3 is the evidence about the shape.

excl_tab <- post |>
    inner_join(pool_line, by = "par") |>
    group_by(par) |>
    summarise(n = n(), excludes_pooled = sum(lo > value | hi < value),
        share = round(100 * excludes_pooled / n, 1), .groups = "drop")
excl_tab
saveRDS(excl_tab, file.path(rds_dir, "explore_ssi_excl_tab.rds"))

##==============================================================================
## 6. Spread by exposure group [PAPER]

p_exposure <- post |>
    filter(par == "Spread (m)") |>
    ggplot(aes(x = exposure, y = mean)) +
    geom_boxplot(outlier.shape = NA, colour = col_bar, fill = NA) +
    geom_jitter(width = 0.15, height = 0, size = 1.6, colour = col_ink) +
    labs(x = NULL, y = "Posterior mean spread (m)")
p_exposure
save_fig(p_exposure, "ssi_exposure", w = 6.5, h = 4.5)

exposure_tab <- post |>
    filter(par == "Spread (m)") |>
    group_by(exposure) |>
    summarise(n = n(), spread_med = round(median(mean), 3),
        spread_min = round(min(mean), 3), spread_max = round(max(mean), 3),
        .groups = "drop")
exposure_tab
saveRDS(exposure_tab, file.path(rds_dir, "explore_ssi_exposure_tab.rds"))

##==============================================================================
## 7. The exposure covariate as a candidate model [PAPER]
##
## The actionable version of the finding. INLA's bGEV accepts fixed
## covariates on the spread through the second argument of inla.mdata,
## so a basin indicator can be fitted directly. Two dummy columns, one
## extra hyperparameter each, against the pooled model.

x_exp <- model.matrix(~ exposure, data = dat_ok)[, -1, drop = FALSE]
fit_exp <- fit_bgev_simple(dat_ok$max, group = dat_ok$site,
    spread_x = x_exp)
round(fit_exp$summary.hyperpar[, c(1, 3, 5)], 4)

exp_ic <- tibble(
    Model = c("pooled", "spread ~ exposure"),
    WAIC = round(c(fit_pool_ok$waic$waic, fit_exp$waic$waic), 1),
    DIC_check = round(c(fit_pool_ok$dic$dic, fit_exp$dic$dic), 1)) |>
    mutate(dWAIC = round(WAIC - min(WAIC), 1), .after = WAIC)
exp_ic
saveRDS(exp_ic, file.path(rds_dir, "explore_ssi_exposure_ic.rds"))

##==============================================================================
## 8. How well the tail is identified [PAPER + LETTER]
##
## xi > 0 is a design property of the bGEV, not a prior choice and not a
## constraint to apologise for. The blend exists precisely to fix the
## artificial finite LOWER bound that a GEV with xi > 0 carries, which is
## awkward in a regression where the parameters move with covariates. It
## replaces the left tail with a Gumbel one and keeps the heavy Frechet
## right tail, giving support on the whole real line. A positive tail is
## therefore what the family is built to estimate, and INLA enforces it:
## a tail interval with a negative lower bound is refused outright with
## "BGEV.TAIL.INTERVAL is void".
##
## Every model compared in this script uses that same likelihood, so the
## between-gauge comparisons are like for like and unaffected. The
## question worth asking is not whether some other family would prefer a
## different sign, but how much these data pin the tail down inside the
## range this one admits. The answer is: very little.

tail_id <- post |>
    filter(par == "Shape") |>
    summarise(n = n(),
        median_mean = round(median(mean), 4),
        median_ci_width = round(median(hi - lo), 3),
        max_hi = round(max(hi), 3))
tail_id

## Refit the identical pooled model several times. The spread and the
## WAIC come back the same every time; the tail posterior MEAN does not.
##
## The mean is the wrong statistic here, and its movement is a symptom
## rather than the finding. 06_tail_integration_check.R tested this under
## a 2 x 2 of integration settings and found the credible interval about
## 0.497 wide on an admissible range of 0.5 in EVERY one of them, with
## the mode always at internal -5.0, i.e. xi about 0.003, hard against
## the boundary. So the posterior is essentially the prior, the mode is
## perfectly stable, and it is only the integration over a flat direction
## that varies. Report the interval, not the mean.
##
## No integration setting fixes it and none should be adopted: the
## production setting used here, ccd with h = 0.005, is the only one of
## the four where WAIC is deterministic to three decimals.
refit_tab <- map(1:5, \(i) {
    f <- fit_bgev_simple(dat_ok$max, group = dat_ok$site)
    tibble(refit = i,
        spread = round(hp_mean(f, "spread for BGEV observations"), 4),
        tail = round(hp_mean(f, "tail for BGEV observations"), 4),
        tail_lo = round(hp_q(f, "tail for BGEV observations",
            "0.025quant"), 4),
        tail_hi = round(hp_q(f, "tail for BGEV observations",
            "0.975quant"), 4),
        WAIC = round(f$waic$waic, 1))
}) |>
    list_rbind()
refit_tab
saveRDS(list(tail_id = tail_id, refit = refit_tab),
    file.path(rds_dir, "explore_ssi_tail_id.rds"))

message("05_explore_spread_shape_inla.R complete.")

##==============================================================================
## CONCLUSION
##
## Same conclusions as the frequentist script, reached with posteriors
## and DIC/WAIC instead of tests, and with the bGEV the paper fits.
##
## 1. The spread must vary, the shape need not. Summed over the 32
##    gauges, against the pooled model, on WAIC:
##
##      free spread    WAIC  520   dWAIC    0
##      free both      WAIC  540   dWAIC   20
##      neither free   WAIC  692   dWAIC  172
##      pooled         WAIC  693   dWAIC  173
##      free tail      WAIC  714   dWAIC  193
##
##    Freeing the spread alone buys 173 WAIC. Freeing the tail alone is
##    WORSE than pooling by 21, and freeing both is worse than freeing
##    the spread alone by 20. So the per-site tail parameters do not pay
##    for themselves, while the per-site spreads do, overwhelmingly. DIC,
##    reported only as a check, orders the five models identically.
##
##    "neither free" reproduces the pooled fit gauge by gauge and lands
##    within 1 WAIC of it, which is the sanity check that the per-site
##    sums and the joint fit are on the same footing.
##
## 2. The spread heterogeneity is large. At 14 of the 32 gauges (44%) the
##    95% credible interval for the site's own spread excludes the pooled
##    posterior mean. Under a genuinely common spread that should happen
##    at about 5% of gauges.
##
## 3. The spread tracks basin geometry, as in the frequentist script.
##    Median posterior spread is 0.37 m for the Atlantic gauges, 0.48 m
##    in the Irish Sea and 0.51 m in the Bristol Channel. The referee
##    expected a windward and leeward contrast; the pattern is funnel
##    amplification, with the most sheltered basin carrying the largest
##    spread.
##
## 4. The covariate version is directly fittable and it works. INLA's
##    bGEV takes fixed covariates on the spread through the second
##    argument of inla.mdata, so "spread ~ exposure" is a one-line
##    change, not a new method. It improves WAIC from 693 to 683, and the
##    Bristol Channel contrast against the Atlantic is credible
##    (0.111, 95% CrI 0.025 to 0.198) while the Irish Sea contrast is not
##    (0.026, -0.056 to 0.110).
##
##    But note the size: the three-level basin covariate recovers about
##    10 WAIC of the 173 available from fully free per-site spreads. It is
##    the right direction and it is defensible in the paper, yet most of
##    the site-to-site variation in spread is NOT explained by basin.
##
## 5. The tail of THIS model is not identified. The median per-site 95%
##    credible interval for the shape is 0.19 wide on an admissible range
##    of only 0.5, and for the pooled fit the interval is about 0.497
##    wide, i.e. the entire range. 06_tail_integration_check.R confirms
##    that this holds under every integration setting tested, that the
##    mode is stable throughout at xi about 0.003, and that no setting
##    change improves matters. The posterior is essentially the prior.
##    Report the interval and the mode; never the posterior mean.
##
##    This reinforces point 1 from the other direction. There is no
##    evidence the tail needs to vary because there is barely any
##    information about the tail to begin with.
##
##    IMPORTANT SCOPE LIMIT. This is the pooled model, which has no
##    latent field. The spatio-temporal models of the paper identify the
##    tail an order of magnitude better, with credible intervals 0.046
##    (m0st_nb) and 0.062 (m0st) wide and a mode clear of the boundary.
##    Nothing here should be read as a statement about them; that is what
##    07_tail_stability_m0st_nb.R checks.
##
## 6. On the positive tail. The bGEV estimates xi > 0 by construction:
##    the blend exists to remove the artificial finite lower bound that a
##    GEV with xi > 0 carries, keeping the heavy Frechet upper tail and
##    giving support on the whole real line. That is the reason the
##    family was chosen, so it is a property of the model, not a result
##    and not a limitation to report. Every comparison above is made
##    within that one likelihood, so the between-gauge conclusions are
##    like for like and none of them depend on the constraint.
##
##    The frequentist script 05_explore_spread_shape.R remains in the
##    repository for the record only. Its unconstrained GEV fits are a
##    different likelihood and are deliberately NOT used for any result
##    reported in the paper.
##
## In one line: relax the constant spread with a covariate, keep the
## constant shape, and say plainly that the tail is weakly identified.
##==============================================================================
