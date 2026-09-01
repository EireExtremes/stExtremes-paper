##==============================================================================
## Empirical trend check for the annual-max tidal residuals.
##
## Ground-truths the fitted bGEV temporal trend (trendc) against the raw data,
## independent of INLA / the bGEV machinery, to see whether the ~10-17 mm/yr
## fitted trend is present in the observations or is inflated by the AR(1)
## field and the changing station network. Self-contained: reads dc_max
## directly, so it can be sourced after 01_ or 02_ but does not depend on them.
##
## Three views, all on the annual maxima (NA years dropped):
##   1. per-site Theil-Sen slope (robust, a trend-in-the-median estimate, the
##      analog of the bGEV location trend) and per-site OLS slope;
##   2. pooled within-site OLS trend (max ~ year + site_name) on all sites;
##   3. the same pooled trend restricted to long-record sites, to expose any
##      network / NA artifact.
## All slopes are in m/yr and directly comparable to trendc (per year_idx =
## per year). Reference values: the fitted bGEV trends and the regional MSL
## band are printed and drawn for eyeballing.
##==============================================================================

here::i_am("03_trend_check.R")

suppressPackageStartupMessages(library(tidyverse))

theme_set(theme_bw(base_size = 14))

figs_dir <- here::here("figs")
if(!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

## Tunables and reference values.
min_years <- 50                                  # long-record threshold
bgev_trend <- c(barrier = 0.0096, stationary = 0.0172)   # fitted trendc means
msl_band <- c(0.002, 0.004)                      # regional MSL rise (m/yr)

dat <- read_csv(here::here("data", "dc_max.csv"),
    show_col_types = FALSE) |>
    filter(!is.na(max)) |>
    select(country, site_name, year, year_idx, max)

##==============================================================================
## Per-site slopes
##==============================================================================
ols_slope <- function(y, x) {
    if(sum(is.finite(x) & is.finite(y)) < 3) return(NA_real_)
    unname(coef(lm(y ~ x))[2])
}

theil_sen <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    if(length(x) < 3) return(NA_real_)
    cb <- combn(length(x), 2)
    dx <- x[cb[2, ]] - x[cb[1, ]]
    dy <- y[cb[2, ]] - y[cb[1, ]]
    median(dy[dx != 0] / dx[dx != 0])
}

site_slopes <- dat |>
    summarise(
        n_obs = n(),
        slope_ols = ols_slope(max, year),
        slope_sen = theil_sen(year, max),
        .by = c(country, site_name))

site_slope_summary <- site_slopes |>
    summarise(
        n_sites = n(),
        med_ols = median(slope_ols, na.rm = TRUE),
        med_sen = median(slope_sen, na.rm = TRUE),
        q25_sen = quantile(slope_sen, 0.25, na.rm = TRUE),
        q75_sen = quantile(slope_sen, 0.75, na.rm = TRUE),
        frac_pos = mean(slope_sen > 0, na.rm = TRUE))
print(site_slope_summary)

##==============================================================================
## Pooled within-site trend (common year slope, site fixed effects)
##==============================================================================
pooled_ols <- function(d) {
    m <- lm(max ~ year + site_name, data = d)
    ci <- confint(m, "year")
    tibble(slope = unname(coef(m)[["year"]]),
        lwr = ci[1], upr = ci[2],
        n_site = n_distinct(d$site_name), n_obs = nrow(d))
}

long_sites <- site_slopes |>
    filter(n_obs >= min_years) |>
    pull(site_name)

pooled_tab <- bind_rows(
    mutate(pooled_ols(dat), sample = "all sites", .before = 1),
    mutate(pooled_ols(filter(dat, site_name %in% long_sites)),
        sample = sprintf(">=%d yrs", min_years), .before = 1))
print(pooled_tab)

reference <- tibble(
    source = c("bGEV m1st (barrier)", "bGEV m1st_nb (stationary)",
        "MSL rise low", "MSL rise high"),
    slope = c(bgev_trend[["barrier"]], bgev_trend[["stationary"]],
        msl_band[1], msl_band[2]))
print(reference)

##==============================================================================
## Figure: per-site Theil-Sen slopes vs the fitted trends and the MSL band
##==============================================================================
ref_df <- tibble(model = names(bgev_trend), trend = unname(bgev_trend))

p_slopes <- ggplot(site_slopes, aes(slope_sen)) +
    annotate("rect", xmin = msl_band[1], xmax = msl_band[2],
        ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.6) +
    geom_histogram(bins = 20, fill = "grey40", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_vline(data = ref_df, aes(xintercept = trend, colour = model),
        linewidth = 0.9) +
    scale_colour_manual(name = "Fitted bGEV trend",
        values = c(barrier = "firebrick", stationary = "dodgerblue3")) +
    labs(
        ## Label kept non-technical on purpose: the estimator is
        ## Theil-Sen (see above), but the paper and the response letter
        ## describe it simply as a trend fitted at each gauge, since the
        ## choice of robust estimator is not the point being made.
        x = "Trend in the annual maxima at each gauge (m/yr)",
        y = "Number of sites",
        title = "Empirical trend vs fitted bGEV trend (grey band = MSL rise)")
ggsave(file.path(figs_dir, "trend_check_site_slopes.png"), p_slopes,
    width = 9, height = 6, dpi = 300, device = png)
ggsave(file.path(figs_dir, "trend_check_site_slopes.pdf"), p_slopes,
    width = 9, height = 6, dpi = 300, device = cairo_pdf)

message("## Done. See site_slope_summary, pooled_tab, reference, figure.")
