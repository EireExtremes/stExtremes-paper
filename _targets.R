##==============================================================================
## targets pipeline for the sea-level extremes analysis
##
## Run it with:
##   targets::tar_make()
##
## and inspect it without running anything with:
##   targets::tar_manifest()
##   targets::tar_visnetwork()   # needs visNetwork
##
## HOW THE ANALYSIS IS WIRED INTO targets
##
## The analysis is a set of scripts that run for their side effects: each
## writes cached fits into rds/ and figures into figs/. Rather than
## rewrite them as functions, each target sources one script and returns
## the outputs that later steps actually consume. Those returns are
## format = "file" targets, so targets tracks the outputs on disk and
## reruns a step when its script, its inputs, or its outputs change.
##
## Each target declares a CONTRACT: the files that step must produce. If a
## script finishes without writing them, the target fails rather than
## letting a later step read a stale file. The contract is deliberately
## the small set of outputs the pipeline depends on, not every file a
## script writes; the rest are side effects in the same two directories.
##
## Two things about the shape of the graph:
##
## - 00_explore_data.R runs ONCE here, after 01, not twice as the old
##   shell driver ran it. Its data figures need only the CSV, but its mesh
##   figure needs rds/pinned_geometry.rds, which 01 writes. Ordering it
##   after 01 gets the whole set in one pass.
## - 04, 08 and 09 source 01 themselves, for the model settings. That is
##   cheap: 01 caches its fits, so the second pass reads them back rather
##   than refitting. They still depend on the models target, so the cache
##   is guaranteed to exist first.
##
## Not in this pipeline, and run on their own when wanted:
## 05_explore_spread_shape.R (the frequentist companion to the INLA
## version), 06_tail_integration_check.R, 07_tail_stability_m0st_nb.R and
## verify_return_levels.R, which is sourced after 02.
##
## Expect about 11 hours for a cold run on 8 cores, nearly all of it in
## 01 (~6.5 h), 09 (~3 h) and 08 (~1.2 h).
##==============================================================================

library(targets)

tar_option_set(
    packages = "here",
    ## Each target is evaluated in the same external R process, one after
    ## another. The INLA fits are already internally parallel and the
    ## scripts pin BLAS to one thread, so there is nothing to gain from
    ## running targets in parallel here, and memory would suffer.
    format = "rds")

##==============================================================================
## Helper

## Source one analysis script and return the files it was supposed to
## write. `depends` exists so a target can name its upstream targets and
## have targets draw the edge; the value itself is not used.
run_script <- function(script, outputs, depends = NULL) {
    path <- here::here(script)
    stopifnot(file.exists(path))
    env <- new.env(parent = globalenv())
    source(path, local = env, echo = FALSE)
    produced <- here::here(outputs)
    missing <- outputs[!file.exists(produced)]
    if(length(missing) > 0) {
        stop(script, " finished without writing: ",
            paste(missing, collapse = ", "), call. = FALSE)
    }
    produced
}

##==============================================================================
## Pipeline

list(
    ## The one input the whole analysis reads.
    tar_target(data_csv, here::here("data", "dc_max.csv"), format = "file"),

    ## Cached days-with-data per site-year. Tracked in the repository, so
    ## the coverage blocks of 00 run without the raw sub-hourly record.
    tar_target(coverage_cache,
        here::here("rds", "explore_data_cov_year.rds"), format = "file"),

    ##--------------------------------------------------------------------
    ## Fitting

    ## The expensive one: fits and caches every retained model, and pins
    ## the mesh so downstream steps use the same geometry.
    tar_target(models,
        run_script("01_explore_models.R",
            "rds/pinned_geometry.rds",
            depends = data_csv),
        format = "file"),

    ##--------------------------------------------------------------------
    ## Data section

    ## After the models, because the mesh figure needs the pinned geometry.
    tar_target(explore_data,
        run_script("00_explore_data.R",
            c("figs/map_study_area.png",
              "figs/data_record_length.png",
              "figs/data_glyph_map.png",
              "figs/map_mesh.png"),
            depends = c(data_csv, models, coverage_cache)),
        format = "file"),

    ##--------------------------------------------------------------------
    ## Results

    ## Barrier vs stationary: fields and return levels.
    tar_target(compare_barrier,
        run_script("02_compare_barrier.R",
            c("figs/compare_qalpha_mean.png",
              "figs/compare_qalpha_diff.png",
              "figs/compare_rlevel_diff.png",
              "figs/compare_field_mean.png",
              "figs/compare_field_sd.png"),
            depends = models),
        format = "file"),

    ## Empirical trend check. Self-contained: reads the CSV only.
    tar_target(trend_check,
        run_script("03_trend_check.R",
            "figs/trend_check_site_slopes.png",
            depends = data_csv),
        format = "file"),

    ## Posterior-predictive return levels, with full uncertainty.
    tar_target(posterior_predictive,
        run_script("04_posterior_predictive.R",
            c("figs/pp_rlevel_median_map.png",
              "figs/pp_rlevel_ciwidth_map.png",
              "figs/pp_rlevel_sites_key.png",
              "figs/pp_rlevel_sites_all.png",
              "figs/pp_rlevel_quantile_map_50yr.png"),
            depends = models),
        format = "file"),

    ## Does spread or shape need to vary by gauge? (Referee 1, point 3.)
    tar_target(spread_shape,
        run_script("05_explore_spread_shape_inla.R",
            "rds/explore_ssi_post.rds",
            depends = data_csv),
        format = "file"),

    ## Tail prior sensitivity (Referee 1, point 3).
    tar_target(tail_prior,
        run_script("08_tail_prior_sensitivity.R",
            c("rds/tailprior_shift.rds",
              "rds/tailprior_summary.rds",
              "figs/tailprior_rlevel_shift.png"),
            depends = models),
        format = "file"),

    ## Mesh convergence (Referee 2, Sec. 4.2).
    tar_target(mesh_convergence,
        run_script("09_mesh_convergence.R",
            c("rds/mesh_conv_hyper.rds",
              "rds/mesh_conv_rlevel.rds",
              "figs/mesh_conv_rlevel.png"),
            depends = models),
        format = "file"),

    ##--------------------------------------------------------------------
    ## Figures

    ## Trims the baked-in white margins. It rewrites the figures in place,
    ## so it returns a manifest of what it touched rather than the paths
    ## the steps above already own.
    tar_target(figures_cropped, {
        deps <- c(explore_data, compare_barrier, trend_check,
            posterior_predictive, spread_shape, tail_prior,
            mesh_convergence)
        run_script("crop_figs.R", character(0))
        pngs <- list.files(here::here("figs"), pattern = "\\.png$",
            full.names = TRUE)
        data.frame(
            file = basename(pngs),
            bytes = file.size(pngs),
            row.names = NULL)
    })
)
