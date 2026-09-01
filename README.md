# stExtremes-paper

Data and analysis code for:

> **A Bayesian Hierarchical Spatiotemporal Model With Physical Barriers
> for Extreme Sea-Level Prediction in Ireland**
>
> Fernando Mayer, Niamh Mimnagh and Niamh Cahill.
> Department of Mathematics and Statistics, Maynooth University.
>
> *Environmetrics*, 2026.
> <https://onlinelibrary.wiley.com/doi/10.1002/env.70135>

The paper estimates extreme sea levels at gauged and ungauged coastal
locations around Ireland and the west coast of Great Britain, modelling
annual maxima of tidal residuals with the blended generalised extreme
value (bGEV) distribution, whose location quantile varies in space and
time through a latent Gaussian field, with coastlines entering as
physical barriers on the spatial correlation. Inference is by INLA with
an SPDE representation.

This repository is the target of the paper's Data Availability Statement.
It holds the modelled dataset and the code that produces every result in
the paper. The manuscript itself is not here.

## What is here

| Path | What it is |
|---|---|
| `data/dc_max.csv` | The modelled dataset: annual maximum tidal residuals at 44 gauges, 1968-2025 |
| `data/dc_max.R` | The filter that produces it from the all-time residual maxima |
| `data/dc_max_all_time_residual.rda` | Annual maxima before the site and year filters |
| `data/dc_max_files/` | The tidal-harmonic scripts upstream of that, and their two annual-maxima inputs |
| `00_`-`09_*.R` | The analysis |
| `verify_return_levels.R` | Checks the bGEV return-level parametrisation. Source it after `02_` |
| `crop_figs.R` | Trims the baked-in white margins from the saved figures |
| `_targets.R` | The pipeline: what runs, in what order, and what each step must produce |
| `renv.lock` | The pinned computational environment |

## The data

`data/dc_max.csv` is the single input every analysis script reads: 2552
site-years over 44 gauges and 1968-2025, of which 1066 have an observed
maximum. Years a gauge did not record are present as `NA` rather than
absent, so the panel is rectangular. `max` is the annual maximum tidal
residual in metres, that is, observed sea level minus the predicted tide
from a harmonic fit, so it is the non-tidal (surge) component.

The raw sub-hourly records the residuals come from are not redistributed
here. Only the within-year coverage blocks of `00_explore_data.R` need
them, and those blocks skip themselves with a message when the file is
absent; their result is cached in `rds/explore_data_cov_year.rds`.

`data/dc_max_files/03b_fit_ftide_models.R` carries three commented-out
`source()` calls for helpers that are not in this repository. That is why
the two annual-maxima tables under `data/dc_max_files/rds/` are tracked:
they are the inputs `05_dc_max.R` needs and they cannot be regenerated
here.

## Setting up

The environment is pinned with [renv](https://rstudio.github.io/renv/).
From a clone:

```r
renv::restore()
```

That installs everything at the recorded versions, into a project-local
library that leaves the rest of your R installation alone. It includes
`INLA` 26.06.08 from <https://inla.r-inla-download.org/R/testing> (not on
CRAN) and the `stExtremes` package from GitHub at the exact commit the
analysis was run against. Opening R in this directory activates renv
automatically through `.Rprofile`; until `renv::restore()` has been run,
the project library is empty and nothing will load.

## Running it

```r
targets::tar_make()
```

`_targets.R` defines the pipeline and its dependencies, so a re-run only
repeats what is out of date. `targets::tar_visnetwork()` draws the graph;
`targets::tar_manifest()` lists the steps without running them.

Each step declares the files it must produce, and fails if it finishes
without writing them, so a later step cannot quietly read a stale result.

Expect about 11 hours for a cold run on 8 cores, nearly all of it in
`01_explore_models.R` (~6.5 h), `09_mesh_convergence.R` (~3 h) and
`08_tail_prior_sensitivity.R` (~1.2 h). Everything else takes minutes.
Fitted models are also cached under `rds/`, so even a forced re-run skips
refitting unless that cache is deleted.

Note that `00_explore_data.R` runs after `01_explore_models.R`, not
before: its data figures need only the CSV, but its mesh figure needs the
pinned geometry that `01` writes.

Four scripts sit outside the pipeline and are run on their own when
wanted: `05_explore_spread_shape.R`, the frequentist companion to the
INLA version; `06_tail_integration_check.R` and
`07_tail_stability_m0st_nb.R`, which diagnose the stability of the tail
posterior; and `verify_return_levels.R`, which is sourced after `02`.

## Reproducibility

Beyond the pinned library, the scripts control what makes INLA runs
comparable across machines: the mesh is built once and reused from
`rds/pinned_geometry.rds`, BLAS and OpenMP are held to one thread, INLA's
integration strategy is fixed, and each machine's R, BLAS/LAPACK,
GEOS/PROJ/GDAL and package versions are recorded in
`rds/env_info_<hostname>.rds` on first run.

`renv.lock` pins **INLA 26.06.08**, the version that produced the results
in the paper, rather than whatever is current. The INLA testing
repository keeps its historical builds, so `renv::restore()` retrieves
that exact version; this has been checked by restoring it. Note that a
newer INLA is released roughly monthly and will not be used unless the
lockfile is deliberately updated.

One caveat remains: results move a little between machines even at
identical package versions, mostly through the geometry libraries
(GEOS/PROJ/GDAL), which renv does not control. That is why the
environment is recorded per run rather than assumed.

## Outputs

Generated files are not tracked: `rds/` (the cached fits, several GB),
`figs/`, and `_targets/`. All are rebuilt by `tar_make()`.

## Licence

GPL-3. See `LICENSE`.
