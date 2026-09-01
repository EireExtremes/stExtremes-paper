safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  sqrt(mean(x^2))
}

safe_mae <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(abs(x))
}

extract_tide_features <- function(fit) {
  out <- tibble::tibble(
    form_factor = if (!is.null(fit$fval)) unname(fit$fval) else NA_real_,
    MLWS = NA_real_,
    MLWN = NA_real_,
    MSL_feature = NA_real_,
    MHWN = NA_real_,
    MHWS = NA_real_,
    MLLW = NA_real_,
    MHLW = NA_real_,
    MLHW = NA_real_,
    MHHW = NA_real_
  )
  
  if (!is.null(fit$features1)) {
    f1 <- fit$features1
    
    if ("MLWS" %in% names(f1)) out$MLWS <- unname(f1["MLWS"])
    if ("MLWN" %in% names(f1)) out$MLWN <- unname(f1["MLWN"])
    if ("MSL"  %in% names(f1)) out$MSL_feature <- unname(f1["MSL"])
    if ("MHWN" %in% names(f1)) out$MHWN <- unname(f1["MHWN"])
    if ("MHWS" %in% names(f1)) out$MHWS <- unname(f1["MHWS"])
  }
  
  if (!is.null(fit$features2)) {
    f2 <- fit$features2
    
    if ("MLLW" %in% names(f2)) out$MLLW <- unname(f2["MLLW"])
    if ("MHLW" %in% names(f2)) out$MHLW <- unname(f2["MHLW"])
    if ("MSL" %in% names(f2) && is.na(out$MSL_feature)) {
      out$MSL_feature <- unname(f2["MSL"])
    }
    if ("MLHW" %in% names(f2)) out$MLHW <- unname(f2["MLHW"])
    if ("MHHW" %in% names(f2)) out$MHHW <- unname(f2["MHHW"])
  }
  
  out
}

infer_prediction_step_seconds <- function(time_vec) {
  time_vec <- sort(unique(time_vec))
  
  if (length(time_vec) < 2) return(3600)
  
  dt <- diff(as.numeric(time_vec))
  dt <- dt[is.finite(dt) & dt > 0]
  
  if (length(dt) == 0) return(3600)
  
  max(1, round(stats::median(dt)))
}

fit_ftide_one_group <- function(df,
                                min_obs_required,
                                hcn_use,
                                use_nodal = TRUE,
                                use_smsl = TRUE) {
  
  df <- df |>
    dplyr::arrange(time) |>
    dplyr::distinct(time, .keep_all = TRUE)
  
  n_rows <- nrow(df)
  n_nonmissing <- sum(!is.na(df$sea_level))
  n_missing <- sum(is.na(df$sea_level))
  
  if (n_rows == 0) {
    return(list(ok = FALSE, reason = "No rows in group"))
  }
  
  if (n_nonmissing < min_obs_required) {
    return(list(
      ok = FALSE,
      reason = paste0(
        "Too few non-missing observations: ",
        n_nonmissing,
        " < ",
        min_obs_required
      )
    ))
  }
  
  time_start <- min(df$time, na.rm = TRUE)
  time_end   <- max(df$time, na.rm = TRUE)
  
  duration_days <- as.numeric(
    difftime(time_end, time_start, units = "days")
  )
  
  fit <- tryCatch(
    TideHarmonics::ftide(
      x     = df$sea_level,
      dto   = df$time,
      hcn   = hcn_use,
      nodal = use_nodal,
      smsl  = use_smsl
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(list(
      ok = FALSE,
      reason = paste("ftide error:", conditionMessage(fit))
    ))
  }
  
  step_seconds <- infer_prediction_step_seconds(df$time)
  step_hours   <- step_seconds / 3600
  
  pred <- tryCatch(
    predict(
      fit,
      from  = time_start,
      to    = time_end,
      by    = step_hours,
      msl   = TRUE,
      split = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(pred, "error")) {
    return(list(
      ok = FALSE,
      reason = paste("Prediction error:", conditionMessage(pred))
    ))
  }
  
  pred_time <- seq(
    from = time_start,
    to   = time_end,
    by   = step_seconds
  )
  
  n_pred <- min(length(pred), length(pred_time))
  pred_time <- pred_time[seq_len(n_pred)]
  pred_vals <- as.numeric(pred)[seq_len(n_pred)]
  
  tide_pred_at_obs <- approx(
    x = as.numeric(pred_time),
    y = pred_vals,
    xout = as.numeric(df$time),
    rule = 2
  )$y
  
  residuals_df <- df |>
    dplyr::mutate(
      tide_pred = tide_pred_at_obs,
      residual = sea_level - tide_pred,
      skew_surge_date = as.Date(time, tz = "UTC")
    ) |>
    dplyr::group_by(skew_surge_date) |>
    dplyr::mutate(
      daily_max_sea_level = if (all(is.na(sea_level))) {
        NA_real_
      } else {
        max(sea_level, na.rm = TRUE)
      },
      daily_max_tide_pred = if (all(is.na(tide_pred))) {
        NA_real_
      } else {
        max(tide_pred, na.rm = TRUE)
      },
      surge_estimate = daily_max_sea_level - daily_max_tide_pred
    ) |>
    dplyr::ungroup()
  
  summary_tbl <- tibble::tibble(
    time_start = time_start,
    time_end = time_end,
    duration_days = duration_days,
    n_rows = n_rows,
    n_nonmissing = n_nonmissing,
    n_missing = n_missing,
    prop_missing = n_missing / n_rows,
    n_predicted = sum(!is.na(residuals_df$tide_pred)),
    prediction_coverage = mean(!is.na(residuals_df$tide_pred)),
    inferred_step_seconds = step_seconds,
    rmse_residual = safe_rmse(residuals_df$residual),
    mae_residual = safe_mae(residuals_df$residual),
    mean_residual = mean(residuals_df$residual, na.rm = TRUE),
    sd_residual = sd(residuals_df$residual, na.rm = TRUE),
    max_abs_residual = if (all(is.na(residuals_df$residual))) {
      NA_real_
    } else {
      max(abs(residuals_df$residual), na.rm = TRUE)
    }
  ) |>
    dplyr::bind_cols(extract_tide_features(fit))
  
  list(
    ok = TRUE,
    fit = fit,
    residuals = residuals_df,
    summary = summary_tbl
  )
}
