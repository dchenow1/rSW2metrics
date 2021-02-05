
# daily time step series:
#   - precipitation,
#   - mean air temperature;
#   - SWA at -3.9 MPa at shallow 0-20 and 0-100 cm;
#   - drainage
# watershed maps of average conditions in June-July-August: same variables

metric_PPT_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  get_variable_in_months(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    sw2_out = "PRECIP",
    sw2_var = "ppt",
    var_label = "PPT_sum_JJA_mm",
    months = 6:8,
    fun_time = sum,
    var_scaler = 10,
    include_year = include_year
  )
}


metric_Tmean_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  get_variable_in_months(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    sw2_out = "TEMP",
    sw2_var = "avg_C",
    var_label = "T_mean_JJA_C",
    months = 6:8,
    fun_time = mean,
    var_scaler = 1,
    include_year = include_year
  )
}


get_SWA_JJA <- function(
  path, name_sw2_run, id_scen_used,
  out_label,
  list_years_scen_used,
  include_year = FALSE,
  soils,
  used_depth_range_cm = NULL,
  SWP_limit_MPa = -Inf,
  ...
) {
  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    swa_daily <- calc_SWA(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      soils = soils,
      SWP_limit_MPa = SWP_limit_MPa,
      used_depth_range_cm = used_depth_range_cm
    )

    # Helper variables
    ts_years <- unique(swa_daily[["time"]][, "Year"])

    # Calculate and format
    x_monthly <- aggregate(
      x = swa_daily[["values"]][[1]],
      by = list(
        Month = swa_daily[["time"]][, "Month"],
        Year = swa_daily[["time"]][, "Year"]
      ),
      FUN = mean
    )

    res[[k1]] <- format_yearly_to_matrix(
      x = list(tapply(
        X = x_monthly[["x"]],
        INDEX = x_monthly["Year"],
        FUN = function(x) mean(x[6:8])
      )),
      years = ts_years,
      out_labels = out_label
    )

    if (include_year) {
      res[[k1]] <- rbind(
        Year = ts_years,
        res[[k1]]
      )
    }
  }

  res
}



metric_SWAat0to020cm39bar_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_JJA(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat0to020cm39bar_mean_JJA_mm",
    soils = soils,
    used_depth_range_cm = c(0, 20),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}

metric_SWAat0to100cm39bar_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_JJA(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat0to100cm39bar_mean_JJA_mm",
    soils = soils,
    used_depth_range_cm = c(0, 100),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}

metric_SWAat20to100cm39bar_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_JJA(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat20to100cm39bar_mean_JJA_mm",
    soils = soils,
    used_depth_range_cm = c(20, 100),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}

metric_SWAat20to040cm39bar_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_JJA(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat20to040cm39bar_mean_JJA_mm",
    soils = soils,
    used_depth_range_cm = c(20, 40),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}


metric_SWAat40to060cm39bar_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_JJA(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat40to060cm39bar_mean_JJA_mm",
    soils = soils,
    used_depth_range_cm = c(40, 60),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}


metric_SWAat60to080cm39bar_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_JJA(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat60to080cm39bar_mean_JJA_mm",
    soils = soils,
    used_depth_range_cm = c(60, 80),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}



metric_SWAat80to100cm39bar_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_JJA(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat80to100cm39bar_mean_JJA_mm",
    soils = soils,
    used_depth_range_cm = c(80, 100),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}

metric_DR_JJA <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  get_variable_in_months(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    sw2_out = "DEEPSWC",
    sw2_var = "lowLayerDrain_cm",
    var_label = "DR_sum_JJA_mm",
    months = 6:8,
    fun_time = sum,
    var_scaler = 10,
    include_year = include_year
  )
}


metric_PPT_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    tmp_daily <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Day",
      sw2_outs = "PRECIP",
      sw2_vars = "ppt",
      varnames_are_fixed = TRUE
    )

    res[[k1]] <- format_daily_to_matrix(
      x = 10 * tmp_daily[["values"]][["ppt"]],
      time = tmp_daily[["time"]],
      out_labels = "PPT_mm",
      include_year = include_year
    )
  }

  res
}


metric_Tmean_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    tmp_daily <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Day",
      sw2_outs = "TEMP",
      sw2_vars = "avg_C",
      varnames_are_fixed = TRUE
    )

    res[[k1]] <- format_daily_to_matrix(
      x = tmp_daily[["values"]][["avg_C"]],
      time = tmp_daily[["time"]],
      out_labels = "Tmean_C",
      include_year = include_year
    )
  }

  res
}

get_SWA_daily <- function(
  path, name_sw2_run, id_scen_used,
  list_years_scen_used,
  out_label,
  include_year = FALSE,
  soils,
  used_depth_range_cm = NULL,
  SWP_limit_MPa = -Inf,
  ...
) {
  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    tmp_daily <- calc_SWA(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      soils = soils,
      SWP_limit_MPa = SWP_limit_MPa,
      used_depth_range_cm = used_depth_range_cm
    )

    res[[k1]] <- format_daily_to_matrix(
      x = 10 * tmp_daily[["values"]][[1]],
      time = tmp_daily[["time"]],
      out_labels = out_label,
      include_year = include_year
    )
  }

  res
}


metric_SWAat0to020cm39bar_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_daily(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat0to020cm39bar_mm",
    soils = soils,
    used_depth_range_cm = c(0, 20),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}

metric_SWAat0to100cm39bar_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_daily(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat0to100cm39bar_mm",
    soils = soils,
    used_depth_range_cm = c(0, 100),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}


metric_SWAat20to100cm39bar_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_daily(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat20to100cm39bar_mm",
    soils = soils,
    used_depth_range_cm = c(20, 100),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}

metric_SWAat20to040cm39bar_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_daily(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat20to040cm39bar_mm",
    soils = soils,
    used_depth_range_cm = c(20, 40),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}


metric_SWAat40to060cm39bar_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_daily(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat40to060cm39bar_mm",
    soils = soils,
    used_depth_range_cm = c(40, 60),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}


metric_SWAat60to080cm39bar_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_daily(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat60to080cm39bar_mm",
    soils = soils,
    used_depth_range_cm = c(60, 80),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}


metric_SWAat80to100cm39bar_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA_daily(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    out_label = "SWAat80to100cm39bar_mm",
    soils = soils,
    used_depth_range_cm = c(80, 100),
    SWP_limit_MPa = -3.9,
    include_year = include_year
  )
}


metric_DR_daily <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    tmp_daily <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Day",
      sw2_outs = "DEEPSWC",
      sw2_vars = "lowLayerDrain_cm",
      varnames_are_fixed = TRUE
    )

    res[[k1]] <- format_daily_to_matrix(
      x = 10 * tmp_daily[["values"]][["lowLayerDrain_cm"]],
      time = tmp_daily[["time"]],
      out_labels = "DR_mm",
      include_year = include_year
    )
  }

  res
}
