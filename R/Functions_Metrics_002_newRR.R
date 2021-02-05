
#------ New metrics ------

get_rh <- function(path, name_sw2_run, id_scen, years) {
  # Extract a variable from outputs as template
  res <- extract_from_sw2(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen = id_scen,
    years = years,
    sw2_tp = "Day",
    sw2_outs = "TEMP",
    sw2_vars = "avg_C",
    varnames_are_fixed = TRUE
  )

  # Provide correct name and initialize
  res[["values"]][[1]][] <- NA
  names(res[["values"]]) <- "rh"

  # Extract RH (from inputs)
  sim_input <- new.env(parent = emptyenv())
  load(
    file = file.path(path, name_sw2_run, "sw_input.RData"),
    envir = sim_input
  )

  # nolint start
  mm <- rSOILWAT2::swCloud_Humidity(
    sim_input[["swRunScenariosData"]][[id_scen]]
  )

  # Interpolate from monthly normals to daily values
  mon <- seq_len(12)
  sp <- splines::periodicSpline(mm ~ mon, period = 12)
  # nolint end

  for (yr in unique(res[["time"]][, "Year"])) {
    ids <- res[["time"]][, "Year"] %in% yr
    doys <- res[["time"]][ids, "Day"]
    res[["values"]][["rh"]][ids] <- predict(
      sp,
      doys * 12 / length(doys)
    )[["y"]]
  }

  res
}


get_vpd <- function(
  path, name_sw2_run,
  id_scen, years, group_by_month, first_month_of_year, ...
) {

  if (!missing(years)) {
    warning("'years' is not implemented but provided as argument!")
  }

  temp_min <- get_values_from_sw2(
    id_scen = id_scen,
    path, name_sw2_run,
    group_by_month = seq_len(12), #rep(0, 12),
    first_month_of_year = first_month_of_year,
    sw2_tp = "Day",
    sw2_out = "TEMP",
    sw2_var = "min_C",
    varnames_are_fixed = TRUE
  )

  temp_max <- get_values_from_sw2(
    id_scen = id_scen,
    path, name_sw2_run,
    group_by_month = seq_len(12), #rep(0, 12),
    first_month_of_year = first_month_of_year,
    sw2_tp = "Day",
    sw2_out = "TEMP",
    sw2_var = "max_C",
    varnames_are_fixed = TRUE
  )

  rh <- get_rh(path, name_sw2_run, id_scen = id_scen)


  cbind(
    rh[, 1:2],
    rSFSW2:::vpd(
      Tmin = temp_min[["vals"]][[1]],
      Tmax = temp_max[["vals"]][[1]],
      RHmean = rh[, "rh"]
    )
  )
}



#--- CWD = climatic water deficit [mm] = PET - ET
metric_CWD <- function(
  path, name_sw2_run,
  id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    # Daily PET and ET
    tmp_daily <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Day",
      sw2_outs = c("PET", "AET"),
      sw2_vars = c("pet_cm", "evapotr_cm"),
      varnames_are_fixed = TRUE
    )

    cwd_daily <- list(
      time = tmp_daily[["time"]],
      values = list(
        cwd = 10 * (
          tmp_daily[["values"]][["pet_cm"]] -
          tmp_daily[["values"]][["evapotr_cm"]]
        )
      )
    )


    res[[k1]] <- t(get_new_yearly_aggregations(
      x_daily = cwd_daily,
      # (Monthly) mean air temperature
      temp_monthly = extract_from_sw2(
        path = path,
        name_sw2_run = name_sw2_run,
        id_scen = id_scen_used[k1],
        years = list_years_scen_used[[k1]],
        sw2_tp = "Month",
        sw2_outs = "TEMP",
        sw2_vars = "avg_C",
        varnames_are_fixed = TRUE
      ),
      fun_time = sum,
      fun_extreme = max,
      output = c(
        "values", "seasonal_variability", "seasonality",
        "extreme_mean010day"
      )
    ))
  }

  res
}


#--- MDD = Soil moisture degree days [C x day] =
# degree days where soil water potential lt or gt limit

# wet based on any(SWC[i] > SWC_limit)
# dry based on all(SWC[i] < SWC_limit)
get_wetdry <- function(
  path, name_sw2_run, id_scen,
  years,
  soils,
  used_depth_range_cm = NULL,
  sm_periods = list(op = `>`, limit = -Inf)
) {

  #--- Soil moisture (SWP)
  sm_daily <- extract_from_sw2(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen = id_scen,
    years = years,
    sw2_tp = "Day",
    sw2_outs = "SWPMATRIC",
    sw2_vars = "Lyr",
    varnames_are_fixed = FALSE
  )

  widths_cm <- calc_soillayer_weights(
    soil_depths_cm = soils[["depth_cm"]],
    used_depth_range_cm = used_depth_range_cm
  )

  id_slyrs <- which(!is.na(widths_cm))
  widths_cm <- widths_cm[id_slyrs]


  # Days that meet soil moisture criterion
  res <- sm_daily
  sm <- do.call(
    what = sm_periods[["op"]],
    args = list(
      - 1 / 10 * sm_daily[["values"]][[1]][, id_slyrs, drop = FALSE],
      sm_periods[["limit"]]
    )
  )

  res[["values"]][[1]] <- if (do.call(sm_periods[["op"]], list(1, 0))) {
    # wet: op = `>` --> wet(profile) = any(wet[i])
    apply(sm, 1, any)
  } else {
    # dry: op = `<` --> dry(profile) = all(dry[i])
    apply(sm, 1, all)
  }

  res
}

# v1b: Temp > limit & sm <> limit & snow == 0
# wet based on any(SWC[i] > SWC_limit)
# dry based on all(SWC[i] < SWC_limit)
get_MDD <- function(
  path, name_sw2_run, id_scen,
  years,
  soils,
  used_depth_range_cm = NULL,
  t_periods = list(op = `>`, limit = 5),
  sm_periods = list(op = `>`, limit = Inf),
  snow_periods = list(op = `<=`, limit = 0)
) {

  # Daily wet/dry conditions
  sm <- get_wetdry(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen = id_scen,
    years = years,
    soils = soils,
    used_depth_range_cm = used_depth_range_cm,
    sm_periods = sm_periods
  )

  #--- (Daily) mean air temperature and snow cover
  tmp_daily <- extract_from_sw2(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen = id_scen,
    years = years,
    sw2_tp = "Day",
    sw2_outs = c("TEMP", "SNOWPACK"),
    sw2_vars = c("avg_C", "snowpackWaterEquivalent_cm"),
    varnames_are_fixed = TRUE
  )

  # Days that meet air temperature criterion
  dg <- do.call(
    what = t_periods[["op"]],
    args = list(
      tmp_daily[["values"]][["avg_C"]],
      t_periods[["limit"]]
    )
  )

  # Days that meet snow criterion
  snw <- do.call(
    what = snow_periods[["op"]],
    args = list(
      tmp_daily[["values"]][["snowpackWaterEquivalent_cm"]],
      snow_periods[["limit"]]
    )
  )

  # Temperature when all criteria are met
  ids <- sm[["values"]][[1]] & dg & snw
  mdd <- rep(0, length(ids))
  mdd[ids] <- tmp_daily[["values"]][["avg_C"]][ids] - t_periods[["limit"]]

  list(
    time = tmp_daily[["time"]],
    values = list(mdd = mdd)
  )
}


#--- TDD = Total degree days [C x day] = MDD(op = `>`, limit_MPa = -Inf)
calc_TDD <- function(
  path, name_sw2_run, id_scen_used,
  list_years_scen_used,
  soils,
  used_depth_range_cm = NULL,
  Temp_limit_C = 5,
  ...
) {
  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    tdd_daily <- get_MDD(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      soils = soils,
      used_depth_range_cm = used_depth_range_cm,
      t_periods = list(op = `>`, limit = Temp_limit_C),
      sm_periods = list(op = `>`, limit = -Inf)
    )

    res[[k1]] <- t(get_new_yearly_aggregations(
      x_daily = tdd_daily,
      fun_time = sum,
      fun_extreme = max,
      periods = list(op = `>`, limit = 0),
      output = c(
        "values", "seasonal_variability",
        "extreme_duration_consecutive_periods_days"
      )
    ))
  }

  res
}

metric_TDDat5C <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils, ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = "depth_cm"
  ))

  calc_TDD(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    soils = soils,
    Temp_limit_C = 5,
    ...
  )
}


#--- WDD = Wet degree days [C x day] = MDD(op = `>`, limit_MPa = -1.5 MPa)
calc_WDD <- function(
  path, name_sw2_run, id_scen_used,
  list_years_scen_used,
  soils,
  used_depth_range_cm = NULL,
  Temp_limit_C = 5,
  SWP_limit_MPa = -1.5, ...
) {

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    wdd_daily <- get_MDD(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      soils = soils,
      used_depth_range_cm = used_depth_range_cm,
      t_periods = list(op = `>`, limit = Temp_limit_C),
      sm_periods = list(op = `>`, limit = SWP_limit_MPa)
    )

    res[[k1]] <- t(get_new_yearly_aggregations(
      x_daily = wdd_daily,
      # (Monthly) mean air temperature
      temp_monthly = extract_from_sw2(
        path = path,
        name_sw2_run = name_sw2_run,
        id_scen = id_scen_used[k1],
        years = list_years_scen_used[[k1]],
        sw2_tp = "Month",
        sw2_outs = "TEMP",
        sw2_vars = "avg_C",
        varnames_are_fixed = TRUE
      ),
      fun_time = sum,
      output = c("values", "seasonality")
    ))
  }

  res
}

metric_WDDat5C0to100cm15bar <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils, ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = "depth_cm"
  ))

  calc_WDD(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    soils = soils,
    Temp_limit_C = 5,
    SWP_limit_MPa = -1.5,
    used_depth_range_cm = c(0, 100),
    ...
  )
}


#--- DDD = Dry degree days [C x day] = MDD(op = `<`, limit_MPa = -3.0 MPa)
calc_DDD <- function(
  path, name_sw2_run, id_scen_used,
  list_years_scen_used,
  soils,
  used_depth_range_cm = NULL,
  Temp_limit_C = 5,
  SWP_limit_MPa = -3,
  output = c("values", "extreme_value_consecutive_periods"),
  ...
) {
  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    ddd_daily <- get_MDD(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      soils = soils,
      used_depth_range_cm = used_depth_range_cm,
      t_periods = list(op = `>`, limit = Temp_limit_C),
      sm_periods = list(op = `<`, limit = SWP_limit_MPa)
    )

    res[[k1]] <- t(get_new_yearly_aggregations(
      x_daily = ddd_daily,
      fun_time = sum,
      fun_extreme = max,
      periods = list(op = `>`, limit = 0),
      output = output
    ))
  }

  res
}

metric_DDDat5C0to030cm30bar <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils,
  ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = "depth_cm"
  ))

  calc_DDD(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    soils = soils,
    Temp_limit_C = 5,
    SWP_limit_MPa = -3,
    used_depth_range_cm = c(0, 30),
    output = "extreme_value_consecutive_periods",
    ...
  )
}

metric_DDDat5C0to100cm30bar <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils, ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = "depth_cm"
  ))

  calc_DDD(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    soils = soils,
    Temp_limit_C = 5,
    SWP_limit_MPa = -3,
    used_depth_range_cm = c(0, 100),
    output = c("values", "extreme_value_consecutive_periods"),
    ...
  )
}


#--- SWA = Soil water availability [mm] (SWA above -3.0 and -3.9 MPa):
# calculate monthly means
calc_SWA <- function(
  path, name_sw2_run, id_scen,
  years,
  soils,
  used_depth_range_cm = NULL,
  SWP_limit_MPa = -Inf
) {

  #--- Soil moisture
  sm_daily <- extract_from_sw2(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen = id_scen,
    years = years,
    sw2_tp = "Day",
    sw2_outs = "VWCMATRIC",
    sw2_vars = "Lyr",
    varnames_are_fixed = FALSE
  )

  widths_cm <- calc_soillayer_weights(
    soil_depths_cm = soils[["depth_cm"]],
    used_depth_range_cm = used_depth_range_cm
  )

  id_slyrs <- which(!is.na(widths_cm))

  if (length(id_slyrs) > 0) {
    widths_cm <- widths_cm[id_slyrs]

    # Determine SWA for each soil layer
    base_SWCmatric_mm <- if (is.finite(SWP_limit_MPa)) {
      10 * widths_cm * rSOILWAT2::SWPtoVWC(
        swp = SWP_limit_MPa,
        sand = soils[["sand_frac"]][id_slyrs],
        clay = soils[["clay_frac"]][id_slyrs]
      )
    } else {
      rep(0, length(id_slyrs))
    }

    swc_matric_mm <- 10 * sweep(
      x = sm_daily[["values"]][[1]][, id_slyrs, drop = FALSE],
      MARGIN = 2,
      STATS = widths_cm,
      FUN = "*"
    )

    swa_by_layer <- sweep(
      x = swc_matric_mm,
      MARGIN = 2,
      STATS = base_SWCmatric_mm,
      FUN = "-"
    )
    swa_by_layer[swa_by_layer < 0] <- 0

    # Sum SWA across soil profile
    swa_daily_values <- apply(
      X = swa_by_layer,
      MARGIN = 1,
      FUN = sum
    )

  } else {
    # No soil layers in the depth range
    swa_daily_values <- rep(NA, nrow(sm_daily[["time"]]))
  }

  res <- sm_daily
  res[["values"]][[1]] <- swa_daily_values
  res
}



get_SWA <- function(
  path, name_sw2_run, id_scen_used,
  list_years_scen_used,
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

    res[[k1]] <- t(get_new_yearly_aggregations(
      x_daily = swa_daily,
      # (Monthly) mean air temperature
      temp_monthly = extract_from_sw2(
        path = path,
        name_sw2_run = name_sw2_run,
        id_scen = id_scen_used[k1],
        years = list_years_scen_used[[k1]],
        sw2_tp = "Month",
        sw2_outs = "TEMP",
        sw2_vars = "avg_C",
        varnames_are_fixed = TRUE
      ),
      fun_time = mean,
      output = c("values", "seasonal_variability", "seasonality")
    ))
  }

  res
}


metric_SWAat0to100cm30bar <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils, ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    soils = soils,
    SWP_limit_MPa = -3,
    used_depth_range_cm = c(0, 100),
    ...
  )
}

metric_SWAat0to100cm39bar <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils, ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = c("depth_cm", "sand_frac", "clay_frac")
  ))

  get_SWA(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    soils = soils,
    SWP_limit_MPa = -3.9,
    used_depth_range_cm = c(0, 100),
    ...
  )
}


#--- DSI = Duration of dry soil intervals at -1.5 and -3.0 SWP thresholds

# dry based on MDD
calc_DSI <- function(
  path, name_sw2_run, id_scen_used,
  list_years_scen_used,
  soils,
  used_depth_range_cm = NULL,
  SWP_limit_MPa = -Inf,
  fun_periods = function(x) c(max = max(x), mean = mean(x), N = length(x)),
  include_year = FALSE,
  ...
) {
  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    dry_daily <- get_wetdry(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      soils = soils,
      used_depth_range_cm = used_depth_range_cm,
      sm_periods = list(op = `<`, limit = SWP_limit_MPa)
    )

    tmp <- lapply(
      X = tapply(
        X = dry_daily[["values"]][[1]],
        INDEX = dry_daily[["time"]][, "Year"],
        FUN = function(x) {
          tmp <- rle(x)
          if (any(tmp[["values"]] == 1)) {
            tmp[["lengths"]][tmp[["values"]]]
          } else {
            0
          }
        }
      ),
      FUN = fun_periods
    )

    tmp2 <- array(
      unlist(tmp),
      dim = c(length(tmp[[1]]), length(tmp)),
      dimnames = list(names(tmp[[1]]), NULL)
    )

    res[[k1]] <- if (include_year) {
      rbind(
        Year = as.integer(names(tmp)),
        tmp2
      )
    } else {
      tmp2
    }
  }

  res
}

metric_DSIat0to100cm15bar <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils, ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = "depth_cm"
  ))

  calc_DSI(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    soils = soils,
    SWP_limit_MPa = -1.5,
    used_depth_range_cm = c(0, 100),
    fun_periods = function(x) c(mean = mean(x), N = length(x))
  )
}

metric_DSIat0to100cm30bar <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils, ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = "depth_cm"
  ))

  calc_DSI(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    soils = soils,
    SWP_limit_MPa = -3.0,
    used_depth_range_cm = c(0, 100),
    fun_periods = function(x) c(max = max(x))
  )
}


#--- Frost: Consistency (across years) of last and first frost events
# (StDev of DOY)
# For the last event before mid-year and first event after mid-year where
# mid-year = July 15 (circa summer solstice + 1 month)
get_frost_doy <- function(
  path, name_sw2_run, id_scen, years,
  Temp_limit_C = -5,
  hemisphere_NS = c("N", "S"),
  include_year = FALSE,
  ...
) {
  hemisphere_NS <- match.arg(hemisphere_NS)
  stopifnot(hemisphere_NS == "N") #TODO: implement for southern hemisphere

  #--- (Daily) minimum air temperature
  temp_daily <- extract_from_sw2(
    path = path,
    name_sw2_run = name_sw2_run,
    years = years,
    id_scen = id_scen,
    sw2_tp = "Day",
    sw2_outs = "TEMP",
    sw2_vars = "min_C",
    varnames_are_fixed = TRUE
  )

  is_frost <- temp_daily[["values"]][[1]] < Temp_limit_C

  # Mid-year: summer solstice + 1 month
  # North: June solstice (Jun 20-22 = 171-173)
  # South: December solstice (Dec 20-23 = 354-357)
  doy_mid <- 196 # July 15 (in non-leap year)

  res <- as.matrix(aggregate(
    x = is_frost,
    by = list(Year = temp_daily[["time"]][, "Year"]),
    function(x) {
      tmp <- which(x)
      c(
        last = max(tmp[tmp <= doy_mid]),
        first = min(tmp[tmp > doy_mid])
      )
    }
  ))
  colnames(res) <- c("Year", "LastFrost", "FirstFrost")

  if (include_year) res else res[, -1]
}

metric_FrostDaysAtNeg5C <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    res[[k1]] <- t(get_frost_doy(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      Temp_limit_C = -5
    ))
  }

  res
}

# CorTempPPT: Correlation between monthly temperature and
# precipitation by year
metric_CorTempPPT <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    vals_monthly <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Month",
      sw2_outs = c("TEMP", "PRECIP"),
      sw2_vars = c("min_C", "ppt"),
      varnames_are_fixed = TRUE
    )

    tmp <- as.vector(by(
      data = cbind(
        vals_monthly[["values"]][["min_C"]],
        vals_monthly[["values"]][["ppt"]]
      ),
      INDICES = vals_monthly[["time"]][, "Year"],
      FUN = function(x) cor(x[, 1], x[, 2])
    ))

    res[[k1]] <- if (include_year) {
      rbind(
        Year = unique(vals_monthly[["time"]][, "Year"]),
        seasonality = tmp
      )
    } else {
      rbind(seasonality = tmp)
    }
  }

  res
}


#--- Plant recruitment index:
calc_RecruitmentIndex_v2 <- function(
  path, name_sw2_run, id_scen,
  years,
  soils,
  hemisphere_NS = c("N", "S"),
  recruitment_depth_range_cm = c(5, 30),
  Temp_limit_C = 5,
  Wet_SWP_limit_MPa = -1.5,
  Dry_SWP_limit_MPa = -3,
  init_WDD = 15,
  init_days = 3,
  init_depth_range_cm = c(0, 5),
  stop_DDD = 0,
  stop_days_DDD = 0,
  stop_depth_range_cm = c(0, 30),
  stop_TDD = 0,
  stop_days_TDD = 0,
  include_year = FALSE,
  ...
) {
  stopifnot(requireNamespace("zoo", quietly = TRUE))

  hemisphere_NS <- match.arg(hemisphere_NS)
  stopifnot(hemisphere_NS == "N") #TODO: implement for southern hemisphere

  # Mid-year: summer solstice + 1 month
  # North: June solstice (Jun 20-22 = 171-173)
  # South: December solstice (Dec 20-23 = 354-357)
  doy_mid <- 196 # July 15 (in non-leap year)

  # WDD that initiates a recruitment period (germination window)
  wdd_start <- get_MDD(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen = id_scen,
    years = years,
    soils = soils,
    used_depth_range_cm = init_depth_range_cm,
    t_periods = list(op = `>`, limit = Temp_limit_C),
    sm_periods = list(op = `>`, limit = Wet_SWP_limit_MPa)
  )

  # WDD for recruitment
  wdd_recruit <- get_MDD(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen = id_scen,
    years = years,
    soils = soils,
    used_depth_range_cm = recruitment_depth_range_cm,
    t_periods = list(op = `>`, limit = Temp_limit_C),
    sm_periods = list(op = `>`, limit = Wet_SWP_limit_MPa)
  )

  N_days <- length(wdd_recruit[["values"]][[1]])

  # DDD that stops a recruitment period
  ddd_stop <- get_MDD(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen = id_scen,
    years = years,
    soils = soils,
    used_depth_range_cm = stop_depth_range_cm,
    t_periods = list(op = `>`, limit = Temp_limit_C),
    sm_periods = list(op = `<`, limit = Dry_SWP_limit_MPa)
  )

  # (Absence of) TDD that stops a recruitment period
  tdd_nostop <- get_MDD(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen = id_scen,
    years = years,
    soils = soils,
    t_periods = list(op = `>`, limit = Temp_limit_C),
    sm_periods = list(op = `>`, limit = -Inf)
  )


  # List of possible start days
  if (init_WDD <= 0) {
    ids_start <- which(wdd_start[["values"]][[1]] > 0)

  } else {
    # (i) after `init_days` with WDD
    tmp1a <- zoo::rollsum(
      wdd_start[["values"]][[1]] > 0,
      k = init_days,
      fill = 0,
      align = "right"
    ) >= init_days

    # (ii) and with a sum of `init_WDD`
    tmp1b <- zoo::rollsum(
      wdd_start[["values"]][[1]],
      k = init_days,
      fill = 0,
      align = "right"
    ) >= init_WDD

    ids_start <- 1 + which(tmp1a & tmp1b)
    tmp <- length(ids_start)
    if (tmp > 0 && ids_start[tmp] > N_days) {
      ids_start[tmp] <- ids_start[tmp] - 1
      ids_start <- unique(ids_start)
    }
  }

  # List of end/stop days due to DDD
  if (stop_DDD <= 0) {
    ids_stop_DDD <- which(ddd_stop[["values"]][[1]] > 0)

  } else {
    # (i) after `stop_days_DDD` with DDD
    tmp2a <- zoo::rollsum(
      ddd_stop[["values"]][[1]] > 0,
      k = stop_days_DDD,
      fill = 0,
      align = "right"
    ) >= stop_days_DDD

    # (ii) and with a sum of `stop_DDD`
    tmp2b <- zoo::rollsum(
      ddd_stop[["values"]][[1]],
      k = stop_days_DDD,
      fill = 0,
      align = "right"
    ) >= stop_DDD

    ids_stop_DDD <- 1 + which(tmp2a & tmp2b)
    tmp <- length(ids_stop_DDD)
    if (tmp > 0 && ids_stop_DDD[tmp] > N_days) {
      ids_stop_DDD[tmp] <- ids_stop_DDD[tmp] - 1
      ids_stop_DDD <- unique(ids_stop_DDD)
    }
  }


  # List of end/stop days due to (absence of) TDD
  if (stop_TDD <= 0 && stop_days_TDD < 1) {
    ids_stop_TDD <- which(tdd_nostop[["values"]][[1]] <= 0)

  } else {
    # (i) after `stop_days_TDD` with TDD
    tmp3a <- zoo::rollsum(
      tdd_nostop[["values"]][[1]] <= 0,
      k = stop_days_TDD,
      fill = 0,
      align = "right"
    ) >= stop_days_TDD

    # (ii) and with a sum of `stop_TDD`
    tmp3b <- zoo::rollsum(
      tdd_nostop[["values"]][[1]],
      k = stop_days_TDD,
      fill = 0,
      align = "right"
    ) <= stop_TDD

    ids_stop_TDD <- 1 + which(tmp3a & tmp3b)
    tmp <- length(ids_stop_TDD)
    if (tmp > 0 && ids_stop_TDD[tmp] > N_days) {
      ids_stop_TDD[tmp] <- ids_stop_TDD[tmp] - 1
      ids_stop_TDD <- unique(ids_stop_TDD)
    }
  }


  # Combine all stopping days and add day after last simulated day as end day
  ids_stop <- unique(sort(c(ids_stop_DDD, ids_stop_TDD)))
  ids_stop[length(ids_stop) + 1] <- 1 + length(ddd_stop[["values"]][[1]])


  # List start/end of all suitable periods
  periods <- list()

  k0 <- 1
  for (k1 in seq_along(ids_start)) {
    # Identify start day and locate earliest stop day
    tmp1 <- ids_start[k1]
    periods[[k0]] <- c(
      start = tmp1,
      end = min(ids_stop[ids_stop >= tmp1])
    )
    k0 <- k0 + 1
  }

  periods <- do.call(rbind, periods)


  # Recruitment potential: sum of WDD within suitable soil depths
  ts_years <- unique(wdd_recruit[["time"]][, "Year"])
  jan0 <- as.Date(paste0(ts_years[1] - 1, "-12-31"))

  res <- array(
    data = 0,
    dim = c(length(ts_years), 6 + as.integer(include_year)),
    dimnames = list(NULL,
      c(
        if (include_year) "Year",
        paste0(
          rep(c("Spring", "Fall"), each = 3),
          "Recruitment_",
          rep(c("maxWDD", "DOY", "DurationDays"), times = 2)
        )
      )
    )
  )

  if (include_year) {
    res[, "Year"] <- ts_years
  }

  # Loop over years
  for (k1 in seq_along(ts_years)) {
    doy_mid_lyr <- doy_mid + rSW2utils::isLeapYear(ts_years[k1])

    # Identify which periods start or end during current year
    ids_yr <- which(wdd_recruit[["time"]][, "Year"] == ts_years[k1])
    ids_periods <- which(
      periods[, "start"] %in% ids_yr | periods[, "end"] %in% ids_yr
    )

    # Loop over periods in current year
    for (k2 in seq_along(ids_periods)) {
      # Identify start/end of current period
      id_mid_yr <- ids_yr[1] + doy_mid_lyr - 1
      lims <- c(
        max(ids_yr[1], periods[ids_periods[k2], "start"]),
        min(ids_yr[length(ids_yr)], periods[ids_periods[k2], "end"])
      )

      # Identify maximum (cumulative) WDD (and starting DOY) of
      # periods in current year
      if (lims[1] < id_mid_yr && lims[2] >= id_mid_yr) {
        # Current period crosses mid-year date
        ids1 <- seq(from = lims[1], to = id_mid_yr - 1)
        ids2 <- seq(from = id_mid_yr, to = lims[2])
        tmp <- c(
          sum(wdd_recruit[["values"]][[1]][ids1]),
          sum(wdd_recruit[["values"]][[1]][ids2])
        )

        if (tmp[1] > res[k1, "SpringRecruitment_maxWDD"]) {
          res[k1, "SpringRecruitment_maxWDD"] <- tmp[1]
          res[k1, "SpringRecruitment_DOY"] <-
            as.POSIXlt(jan0 + lims[1])$yday + 1
          res[k1, "SpringRecruitment_DurationDays"] <- length(ids1)
        }

        if (tmp[2] > res[k1, "FallRecruitment_maxWDD"]) {
          res[k1, "FallRecruitment_maxWDD"] <- tmp[2]
          res[k1, "FallRecruitment_DOY"] <-
            as.POSIXlt(jan0 + id_mid_yr)$yday + 1
          res[k1, "FallRecruitment_DurationDays"] <- length(ids2)
        }

      } else {
        ids <- seq(from = lims[1], to = lims[2])
        tmp <- sum(wdd_recruit[["values"]][[1]][ids])

        if (all(lims < id_mid_yr)) {
          # Current period is completely before mid-year date
          if (tmp > res[k1, "SpringRecruitment_maxWDD"]) {
            # Current spring period is larger than previous ones -> replace
            res[k1, "SpringRecruitment_maxWDD"] <- tmp
            res[k1, "SpringRecruitment_DOY"] <-
              as.POSIXlt(jan0 + lims[1])$yday + 1
            res[k1, "SpringRecruitment_DurationDays"] <- length(ids)
          }

        } else {
          # Current period is completely after mid-year date
          if (tmp > res[k1, "FallRecruitment_maxWDD"]) {
            # Current fall period is larger than previous ones -> replace
            res[k1, "FallRecruitment_maxWDD"] <- tmp
            res[k1, "FallRecruitment_DOY"] <-
              as.POSIXlt(jan0 + lims[1])$yday + 1
            res[k1, "FallRecruitment_DurationDays"] <- length(ids)
          }
        }
      }
    }
  }

  res[res == 0] <- NA

  res
}


# max WDD across recruitment events before/after mid-summer (July 15): where
# recruitment potential is accumulated WDD at 5-20 cm during intervals which
# start after 3-day periods with WDD > 0 that sum to >= 15 WDD in 0-5 cm and
# end either after 3-day periods with DDD > 0 that sum to >= 15 DDD in 0-20 cm
# or after 3-day periods with TDD == 0
metric_RecruitmentIndex_v4 <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  soils, ...
) {
  stopifnot(check_metric_arguments(
    out = match.arg(out),
    req_soil_vars = "depth_cm"
  ))

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    res[[k1]] <- t(calc_RecruitmentIndex_v2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      soils = soils,
      hemisphere_NS = "N",
      recruitment_depth_range_cm = c(5, 20),
      Temp_limit_C = 5,
      Wet_SWP_limit_MPa = -1.5,
      Dry_SWP_limit_MPa = -3,
      init_WDD = 15,
      init_days = 3,
      init_depth_range_cm = c(0, 5),
      stop_DDD = 15,
      stop_days_DDD = 3,
      stop_depth_range_cm = c(0, 20),
      stop_TDD = 0,
      stop_days_TDD = 3
    ))
  }

  res
}



calc_peaks_v1 <- function(
  x,
  bufferx_size = 0,
  peak_type = c("value", "volume"),
  days_window = 91,
  peak2_drop_factor = 2 / 3,
  peak2_increase_factor = 1 / 6,
  peak2_min_factor = 1 / 3
) {
  # Identify DOY of candidate peaks
  # (peaks are at least a half-window apart)
  pids <- identify_peaks(x, days_window)

  if (length(pids) > 0) {
    # Calculate peak size
    psize <- peak_size_v1(pids, peak_type, values = x)

    pvals <- switch(
      EXPR = peak_type,
      value = psize,
      volume = x[pids]
    )

    # Primary peak: largest size
    tmp <- which.max(psize)
    p1 <- c(pids[tmp], pvals[tmp])
    # Secondary peaks are at least `peak2_min_factor * p1[2]` tall
    is2 <- pids != p1[1] & pvals >= peak2_min_factor * p1[2]

    # peak 2
    p2 <- rep(NA, 2)

    if (sum(is2) > 0) {
      pids2 <- pids[is2]
      ppos <- order(psize[is2], decreasing = TRUE)

      # Loop over secondary peak candidates (descending values)
      k <- 1
      while (k <= length(pids2)) {
        tmp_p2 <- c(pids2[ppos[k]], pvals[is2][ppos[k]])

        # Days in-between peak 1 and candidate peak 2
        ids_inbetween <- sort(tmp_p2[1]:p1[1])

        if (length(ids_inbetween) > 0) {
          # Sufficient drop?
          ids_dropped <- x[ids_inbetween] < peak2_drop_factor * p1[2]
          has_dropped <- any(ids_dropped)

          if (has_dropped) {
            # Sufficient increase?
            increase <- tmp_p2[2] - min(x[ids_inbetween[ids_dropped]])
            has_increased <- increase >= peak2_increase_factor * p1[2]

            # Is candidate our peak 2?
            if (has_increased) {
              p2 <- tmp_p2
              break
            }
          }
        }

        k <- k + 1
      }
    }
    c(p1, p2)
  } else {
    rep(NA, 4)
  }
}


calc_peaks_v2 <- function(
  x,
  bufferx_size = 0,
  peak_type = c("value", "volume"),
  days_window = 91,
  peaks_drop_factor = 2 / 3,
  peaks_increase_factor = 1 / 6,
  peaksize_min_factor = 1 / 3,
  N_peaks_max_reported = 4
) {
  vars <- c(
    "Peak_N",
    paste0(
      "Peak",
      rep(seq_len(4), each = 2),
      rep(c("_DOY", "_rvol"), times = N_peaks_max_reported)
    )
  )
  res <- rep(NA, length(vars))
  names(res) <- vars

  #--- Identify DOY of candidate peaks
  # (peaks are at least a half-window apart)
  pids <- identify_peaks(x, days_window)
  res["Peak_N"] <- length(pids)

  if (res["Peak_N"] > 0) {
    #--- Peak value, relative volume, and size (as determined by `peak_type`)
    tmp <- cbind(
      ID = seq_along(pids),
      DOY = pids,
      value = x[pids],
      volume = peak_size_v2(
        pids,
        peak_type = "volume",
        values = x,
        window = days_window
      ) / sum(x, na.rm = TRUE)
    )
    peaks <- cbind(
      tmp,
      meas = switch(
        EXPR = peak_type,
        value = tmp[, "value"],
        volume = tmp[, "volume"]
      )
    )

    # ID of largest peak by size
    idmax <- peaks[which.max(peaks[, "meas"]), "ID"]
    tmp <- peaks[, "ID"] == idmax
    maxpeak_size <- peaks[tmp, "meas"]
    maxpeak_value <- peaks[tmp, "value"]

    # Retain sufficiently sized candidate peaks
    tmp <- peaks[, "meas"] >= peaksize_min_factor * maxpeak_size
    peaks <- peaks[tmp, , drop = FALSE]

    res["Peak_N"] <- nrow(peaks)

    if (res["Peak_N"] > 0) {
      # Growing list of identified/confirmed peak IDs
      ids_good <- idmax

      #--- Peaks need sufficient drop and increase from neighboring peaks
      # (loop over peaks in descending size)
      while (nrow(peaks) > 0 && length(ids_good) < nrow(peaks)) {
        # Largest sized peak among candidate peaks
        cid <- peaks[-ids_good, "ID"][which.max(peaks[-ids_good, "meas"])]

        tmp <- peaks[, "ID"] == cid
        cdoy <- peaks[tmp, "DOY"]
        cvalue <- peaks[tmp, "value"]

        tmp <- cid + c(-1, 1)
        ids_neighs <- tmp[tmp %in% peaks[, "ID"]]

        # Sufficiently large drop and increase?
        is_good <- rep(FALSE, length(ids_neighs))

        # Loop over neighboring peaks and check if candidate is a good peak
        for (k in seq_along(ids_neighs)) {

          # Days in-between current candidate peak cid and neighboring peak
          ids_inbetween <- sort(seq.int(
            from = cdoy,
            to = peaks[peaks[, "ID"] == ids_neighs[k], "DOY"]
          ))

          if (length(ids_inbetween) > 0) {
            # Sufficient drop?
            ids_dropped <-
              x[ids_inbetween] < peaks_drop_factor * maxpeak_value
            has_dropped <- any(ids_dropped)

            if (has_dropped) {
              # And sufficient increase?
              increase <- cvalue - min(x[ids_inbetween[ids_dropped]])

              # Is candidate a good peak in relation to neighbor k?
              is_good[k] <-
                increase >= peaks_increase_factor * maxpeak_value
            }
          }
        }

        if (all(is_good)) {
          # Candidate is a good peak -> add to list of good peaks
          ids_good <- c(ids_good, cid)
        } else {
          # Candidate is not a suitable peak -> remove from peak data.frame
          peaks <- peaks[peaks[, "ID"] != cid, , drop = FALSE]
        }
      }


      #--- Update peak volume
      # (volume may have changed due to removal of unsuitable candidate peaks)
      peaks[, "volume"] <- peak_size_v2(
        peaks[, "DOY"],
        peak_type = "volume",
        values = x,
        window = days_window
      ) / sum(x, na.rm = TRUE)

      #--- Copy peak information to results
      res["Peak_N"] <- nrow(peaks)

      tmp <- seq_len(min(N_peaks_max_reported, res["Peak_N"]))
      res[paste0("Peak", tmp, "_DOY")] <- peaks[tmp, "DOY"]
      res[paste0("Peak", tmp, "_rvol")] <- peaks[tmp, "volume"]
    }
  }

  res
}


#--- DOY of two maxima of a smoothed transpiration
calc_transp_peaks_v1 <- function(
  x,
  time,
  # Identify peak size by volume (valley - peak - valley) or by value at peak
  peak_type = c("value", "volume"),
  # Smoothing width of daily transpiration (should be odd)
  days_smoothing = 15,
  # Window for determining secondary candidate peaks (should be odd)
  days_window = 91,
  # Secondary peak should be separated by a valley that
  # dropped below `peak2_drop_factor` * main peak size
  peak2_drop_factor = 2 / 3,
  # Secondary peak should rise above the valley
  # by at least `peak2_drop_factor` * main peak size
  peak2_increase_factor = 1 / 6,
  # Secondary peak should be at least `peak2_min_factor` * main peak size
  peak2_min_factor = 1 / 3
) {
  stopifnot(requireNamespace("zoo", quietly = TRUE))
  stopifnot(days_window %% 2L == 1L) # check that window is odd

  tsmoothed <- if (days_smoothing > 1) {
    zoo::rollapply(
      x,
      width = days_smoothing,
      FUN = mean,
      na.rm = TRUE,
      align = "center",
      partial = TRUE,
      fill = NA
    )
  } else {
    x
  }

  # Peak 1 = annual maximum peak (doy, value)
  # Peak 2 = largest local maximum if exists, i.e., separated from peak 1
  #   - by drop of >= 1/3 of max and
  #   - increase >= 1/6 of max
  ts_years <- unique(time)
  peaks <- matrix(NA, nrow = length(ts_years), ncol = 4)
  bufferx_size <- (days_window + 1) / 2

  for (k1 in seq_along(ts_years)) {
    # Add NA-buffer of size `bufferx_size` around values for current year
    ids_vals <- which(time == ts_years[k1])
    lims_bufft <- bufferx_size + c(1, length(ids_vals))
    tmpx <- rep(NA, length(ids_vals) + 2 * bufferx_size)
    tmpx[lims_bufft[1]:lims_bufft[2]] <- tsmoothed[ids_vals]

    tmp <- calc_peaks_v1(
      x = tmpx,
      bufferx_size = bufferx_size,
      peak_type = peak_type,
      days_window = days_window,
      peak2_drop_factor = peak2_drop_factor,
      peak2_increase_factor = peak2_increase_factor,
      peak2_min_factor = peak2_min_factor
    )

    # Subtract buffer from DOYs
    tmp[c(1, 3)] <- tmp[c(1, 3)] - bufferx_size

    peaks[k1, ] <- tmp
  }

  # Create data.frame formatted as if returned by `aggregate`
  structure(
    list(ts_years, peaks),
    class = "data.frame",
    .Names = c("Year", "x"),
    row.names = seq_along(ts_years)
  )
}


#--- Count, DOY, and relative volume of smoothed transpiration peaks
calc_transp_peaks_v2 <- function(
  x,
  time,
  # Maximum number of peaks for which to return values
  N_peaks_max_reported = 4,
  # Identify peak size by volume (valley - peak - valley) or by value at peak
  peak_type = c("volume", "value"),
  # Smoothing width of daily transpiration (should be odd)
  days_smoothing = 15,
  # Window for determining candidate peaks (should be odd)
  days_window = 91,
  # Peak values should be separated by a valley that
  # dropped below `peaks_drop_factor` * maximum peak value
  peaks_drop_factor = 2 / 3,
  # Peak values should rise above the valley
  # by at least `peaks_drop_factor` * maximum peak value
  peaks_increase_factor = 1 / 6,
  # Peaks should be at least `peaksize_min_factor` * maximum peak size
  peaksize_min_factor = 1 / 4
) {
  stopifnot(requireNamespace("zoo", quietly = TRUE))
  stopifnot(days_window %% 2L == 1L) # check that window is odd

  tsmoothed <- if (days_smoothing > 1) {
    zoo::rollapply(
      x,
      width = days_smoothing,
      FUN = mean,
      na.rm = TRUE,
      align = "center",
      partial = TRUE,
      fill = NA
    )
  } else {
    x
  }

  ts_years <- unique(time)
  peaks <- matrix(
    NA,
    nrow = length(ts_years),
    ncol = 1 + 2 * N_peaks_max_reported
  )
  bufferx_size <- (days_window + 1) / 2

  for (k1 in seq_along(ts_years)) {
    # Add NA-buffer of size `bufferx_size` around values for current year
    ids_vals <- which(time == ts_years[k1])
    lims_bufft <- bufferx_size + c(1, length(ids_vals))
    tmpx <- rep(NA, length(ids_vals) + 2 * bufferx_size)
    tmpx[lims_bufft[1]:lims_bufft[2]] <- tsmoothed[ids_vals]

    tmp <- calc_peaks_v2(
      x = tmpx,
      bufferx_size = bufferx_size,
      peak_type = peak_type,
      days_window = days_window,
      peaks_drop_factor = peaks_drop_factor,
      peaks_increase_factor = peaks_increase_factor,
      peaksize_min_factor = peaksize_min_factor,
      N_peaks_max_reported = N_peaks_max_reported
    )

    # Subtract buffer from DOYs
    ids <- grep("_DOY", names(tmp))
    tmp[ids] <- tmp[ids] - bufferx_size

    peaks[k1, ] <- tmp
  }

  colnames(peaks) <- names(tmp)

  # Create data.frame formatted as if returned by `aggregate`
  structure(
    list(ts_years, peaks),
    class = "data.frame",
    .Names = c("Year", "x"),
    row.names = seq_along(ts_years)
  )
}



calc_TranspirationSeasonality_v1 <- function(
  path, name_sw2_run, id_scen_used,
  list_years_scen_used,
  include_year = FALSE,
  # Identify peak size by volume (valley - peak - valley) or by value at peak
  peak_type = c("value", "volume"),
  # DOY of `transp_timing_probs` quantiles of cumulative transpiration
  transp_timing_probs = c(0.1, 0.5, 0.9),
  # Smoothing width of daily transpiration (should be odd)
  days_smoothing = 15,
  # Window for determining secondary candidate peaks (should be odd)
  days_window = 91,
  # Secondary peak should be separated by a valley that
  # dropped below `peak2_drop_factor` * main peak size
  peak2_drop_factor = 2 / 3,
  # Secondary peak should rise above the valley
  # by at least `peak2_drop_factor` * main peak size
  peak2_increase_factor = 1 / 6,
  # Secondary peak should be at least `peak2_min_factor` * main peak size
  peak2_min_factor = 1 / 3,
  ...
) {
  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    T_by_lyr <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Day",
      sw2_outs = "TRANSP",
      sw2_vars = "transp_total_Lyr",
      varnames_are_fixed = FALSE
    )

    transp <- 10 * rowSums(T_by_lyr[["values"]][[1]])

    # DOY of 10%, 50%, and 90% percentile of cumulative transpiration
    # Annual total transpiration (mm)
    xtmp1 <- calc_transp_seasonality(
      x = transp,
      time = T_by_lyr[["time"]][, "Year"],
      probs = transp_timing_probs
    )

    tmp1 <- t(xtmp1[["x"]])
    rownames(tmp1) <- c(
      "Transpiration_mm",
      paste0("Tcumulative_", round(100 * transp_timing_probs), "perc_DOY")
    )

    # DOY of two maxima of a smoothed transpiration
    xtmp2 <- calc_transp_peaks_v1(
      x = transp,
      time = T_by_lyr[["time"]][, "Year"],
      peak_type = peak_type,
      days_smoothing = days_smoothing,
      days_window = days_window,
      peak2_drop_factor = peak2_drop_factor,
      peak2_increase_factor = peak2_increase_factor,
      peak2_min_factor = peak2_min_factor
    )

    tmp2 <- t(xtmp2[["x"]][, c(1, 3), drop = FALSE])
    rownames(tmp2) <- paste0("Transpiration_peak", 1:2, "_DOY")


    res[[k1]] <- if (include_year) {
      rbind(Year = xtmp1[["Year"]], tmp1, tmp2)
    } else {
      rbind(tmp1, tmp2)
    }
  }

  res
}



calc_TranspirationSeasonality_v2 <- function(
  path, name_sw2_run, id_scen_used,
  list_years_scen_used,
  include_year = FALSE,
  # DOY of `transp_timing_probs` quantiles of cumulative transpiration
  transp_timing_probs = c(0.1, 0.5, 0.9),
  # Maximum number of peaks for which to return values
  N_peaks_max_reported = 4,
  # Identify peak size by volume (valley - peak - valley) or by value at peak
  peak_type = c("volume", "value"),
  # Smoothing width of daily transpiration (should be odd)
  days_smoothing = 15,
  # Window for determining candidate peaks (should be odd)
  days_window = 91,
  # Peak values should be separated by a valley that
  # dropped below `peaks_drop_factor` * maximum peak value
  peaks_drop_factor = 2 / 3,
  # Peak values should rise above the valley
  # by at least `peaks_drop_factor` * maximum peak value
  peaks_increase_factor = 1 / 6,
  # Peaks should be at least `peaksize_min_factor` * maximum peak size
  peaksize_min_factor = 1 / 4,
  ...
) {
  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    T_by_lyr <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Day",
      sw2_outs = "TRANSP",
      sw2_vars = "transp_total_Lyr",
      varnames_are_fixed = FALSE
    )

    transp <- 10 * rowSums(T_by_lyr[["values"]][[1]])

    # DOY of 10%, 50%, and 90% percentile of cumulative transpiration
    # Annual total transpiration (mm)
    xtmp1 <- calc_transp_seasonality(
      x = transp,
      time = T_by_lyr[["time"]][, "Year"],
      probs = transp_timing_probs
    )

    tmp1 <- t(xtmp1[["x"]])
    rownames(tmp1) <- c(
      "Transpiration_mm",
      paste0("Tcumulative_", round(100 * transp_timing_probs), "perc_DOY")
    )

    # Count, DOY, and relative volume of smoothed transpiration peaks
    xtmp2 <- calc_transp_peaks_v2(
      x = transp,
      time = T_by_lyr[["time"]][, "Year"],
      N_peaks_max_reported = N_peaks_max_reported,
      peak_type = peak_type,
      days_smoothing = days_smoothing,
      days_window = days_window,
      peaks_drop_factor = peaks_drop_factor,
      peaks_increase_factor = peaks_increase_factor,
      peaksize_min_factor = peaksize_min_factor
    )

    tmp2 <- t(xtmp2[["x"]])
    rownames(tmp2) <- paste0("Transpiration_", rownames(tmp2))


    res[[k1]] <- if (include_year) {
      rbind(Year = xtmp1[["Year"]], tmp1, tmp2)
    } else {
      rbind(tmp1, tmp2)
    }
  }

  res
}


metric_TranspirationSeasonality_v1 <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  calc_TranspirationSeasonality_v1(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    include_year = FALSE,
    peak_type = "volume",
    transp_timing_probs = c(0.1, 0.5, 0.9),
    days_smoothing = 15,
    days_window = 91,
    peak2_drop_factor = 2 / 3,
    peak2_increase_factor = 1 / 6,
    peak2_min_factor = 1 / 3
  )
}

metric_TranspirationSeasonality_v2 <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  calc_TranspirationSeasonality_v1(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    include_year = FALSE,
    peak_type = "value",
    transp_timing_probs = c(0.1, 0.5, 0.9),
    days_smoothing = 15,
    days_window = 91,
    peak2_drop_factor = 2 / 3,
    peak2_increase_factor = 1 / 6,
    peak2_min_factor = 1 / 3
  )
}


metric_TranspirationSeasonality_v3 <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  calc_TranspirationSeasonality_v1(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    include_year = FALSE,
    peak_type = "value",
    transp_timing_probs = c(0.1, 0.5, 0.9),
    days_smoothing = 31,
    days_window = 91,
    peak2_drop_factor = 2 / 3,
    peak2_increase_factor = 1 / 6,
    peak2_min_factor = 1 / 3
  )
}



metric_TranspirationSeasonality_v4 <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  calc_TranspirationSeasonality_v2(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    include_year = FALSE,
    transp_timing_probs = c(0.1, 0.5, 0.9),
    N_peaks_max_reported = 4L,
    peak_type = "value",
    days_smoothing = 15,
    days_window = 91,
    peaks_drop_factor = 2 / 3,
    peaks_increase_factor = 1 / 6,
    peaksize_min_factor = 1 / 3
  )
}

metric_TranspirationSeasonality_v5 <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  calc_TranspirationSeasonality_v2(
    path, name_sw2_run, id_scen_used,
    list_years_scen_used = list_years_scen_used,
    include_year = FALSE,
    transp_timing_probs = c(0.1, 0.5, 0.9),
    N_peaks_max_reported = 4L,
    peak_type = "volume",
    days_smoothing = 15,
    days_window = 91,
    peaks_drop_factor = 2 / 3,
    peaks_increase_factor = 1 / 6,
    peaksize_min_factor = 1 / 4
  )
}



get_SW2flux <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  sw2_out, sw2_var,
  transform = function(x) x,
  out_labels,
  include_year = FALSE,
  ...
) {
  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    tmp_mon <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Month",
      sw2_outs = sw2_out,
      sw2_vars = sw2_var,
      varnames_are_fixed = TRUE
    )

    tmp_yr <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Year",
      sw2_outs = sw2_out,
      sw2_vars = sw2_var,
      varnames_are_fixed = TRUE
    )

    ts_years <- tmp_yr[["time"]][, "Year"]

    res[[k1]] <- rbind(
      format_yearly_to_matrix(
        x = lapply(tmp_yr[["values"]], transform),
        years = ts_years,
        out_labels = out_labels
      ),
      format_monthly_to_matrix(
        x = lapply(tmp_mon[["values"]], transform),
        years = ts_years,
        out_labels = out_labels
      )
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

#--- Annual and monthly fluxes:
#' Evapotranspiration (ET) [mm]
metric_ET <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  get_SW2flux(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    sw2_out = "AET",
    sw2_var = "evapotr_cm",
    transform = function(x) 10 * x,
    out_labels = "ET_mm"
  )
}

#' Diffuse Recharge (DR) = Deep Drainage [mm]
metric_DR <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  get_SW2flux(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    sw2_out = "DEEPSWC",
    sw2_var = "lowLayerDrain_cm",
    transform = function(x) 10 * x,
    out_labels = "DeepDrainage_mm"
  )
}


#' Radiation (H) [MJ/m2]
metric_Radiation <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  get_SW2flux(
    path = path,
    name_sw2_run = name_sw2_run,
    id_scen_used = id_scen_used,
    list_years_scen_used = list_years_scen_used,
    sw2_out = "PET",
    sw2_var = c("H_oh_MJm-2", "H_gt_MJm-2"),
    out_labels = c("H_oh_MJm-2", "H_gt_MJm-2")
  )
}



metric_Climate_Annual <- function(
  path, name_sw2_run,
  id_scen_used, list_years_scen_used,
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
      sw2_outs = c("TEMP", "TEMP", "TEMP"),
      sw2_vars = c("max_C", "min_C", "avg_C"),
      varnames_are_fixed = TRUE
    )

    tmp_mon <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Month",
      sw2_outs = c(
        "PRECIP", "PRECIP", "PRECIP", "PRECIP",
        "TEMP", "TEMP", "TEMP",
        "SNOWPACK",
        "PET", "AET"
      ),
      sw2_vars = c(
        "ppt", "rain", "snow_fall", "snowmelt",
        "max_C", "avg_C", "min_C",
        "snowpackWaterEquivalent_cm",
        "pet_cm", "evapotr_cm"
      ),
      varnames_are_fixed = TRUE
    )

    tmp_yr <- extract_from_sw2(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      sw2_tp = "Year",
      sw2_outs = c(
        "PRECIP", "PRECIP", "PRECIP",
        "TEMP", "TEMP", "TEMP",
        "SNOWPACK",
        "PET", "AET"
      ),
      sw2_vars = c(
        "ppt", "rain", "snow_fall",
        "max_C", "avg_C", "min_C",
        "snowpackWaterEquivalent_cm",
        "pet_cm", "evapotr_cm"
      ),
      varnames_are_fixed = TRUE
    )

    # Helper variables
    ts_years <- tmp_yr[["time"]][, "Year"]

    PPTinJAS <- tapply(
      X = tmp_mon[["values"]][["ppt"]],
      INDEX = tmp_mon[["time"]][, "Year"],
      FUN = function(x) sum(x[7:9])
    )


    # Output
    res[[k1]] <- rbind(
      # Potential evapotranspiration
      format_yearly_to_matrix(
        x = 10 * tmp_yr[["values"]][["pet_cm"]],
        years = ts_years,
        out_labels = "PET_mm"
      ),
      format_monthly_to_matrix(
        x = 10 * tmp_mon[["values"]][["pet_cm"]],
        years = ts_years,
        out_labels = "PET_mm"
      ),

      # Climatic water deficit
      format_yearly_to_matrix(
        x = 10 * (
          tmp_yr[["values"]][["pet_cm"]] -
          tmp_yr[["values"]][["evapotr_cm"]]
        ),
        years = ts_years,
        out_labels = "CWD_mm"
      ),
      format_monthly_to_matrix(
        x = 10 * (
          tmp_mon[["values"]][["pet_cm"]] -
          tmp_mon[["values"]][["evapotr_cm"]]
        ),
        years = ts_years,
        out_labels = "CWD_mm"
      ),

      # Maximum temperature (mean across year)
      format_yearly_to_matrix(
        x = tmp_yr[["values"]][["max_C"]],
        years = ts_years,
        out_labels = "Tmax_mean_C"
      ),
      format_monthly_to_matrix(
        x = tmp_mon[["values"]][["max_C"]],
        years = ts_years,
        out_labels = "Tmax_mean_C"
      ),

      # Mean temperature
      format_yearly_to_matrix(
        x = tmp_yr[["values"]][["avg_C"]],
        years = ts_years,
        out_labels = "Tmean_mean_C"
      ),
      format_monthly_to_matrix(
        x = tmp_mon[["values"]][["avg_C"]],
        years = ts_years,
        out_labels = "Tmean_mean_C"
      ),

      # Minimum temperature (mean across year)
      format_yearly_to_matrix(
        x = tmp_yr[["values"]][["min_C"]],
        years = ts_years,
        out_labels = "Tmin_mean_C"
      ),
      format_monthly_to_matrix(
        x = tmp_mon[["values"]][["min_C"]],
        years = ts_years,
        out_labels = "Tmin_mean_C"
      ),

      # Temperature range (difference between tmax and tmin)
      Trange_diurnal_C = tapply(
        X =
          tmp_daily[["values"]][["max_C"]] - tmp_daily[["values"]][["min_C"]],
        INDEX = tmp_daily[["time"]][, "Year"],
        FUN = mean
      ),
      format_yearly_to_matrix(
        x = tapply(
          X =
            tmp_daily[["values"]][["max_C"]] -
            tmp_daily[["values"]][["min_C"]],
          INDEX = tmp_daily[["time"]][, "Year"],
          FUN = mean
        ),
        years = ts_years,
        out_labels = "Trange_diurnal_C"
      ),
      format_monthly_to_matrix(
        x = t(tapply(
          X =
            tmp_daily[["values"]][["max_C"]] -
            tmp_daily[["values"]][["min_C"]],
          INDEX = list(
            tmp_daily[["time"]][, "Year"],
            tmp_daily[["time"]][, "Month"]
          ),
          FUN = mean
        )),
        years = ts_years,
        out_labels = "Trange_diurnal_C"
      ),

      # SD of mean temperature
      format_yearly_to_matrix(
        x = tapply(
          X = tmp_daily[["values"]][["avg_C"]],
          INDEX = tmp_daily[["time"]][, "Year"],
          FUN = sd
        ),
        years = ts_years,
        out_labels = "Tmean_C_SD"
      ),
      format_monthly_to_matrix(
        x = t(tapply(
          X = tmp_daily[["values"]][["avg_C"]],
          INDEX = list(
            tmp_daily[["time"]][, "Year"],
            tmp_daily[["time"]][, "Month"]
          ),
          FUN = sd
        )),
        years = ts_years,
        out_labels = "Tmean_C_SD"
      ),

      # Temperature of hottest month
      format_yearly_to_matrix(
        x = tapply(
          X = tmp_mon[["values"]][["avg_C"]],
          INDEX = tmp_mon[["time"]][, "Year"],
          FUN = max
        ),
        years = ts_years,
        out_labels = "Tmean_hottestmonth_C"
      ),

      # Temperature of coldest month
      format_yearly_to_matrix(
        x = tapply(
          X = tmp_mon[["values"]][["avg_C"]],
          INDEX = tmp_mon[["time"]][, "Year"],
          FUN = min
        ),
        years = ts_years,
        out_labels = "Tmean_coldestmonth_C"
      ),

      # Precipitation
      format_yearly_to_matrix(
        x = 10 * tmp_yr[["values"]][["ppt"]],
        years = ts_years,
        out_labels = "PPT_mm"
      ),
      format_monthly_to_matrix(
        x = 10 * tmp_mon[["values"]][["ppt"]],
        years = ts_years,
        out_labels = "PPT_mm"
      ),

      format_yearly_to_matrix(
        x = 10 * tmp_yr[["values"]][["rain"]],
        years = ts_years,
        out_labels = "Rain_mm"
      ),
      format_monthly_to_matrix(
        x = 10 * tmp_mon[["values"]][["rain"]],
        years = ts_years,
        out_labels = "Rain_mm"
      ),

      format_yearly_to_matrix(
        x = 10 * tmp_yr[["values"]][["snow_fall"]],
        years = ts_years,
        out_labels = "Snowfall_mm"
      ),
      format_monthly_to_matrix(
        x = 10 * tmp_mon[["values"]][["snow_fall"]],
        years = ts_years,
        out_labels = "Snowfall_mm"
      ),

      format_yearly_to_matrix(
        x = 10 * tmp_yr[["values"]][["snowpackWaterEquivalent_cm"]],
        years = ts_years,
        out_labels = "Snowpack_SWE_mm"
      ),
      format_monthly_to_matrix(
        x = 10 * tmp_mon[["values"]][["snowpackWaterEquivalent_cm"]],
        years = ts_years,
        out_labels = "Snowpack_SWE_mm"
      ),

      # Monthly snowmelt
      format_monthly_to_matrix(
        x = 10 * tmp_mon[["values"]][["snowmelt"]],
        years = ts_years,
        out_labels = "Snowmelt_mm"
      ),

      # Ratio of rain to all precipitation
      format_yearly_to_matrix(
        x = tmp_yr[["values"]][["rain"]] / tmp_yr[["values"]][["ppt"]],
        years = ts_years,
        out_labels = "Rain_to_PPT"
      ),

      # Ratio of July, August, and September precipitation to annual
      format_yearly_to_matrix(
        x = PPTinJAS / tmp_yr[["values"]][["ppt"]],
        years = ts_years,
        out_labels = "PPTinJAS_to_PPT"
      ),

      # Sum of precipitation in months of July, August, and September
      format_yearly_to_matrix(
        x = 10 * PPTinJAS,
        years = ts_years,
        out_labels = "PPTinJAS_mm"
      ),

      # Precipitation of the wettest month
      format_yearly_to_matrix(
        x = 10 * tapply(
          X = tmp_mon[["values"]][["ppt"]],
          INDEX = tmp_mon[["time"]][, "Year"],
          FUN = max
        ),
        years = ts_years,
        out_labels = "PPT_wettestmonth_mm"
      ),


      # Precipitation of the driest month
      format_yearly_to_matrix(
        x = 10 * tapply(
          X = tmp_mon[["values"]][["ppt"]],
          INDEX = tmp_mon[["time"]][, "Year"],
          FUN = min
        ),
        years = ts_years,
        out_labels = "PPT_driestmonth_mm"
      )
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


#--- Soil moisture regimes / soil temperature regimes
metric_SMTRs <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "across_years",
  ...
) {
  stopifnot(requireNamespace("rSW2funs", quietly = TRUE))
  stopifnot(check_metric_arguments(out = match.arg(out)))

  res <- list()

  sim_input <- new.env(parent = emptyenv())
  load(
    file = file.path(path, name_sw2_run, "sw_input.RData"),
    envir = sim_input
  )

  for (k1 in seq_along(id_scen_used)) {
    sim_data <- new.env(parent = emptyenv())
    load(
      file = file.path(
        path,
        name_sw2_run,
        paste0("sw_output_sc", id_scen_used[k1], ".RData")
      ),
      envir = sim_data
    )

    tmp <- lapply(
      list_years_scen_used[[k1]],
      function(yrs) {
        # nolint start
        rSOILWAT2::swYears_EndYear(sim_input[["swRunScenariosData"]][[k1]]) <- 9999
        rSOILWAT2::swYears_StartYear(sim_input[["swRunScenariosData"]][[k1]]) <-
          yrs[1]
        rSOILWAT2::swYears_EndYear(sim_input[["swRunScenariosData"]][[k1]]) <-
          yrs[length(yrs)]
        # nolint end

        tmp <- rSW2funs::calc_SMTRs(
          sim_in = sim_input[["swRunScenariosData"]][[k1]],
          sim_out = sim_data[["runDataSC"]]
        )

        cbind(tmp[["STR"]], tmp[["SMR"]])
      }
    )

    res[[k1]] <- t(do.call(rbind, tmp))
  }

  res
}
