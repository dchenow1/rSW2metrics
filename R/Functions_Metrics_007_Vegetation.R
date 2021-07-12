

#--- Fractional land cover (available with rSOILWAT2 v3.1.2/SOILWAT v5.2.0)
metric_land_cover_v2 <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    sim_data <- collect_sw2_sim_data(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      output_sets = list(
        cover = list(
          sw2_tp = "Year",
          sw2_outs = "BIOMASS",
          sw2_vars = c(
            "fCover_BareGround",
            "fCover_tree", "fCover_shrub", "fCover_forbs", "fCover_grass"
          ),
          varnames_are_fixed = TRUE
        )
      )
    )

    res[[k1]] <- t(do.call(cbind, sim_data[["cover"]][["values"]]))
  }

  res
}


#--- Monthly and annual vegetation biomass
# available with rSOILWAT2 v3.1.2/SOILWAT v5.2.0
metric_veg_biomass_v2 <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  tmp_var <- c("total", "tree", "shrub", "forbs", "grass")
  tmp_veg <- c(
    paste0("Biomass_", tmp_var),
    "Biomass_litter",
    paste0("Biolive_", tmp_var)
  )

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    sim_data <- collect_sw2_sim_data(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      output_sets = list(
        yr = list(
          sw2_tp = "Year",
          sw2_outs = "BIOMASS",
          sw2_vars = tmp_veg,
          varnames_are_fixed = TRUE
        ),
        mon = list(
          sw2_tp = "Month",
          sw2_outs = "BIOMASS",
          sw2_vars = tmp_veg,
          varnames_are_fixed = TRUE
        )
      )
    )

    res[[k1]] <- rbind(
      format_yearly_to_matrix(
        x = sim_data[["yr"]][["values"]],
        years = sim_data[["yr"]][["time"]][, "Year"],
        out_labels = tmp_veg
      ),
      format_monthly_to_matrix(
        x = sim_data[["mon"]][["values"]],
        years = sim_data[["yr"]][["time"]][, "Year"],
        out_labels = tmp_veg
      )
    )
  }

  res
}


#--- Monthly and annual vegetation biomass
# available with rSOILWAT2 v1.7.0/SOILWAT v3.5.0
# and before rSOILWAT2 v3.1.2/SOILWAT v5.2.0)
metric_veg_biomass_v1 <- function(
  path, name_sw2_run, id_scen_used, list_years_scen_used,
  out = "ts_years",
  include_year = FALSE,
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  tmp_var <- c("total", "tree", "shrub", "forbs", "grass")
  tmp_veg_metrics <- c(
    paste0("Biomass_", tmp_var),
    "Biomass_litter",
    paste0("Biolive_", tmp_var)
  )

  tmp_var2 <- c("Total", "Tree", "Shrub", "Forb", "Grass")
  tmp_veg_tbm <- paste0(tmp_var2, "Biomass")
  tmp_veg_lbm <- paste0(tmp_var2, "Biolive")
  tmp_veg_sw2old <- c(tmp_veg_tbm, tmp_veg_lbm)

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    sim_data <- collect_sw2_sim_data(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      output_sets = list(
        yr = list(
          sw2_tp = "Year",
          sw2_outs = "CO2EFFECTS",
          sw2_vars = tmp_veg_sw2old,
          varnames_are_fixed = TRUE
        ),
        mon = list(
          sw2_tp = "Month",
          sw2_outs = "CO2EFFECTS",
          sw2_vars = tmp_veg_sw2old,
          varnames_are_fixed = TRUE
        )
      ),
      fail = FALSE
    )


    res[[k1]] <- rbind(
      format_yearly_to_matrix(
        x = c(
          sim_data[["yr"]][["values"]][tmp_veg_tbm],
          Biomass_litter = list(rep(NA, nrow(sim_data[["yr"]][["time"]]))),
          sim_data[["yr"]][["values"]][tmp_veg_lbm]
        ),
        years = sim_data[["yr"]][["time"]][, "Year"],
        out_labels = tmp_veg_metrics
      ),
      format_monthly_to_matrix(
        x = c(
          sim_data[["mon"]][["values"]][tmp_veg_tbm],
          Biomass_litter = list(rep(NA, nrow(sim_data[["mon"]][["time"]]))),
          sim_data[["mon"]][["values"]][tmp_veg_lbm]
        ),
        years = sim_data[["yr"]][["time"]][, "Year"],
        out_labels = tmp_veg_metrics
      )
    )
  }

  res
}
