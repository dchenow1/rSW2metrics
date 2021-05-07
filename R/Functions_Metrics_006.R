
#--- AI = aridity index [mm/mm] = PPT/PET
metric_AI <- function(
  path, name_sw2_run,
  id_scen_used, list_years_scen_used,
  out = "ts_years",
  ...
) {
  stopifnot(check_metric_arguments(out = match.arg(out)))

  res <- list()

  for (k1 in seq_along(id_scen_used)) {
    # Annual PET and PPT
    sim_data <- collect_sw2_sim_data(
      path = path,
      name_sw2_run = name_sw2_run,
      id_scen = id_scen_used[k1],
      years = list_years_scen_used[[k1]],
      output_sets = list(
        yr = list(
          sw2_tp = "Year",
          sw2_outs = c("PET", "PRECIP"),
          sw2_vars = c(pet = "pet_cm", ppt = "ppt"),
          varnames_are_fixed = TRUE
        )
      )
    )

    res[[k1]] <- matrix(
      data =
        sim_data[["yr"]][["values"]][["ppt"]] /
        sim_data[["yr"]][["values"]][["pet"]],
      nrow = 1,
      dimnames = list("AI", NULL)
    )
  }

  res
}
