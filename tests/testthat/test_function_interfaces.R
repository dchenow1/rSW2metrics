
required_metrics_arguments <- c(
  "path",
  "name_sw2_run",
  "id_scen_used",
  "list_years_scen_used",
  "out",
  "..."
)

test_that("Check parameters of outside-facing metrics", {
   funs_check <- list_all_metrics()

  for (fun in funs_check) {
    ff <- formals(fun)

    tmp <- !(required_metrics_arguments %in% names(ff))
    tmp_args <- required_metrics_arguments[tmp]
    expect_equal(
      tmp_args, character(0),
      label = paste(shQuote(fun), "is missing required arguments:")
    )

    expect_true(
      all(eval(ff[["out"]]) %in% c("across_years", "ts_years", "raw"))
    )
  }
})


required_inputcollectors_arguments <- c(
  "path",
  "name_sw2_run",
  "..."
)

test_that("Check parameters of outside-facing metrics", {
  funs_check <- list_all_input_collectors()

  for (fun in funs_check) {
    ff <- formals(fun)

    tmp <- !(required_inputcollectors_arguments %in% names(ff))
    tmp_args <- required_inputcollectors_arguments[tmp]
    expect_equal(
      tmp_args, character(0),
      label = paste(shQuote(fun), "is missing required arguments:")
    )
  }
})
