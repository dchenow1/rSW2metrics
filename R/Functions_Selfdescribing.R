
#' List all metric functions
#' @export
list_all_metrics <- function() {
  tmp <- ls(getNamespace("rSW2metrics"))
  tmp <- tmp[grepl("^metric_", tmp)]
  tmp[sapply(tmp, function(x) is.function(get0(x)))]
}


has_fun_soils_as_arg <- function(fun) {
  isTRUE("soils" %in% names(formals(fun)))
}

has_fun_ts_as_output <- function(fun) {
  tmp <- eval(formals(fun)[["out"]])
  isTRUE(tmp == "ts_years" || tmp == c("ts_years", "raw"))
}

is_fun_collecting_inputs <- function(fun) {
  isTRUE(grepl("collect_input_", fun, fixed = TRUE))
}




check_metric_arguments <- function(out, req_soil_vars) {
  fun_args <- as.list(match.call(
    definition = sys.function(sys.parent()),
    call = sys.call(sys.parent()),
    envir = parent.frame(2L)
  ))

  # Evaluate arguments, but not function name (first element)
  fun_name <- if (is.function(fun_args[[1]])) {
    "anonymous"
  } else {
    as.character(fun_args[[1]])
  }

  fun_args <- lapply(fun_args[-1], eval)

  if (missing(out) || is.null(out)) {
    stop("'out' is a required argument to function ", shQuote(fun_name))
  }

  if ("out" %in% names(fun_args) && !(fun_args[["out"]] %in% c(out, "raw"))) {
    stop(
      "Inconsistency in 'out': ",
      shQuote(out), " versus ", shQuote(fun_args[["out"]])
    )
  }

  if (out %in% "ts_years") {
    if (
      !(
        isTRUE(inherits(fun_args[["list_years_scen_used"]], "list")) &&
          all(sapply(fun_args[["list_years_scen_used"]], is.numeric))
      )
    ) {
      stop(
        "'list_years_scen_used' ",
        "should be a list with integer vectors ",
        "for function ", shQuote(fun_name)
      )
    }

  } else if (out %in% "across_years") {
    if (
      !(
        isTRUE(inherits(fun_args[["list_years_scen_used"]], "list")) &&
          all(
            sapply(fun_args[["list_years_scen_used"]], inherits, what = "list")
          ) &&
          all(
            unlist(lapply(
              fun_args[["list_years_scen_used"]],
              function(x) sapply(x, is.numeric)
            ))
          )
      )
    ) {
      stop(
        "'list_years_scen_used' ",
        "should be a list of lists with integer vectors ",
        "for function ", shQuote(fun_name)
      )
    }
  } else if (out %in% "raw") {
  }

  if (!missing(req_soil_vars)) {
    if ("soils" %in% names(fun_args)) {
      if (any(tmp <- !(req_soil_vars %in% names(list_soil_variables())))) {
        stop(
          "Requested soil variable(s) ",
          paste(shQuote(req_soil_vars[tmp]), collapse = ", "),
          " are not implemented, see `list_soil_variables()`."
        )
      }

      has_vars <-
        req_soil_vars %in% names(fun_args[["soils"]]) &
        sapply(
          req_soil_vars,
          function(x) !is.null(fun_args[["soils"]][[x]])
        )

      if (any(!has_vars)) {
        stop(
          "Requested soil variable(s) ",
          paste(shQuote(req_soil_vars[!has_vars]), collapse = ", "),
          " are missing from 'soils' object."
        )
      }

    } else {
      stop(
        "Soil variables are required for function ", shQuote(fun_name),
        " but there is no 'soils' object."
      )
    }
  }

  invisible(TRUE)
}



#' List all input collecting functions
#' @export
list_all_input_collectors <- function() {
  tmp <- ls(getNamespace("rSW2metrics"))
  tmp <- tmp[grepl("^collect_input_", tmp)]
  tmp[sapply(tmp, function(x) is.function(get0(x)))]
}


#' List of implemented soil properties
#'
#' @export
list_soil_variables <- function() {
  c(
    depth_cm = "depth",
    sand_frac = "sand",
    clay_frac = "clay",
    gravel_content = "gravel"
  )
}
