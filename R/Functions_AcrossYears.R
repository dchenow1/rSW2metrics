
add_aggs_across_years <- function(
  x,
  fun,
  id_scen_used,
  list_years_scen_used
) {
  fun_across_years <- match.fun(fun)
  stopifnot(
    inherits(list_years_scen_used, "list"),
    sapply(list_years_scen_used, is.numeric)
  )

  # Prepare container
  tmp1 <- names(fun_across_years(1))
  tmp2 <- sapply(
    seq_along(id_scen_used),
    function(k1) {
      ntmp <- length(list_years_scen_used[[k1]])
      paste0(
        "sc", id_scen_used[k1], "_",
        list_years_scen_used[[k1]][1], "-",
        list_years_scen_used[[k1]][ntmp]
      )
    }
  )

  rn_aggs <- paste0(tmp1, "_", rep(tmp2, each = length(tmp1)))

  x_aggs <- array(
    data = NA,
    dim = c(length(rn_aggs), nrow(x)),
    dimnames = list(rn_aggs, NULL)
  )

  cn_vals <- colnames(x)

  tmp <- unique(nchar(as.character(range(unlist(list_years_scen_used)))))
  tag_cn_yrs <- paste0("_[[:digit:]{", paste0(tmp, collapse = ","), "}]")

  # Calculate aggregations for each scenario
  for (k1 in seq_along(id_scen_used)) {
    tmp <- paste0("sc", id_scen_used[k1], tag_cn_yrs)
    tmp1 <- grep(tmp, rn_aggs)
    tmp2 <- grep(tmp, cn_vals)

    x_aggs[tmp1, ] <- apply(
      X = x[, tmp2],
      MARGIN = 1,
      FUN = fun_across_years,
      na.rm = TRUE
    )
  }

  # Combine with data
  cn_header <- c("site", "group")
  cbind(
    x[, cn_header],
    t(x_aggs),
    x[, - sapply(cn_header, grep, x = cn_vals)]
  )
}



#' Coefficient of variation
cv <- function(x, ...) {
  tmp <- mean(x)
  if (isTRUE(abs(tmp) > sqrt(.Machine$double.eps))) {
    sd(x) / tmp
  } else {
    NA
  }
}

#' Sen Slope of Mann-Kendall Trend Test
sen_slope <- function(x, ...) {
  stopifnot(requireNamespace("modifiedmk", quietly = TRUE))
  if (sum(is.finite(x)) > 3) {
    unname(modifiedmk::mkttest(x)["Sen's slope"])
  } else {
    NA
  }
}


#' Calculate aggregate statistics of values across years
#'
#' @param x A numeric vector.
#' @param na.rm A logical value
#'
#' @name aggs_across_years
NULL

#' @rdname aggs_across_years
#' @export
mean_cv <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[is.finite(x)]

  c(mean = mean(x), cv = cv(x))
}

#' @rdname aggs_across_years
#' @export
mean_cv_trend <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[is.finite(x)]

  c(mean = mean(x), cv = cv(x), trend = sen_slope(x))
}
