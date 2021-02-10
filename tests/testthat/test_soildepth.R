

test_that("Soil layer widths", {
  soil_depths_cm <- list(
    NULL,
    5,
    10,
    c(5, 10, 30),
    c(5, 10, 30, 50, 100, 200),
    c(30, 50, 100, 200)
  )


  #--- Use complete depth range
  # Including: no soils --> `numeric(0)`
  for (k in seq_along(soil_depths_cm)) {
    expect_length(
      calc_soillayer_weights(
        soil_depths_cm = soil_depths_cm[[k]],
        used_depth_range_cm = NULL
      ),
      length(soil_depths_cm[[k]])
    )
  }

  expect_equal(
    calc_soillayer_weights(
      soil_depths_cm = c(10, 20, 30),
      used_depth_range_cm = c(5, 30)
    ),
    c(10, 10, 10)
  )

  expect_equal(
    calc_soillayer_weights(
      soil_depths_cm = c(10, 20, 30),
      used_depth_range_cm = c(5, 100)
    ),
    c(10, 10, 10)
  )

  expect_equal(
    calc_soillayer_weights(
      soil_depths_cm = c(5, 20, 30),
      used_depth_range_cm = c(5, 100)
    ),
    c(NA, 15, 10)
  )

  expect_equal(
    calc_soillayer_weights(
      soil_depths_cm = c(5, 20, 30, 40, 50),
      used_depth_range_cm = c(5, 30)
    ),
    c(NA, 15, 10, NA, NA)
  )
})


test_that("Soil layer availability", {

  expect_warning(
    check_soillayer_availability(
      soil_depths_cm = c(10, 20, 30),
      used_depth_range_cm = c(5, 30),
      strict = c(TRUE, FALSE),
      type = "warn"
    )
  )

  expect_warning(
    check_soillayer_availability(
      soil_depths_cm = c(10, 20, 30),
      used_depth_range_cm = c(5, 30),
      strict = TRUE,
      type = "warn"
    )
  )

  expect_error(
    check_soillayer_availability(
      soil_depths_cm = c(10, 20, 30),
      used_depth_range_cm = c(5, 30),
      strict = TRUE,
      type = "error"
    )
  )


  expect_silent(
    check_soillayer_availability(
      soil_depths_cm = c(10, 20, 30),
      used_depth_range_cm = c(5, 30),
      strict = c(FALSE, TRUE),
      type = "warn"
    )
  )

  expect_silent(
    check_soillayer_availability(
      soil_depths_cm = c(10, 20, 30),
      used_depth_range_cm = c(5, 30),
      strict = FALSE,
      type = "warn"
    )
  )
})
