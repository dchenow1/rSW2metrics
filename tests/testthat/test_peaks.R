

test_that("Peaks", {
  ## No or less than two plateaus -> central candidate from among all values
  expect_equal(central_candidate(1), 1)
  expect_equal(central_candidate(NULL), NULL)
  expect_equal(central_candidate(1:2), 1)
  expect_equal(central_candidate(1:3), 2)
  expect_equal(central_candidate(1:4), 2)
  expect_equal(central_candidate(c(1, 3, 5)), 3)
  expect_equal(central_candidate(c(1, 3, 5, 7, 9, 17:19)), 7)
  expect_equal(central_candidate(c(1, 20:29)), 24)

  ## More than one plateau -> central candidate from longest plateau
  expect_equal(central_candidate(c(1, 20:30, 40:60)), 50)
  expect_equal(central_candidate(c(1, 20:30, 40:61)), 50)
  expect_equal(central_candidate(40:60), 50)
  expect_equal(central_candidate(c(1, 20:30, 40:60, 70:80, 100:110)), 50)
  expect_equal(central_candidate(c(1, 20:30, 40:60, 70:90)), 80)
})
