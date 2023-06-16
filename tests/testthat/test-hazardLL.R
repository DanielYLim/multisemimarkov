test_that("Is this double data?", {
  expect_type(hazardLL(2,4,3), "double")
})


test_that("Is this greater than 0?", {
  expect_gt(hazardLL(2,4,3), 0)
})
