test_that("Is this double data?", {
  expect_type(hazardWei(2,4,3), "double")
})


test_that("Is this greater than 0?", {
  expect_gt(hazardWei(2,4,3), 0)
})
