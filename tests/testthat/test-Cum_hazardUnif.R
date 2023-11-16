test_that("CDF by Cum_hazardUnif is equal to the area of rectangle?", {

  num1 <- Cum_hazardUnif(3,6,5)
  F_t = 1 - exp(-num1) # Cumulative Density Function by Cumulative Hazard Fn

  num2 <- (1/(6-3))*(5-3)  # area of rectangle

  expect_equal(F_t, num2)
})
