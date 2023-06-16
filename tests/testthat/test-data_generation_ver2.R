test_that("Cause 1->2 and Cause 1->3 were properly generated?", {

  temp <- data_generation_ver2(1000)
  temp <- temp$status12 + temp$status13
  temp2 <- sum(temp>1)



  expect_equal(temp2, 0)
})




test_that("Cause 1->2 and Cause 2->3 were properly generated?", {

  temp <- data_generation_ver2(1000)
  temp < - temp$status23 + temp$status12
  temp2 <- sum(temp==-1)



  expect_equal(temp2, 0)
})
