test_that("external_flux works on tiny example", {
  Ft <- Matrix::Matrix(0, 3, 3, sparse = TRUE)
  Ft[1, 2] <- 4
  Ft[2, 3] <- 2
  a <- c(6, 3, 0)
  b <- c(0, 4, 5)
  ef <- external_flux(Ft, a, b)
  expect_true(ef$ext_out >= 0 && ef$ext_in >= 0)
})
