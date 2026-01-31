test_that("flows_summary returns expected components", {
  F1 <- Matrix::Matrix(0, 3, 3, sparse = TRUE); F1[1, 2] <- 1
  F2 <- Matrix::Matrix(0, 3, 3, sparse = TRUE); F2[2, 3] <- 2
  out <- flows_summary(list(F1, F2))
  expect_true(all(c("Psum", "inflow", "outflow", "hub_score", "mean_distance") %in% names(out)))
  expect_equal(length(out$inflow), 3)
})
