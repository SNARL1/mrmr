context("getZ")

test_that("Dimensionality of Z is correct", {
  d <- load_mr("ex.xlsx")
  Y <- getY(d, M = 100)
  Z <- getZ(Y)
  expect_equal(ncol(Z), 9)
  expect_equal(nrow(Z), 136)
})

