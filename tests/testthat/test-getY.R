context("getY")

test_that("Dimensionality of Y is correct", {
  d <- load_mr("ex.xlsx")
  Y <- getY(d, M = 100)
  expect_equal(ncol(Y), 8)
  expect_equal(nrow(Y), 136)
})

