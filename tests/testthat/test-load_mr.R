context("load_mr")

test_that("valid files are validated", {
  d <- load_mr("ex.xlsx")
  expect_equal(ncol(d$surveys), 16)
  expect_equal(ncol(d$translocations), 3)
})
