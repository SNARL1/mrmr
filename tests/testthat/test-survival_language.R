
test_that("plot_model errors when what='survival' but no translocations", {
  mock_model_obj <- list(data=list(translocations=NA))
  expect_error(survival_table(mock_model_obj),
               regexp = "No translocation data are present")
})