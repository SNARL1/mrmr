
test_that("plot_model errors when what='survival' but no translocations", {
  mock_model_obj <- list(data=list(translocations=NA))
  expect_error(plot_model(mock_model_obj, what='survival'),
               regexp = "No translocation data are present")
})