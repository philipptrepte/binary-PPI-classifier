data("example_ppi_prediction")

testthat::test_that("Test plot_probDis helper function", {
  test_result <- probDis.plot(example_ppi_prediction)
  testthat::expect_true(ggplot2::is.ggplot(test_result))
})
