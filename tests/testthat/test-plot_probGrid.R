library(testthat)
data("example_ppi_prediction")

testthat::test_that("Test plot_probGrid helper function", {
  test_result <- probGrid.plot(example_ppi_prediction, ylim = c(-0.1, 0.8))
  testthat::expect_true(ggplot2::is.ggplot(test_result))
})
