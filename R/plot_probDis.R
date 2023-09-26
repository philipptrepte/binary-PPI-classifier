#' Plot a scatter plot of the probability distribution against the primary assay feature
#'
#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#'
#' @param ppi_prediction_result: result object from the ppi.prediction() function.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
probDis.plot <- function(ppi_prediction_result, set="test") {
  assertthat::assert_that(set %in% c("test", "train"), msg = "'set=' must be 'test' or 'train'.")
  if(set == "test") {
    plot <- ppi_prediction_result$predDf %>%
      ggplot2::ggplot(aes(x = mean_cBRET, y = predMat, fill = complex)) +
      ggplot2::geom_point(shape = 21, alpha = 0.7, size = 1.5) +
      ggpubr::theme_pubr() +
      ggplot2::theme(legend.position = "none") +
      viridis::scale_fill_viridis(discrete = TRUE, option = "D") +
      ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
      ggplot2::labs(x = "cBRET ratio", y = "probability") +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12),
                     plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     legend.position = "right")
  } else if(set == "train") {
    plot <- ppi_prediction_result$predTrainDf %>%
      ggplot2::ggplot(aes(x = mean_cBRET, y = predTrainMat, fill = reference)) +
      ggplot2::geom_point(shape = 21, alpha = 0.7, size = 1.5) +
      ggpubr::theme_pubr() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_fill_manual(values = c("#6CA6C1", "#D930C5")) +
      ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
      ggplot2::labs(x = "cBRET ratio", y = "probability") +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12),
                     plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     legend.position = "right")
  }
  return(plot)
}
