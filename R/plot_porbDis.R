#' Plot a scatter plot of the probability distribution against the primary assay feature
#'
#' @param ppi_prediction_result
#'
#' @return
#' @export
#'
#' @examples
plot.probDis <- function(ppi_prediction_result) {
  ppi_prediction_result$predDf %>%
    ggplot(aes(x = mean_cBRET, y = predMat, fill = complex)) +
    geom_point(shape = 21, alpha = 0.7, size = 1.5) +
    theme_pubr() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#6CA6C1", "#D930C5")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
    labs(x = "cBRET ratio", y = "probability") +
    theme_pubr() +
    ggplot2::theme(text = element_text(family = "Avenir"),
                   plot.title = element_text(size = 12),
                   plot.subtitle = element_text(size = 8),
                   axis.title = element_text(family = "Avenir Medium"),
                   legend.position = "right")
}
