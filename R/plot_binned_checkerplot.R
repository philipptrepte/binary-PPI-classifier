#' Binned checker plot to display results from all-by-all screens
#'
#' @param ppi_prediction_result: result object from the ppi.prediction() function.
#' @param fraction_na: fraction of interactions allowed to be *NA* to still draw the checker plot
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data("example_ppi_prediction")
#' binnedhm.plot(example_ppi_prediction)
binnedcp.plot <- function(ppi_prediction_result, fraction_na = 0.05) {
  palette <- colorRampPalette(c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"))(100)

  na_values <- ppi_prediction_result$predDf %>%
    dplyr::mutate(pos = ifelse(predMat >= 0.95, 3,
                               ifelse(predMat >= 0.75, 2, 0))) %>%
    group_by(Donor_protein, Acceptor_protein, interaction) %>%
    dplyr::summarise(pos = max(pos, na.rm = TRUE)) %>%
    dplyr::select(Donor_protein, Acceptor_protein, pos) %>%
    pivot_wider(names_from = Acceptor_protein, values_from = pos)
  fraction_na_values <- sum(is.na(na_values)) / (dim(na_values)[1] * dim(na_values)[2])
  assertthat::assert_that(fraction_na_values <= fraction_na, msg = paste0("The fraction of NA values is ", round(fraction_na_values, 3), " and to high (>", fraction_na,") to plot a meaningful checker plot."))

  ppi_prediction_result$predDf %>%
    dplyr::mutate(pos = ifelse(predMat >= 0.95, 3,
                               ifelse(predMat >= 0.75, 2, 0))) %>%
    group_by(Donor_protein, Acceptor_protein, interaction) %>%
    dplyr::summarise(pos = max(pos, na.rm = TRUE)) %>%
    dplyr::select(Donor_protein, Acceptor_protein, pos) %>%
    ggplot(mapping = aes(x = Donor_protein, y = Acceptor_protein, fill = factor(pos, levels = c(3,2,0)))) +
    geom_tile() +
    scale_fill_manual(name = "probability bins", labels = c(expression("\u2265" * "95%"), expression("\u2265" * "75%"), ""), values = c(palette[95], palette[75], "white")) +
    theme(strip.text = element_text(colour = "white"),
          strip.background = element_rect(colour = "black", fill = "black"),
          plot.background = element_rect(fill = "white", colour = "white"),
          axis.line = element_line(colour = "white"),
          panel.background = element_rect(fill = "white", color = "white"),
          text = element_text(family = "Avenir"),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8),
          axis.text = element_text(size = 12),
          axis.title = element_text(family = "Avenir Medium"),
          legend.position = "right") +
    rotate_x_text(angle = 90, hjust = 1, vjust = 0.5) +
    ggplot2::labs(x = "Protein X", y = "Protein Y", title = paste0("Checker plot for the ", ppi_prediction_result$ensembleSize, " ", ppi_prediction_result$model.type, " models"),
                  subtitle = paste0("sampling: ", base::ifelse(ppi_prediction_result$sampling == "unweighted", "unweighted", ppi_prediction_result$sampling),
                                    if(ppi_prediction_result$model.type == "svm") paste0(" | kernel type: ", ppi_prediction_result$kernelType),
                                    " | cutoff: ", ppi_prediction_result$cutoff,
                                    " | iterations: ", ppi_prediction_result$iter))
}
