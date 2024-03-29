#' Plot a bar diagram of the recovery rates at 50%, 75% and 95% probability
#'
#' @import dplyr
#' @import ggplot2
#' @import stringr
#' @import ggpubr
#'
#' @param ppi_prediction_result: result object from the ppi.prediction() function.
#' @param set: PPI set to generate the plot for: "test" or "train"
#' @param y: specify to plot number ('count') or fraction ('frac') of positive interactions
#' @param groupby: for the training set group by reference or complex
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' data("example_ppi_prediction")
#' recovery.plot(example_ppi_prediction)
recovery.plot <- function(ppi_prediction_result, set="train", y='count', groupby = 'reference') {
  assertthat::assert_that(set %in% c("test", "train"), msg = "'set=' must be 'test' or 'train'.")
  assertthat::assert_that(y %in% c("count", "frac"), msg = "'y=' must be 'count' or 'frac'.")
  assertthat::assert_that(groupby %in% c("reference", "complex"), msg = "'groupby=' must be 'reference' or 'complex'.")

  if(set == "test") {
    Df <- ppi_prediction_result$predDf %>%
      dplyr::mutate(reference = complex) %>%
      dplyr::group_by(interaction, reference) %>%
      dplyr::summarise(score = max(predMat, na.rm = TRUE)) %>%
      dplyr::mutate(pos = base::ifelse(score > 0.5, 1, 0),
                    probability = 0.5) %>%
      dplyr::group_by(reference, probability) %>%
      dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                       total = n(),
                       recovery = pos / total * 100,
                       sep = sqrt(recovery*(100-recovery)/total)) %>%
      rbind(
        ppi_prediction_result$predDf %>%
          dplyr::mutate(reference = complex) %>%
          dplyr::group_by(interaction, reference) %>%
          dplyr::summarise(score = max(predMat, na.rm = TRUE)) %>%
          dplyr::mutate(pos = base::ifelse(score > 0.75, 1, 0),
                        probability = 0.75) %>%
          dplyr::group_by(reference, probability) %>%
          dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                           total = n(),
                           recovery = pos / total * 100,
                           sep = sqrt(recovery*(100-recovery)/total))
      ) %>%
      rbind(
        ppi_prediction_result$predDf %>%
          dplyr::mutate(reference = complex) %>%
          dplyr::group_by(interaction, reference) %>%
          dplyr::summarise(score = max(predMat, na.rm = TRUE)) %>%
          dplyr::mutate(pos = base::ifelse(score > 0.95, 1, 0),
                        probability = 0.95) %>%
          dplyr::group_by(reference, probability) %>%
          dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                           total = n(),
                           recovery = pos / total * 100,
                           sep = sqrt(recovery*(100-recovery)/total)))
  } else if(set == "train") {
    if(groupby == "reference") {
      Df <- ppi_prediction_result$predTrainDf %>%
        dplyr::mutate(reference = base::ifelse(stringr::str_detect(complex, paste(ppi_prediction_result$negative.reference, collapse = "|")), "RRS", "PRS")) %>%
        dplyr::group_by(interaction, reference) %>%
        dplyr::summarise(score = max(predTrainMat, na.rm = TRUE)) %>%
        dplyr::mutate(pos = base::ifelse(score > 0.5, 1, 0),
                      probability = 0.5) %>%
        dplyr::group_by(reference, probability) %>%
        dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                         total = n(),
                         recovery = pos / total * 100,
                         sep = sqrt(recovery*(100-recovery)/total)) %>%
        rbind(
          ppi_prediction_result$predTrainDf %>%
            dplyr::mutate(reference = base::ifelse(stringr::str_detect(complex, paste(ppi_prediction_result$negative.reference, collapse = "|")), "RRS", "PRS")) %>%
            dplyr::group_by(interaction, reference) %>%
            dplyr::summarise(score = max(predTrainMat, na.rm = TRUE)) %>%
            dplyr::mutate(pos = base::ifelse(score > 0.75, 1, 0),
                          probability = 0.75) %>%
            dplyr::group_by(reference, probability) %>%
            dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                             total = n(),
                             recovery = pos / total * 100,
                             sep = sqrt(recovery*(100-recovery)/total))
        ) %>%
        rbind(
          ppi_prediction_result$predTrainDf %>%
            dplyr::mutate(reference = base::ifelse(stringr::str_detect(complex, paste(ppi_prediction_result$negative.reference, collapse = "|")), "RRS", "PRS")) %>%
            dplyr::group_by(interaction, reference) %>%
            dplyr::summarise(score = max(predTrainMat, na.rm = TRUE)) %>%
            dplyr::mutate(pos = base::ifelse(score > 0.95, 1, 0),
                          probability = 0.95) %>%
            dplyr::group_by(reference, probability) %>%
            dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                             total = n(),
                             recovery = pos / total * 100,
                             sep = sqrt(recovery*(100-recovery)/total)))
    } else if(groupby == "complex") {
      Df <- ppi_prediction_result$predTrainDf %>%
        dplyr::mutate(reference = complex) %>%
        dplyr::group_by(interaction, reference) %>%
        dplyr::summarise(score = max(predTrainMat, na.rm = TRUE)) %>%
        dplyr::mutate(pos = base::ifelse(score > 0.5, 1, 0),
                      probability = 0.5) %>%
        dplyr::group_by(reference, probability) %>%
        dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                         total = n(),
                         recovery = pos / total * 100,
                         sep = sqrt(recovery*(100-recovery)/total)) %>%
        rbind(
          ppi_prediction_result$predTrainDf %>%
            dplyr::mutate(reference = complex) %>%
            dplyr::group_by(interaction, reference) %>%
            dplyr::summarise(score = max(predTrainMat, na.rm = TRUE)) %>%
            dplyr::mutate(pos = base::ifelse(score > 0.75, 1, 0),
                          probability = 0.75) %>%
            dplyr::group_by(reference, probability) %>%
            dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                             total = n(),
                             recovery = pos / total * 100,
                             sep = sqrt(recovery*(100-recovery)/total))
        ) %>%
        rbind(
          ppi_prediction_result$predTrainDf %>%
            dplyr::mutate(reference = complex) %>%
            dplyr::group_by(interaction, reference) %>%
            dplyr::summarise(score = max(predTrainMat, na.rm = TRUE)) %>%
            dplyr::mutate(pos = base::ifelse(score > 0.95, 1, 0),
                          probability = 0.95) %>%
            dplyr::group_by(reference, probability) %>%
            dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                             total = n(),
                             recovery = pos / total * 100,
                             sep = sqrt(recovery*(100-recovery)/total)))
    }
  }


  if(y == 'count') {
    plot <- Df %>%
      ggplot2::ggplot(aes(x = reference, y = pos, fill = reference)) +
      ggplot2::geom_bar(stat = "identity", col = "black") +
      ggplot2::geom_text(mapping = aes(label = round(pos, 0)), vjust = -3, color = "black", size = 2) +
      ggplot2::ylim(0, NA) +
      ggplot2::scale_fill_manual(values = c("#6CA6C1", "#D930C5", "#FFE66D", "#8ED2C6", "#01A79D", "#FFB400", "#F6511D", "#00a6ed", "#0a2463")) +
      ggplot2::facet_wrap(~ probability, ncol = 3) +
      ggplot2::labs(x = "", y = "Value", title = paste0("Mean + IQR Learning Curves for the ", ppi_prediction_result$ensembleSize, " ", ppi_prediction_result$model.type, " models"),
                    subtitle = paste0("sampling: ", base::ifelse(ppi_prediction_result$sampling == "unweighted", "unweighted", ppi_prediction_result$sampling),
                                      if(ppi_prediction_result$model.type == "svm") paste0(" | kernel type: ", ppi_prediction_result$kernelType),
                                      " | cutoff: ", ppi_prediction_result$cutoff,
                                      " | iterations: ", ppi_prediction_result$iter)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12),
                     plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                     legend.position = "right")
  } else if (y == 'frac') {
    plot <- Df %>%
      ggplot2::ggplot(aes(x = reference, y = recovery, fill = reference, ymin = recovery, ymax = recovery + sep)) +
      ggplot2::geom_bar(stat = "identity", col = "black") +
      ggplot2::geom_errorbar(width = 0.2) +
      ggplot2::geom_text(mapping = aes(label = round(recovery, 1)), vjust = -3, color = "black", size = 2) +
      ggplot2::ylim(0, 100) +
      ggplot2::scale_fill_manual(values = c("#6CA6C1", "#D930C5", "#FFE66D", "#8ED2C6", "#01A79D", "#FFB400", "#F6511D", "#00a6ed", "#0a2463")) +
      ggplot2::facet_wrap(~ probability, ncol = 3) +
      ggplot2::labs(x = "", y = "Value", title = paste0("Mean + IQR Learning Curves for the ", ppi_prediction_result$ensembleSize, " ", ppi_prediction_result$model.type, " models"),
                    subtitle = paste0("sampling: ", base::ifelse(ppi_prediction_result$sampling == "unweighted", "unweighted", ppi_prediction_result$sampling),
                                      if(ppi_prediction_result$model.type == "svm") paste0(" | kernel type: ", ppi_prediction_result$kernelType),
                                      " | cutoff: ", ppi_prediction_result$cutoff,
                                      " | iterations: ", ppi_prediction_result$iter)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12),
                     plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                     legend.position = "right")
  }
  return(plot)
}
