#' Plot a bar diagram of the recovery rates at 50%, 75% and 95% probability
#'
#' @import dplyr
#' @import ggplot2
#' @import stringr
#' @import ggpubr
#'
#' @param ppi_prediction_result: result object from the ppi.prediction() function.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
recovery.plot <- function(ppi_prediction_result) {
  ppi_prediction_result$predDf %>%
    dplyr::mutate(reference = base::ifelse(stringr::str_detect(complex, "inter-complex|RRS"), "RRS", "PRS")) %>%
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
      example_ppi_prediction$predDf %>%
        dplyr::mutate(reference = base::ifelse(stringr::str_detect(complex, "inter-complex|RRS"), "RRS", "PRS")) %>%
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
      example_ppi_prediction$predDf %>%
        dplyr::mutate(reference = base::ifelse(stringr::str_detect(complex, "inter-complex|RRS"), "RRS", "PRS")) %>%
        dplyr::group_by(interaction, reference) %>%
        dplyr::summarise(score = max(predMat, na.rm = TRUE)) %>%
        dplyr::mutate(pos = base::ifelse(score > 0.95, 1, 0),
                      probability = 0.95) %>%
        dplyr::group_by(reference, probability) %>%
        dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                         total = n(),
                         recovery = pos / total * 100,
                         sep = sqrt(recovery*(100-recovery)/total))
    ) %>%
    ggplot2::ggplot(aes(x = reference, y = recovery, fill = reference, ymin = recovery, ymax = recovery + sep)) +
    ggplot2::geom_bar(stat = "identity", col = "black") +
    ggplot2::geom_errorbar(width = 0.2) +
    ggplot2::geom_text(mapping = aes(label = round(recovery, 1)), vjust = -3, color = "black", size = 2) +
    ggplot2::ylim(0, 100) +
    ggplot2::scale_fill_manual(values = c("#6CA6C1", "#D930C5")) +
    ggplot2::facet_wrap(~ probability, ncol = 3) +
    ggpubr::theme_pubr() +
    ggplot2::theme(text = element_text(family = "Avenir"),
                   plot.title = element_text(size = 12),
                   plot.subtitle = element_text(size = 8),
                   axis.title = element_text(family = "Avenir Medium"),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                   legend.position = "right")
}
