#' Plot a bar diagram of the recovery rates at 50%, 75% and 95% probability
#'
#' @param ppi_prediction_result
#'
#' @return
#' @export
#'
#' @examples
recovery.plot <- function(ppi_prediction_result) {
  ppi_prediction_result$predDf %>%
    dplyr::mutate(reference = ifelse(str_detect(complex, "inter-complex|RRS"), "RRS", "PRS")) %>%
    group_by(interaction, reference) %>%
    dplyr::summarise(score = max(predMat, na.rm = TRUE)) %>%
    dplyr::mutate(pos = ifelse(score > 0.5, 1, 0),
                  probability = 0.5) %>%
    group_by(reference, probability) %>%
    dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                     total = n(),
                     recovery = pos / total * 100,
                     sep = sqrt(recovery*(100-recovery)/total)) %>%
    rbind(
      example_ppi_prediction$predDf %>%
        dplyr::mutate(reference = ifelse(str_detect(complex, "inter-complex|RRS"), "RRS", "PRS")) %>%
        group_by(interaction, reference) %>%
        dplyr::summarise(score = max(predMat, na.rm = TRUE)) %>%
        dplyr::mutate(pos = ifelse(score > 0.75, 1, 0),
                      probability = 0.75) %>%
        group_by(reference, probability) %>%
        dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                         total = n(),
                         recovery = pos / total * 100,
                         sep = sqrt(recovery*(100-recovery)/total))
    ) %>%
    rbind(
      example_ppi_prediction$predDf %>%
        dplyr::mutate(reference = ifelse(str_detect(complex, "inter-complex|RRS"), "RRS", "PRS")) %>%
        group_by(interaction, reference) %>%
        dplyr::summarise(score = max(predMat, na.rm = TRUE)) %>%
        dplyr::mutate(pos = ifelse(score > 0.95, 1, 0),
                      probability = 0.95) %>%
        group_by(reference, probability) %>%
        dplyr::summarise(pos = sum(pos, na.rm = TRUE),
                         total = n(),
                         recovery = pos / total * 100,
                         sep = sqrt(recovery*(100-recovery)/total))
    ) %>%
    ggplot(aes(x = reference, y = recovery, fill = reference, ymin = recovery, ymax = recovery + sep)) +
    geom_bar(stat = "identity", col = "black") +
    geom_errorbar(width = 0.2) +
    geom_text(mapping = aes(label = round(recovery, 1)), vjust = -3, color = "black", size = 2) +
    ylim(0, 100) +
    scale_fill_manual(values = c("#6CA6C1", "#D930C5")) +
    facet_wrap(~ probability, ncol = 3) +
    theme_pubr() +
    theme(legend.position = "none") +
    theme_pubr() +
    ggplot2::theme(text = element_text(family = "Avenir"),
                   plot.title = element_text(size = 12),
                   plot.subtitle = element_text(size = 8),
                   axis.title = element_text(family = "Avenir Medium"),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                   legend.position = "right")
}
