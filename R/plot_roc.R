#' Plot ROC curve of the predicted probabilities against the assay features the ML was trained on
#'
#' @import dplyr
#' @import ggplot2
#' @import stringr
#' @import ggpubr
#' @import plotROC
#' @import DescTools
#' @import viridis
#' @import purrr
#'
#' @param ppi_prediction_result: result object from the ppi.prediction() function.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
roc.plot <- function(ppi_prediction_result) {
  n.assays <- length(ppi_prediction_result$assay)
  prob.aggregated <- ppi_prediction_result$predDf %>%
    dplyr::mutate(D = ifelse(stringr::str_detect(complex, "inter-complex|RRS"), "0", "1")) %>%
    dplyr::group_by(interaction, D) %>%
    dplyr::summarise(M = max(predMat, na.rm = TRUE))

  assay.aggregated <- list()
  for(a in seq(n.assays)) {
    assay.aggregated[[a]] <- ppi_prediction_result$predDf %>%
      dplyr::mutate(D = ifelse(stringr::str_detect(complex, "inter-complex|RRS"), "0", "1")) %>%
      dplyr::group_by(interaction, D) %>%
      dplyr::summarise(M = max(!!(rlang::sym(ppi_prediction_result$assay[a])), na.rm = TRUE))
  }

  prob.aggregated.roc <- plotROC::calculate_roc(M = prob.aggregated$M, D = prob.aggregated$D) %>% #, cutoffs = seq(from = 0, to = 1, by = 1/(nrow(referenceSet.BRET.prob.aggregated)-1))) %>%
    dplyr::group_by(TPF, FPF) %>%
    dplyr::summarise(c = max(c, na.rm = TRUE)) %>%
    dplyr::mutate(data = "probability")

  assay.aggregated.roc <- list()
  for(r in seq(n.assays)) {
    assay.aggregated.roc[[r]] <- plotROC::calculate_roc(M = assay.aggregated[[r]]$M, D = assay.aggregated[[r]]$D) %>%
      dplyr::group_by(TPF, FPF) %>%
      dplyr::summarise(c = max(c, na.rm = TRUE)) %>%
      dplyr::mutate(data = ppi_prediction_result$assay[r])
  }

  annotations <- list()
  for(t in seq(n.assays)) {
    if(t == 1) {
      text <- base::paste0("AUC = ", round(DescTools::AUC(x = prob.aggregated.roc$FPF, y = prob.aggregated.roc$TPF), 2))
      x_pos = 0.8
      y_pos = 0.0 + 0.1*(t-1)
      annotations[[t]] <- ggplot2::annotate(geom = "text", x = x_pos, y = y_pos, color = viridis::viridis_pal()(n.assays+1)[t],
                                            label = text)
    }
    text <- paste0("AUC = ", round(DescTools::AUC(x = assay.aggregated.roc[[t]]$FPF, y = assay.aggregated.roc[[t]]$TPF), 2))
    x_pos = 0.8
    y_pos = 0.0 + 0.1*(t)
    annotations[[t+1]] <- ggplot2::annotate(geom = "text", x = x_pos, y = y_pos, color = viridis::viridis_pal()(n.assays+1)[t+1],
                                            label = text)
  }

  prob.aggregated.roc %>%
    base::rbind(purrr::list_rbind(assay.aggregated.roc)) %>%
    dplyr::mutate(data = base::factor(data, levels = c("probability", ppi_prediction_result$assay))) %>%
    ggplot2::ggplot(aes(x = FPF, y = TPF, color = data)) +
    ggplot2::geom_line() +
    viridis::scale_color_viridis(discrete = TRUE) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.1)) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.1)) +
    annotations +
    ggpubr::theme_pubr() +
    ggplot2::theme(text = element_text(family = "Avenir"),
                   plot.title = element_text(size = 12),
                   plot.subtitle = element_text(size = 8),
                   axis.title = element_text(family = "Avenir Medium"),
                   legend.position = "right") +
    ggplot2::labs(color = "")
}
