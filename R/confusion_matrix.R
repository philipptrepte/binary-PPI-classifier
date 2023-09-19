#' Function to calculate a confusion matrix for
#'
#' @import dplyr
#' @import caret
#'
#' @param ppi_prediction_result: result object from the function ppi.prediction
#'
#' @return a list with elements
#' @export
#'
#' @examples
#' data("example_ppi_prediction")
#' example_confusion_matrix <- confusion.matrix(example_ppi_prediction)
#'
confusion.matrix <- function(ppi_prediction_result) {
  data <- ppi_prediction_result$predDf %>%
    dplyr::filter(!is.na(predMat)) %>%
    dplyr::mutate(complex = ifelse(stringr::str_detect(complex, paste(ppi_prediction_result$negative.reference, collapse = "|")),
                                   "RRS", "PRS")) %>%
    dplyr::group_by(interaction, complex) %>%
    dplyr::summarise(prob = max(predMat, na.rm = TRUE)) %>%
    ungroup()
  confusion_mat <- list()
  confusion_df <- base::data.frame()
  for (i in 1:3) {
    j = c(0.5, 0.75, 0.95)[i]
    confusion_mat[[i]] <- caret::confusionMatrix(base::factor(base::ifelse(data$prob >= j, "PRS", "RRS")),
                                                 base::factor(data$complex))
    confusion_df <- base::rbind(confusion_df, confusion_mat[[i]]$table %>%
                                  base::cbind(probability = as.numeric(j) * 100))
    if (i == 3) {
      names(confusion_mat) <- c("50%", "75%", "95%")
      rownames(confusion_df) <- c("PRS: 50%", "RRS: 50%",
                                  "PRS: 75%", "RRS: 75%", "PRS: 95%", "RRS: 95%")
    }
  }
  return(list(confusionMatrixList = confusion_mat, confusionMatrixDf = confusion_df))
}
