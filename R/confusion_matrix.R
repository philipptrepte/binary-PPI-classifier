
#' Function to calculate a confusion matrix for
#'
#' @param ppi_prediction_result
#'
#' @return
#' @export
#'
#' @examples
#' data("example_ppi_prediction")
#' confusion_matrix <- confusion.matrix(example_ppi_prediction)
#'
confusion.matrix <- function(ppi_prediction_result) {
  data <- ppi_prediction_result$predDf %>%
    group_by(interaction, complex) %>%
    summarise(prob = max(predMat, na.rm = TRUE))
  confusion_mat <- list()
  confusion_df <- data.frame()
  for(i in 1:3) {
    j = c(0.5, 0.75, 0.95)[i]
    confusion_mat[[i]] <- caret::confusionMatrix(factor(ifelse(data$prob >= j, "PRS", "RRS")),
                                          factor(data$complex))
    confusion_df <- rbind(confusion_df, confusion_mat[[i]]$table %>%
                            cbind(probability = as.numeric(j)*100))
    if(i == 3) {
      names(confusion_mat) <- c("50%", "75%", "95%")
      rownames(confusion_df) <- c("PRS: 50%", "RRS: 50%", "PRS: 75%", "RRS: 75%", "PRS: 95%", "RRS: 95%")
    }
  }
  return(list(confusionMatrixList = confusion_mat,
              confusionMatrixDf = confusion_df))
}


