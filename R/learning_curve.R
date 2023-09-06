usethis::use_package('dplyr')
usethis::use_package('tibble')
usethis::use_package('tidyr')
usethis::use_package('ggplot2')
usethis::use_package('ggpubr')

#' Function to calculate and plot the accuracy, hinge loss and binary cross-entropy loss of a machine learning prediction using the ppi.prediction function.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import ggplot2
#' @import ggpubr
#'
#' @param ppi_prediction_result: result from the function ppi.prediction
#' @param train_sizes: base::sequence of fraction of training sizes to be used for calculation from >0 to 1
#' @param models: Integer of models used to calculate the loss functions. If "all" then all models as specified by the ensembleSize in ppi.prediction will be used.
#'
#' @return
#' @export
#'
#' @examples
#'
learning.curve <- function(ppi_prediction_result, train_sizes = base::seq(0.1, 1, by = 0.1), models = "all") {
  #extract results from ppi.prediction function
  train_data = ppi_prediction_result$training.sets
  test_data = ppi_prediction_result$testMat
  test_data_labels <- base::as.data.frame(test_data) %>% tibble::rownames_to_column("sample") %>%
    tidyr::separate(col = "sample", into = c("complex", "interaction", "sample", "orientation"), sep = ";") %>%
    dplyr::mutate(reference = complex) %>%
    dplyr::pull(reference)
  if(!any(test_data_labels %in% ppi_prediction_result$negative.reference)) {
    stop("The test data must contain positive and negative reference interactions.")
  }
  negative_reference = ppi_prediction_result$negative.reference
  if(models == "all") {
    models = base::seq(ppi_prediction_result$model.e)
  }
  else if(is.double(models)) {
    models = base::seq(round(models, 0))
  }
  else if(is.double(models) | models == "all") {
    stop("'models' has to be an integer or 'all', referring to all models in the results from the ppi.prediction$model.e")
  }

  #define custom functions
  `%ni%` <- Negate(`%in%`)
  binary_cross_entropy_loss <- function(actual_labels, predicted_probs) {
    -base::mean(actual_labels * log(predicted_probs) + (1 - actual_labels) * log(1 - predicted_probs))
  }
  calculate_hinge_loss <- function(actual_labels, predicted_labels) {
    loss <- 1 - actual_labels * predicted_labels
    hinge_loss <- base::ifelse(loss < 0, 0, loss)
    return(hinge_loss)
  }
  multiAdaSampling <- function(train.mat, test.mat, label, model_type, kernelType = kernelType, iter = iter, cost = C, gamma = gamma, degree = degree, coef0 = coef0) {

    X <- train.mat
    Y <- label

    model <- NULL
    pred.mat <- NULL
    accuracy <- rep(NA, iter)

    for (i in seq_len(iter)) {
      tmp <- X
      base::rownames(tmp) <- NULL

      if (is.null(gamma)) {
        gamma <- 0.1 / length(tmp)
      }

      if (model_type == "svm") {
        model <- e1071::svm(Y ~ ., data = data.frame(tmp, Y), type = "C-classification", kernel = kernelType, probability = TRUE, cost = cost, gamma = gamma, degree = degree, coef0 = coef0)
        pred.train <- stats::predict(model, train.mat, decision.values = TRUE, probability = TRUE)
        pred.mat <- attr(pred.train, "probabilities")
      } else if (model_type == "randomforest") {
        model <- randomForest::randomForest(X, base::factor(Y))
        pred.train <- stats::predict(model, train.mat, type = "prob")
        pred.mat <- pred.train[,c(2,1)]
      } else {
        stop("Invalid model_type. Please choose 'svm' or 'randomforest'.")
      }

      X <- c()
      Y <- c()

      for (j in seq_len(ncol(pred.mat))) {
        voteClass <- pred.mat[label == base::colnames(pred.mat)[j], ]
        idx <- base::sample(seq_len(base::nrow(voteClass)),
                      size = base::nrow(voteClass), replace = TRUE,
                      prob = voteClass[, j])
        X <- base::rbind(X, base::as.matrix(train.mat[base::rownames(voteClass)[idx], ]))
        Y <- c(Y, label[base::rownames(voteClass)[idx]])
      }
      base::colnames(X) <- base::colnames(train.mat)
    }

    if(model_type == "svm") {
      pred.values <- stats::predict(model, newdata = test.mat, decision.values = TRUE, probability = TRUE)
      pred <- attr(pred.values, "probabilities")
    } else if (model_type == "randomforest") {
      pred.values <- stats::predict(model, newdata = test.mat, type = "prob")
      pred <- pred.values[,c(2,1)]
    } else {
      stop("Invalid model_type. Please choose 'svm' or 'randomforest'.")
    }


    return(list(pred, pred.values, model))
  }

  train_performance <- base::matrix(0, nrow = length(train_sizes), ncol = length(models))  # Store training performance
  test_performance <- base::matrix(0, nrow = length(train_sizes), ncol = length(models))  # Store test performance
  train_loss <- base::matrix(0, nrow = length(train_sizes), ncol = length(models))  # Store training loss
  test_loss <- base::matrix(0, nrow = length(train_sizes), ncol = length(models))  # Store test loss
  train_hinge <- base::matrix(0, nrow = length(train_sizes), ncol = length(models))  # Store training hinge
  test_hinge <- base::matrix(0, nrow = length(train_sizes), ncol = length(models))  # Store test hinge

  test_data <- base::as.data.frame(test_data) %>% tibble::rownames_to_column("sample") %>%
    tidyr::separate(col = "sample", into = c("complex", "interaction", "sample", "orientation"), sep = ";") %>%
    dplyr::mutate(reference = complex) %>%
    tidyr::unite(complex, reference, interaction, sample, orientation, col = "sample", sep = ";") %>%
    tibble::column_to_rownames("sample")

  for (i in seq_along(models)) {
    train_labels <- train_data[[i]] %>% tibble::as_tibble() %>% mutate("id" = base::rownames(train_data[[i]])) %>% tidyr::separate(col = "id", into = c("reference"), extra = "drop", sep = ";") %>% dplyr::pull(reference)
    train_labels <- as.integer(!base::grepl(paste(negative_reference, collapse = "|"), train_labels))

    for (j in seq_along(train_sizes)) {
      if(j %% 10 == 0) {
        message(train_sizes[j]*100, "% training size for model ", i, " completed.")
      }
      size <- train_sizes[j]

      # Select a subset of the training data
      subset_indices_prs <- base::sample(which(train_labels == 1), size = round(size * length(train_labels[train_labels == 1])))
      subset_indices_rrs <- base::sample(which(train_labels == 0), size = round(size * length(train_labels[train_labels == 0])))
      subset_indices <- c(subset_indices_prs, subset_indices_rrs)
      subset_train_data <- train_data[[i]][subset_indices, , drop = FALSE] %>% base::as.matrix()
      subset_train_labels <- train_labels[subset_indices]
      test_data_subset <- subset(test_data, base::rownames(test_data) %ni% base::rownames(train_data[[i]])) %>% base::as.matrix()
      test_labels <- test_data_subset %>% tibble::as_tibble() %>% mutate("id" = base::rownames(test_data_subset)) %>% tidyr::separate(col = "id", into = c("reference"), extra = "drop", sep = ";") %>% dplyr::pull(reference)
      test_labels <- base::ifelse(test_labels == "PRS", 1, 0)

      # Make predictions on the training and test data using the trained model
      cls = base::as.factor(subset_train_labels+1)
      names(cls) <- base::rownames(subset_train_data)
      model <- multiAdaSampling(train.mat = subset_train_data, test.mat = test_data_subset, model_type = ppi_prediction_result$model.type, label = cls, iter = ppi_prediction_result$iter,
                                kernelType = ppi_prediction_result$kernelType, cost = ppi_prediction_result$C, degree = ppi_prediction_result$degree, gamma = ppi_prediction_result$gamma, coef0 = ppi_prediction_result$coef0)
      train_predictions <- attr(stats::predict(model[[3]], subset_train_data, probability = TRUE), "probabilities")[,1]
      test_predictions <- attr(stats::predict(model[[3]], test_data_subset, probability = TRUE), "probabilities")[,1]
      train_predictions_lables <- base::ifelse(train_predictions>0.5, 1, 0)
      test_predictions_lables <- base::ifelse(test_predictions>0.5, 1, 0)

      # Calculate performance metrics accuracy
      train_perf <- base::mean(train_predictions_lables == subset_train_labels)
      test_perf <- base::mean(test_predictions_lables == test_labels)

      # Calculate 'Binary Cross-Entropy Loss' and 'Hinge loss'
      train_loss_fold <- binary_cross_entropy_loss(train_predictions_lables, train_predictions) #-sum(fold_train_labels * log(train_predictions) + (1 - fold_train_labels) * log(1 - train_predictions))
      test_loss_fold <- binary_cross_entropy_loss(test_predictions_lables, test_predictions) #-sum(fold_val_labels * log(val_predictions) + (1 - fold_val_labels) * log(1 - val_predictions))
      train_hinge_loss <- calculate_hinge_loss(actual_labels = base::ifelse(subset_train_labels == 0, -1, subset_train_labels), predicted_labels = base::ifelse(train_predictions_lables == 0, -1, train_predictions_lables))
      test_hinge_loss <- calculate_hinge_loss(actual_labels = base::ifelse(test_labels == 0, -1, test_labels), predicted_labels = base::ifelse(test_predictions_lables == 0, -1, test_predictions_lables))

      # Store the performance metrics
      train_performance[j, i] <- train_perf
      test_performance[j, i] <- test_perf
      train_loss[j, i] <- train_loss_fold
      test_loss[j, i] <- test_loss_fold
      train_hinge[j, i] <- base::mean(train_hinge_loss)
      test_hinge[j, i] <- base::mean(test_hinge_loss)
    }
  }

  #Store results in a data frame
  df <- data.frame(Training_Size = rep(train_sizes, times = length(ppi_prediction_result$model.e)*2),
                   Label = c(rep("Training", length(train_performance)), rep("Test", length(test_performance))),
                   Accuracy = c(train_performance, test_performance),
                   Loss = c(train_loss, test_loss),
                   Hinge = c(train_hinge, test_hinge),
                   Model = rep(1:length(ppi_prediction_result$model.e), each = length(train_sizes)))
  df$Label <- base::factor(df$Label, levels = c("Training", "Test"))
  #Summarise the results from the 50 models: calculate mean, sd, 5%, 25%, 75%, 95% CI
  df2 <- tidyr::gather(df, key = "Performance_Type", value = "Value",
                       Accuracy, Loss, Hinge) %>%
    dplyr::group_by(Training_Size, Label, Performance_Type) %>%
    dplyr::summarise(Mean = base::mean(Value),
                     Sd = sd(Value),
                     Ci5 = quantile(Value, probs = 0.05),
                     Ci25 = quantile(Value, probs = 0.25),
                     Ci75 = quantile(Value, probs = 0.75),
                     Ci95 = quantile(Value, probs = 0.95)) %>%
    dplyr::ungroup()

  #Plot the accuracy, hinge loss and binary cross-entropy loss
  plot <- ggplot(df2, aes(x = Training_Size, y = Mean, color = Label, fill = Label)) +
    ggplot2::geom_ribbon(aes(ymin = Ci25, ymax = Ci75), alpha = 0.15, linetype = "dotdash", size = 0.5) +
    ggplot2::geom_line(size = 1) +
    ggplot2::facet_wrap(. ~ Performance_Type, scales = "free", labeller = labeller(Performance_Type = c(Loss = "Binary Cross-Entropy Loss", Accuracy = "Accuracy", Hinge = "Hinge Loss"))) +
    ggplot2::labs(x = "Fraction of Training Set Size", y = "Value", title = paste0("Mean + IQR Learning Curves for the ", ensembleSize," SVM models (", base::ifelse(ppi_prediction_result$sampling == "unweighted", "unweighted", ppi_prediction_result$sampling), "sampling)")) +
    ggplot2::scale_color_manual(values = c("Training" = "#6CA6C1", "Test" = "#D930C5")) +
    ggplot2::scale_fill_manual(values = c("Training" = "#6CA6C1", "Test" = "#D930C5")) +
    ggpubr::theme_pubr() +
    ggplot2::theme(text = element_text(family = "Avenir"),
          plot.title = element_text(size = 12),
          axis.title = element_text(family = "Avenir Medium")) +
    ggplot2::scale_x_continuous(limits = c(0,NA), breaks = base::seq(0, 1, by=0.2))

  # Return the learning curves with performance and loss
  return(list(Training_Sizes = train_sizes,
              Training_Performance = train_performance,
              Test_Performance = test_performance,
              Training_Loss = train_loss,
              Test_Loss = test_loss,
              Train_Hinge = train_hinge,
              Test_Hinge = test_hinge,
              Df = df,
              Df_sumstat = df2,
              learning_plot = plot))
}
