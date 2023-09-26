#' Function to calculate and plot the accuracy, hinge loss and binary cross-entropy loss of a machine learning prediction using the ppi.prediction function.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import ggplot2
#' @import ggpubr
#' @import cowplot
#'
#' @param ppi_prediction_result: result from the function ppi.prediction
#' @param train_sizes: sequence of fraction of training sizes to be used for calculation from >0 to 1
#' @param models: Integer of models used to calculate the loss functions. If "all" then all models as specified by the ensembleSize in ppi.prediction will be used.
#' @param verbose: boolean, prints detailed informations
#' @param prob: probability cutoff to calculate accuracy
#'
#' @return a list with elements
#' @export
#'
#' @examples
learning.curve <- function(ppi_prediction_result, train_sizes = base::seq(0.1, 1, by = 0.1), models = "all", verbose = TRUE, prob = 0.5) {
  set.seed(ppi_prediction_result$seed)
  #extract results from ppi.prediction function
  train_data = ppi_prediction_result$training.sets
  test_data = ppi_prediction_result$predTrainDf %>%
    dplyr::filter(!is.na(predTrainMat))
  test_data_labels <- test_data %>%
    dplyr::pull(reference)

  if(!any(sapply(ppi_prediction_result$negative.reference, function(x) any(str_detect(test_data_labels, x))))) {
    stop("The test data must contain positive and negative reference interactions.")
  }
  negative_reference = ppi_prediction_result$negative.reference
  if(models == "all") {
    n.models = base::seq(ppi_prediction_result$model.e)
  } else if(is.double(models)) {
    n.models = base::seq(round(models, 0))
  } else if(is.double(models) | models == "all") {
    stop("'models' has to be an integer or 'all', referring to all models in the results from the ppi.prediction$model.e")
  }

  #define custom functions
  `%ni%` <- Negate(`%in%`)

  binary_cross_entropy_loss <- function(actual_labels, predicted_probs) {
    entropy_loss <- -median(actual_labels * log(predicted_probs) + (1 - actual_labels) * log(1 - predicted_probs))
    return(entropy_loss)
  }

  calculate_hinge_loss <- function(actual_labels, predicted_labels) {
    loss <- 1 - actual_labels * predicted_labels
    hinge_loss <- base::ifelse(loss < 0, 0, loss)
    return(hinge_loss)
  }

  multiAdaSampling <- function(train.mat, label, model_type, kernelType = kernelType, iter = iter, cost = C, gamma = gamma, degree = degree, coef0 = coef0) {

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

    return(model)
  }

  train_performance <- base::matrix(0, nrow = length(train_sizes), ncol = length(n.models))  # Store training performance
  test_performance <- base::matrix(0, nrow = length(train_sizes), ncol = length(n.models))  # Store test performance
  train_loss <- base::matrix(0, nrow = length(train_sizes), ncol = length(n.models))  # Store training loss
  test_loss <- base::matrix(0, nrow = length(train_sizes), ncol = length(n.models))  # Store test loss
  train_hinge <- base::matrix(0, nrow = length(train_sizes), ncol = length(n.models))  # Store training hinge
  test_hinge <- base::matrix(0, nrow = length(train_sizes), ncol = length(n.models))  # Store test hinge

  test_data_mat <- test_data %>%
    tidyr::unite(complex, reference, interaction, sample, orientation, col = "sample", sep = ";") %>%
    dplyr::select(sample, ppi_prediction_result$assay) %>%
    tibble::column_to_rownames("sample") %>% as.matrix()

  for (i in seq_along(n.models)) {
    train_labels <- train_data[[i]] %>% tibble::as_tibble() %>% mutate("id" = base::rownames(train_data[[i]])) %>% tidyr::separate(col = "id", into = c("complex", "reference"), extra = "drop", sep = ";") %>% dplyr::pull(reference)
    train_labels <- as.integer(ifelse(train_labels == "RRS", 0, 1))
    names(train_labels) <- base::rownames(train_data[[i]])

    for (j in seq_along(train_sizes)) {
      if(j == max(seq_along(train_sizes)) & verbose == TRUE) {
        message("learning curve for model ", i, " completed.")
      }
      size <- train_sizes[j]
      max_size <- train_sizes[max(seq_along(train_sizes))]
      if(round(size * length(train_labels[train_labels == 1])) == 1){
        if(verbose == TRUE) {
          message(paste0("at ", train_sizes[j]*100, "% of the training set, the training subset with only one interaction is to small. Training will continue with ", train_sizes[j+1]*100, "% of training set."))
        }
        next()
      }

      # Select a subset of the training data
      subset_indices_prs <- base::sample(which(train_labels == 1), size = round(size * length(train_labels[train_labels == 1])))
      subset_indices_rrs <- base::sample(which(train_labels == 0), size = round(size * length(train_labels[train_labels == 0])))
      subset_indices <- c(subset_indices_prs, subset_indices_rrs)
      subset_train_data <- train_data[[i]][subset_indices, , drop = FALSE] %>% base::as.matrix()
      subset_train_labels <- train_labels[subset_indices]
      assertthat::assert_that(all(str_detect(names(subset_train_labels[subset_train_labels == 0]), paste(negative_reference, collapse = "|"))),
                              msg = 'Train labels do not contain negative reference interactions')

      test_data_subset <- subset(test_data_mat, base::rownames(test_data_mat) %ni% base::rownames(subset_train_data)) %>% base::as.matrix()
      test_data_subset_labels <- test_data_subset %>% tibble::as_tibble() %>% mutate("id" = base::rownames(test_data_subset)) %>% tidyr::separate(col = "id", into = c("complex", "reference"), extra = "drop", sep = ";") %>% dplyr::pull(reference)
      test_data_subset_labels <- as.integer(ifelse(test_data_subset_labels == "RRS", 0, 1))
      n_prs <- length(test_data_subset_labels[test_data_subset_labels == 1])
      n_rrs <- length(test_data_subset_labels[test_data_subset_labels == 0])
      n_subset <- min(n_prs, n_rrs, round(max_size * length(train_labels)))
      test_indices_prs <- sample(rownames(test_data_subset)[test_data_subset_labels == 1], size = n_subset)
      test_indices_rrs <- sample(rownames(test_data_subset)[test_data_subset_labels == 0], size = n_subset)
      test_subset_indices <- c(test_indices_prs, test_indices_rrs)
      test_data_subset <- subset(test_data_subset, rownames(test_data_subset) %in% test_subset_indices)
      assertthat::assert_that(any(base::rownames(test_data_subset) %ni% base::rownames(subset_train_data)),
                              msg="Test set contains training set interactions")
      test_labels <- test_data_subset %>% tibble::as_tibble() %>% mutate("id" = base::rownames(test_data_subset)) %>% tidyr::separate(col = "id", into = c("complex", "reference"), extra = "drop", sep = ";") %>% dplyr::pull(reference)
      test_labels <- as.integer(ifelse(test_labels == "RRS", 0, 1))
      names(test_labels) <- base::rownames(test_data_subset)
      assertthat::assert_that(all(str_detect(names(test_labels[test_labels == 0]), paste(negative_reference, collapse = "|"))),
                              msg = 'Test labels do not contain negative reference interactions')

      # Make predictions on the training and test data using the trained model
      cls = base::as.factor(subset_train_labels+1)
      names(cls) <- base::rownames(subset_train_data)
      model <- multiAdaSampling(train.mat = subset_train_data, model_type = ppi_prediction_result$model.type, label = cls, iter = ppi_prediction_result$iter,
                                kernelType = ppi_prediction_result$kernelType, cost = ppi_prediction_result$C, degree = ppi_prediction_result$degree, gamma = ppi_prediction_result$gamma, coef0 = ppi_prediction_result$coef0)

      if(ppi_prediction_result$model.type == "svm") {
        train_predictions <- stats::predict(model, newdata = subset_train_data, decision.values = TRUE, probability = TRUE)
        train_predictions <- attr(train_predictions, "probabilities")
        train_predictions <- ifelse(train_predictions[,1] == 0, 1e-4, train_predictions[,1])
        train_predictions <- ifelse(train_predictions == 1, 0.999, train_predictions)
        train_predictions_lables <- ifelse(train_predictions > prob, 1, 0)

        test_predictions <- stats::predict(model, newdata = test_data_subset, decision.values = TRUE, probability = TRUE)
        test_predictions <- attr(test_predictions, "probabilities")
        test_predictions <- ifelse(test_predictions[,1] == 0, 1e-4, test_predictions[,1])
        test_predictions <- ifelse(test_predictions == 1, 0.999, test_predictions)
        test_predictions_lables <- ifelse(test_predictions > prob, 1, 0)
      } else if (ppi_prediction_result$model.type == "randomforest") {
        train_predictions <- stats::predict(model, newdata = subset_train_data, type = "prob")
        train_predictions <- ifelse(train_predictions[,2] == 0, 1e-4, train_predictions[,2])
        train_predictions <- ifelse(train_predictions == 1, 0.999, train_predictions)
        train_predictions_lables <- ifelse(train_predictions > prob, 1, 0)

        test_predictions <- stats::predict(model, newdata = test_data_subset, type = "prob")
        test_predictions <- ifelse(test_predictions[,2] == 0, 1e-4, test_predictions[,2])
        test_predictions <- ifelse(test_predictions == 1, 0.999, test_predictions)
        test_predictions_lables <- ifelse(test_predictions > prob, 1, 0)
      } else {
        stop("Invalid model.type. Please choose 'svm' or 'randomforest'.")
      }

      # Calculate performance metrics accuracy
      assertthat::assert_that(all(names(train_predictions_lables) == rownames(train_data[[i]][subset_indices,])),
                              msg = "Predicted training set not equal to input training set")
      train_perf <- base::mean(train_predictions_lables == subset_train_labels)
      assertthat::assert_that(all(names(test_predictions_lables) == rownames(test_data_subset)),
                              msg = "Predicted test set not equal to input test set")
      test_perf <- base::mean(test_predictions_lables == test_labels)

      # Calculate 'Binary Cross-Entropy Loss' and 'Hinge loss'
      assertthat::assert_that(all(names(subset_train_labels) == names(train_predictions)),
                              msg = "Predicted training set not equal to input training set")
      train_loss_fold <- binary_cross_entropy_loss(subset_train_labels, train_predictions)
      assertthat::assert_that(all(names(test_labels) == names(test_predictions)),
                              msg = "Predicted test set not equal to input test set")
      test_loss_fold <- binary_cross_entropy_loss(actual_labels = test_labels, predicted_probs = test_predictions)

      assertthat::assert_that(all(names(subset_train_labels) == names(train_predictions_lables)),
                              msg = "Predicted training set not equal to input training set")
      train_hinge_loss <- calculate_hinge_loss(actual_labels = base::ifelse(subset_train_labels == 0, -1, subset_train_labels),
                                               predicted_labels = base::ifelse(train_predictions_lables == 0, -1, train_predictions_lables))
      assertthat::assert_that(all(names(test_labels) == names(test_predictions_lables)),
                              msg = "Predicted test set not equal to input test set")
      test_hinge_loss <- calculate_hinge_loss(actual_labels = base::ifelse(test_labels == 0, -1, test_labels),
                                              predicted_labels = base::ifelse(test_predictions_lables == 0, -1, test_predictions_lables))

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
    ggplot2::labs(x = "Fraction of Training Set Size", y = "Value", title = paste0("Mean + IQR Learning Curves for the ", ppi_prediction_result$ensembleSize, " ", ppi_prediction_result$model.type, " models"),
                  subtitle = paste0("sampling: ", base::ifelse(ppi_prediction_result$sampling == "unweighted", "unweighted", ppi_prediction_result$sampling),
                                    if(ppi_prediction_result$model.type == "svm") paste0(" | kernel type: ", ppi_prediction_result$kernelType),
                                    " | cutoff: ", ppi_prediction_result$cutoff,
                                    " | iterations: ", ppi_prediction_result$iter)) +
    ggplot2::scale_color_manual(values = c("Training" = "#6CA6C1", "Test" = "#D930C5")) +
    ggplot2::scale_fill_manual(values = c("Training" = "#6CA6C1", "Test" = "#D930C5")) +
    ggpubr::theme_pubr() +
    ggplot2::theme(text = element_text(family = "Avenir"),
                   plot.title = element_text(size = 12),
                   plot.subtitle = element_text(size = 8),
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

