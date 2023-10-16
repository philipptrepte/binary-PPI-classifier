#' Plot a probability grid from the mean probabilities from the 'ensembleSize' number of models
#'
#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#' @import rlang
#' @import stats
#' @import varhandle
#' @import ggnewscale
#'
#' @param ppi_prediction_result: result object from the ppi.prediction() function.
#' @param n: grid size. For 3 assays, will be limited to grid size of 40 to reduce computing time.
#' @param x.log.scale: logical to log-scale x-axis values
#' @param xlim: Numeric vector of two values specifying the left and right limit of the scale
#' @param ylim: Numeric vector of two values specifying the bottom and top limit of the scale
#' @param set: Character. PPI set to generate the plot for: "test" or "train"
#' @param model: Integer (1L) or "all". Plots the decision boundaries for a specific model (e.g. 1L for model 1) or the mean of all models.
#' @param x.nudge: Numerical. Which value to add to log transformation of x-axis values, in case of negative x values.
#' @param type: Character. Specify to plot as "2D" or "3D" plot for trainings with 3 features.
#' @param assay: Character. Specifies which assays to plot against each other. Must be one of the training features.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' data("example_ppi_prediction")
#' probGrid.plot(example_ppi_prediction)
probGrid.plot <- function(ppi_prediction_result, n=100, x.log.scale = TRUE, xlim = c(NA,NA), ylim = NULL, set="train", model="all", training_set = 'all', x.nudge = 1, type = "2D", assay = ppi_prediction_result$assay) {
  assertthat::assert_that(set %in% c("test", "train"), msg = "'set=' must be 'test' or 'train'.")
  assertthat::assert_that(model == "all" | is.numeric(model) & model %in% seq_along(ppi_prediction_result$model.e), msg = "model must be 'all' or integer")
  assertthat::assert_that(training_set == "all" | is.numeric(training_set) & model %in% seq_along(ppi_prediction_result$ensembleSize), msg = "training set must be 'all' or integer")
  assertthat::assert_that(length(ppi_prediction_result$assay) <= 3, msg = "Plots can only be generated for predictions with a maximum of 3 training features.")
  assertthat::assert_that(is.numeric(x.nudge), msg = "x.nudge must be numerical.")
  assertthat::assert_that(is.character(assay) & all(assay %in% ppi_prediction_result$assay), msg = "assay must be character and used as training features in ppi_prediction_result.")
  assertthat::assert_that(type %in% c("2D", "3D") & length(ppi_prediction_result$assay) < 3, msg = "`type =` must be either 2D or 3D and can only be applied when trained on 3 features.")

  if (length(ppi_prediction_result$assay) == 1) {
    a <- assay[1]
    x.log.scale <- FALSE
    base::message("ML algorithm was trained on one feature. x-scale will be linear showing the number of interactions")
  } else if (length(ppi_prediction_result$assay) == 2) {
    a <- assay[1]
    b <- assay[2]
  } else if (length(ppi_prediction_result$assay) == 3) {
    a <- assay[1]
    b <- assay[2]
    c <- ppi_prediction_result$assay[!ppi_prediction_result$assay %in% c(a,b)]
    x.log.scale <- FALSE
    if(n > 40) {
      message("n was set to 40 to reduce calculation time")
      n <- 40
    }
  }
  if (set == "test" & !any(colnames(ppi_prediction_result$predDf) %in% 'reference')) {
    Df <- ppi_prediction_result$predDf %>% dplyr::rename(reference = complex)
  } else if (set == "train" & !any(colnames(ppi_prediction_result$predDf) %in% 'reference')) {
    Df <- ppi_prediction_result$predTrainDf %>%
      dplyr::rename(predMat = predTrainMat)
  } else if (set == 'test' & any(colnames(ppi_prediction_result$predDf) %in% 'reference')) {
    Df <- ppi_prediction_result$predDf
  } else if (set == 'train' & any(colnames(ppi_prediction_result$predDf) %in% 'reference')) {
    Df <- ppi_prediction_result$predTrainDf
  }
  if(is.numeric(training_set)) {
    training_set_interactions <- str_split(rownames(ppi_prediction_result$training.sets[[training_set]]), pattern = ";", simplify = TRUE)[,3:5] %>% as_tibble()
    colnames(training_set_interactions) <- c('interaction', 'sample', 'orientation')
    Df <- ppi_prediction_result$predTrainDf %>%
      dplyr::rename(predMat = predTrainMat) %>%
      inner_join(training_set_interactions)
  }
  if (ppi_prediction_result$model.type == "randomforest") {
    stop("Only available for ppi.prediction() results with 'svm' model")
  }
  lseq <- function(from, to, length.out) {
    exp(base::seq(log(from), log(to), length.out = length.out))
  }
  make.grid <- function(ppi_prediction_result) {
    if (length(ppi_prediction_result$assay) == 1) {
      x <- base::cbind(ppi_prediction_result$testMat,
                       row_number = 1:nrow(ppi_prediction_result$testMat)) %>%
        base::rbind(base::cbind(ppi_prediction_result$predTrainDf %>%
                                  dplyr::select(ppi_prediction_result$assay),
                                row_number = nrow(ppi_prediction_result$testMat)+1:nrow(ppi_prediction_result$predTrainDf)))
    } else if (length(ppi_prediction_result$assay) == 2) {
      x <- rbind(ppi_prediction_result$testMat[,c(a,b)],
                 ppi_prediction_result$predTrainDf[,c(a,b)])
    } else if (length(ppi_prediction_result$assay) == 3) {
      x <- rbind(ppi_prediction_result$testMat[,c(a,b,c)],
                 ppi_prediction_result$predTrainDf[,c(a,b,c)])
    }
    grange <- base::apply(x, 2, FUN = range, na.rm = TRUE)
    x1 <- base::seq(from = grange[1, 1], to = grange[2, 1], length = n)
    if (x.log.scale == TRUE) {
      if (any(grange < 0)) {
        grange <- base::apply(grange, 2, function(x) x + abs(min(x)) + x.nudge)
      }
      x2 <- lseq(from = grange[1, 2], to = grange[2, 2], length = n)
    } else if (x.log.scale == FALSE) {
      x2 <- base::seq(from = grange[1, 2], to = grange[2, 2], length = n)
    }
    if(length(ppi_prediction_result$assay) == 1 | length(ppi_prediction_result$assay) == 2) {
      grid <- base::expand.grid(assay1 = x1, assay2 = x2)
    }
    if(length(ppi_prediction_result$assay) == 3) {
      x3 <- base::seq(from = grange[1, 3], to = grange[2, 3], length = n)
      grid <- base::expand.grid(assay1 = x1, assay2 = x2, assay3 = x3)
    }
    return(grid)
  }
  xgrid = make.grid(ppi_prediction_result)
  if (length(ppi_prediction_result$assay) == 1) {
    x.log.scale <- FALSE
  }
  if (length(ppi_prediction_result$assay) == 1) {
    base::colnames(xgrid) <- c(as.character(rlang::sym(ppi_prediction_result$assay[1])), "id")
  } else if (length(ppi_prediction_result$assay) == 2) {
    base::colnames(xgrid) <- c(as.character(a), as.character(b))
  } else if (length(ppi_prediction_result$assay) == 3) {
    base::message("ML algorithm was trained on more 3 features. Only 2 features are plotted. Please specify under assay")
    base::colnames(xgrid) <- c(as.character(a), as.character(b), as.character(c))
  }
  #loop over all models and calculate mean ygrid
  if (model == "all") {
    for (i in base::seq(1, length(ppi_prediction_result$model.e))) {
      if (i == 1) {
        ygrid <- base::data.frame(base::matrix(nrow = base::nrow(xgrid),
                                               ncol = 0))
        ygrid.prob <- base::data.frame(base::matrix(nrow = base::nrow(xgrid),
                                                    ncol = 0))
      }
      ygrid = base::cbind(ygrid, data.frame(stats::predict(ppi_prediction_result$model.e[[i]], xgrid)))
      ygrid.prob = base::cbind(ygrid.prob, attr(stats::predict(ppi_prediction_result$model.e[[i]], xgrid, probability = TRUE), "probabilities")[, 1])

      if(i == length(ppi_prediction_result$model.e)) {
        colnames(ygrid) <- base::seq(1:length(ppi_prediction_result$model.e))
        colnames(ygrid.prob) <- base::seq(1:length(ppi_prediction_result$model.e))
        ygrid <- base::rowMeans(varhandle::unfactor(ygrid))
        ygrid.prob <- base::rowMeans(ygrid.prob)
      }
    }
  }

  #calculate ygrid for a specific model
  if (model %in% seq_along(ppi_prediction_result$model.e)) {
    ygrid <- base::data.frame(base::matrix(nrow = base::nrow(xgrid),ncol = 0))
    ygrid.prob <- base::data.frame(base::matrix(nrow = base::nrow(xgrid),ncol = 0))
    ygrid = data.frame(ygrid = stats::predict(ppi_prediction_result$model.e[[model]], xgrid))
    ygrid.prob = base::cbind(ygrid.prob, attr(stats::predict(ppi_prediction_result$model.e[[model]],
                                                             xgrid, probability = TRUE), "probabilities")[, 1])
  }
  if (x.log.scale == TRUE) {
    xgrid[,2] <- log(xgrid[,2] + abs(min(xgrid[,2])) + x.nudge)

    min.x <- min(ppi_prediction_result$predDf %>% dplyr::select(b),
                 ppi_prediction_result$predTrainDf %>% dplyr::select(b),
                 na.rm = TRUE)
    if (min.x < 0) {
      ppi_prediction_result$predDf[,b] <- ppi_prediction_result$predDf[,b] + abs(min.x) + x.nudge
      ppi_prediction_result$predTrainDf[,b] <- ppi_prediction_result$predTrainDf[,b] + abs(min.x) + x.nudge
    }
  }
  if (length(ppi_prediction_result$assay) == 1 & model == "all" & type == "2D") {
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = xgrid %>%
                                                   base::cbind(ygrid) %>% base::cbind(ygrid.prob) %>%
                                                   base::as.data.frame(), mapping = ggplot2::aes(x = id,
                                                                                                 y = !!rlang::sym(ppi_prediction_result$assay[1]),
                                                                                                 fill = ygrid.prob), shape = 22, size = 1, stroke = 0) +
      ggplot2::scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), limits = c(0, 1)) +
      ggplot2::labs(fill = "probability", size = "probability") +
      ggnewscale::new_scale("fill") + ggplot2::geom_point(data = Df %>%
                                                            dplyr::mutate(id = dplyr::row_number()), mapping = ggplot2::aes(x = id,
                                                                                                                            y = !!rlang::sym(ppi_prediction_result$assay[1]),
                                                                                                                            size = predMat, shape = reference, col = reference),
                                                          fill = "white", stroke = 0.5, alpha = 1) +
      ggplot2::scale_size_binned(range = c(0.01, 3.5), n.breaks = 4, limits = c(0.75, 1)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                hjust = 1), legend.position = "right")
  } else if (length(ppi_prediction_result$assay) >= 2 & model == "all" & type == "2D") {
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = xgrid %>%
                                                   base::cbind(ygrid) %>% base::cbind(ygrid.prob) %>%
                                                   base::as.data.frame(), mapping = ggplot2::aes(x = !!rlang::sym(b),
                                                                                                 y = !!rlang::sym(a),
                                                                                                 fill = ygrid.prob), shape = 22, size = 1, stroke = 0) +
      ggplot2::scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), limits = c(0, 1)) +
      ggplot2::labs(fill = "probability", size = "probability") +
      ggnewscale::new_scale("fill") + ggplot2::geom_point(data = Df %>%
                                                            dplyr::filter(!is.na(!!rlang::sym(b))) %>%
                                                            dplyr::rowwise() %>% dplyr::mutate(`:=`((!!rlang::sym(b)),
                                                                                                    base::ifelse(x.log.scale == TRUE, log(!!rlang::sym(b)),
                                                                                                                 (!!rlang::sym(b))))),
                                                          mapping = ggplot2::aes(x = !!rlang::sym(b),
                                                                                 y = !!rlang::sym(a),
                                                                                 size = predMat, shape = reference, col = reference),
                                                          fill = "white", stroke = 0.5, alpha = 1) +
      ggplot2::scale_size_binned(range = c(0.01, 3.5), n.breaks = 4, limits = c(0.75, 1)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                hjust = 1), legend.position = "right")
    if (set == "train") {
      p <- p +
        ggplot2::scale_shape_manual(values = c(RRS = 23, PRS = 21)) +
        ggplot2::scale_color_manual(values = c(RRS = "#EC008C", PRS = "#00AEEF"))
    }
    if (set == "test") {
      shapes <- c(23, 21, 22, 24, 25, 0, 1, 2, 5, 6, 15:20, 3, 4, 7:14)
      p <- p +
        ggplot2::scale_shape_manual(values = shapes) +
        viridis::scale_color_viridis(discrete = TRUE, option = "D")
    }
  }
  if (length(ppi_prediction_result$assay) == 1 & model %in% seq_along(ppi_prediction_result$model.e) & type == "2D") {
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = xgrid %>%
                                                   base::cbind(ygrid) %>% base::cbind(ygrid.prob) %>%
                                                   base::as.data.frame(), mapping = ggplot2::aes(x = id,
                                                                                                 y = !!rlang::sym(a),
                                                                                                 fill = ygrid), shape = 22, size = 1, stroke = 0) +
      ggplot2::scale_fill_manual(values = c("#FBF49C", "#48C2C5")) +
      ggplot2::labs(fill = "probability", size = "probability") +
      ggnewscale::new_scale("fill") + ggplot2::geom_point(data = Df %>%
                                                            dplyr::mutate(id = dplyr::row_number()), mapping = ggplot2::aes(x = id,
                                                                                                                            y = !!rlang::sym(a),
                                                                                                                            size = predMat, shape = reference, col = reference),
                                                          fill = "white", stroke = 0.5, alpha = 1) +
      ggplot2::scale_size_binned(range = c(0.01, 3.5), n.breaks = 4, limits = c(0.75, 1)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                hjust = 1), legend.position = "right")
  } else if (length(ppi_prediction_result$assay) >= 2 & model %in% seq_along(ppi_prediction_result$model.e) & type == "2D") {
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = xgrid %>%
                                                   base::cbind(ygrid) %>% base::cbind(ygrid.prob) %>%
                                                   base::as.data.frame(), mapping = ggplot2::aes(x = !!rlang::sym(b),
                                                                                                 y = !!rlang::sym(a),
                                                                                                 fill = ygrid), shape = 22, size = 1, stroke = 0) +
      ggplot2::scale_fill_manual(values = c("#FBF49C", "#48C2C5")) +
      ggplot2::labs(fill = "probability", size = "probability") +
      ggnewscale::new_scale("fill") + ggplot2::geom_point(data = Df %>%
                                                            dplyr::filter(!is.na(!!rlang::sym(b))) %>%
                                                            dplyr::rowwise() %>% dplyr::mutate(`:=`((!!rlang::sym(b)),
                                                                                                    base::ifelse(x.log.scale == TRUE, log(!!rlang::sym(b)),
                                                                                                                 (!!rlang::sym(b))))),
                                                          mapping = ggplot2::aes(x = !!rlang::sym(b),
                                                                                 y = !!rlang::sym(a),
                                                                                 size = predMat, shape = reference, col = reference),
                                                          fill = "white", stroke = 0.5, alpha = 1) +
      ggplot2::scale_size_binned(range = c(0.01, 3.5), n.breaks = 4, limits = c(0.75, 1)) +
      ggpubr::theme_pubr() +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                hjust = 1), legend.position = "right")
    if (set == "train") {
      p <- p +
        ggplot2::scale_shape_manual(values = c(RRS = 23, PRS = 21)) +
        ggplot2::scale_color_manual(values = c(RRS = "#EC008C", PRS = "#00AEEF"))
    }
    if (set == "test") {
      shapes <- c(23, 21, 22, 24, 25, 0, 1, 2, 5, 6, 15:20, 3, 4, 7:14)
      p <- p +
        ggplot2::scale_shape_manual(values = shapes) +
        viridis::scale_color_viridis(discrete = TRUE, option = "D")
    }
  }

  if (!is.null(xlim)) {
    p <- p + ggplot2::scale_x_continuous(limits = xlim)
  }
  if (!is.null(ylim)) {
    p <- p + ggplot2::scale_y_continuous(limits = ylim)
  }

  if(length(ppi_prediction_result$assay) == 3 & model == "all" & type == "3D") {
    grid <- xgrid %>%
      base::cbind(ygrid) %>% base::cbind(ygrid.prob) %>%
      base::as.data.frame()

    p <- plotly::plot_ly(x = grid %>% dplyr::pull(rlang::sym(b)), y = grid %>% dplyr::pull(rlang::sym(a)), z = grid %>% dplyr::pull(rlang::sym(c)),
                          type = "isosurface", value = grid$ygrid.prob, colors = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), alpha = 0.5) %>%
      layout(scene = list(xaxis = list(title = b),
                          yaxis = list(title = a),
                          zaxis = list(title = c)))


    p <- plotly::plot_ly(x = Df %>% dplyr::pull(rlang::sym(b)), y = Df %>% dplyr::pull(rlang::sym(a)), z = Df %>% dplyr::pull(rlang::sym(c)), color = Df$reference, colors = c(RRS = "#EC008C", PRS = "#00AEEF"),
                    type = "scatter3d", mode = "markers", size = Df$predMat, text = paste("prob:", round(Df$predMat*100,1), "%"),
                    symbol = 21) %>%
      layout(scene = list(xaxis = list(title = b),
                          yaxis = list(title = a),
                          zaxis = list(title = c)))
  }
  return(p)
}

