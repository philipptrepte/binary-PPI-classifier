#' Plot a probability grid from the mean probabilities from the 'ensembleSize' number of models
#'
#' @param ppi_prediction_result
#' @param n
#' @param x.log.scale
#' @param xlim
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
plot.probGrid <- function(ppi_prediction_result, n=100, x.log.scale = TRUE, xlim = NULL, ylim = NULL) {
  if(ppi_prediction_result$model.type == "randomforest") {
    stop("Only available for ppi.prediction() results with 'svm' model")
  }
  make.grid <- function(ppi_prediction_result) {
    if(length(ppi_prediction_result$assay) > 2) {
      message("ML algorithm was trained on more than 2 features. The grid will be generated for the first two features as specified under assay.")
    }

    if(length(ppi_prediction_result$assay) == 1) {
      x <- cbind(ppi_prediction_result$testMat, row_number = 1:length(ppi_prediction_result$testMat))
      x.log.scale <- FALSE
      message("ML algorithm was trained on one feature. x-scale will be linear showing the number of interactions")
    } else if(length(ppi_prediction_result$assay) >= 2) {
      x <- ppi_prediction_result$testMat
    }

    grange=apply(x,2,range)

    x1=seq(from=grange[1,1],to=grange[2,1],length=n)
    if(x.log.scale == TRUE) {
      x2=lseq(from=grange[1,2],to=grange[2,2],length=n)
    } else if(x.log.scale == FALSE) {
      x2=seq(from=grange[1,2],to=grange[2,2],length=n)
    }
    expand.grid(assay1 = x1, assay2 = x2)
  }

  xgrid = make.grid(ppi_prediction_result)
  if(length(ppi_prediction_result$assay) == 1) {
    x.log.scale <- FALSE
  }
  if(x.log.scale == TRUE) {
    xgrid=cbind(xgrid[,1], log10(xgrid[,2]))
  }
  if(length(ppi_prediction_result$assay) == 1) {
    colnames(xgrid) <- c((rlang::sym(ppi_prediction_result$assay[1])), "id")
  } else if(length(ppi_prediction_result$assay) >= 1){
    colnames(xgrid) <- c((rlang::sym(ppi_prediction_result$assay[1])), (rlang::sym(ppi_prediction_result$assay[2])))
  }

  for(i in seq(1, length(ppi_prediction_result$model.e))) {
    if(i == 1) {
      ygrid <- data.frame(matrix(nrow = nrow(xgrid), ncol = 0))
      ygrid.prob <- data.frame(matrix(nrow = nrow(xgrid), ncol = 0))
    }

    ygrid = cbind(ygrid, data.frame(predict(ppi_prediction_result$model.e[[i]], xgrid)))
    ygrid.prob = cbind(ygrid.prob, attr(predict(ppi_prediction_result$model.e[[i]], xgrid, probability = TRUE), "probabilities")[,1])

    if(i == length(ppi_prediction_result$model.e)) {
      colnames(ygrid) <- seq(1:length(ppi_prediction_result$model.e))
      colnames(ygrid.prob) <- seq(1:length(ppi_prediction_result$model.e))
      ygrid <- rowMeans(varhandle::unfactor(ygrid))
      ygrid.prob <- rowMeans(ygrid.prob)
    }
  }

  if(length(ppi_prediction_result$assay == 1)) {
    p <- ggplot() +
      geom_point(data = xgrid %>% cbind(ygrid) %>% cbind(ygrid.prob) %>% as.data.frame(),
                 mapping = aes(x = id,
                               y = !!rlang::sym(ppi_prediction_result$assay[1]), fill = ygrid.prob),
                 shape = 22, size = 1, stroke = 0) +
      scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), limits = c(0,1)) +
      labs(fill = "probability", size = "probability") +
      new_scale("fill") +
      geom_point(data = ppi_prediction_result$predDf %>%
                   dplyr::mutate(id = row_number()),
                 mapping = aes(x = id,
                               y = !!rlang::sym(ppi_prediction_result$assay[1]),
                               size = predMat, shape = complex, col = complex),
                 fill = "white", stroke = 0.5, alpha = 1) +
      scale_size_binned(range = c(0.01, 3.5), n.breaks = 4, limits = c(0.75, 1)) +
      scale_shape_manual(values = c("RRS" = 23, "PRS" = 21)) +
      scale_color_manual(values = c("RRS" = "#EC008C", "PRS" = "#00AEEF")) +
      theme_pubr() +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12),
                     plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                     legend.position = "right")
  } else if(length(ppi_prediction_result$assay) >= 2) {
    p <- ggplot() +
      geom_point(data = xgrid %>% cbind(ygrid) %>% cbind(ygrid.prob) %>% as.data.frame(),
                 mapping = aes(x = !!rlang::sym(ppi_prediction_result$assay[2]),
                               y = !!rlang::sym(ppi_prediction_result$assay[1]), fill = ygrid.prob),
                 shape = 22, size = 1, stroke = 0) +
      scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), limits = c(0,1)) +
      labs(fill = "probability", size = "probability") +
      new_scale("fill") +
      geom_point(data = ppi_prediction_result$predDf %>%
                   filter(!is.na(!!rlang::sym(ppi_prediction_result$assay[2]))) %>%
                   rowwise() %>%
                   dplyr::mutate(!!rlang::sym(ppi_prediction_result$assay[2]) :=
                                   ifelse(x.log.scale == TRUE,
                                          log10(!!rlang::sym(ppi_prediction_result$assay[2])),
                                          !!rlang::sym(ppi_prediction_result$assay[2]))),
                 mapping = aes(x = !!rlang::sym(ppi_prediction_result$assay[2]),
                               y = !!rlang::sym(ppi_prediction_result$assay[1]),
                               size = predMat, shape = complex, col = complex),
                 fill = "white", stroke = 0.5, alpha = 1) +
      scale_size_binned(range = c(0.01, 3.5), n.breaks = 4, limits = c(0.75, 1)) +
      scale_shape_manual(values = c("RRS" = 23, "PRS" = 21)) +
      scale_color_manual(values = c("RRS" = "#EC008C", "PRS" = "#00AEEF")) +
      theme_pubr() +
      ggplot2::theme(text = element_text(family = "Avenir"),
                     plot.title = element_text(size = 12),
                     plot.subtitle = element_text(size = 8),
                     axis.title = element_text(family = "Avenir Medium"),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                     legend.position = "right")
  }
  if(!is.null(xlim)) {
    p <- p +
      scale_x_continuous(limits = xlim)
  }
  if(!is.null(ylim)) {
    p <- p +
      scale_y_continuous(limits = ylim)
  }
  p
}


