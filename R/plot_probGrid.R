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
#' @param n: grid size
#' @param x.log.scale: boolean to log-scale x-axis values
#' @param xlim: numeric values, specifying the left and right limit of the scale
#' @param ylim: numeric values, specifying the bottom and top limit of the scale
#' @param set: PPI set to generate the plot for: "test" or "train"
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' data("example_ppi_prediction")
#' probGrid.plot(example_ppi_prediction)
probGrid.plot <- function(ppi_prediction_result, n=100, x.log.scale = TRUE, xlim = c(0,NA), ylim = NULL, set="train") {
  assertthat::assert_that(set %in% c("test", "train"), msg = "'set=' must be 'test' or 'train'.")
  if(set == "test") {
    Df <- ppi_prediction_result$predDf %>% dplyr::rename(reference = complex)
  } else if(set == "train") {
    Df <- ppi_prediction_result$predTrainDf %>%
      dplyr::rename(predMat = predTrainMat)
  }
  if (ppi_prediction_result$model.type == "randomforest") {
    stop("Only available for ppi.prediction() results with 'svm' model")
  }
  lseq <- function(from = 1, to = 1e+05, length.out = 6) {
    exp(base::seq(log(from), log(to), length.out = length.out))
  }
  make.grid <- function(ppi_prediction_result) {
    if (length(ppi_prediction_result$assay) > 2) {
      base::message("ML algorithm was trained on more than 2 features. The grid will be generated for the first two features as specified under assay.")
    }
    if (length(ppi_prediction_result$assay) == 1) {
      x <- base::cbind(ppi_prediction_result$testMat,
                       row_number = 1:nrow(ppi_prediction_result$testMat)) %>%
        base::rbind(base::cbind(ppi_prediction_result$trainMat,
                                row_number = nrow(ppi_prediction_result$testMat)+1:nrow(ppi_prediction_result$trainMat)))
      x.log.scale <- FALSE
      base::message("ML algorithm was trained on one feature. x-scale will be linear showing the number of interactions")
    } else if (length(ppi_prediction_result$assay) >= 2) {
      x <- rbind(ppi_prediction_result$testMat, ppi_prediction_result$trainMat)
    }
    grange <- base::apply(x, 2, range)
    x1 <- base::seq(from = grange[1, 1], to = grange[2, 1], length = n)
    if (x.log.scale == TRUE) {
      if(any(grange < 0)) {
        grange <- base::apply(grange, 2, function(x) x + abs(min(x)) + 1e-4)
      }
      x2 <- lseq(from = grange[1, 2], to = grange[2, 2],
                 length = n)
    }
    else if (x.log.scale == FALSE) {
      x2 <- base::seq(from = grange[1, 2], to = grange[2,
                                                       2], length = n)
    }
    base::expand.grid(assay1 = x1, assay2 = x2)
  }
  xgrid = make.grid(ppi_prediction_result)
  if (length(ppi_prediction_result$assay) == 1) {
    x.log.scale <- FALSE
  }
  if (x.log.scale == TRUE) {
    xgrid <- cbind(xgrid[, 1], log(xgrid[, 2]))
  }
  if (length(ppi_prediction_result$assay) == 1) {
    base::colnames(xgrid) <- c(as.character(rlang::sym(ppi_prediction_result$assay[1])),
                               "id")
  } else if (length(ppi_prediction_result$assay) >= 1) {
    base::colnames(xgrid) <- c(as.character(rlang::sym(ppi_prediction_result$assay[1])),
                               as.character(rlang::sym(ppi_prediction_result$assay[2])))
  }
  for (i in base::seq(1, length(ppi_prediction_result$model.e))) {
    if (i == 1) {
      ygrid <- base::data.frame(base::matrix(nrow = base::nrow(xgrid),
                                             ncol = 0))
      ygrid.prob <- base::data.frame(base::matrix(nrow = base::nrow(xgrid),
                                                  ncol = 0))
    }
    ygrid = base::cbind(ygrid, data.frame(stats::predict(ppi_prediction_result$model.e[[i]],
                                                         xgrid)))
    ygrid.prob = base::cbind(ygrid.prob, attr(stats::predict(ppi_prediction_result$model.e[[i]],
                                                             xgrid, probability = TRUE), "probabilities")[, 1])
    if (i == length(ppi_prediction_result$model.e)) {
      colnames(ygrid) <- base::seq(1:length(ppi_prediction_result$model.e))
      colnames(ygrid.prob) <- base::seq(1:length(ppi_prediction_result$model.e))
      ygrid <- base::rowMeans(varhandle::unfactor(ygrid))
      ygrid.prob <- base::rowMeans(ygrid.prob)
    }
  }
  if (x.log.scale == TRUE) {
    min.x <- min(ppi_prediction_result$predDf %>% dplyr::select(!!rlang::sym(ppi_prediction_result$assay[2])),
                 ppi_prediction_result$predTrainDf %>% dplyr::select(!!rlang::sym(ppi_prediction_result$assay[2])),
                 na.rm = TRUE)
    if(min.x < 0) {
      ppi_prediction_result$predDf[,ppi_prediction_result$assay[2]] <- ppi_prediction_result$predDf[,ppi_prediction_result$assay[2]] + abs(min.x) + 1e-4
      ppi_prediction_result$predTrainDf[,ppi_prediction_result$assay[2]] <- ppi_prediction_result$predTrainDf[,ppi_prediction_result$assay[2]] + abs(min.x) + 1e-4
    }
  }
  if (length(ppi_prediction_result$assay) == 1) {
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = xgrid %>%
                                                   base::cbind(ygrid) %>% base::cbind(ygrid.prob) %>%
                                                   base::as.data.frame(), mapping = ggplot2::aes(x = id,
                                                                                                 y = !!rlang::sym(ppi_prediction_result$assay[1]),
                                                                                                 fill = ygrid.prob), shape = 22, size = 1, stroke = 0) +
      ggplot2::scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), limits = c(0, 1)) + ggplot2::labs(fill = "probability", size = "probability") +
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
  } else if (length(ppi_prediction_result$assay) >= 2) {
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = xgrid %>%
                                                   base::cbind(ygrid) %>% base::cbind(ygrid.prob) %>%
                                                   base::as.data.frame(), mapping = ggplot2::aes(x = !!rlang::sym(ppi_prediction_result$assay[2]),
                                                                                                 y = !!rlang::sym(ppi_prediction_result$assay[1]),
                                                                                                 fill = ygrid.prob), shape = 22, size = 1, stroke = 0) +
      ggplot2::scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), limits = c(0, 1)) +
      ggplot2::labs(fill = "probability", size = "probability") +
      ggnewscale::new_scale("fill") + ggplot2::geom_point(data = Df %>%
                                                            dplyr::filter(!is.na(!!rlang::sym(ppi_prediction_result$assay[2]))) %>%
                                                            dplyr::rowwise() %>% dplyr::mutate(`:=`((!!rlang::sym(ppi_prediction_result$assay[2])),
                                                                                                    base::ifelse(x.log.scale == TRUE, log(!!rlang::sym(ppi_prediction_result$assay[2])),
                                                                                                                 (!!rlang::sym(ppi_prediction_result$assay[2]))))),
                                                          mapping = ggplot2::aes(x = !!rlang::sym(ppi_prediction_result$assay[2]),
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
    if(set == "train") {
      p <- p +
        ggplot2::scale_shape_manual(values = c(RRS = 23, PRS = 21)) +
        ggplot2::scale_color_manual(values = c(RRS = "#EC008C", PRS = "#00AEEF"))
    }
    if(set == "test") {
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
  return(p)
}

