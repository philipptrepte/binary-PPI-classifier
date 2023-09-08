#' Generate a logarithmic sequence
#'
#' @param from
#' @param to
#' @param length.out
#'
#' @return
#' @export
#'
#' @examples
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

#' Generate a grid to plot for BRET predictions
#'
#' @param x: takes as input the from result of ppi.prediction()$testMat
#' @param n: size of the grid
#'
#' @return
#' @export
#'
#' @examples
make.grid.BRET <- function(x,n=100) {
  grange=apply(x,2,range)
  x1=seq(from=grange[1,1],to=grange[2,1],length=n)
  x2=lseq(from=grange[1,2],to=grange[2,2],length=n)
  expand.grid(mean_cBRET = x1, mean_mCit = x2)
}

#' Plot a bar diagram of the recovery rates at 50%, 75% and 95% probability
#'
#' @param ppi_prediction_result
#'
#' @return
#' @export
#'
#' @examples
plot.recovery <- function(ppi_prediction_result) {
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

#' Plot a scatter plot of the probability distribution against the primary assay feature
#'
#' @param ppi_prediction_result
#'
#' @return
#' @export
#'
#' @examples
plot.probDis <- function(ppi_prediction_result) {
  ppi_prediction_result$predDf %>%
    ggplot(aes(x = mean_cBRET, y = predMat, fill = complex)) +
    geom_point(shape = 21, alpha = 0.7, size = 1.5) +
    theme_pubr() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#6CA6C1", "#D930C5")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
    labs(x = "cBRET ratio", y = "probability") +
    theme_pubr() +
    ggplot2::theme(text = element_text(family = "Avenir"),
                   plot.title = element_text(size = 12),
                   plot.subtitle = element_text(size = 8),
                   axis.title = element_text(family = "Avenir Medium"),
                   legend.position = "right")
}

#' Plot a probability grid from the mean probabilities from the 'ensembleSize' number of models for the LuTHy-BRET assay
#'
#' @param ppi_prediction_result
#'
#' @return
#' @export
#'
#' @examples
plot.probGrid.BRET <- function(ppi_prediction_result) {
  xgrid.BRET = make.grid.BRET(data.frame(ppi_prediction_result$testMat))
  for(i in seq(1, length(ppi_prediction_result$model.e))) {
    if(i == 1) {
      ygrid.BRET <- data.frame(matrix(nrow = nrow(xgrid.BRET), ncol = 0))
      ygrid.BRET.prob <- data.frame(matrix(nrow = nrow(xgrid.BRET), ncol = 0))
    }
    ygrid.BRET = cbind(ygrid.BRET, data.frame(predict(ppi_prediction_result$model.e[[i]], xgrid.BRET)))
    ygrid.BRET.prob = cbind(ygrid.BRET.prob, attr(predict(ppi_prediction_result$model.e[[i]], xgrid.BRET, probability = TRUE), "probabilities")[,1])
    if(i == length(ppi_prediction_result$model.e)) {
      colnames(ygrid.BRET) <- seq(1:length(ppi_prediction_result$model.e))
      colnames(ygrid.BRET.prob) <- seq(1:length(ppi_prediction_result$model.e))
      ygrid.BRET <- rowMeans(varhandle::unfactor(ygrid.BRET))
      ygrid.BRET.prob <- rowMeans(ygrid.BRET.prob)
    }
  }

  ggplot() +
    geom_point(data = xgrid.BRET %>% cbind(ygrid.BRET) %>% cbind(ygrid.BRET.prob),
               mapping = aes(x = log10(mean_mCit), y = mean_cBRET, fill = ygrid.BRET.prob),
               shape = 22, size = 1, stroke = 0) +
    scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), limits = c(0,1)) +
    labs(fill = "probability", size = "probability") +
    new_scale("fill") +
    geom_point(data = ppi_prediction_result$predDf %>%
                 filter(!is.na(mean_cBRET)),
               mapping = aes(x = log10(mean_mCit), y = mean_cBRET, size = predMat, shape = complex, col = complex),
               fill = "white", stroke = 0.5, alpha = 1) +
    scale_size_binned(range = c(0.01, 3.5), n.breaks = 4, limits = c(0.75, 1)) +
    scale_shape_manual(values = c("RRS" = 23, "PRS" = 21)) +
    scale_color_manual(values = c("RRS" = "#EC008C", "PRS" = "#00AEEF")) +
    scale_y_continuous(limits = c(-0.1, 0.8), breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8)) +
    theme_pubr() +
    ggplot2::theme(text = element_text(family = "Avenir"),
                   plot.title = element_text(size = 12),
                   plot.subtitle = element_text(size = 8),
                   axis.title = element_text(family = "Avenir Medium"),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                   legend.position = "right")
}

#' Plot a probability grid from the mean probabilities from the 'ensembleSize' number of models for the LuTHy-LuC assay
#'
#' @param ppi_prediction_result
#'
#' @return
#' @export
#'
#' @examples
plot.probGrid.LuC <- function(ppi_prediction_result) {
  make.grid.LuC <- function(x,n=100) {
    grange=apply(x,2,range)
    x1=seq(from=grange[1,1],to=grange[2,1],length=n)
    x2=lseq(from=grange[1,2],to=grange[2,2],length=n)
    expand.grid(mean_cLuC = x1, mean_LumiOUT = x2)
  }
  xgrid.LuC = make.grid.LuC(data.frame(ppi_prediction_result$testMat))

  for(i in seq(1, length(ppi_prediction_result$model.e))) {
    if(i == 1) {
      ygrid.LuC <- data.frame(matrix(nrow = nrow(xgrid.LuC), ncol = 0))
      ygrid.LuC.prob <- data.frame(matrix(nrow = nrow(xgrid.LuC), ncol = 0))
    }
    ygrid.LuC = cbind(ygrid.LuC, data.frame(predict(ppi_prediction_result$model.e[[i]], xgrid.LuC)))
    ygrid.LuC.prob = cbind(ygrid.LuC.prob, attr(predict(ppi_prediction_result$model.e[[i]], xgrid.LuC, probability = TRUE), "probabilities")[,1])
    if(i == length(ppi_prediction_result$model.e)) {
      colnames(ygrid.LuC) <- seq(1:length(ppi_prediction_result$model.e))
      colnames(ygrid.LuC.prob) <- seq(1:length(ppi_prediction_result$model.e))
      ygrid.LuC <- rowMeans(unfactor(ygrid.LuC))
      ygrid.LuC.prob <- rowMeans(ygrid.LuC.prob)
    }
  }

  ggplot() +
    geom_point(data = xgrid.LuC %>% cbind(ygrid.LuC) %>% cbind(ygrid.LuC.prob),
               mapping = aes(x = log10(mean_LumiOUT), y = mean_cLuC, fill = ygrid.LuC.prob),
               shape = 22, size = 1, stroke = 0) +
    scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), limits = c(0,1)) +
    labs(fill = "probability", size = "probability") +

    new_scale("fill") +
    geom_point(data = ppi_prediction_result$predDf %>%
                 filter(!is.na(mean_cLuC)),
               mapping = aes(x = log10(mean_LumiOUT), y = mean_cLuC, size = predMat, shape = complex, col = complex),
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

#' Plot ROC curve of the predicted probabilities against the assay features the ML was trained on
#'
#' @param ppi_prediction_result
#'
#' @return
#' @export
#'
#' @examples
plot.roc <- function(ppi_prediction_result) {
  n.assays <- length(ppi_prediction_result$assay)
  prob.aggregated <- ppi_prediction_result$predDf %>%
    dplyr::mutate(D = ifelse(str_detect(complex, "inter-complex|RRS"), "0", "1")) %>%
    group_by(interaction, D) %>%
    dplyr::summarise(M = max(predMat, na.rm = TRUE))

  assay.aggregated <- list()
  for(a in seq(n.assays)) {
    assay.aggregated[[a]] <- ppi_prediction_result$predDf %>%
      dplyr::mutate(D = ifelse(str_detect(complex, "inter-complex|RRS"), "0", "1")) %>%
      group_by(interaction, D) %>%
      dplyr::summarise(M = max(!!(rlang::sym(ppi_prediction_result$assay[a])), na.rm = TRUE))
  }

  prob.aggregated.roc <- plotROC::calculate_roc(M = prob.aggregated$M, D = prob.aggregated$D) %>% #, cutoffs = seq(from = 0, to = 1, by = 1/(nrow(referenceSet.BRET.prob.aggregated)-1))) %>%
    group_by(TPF, FPF) %>%
    dplyr::summarise(c = max(c, na.rm = TRUE)) %>%
    dplyr::mutate(data = "probability")

  assay.aggregated.roc <- list()
  for(r in seq(n.assays)) {
    assay.aggregated.roc[[r]] <- plotROC::calculate_roc(M = assay.aggregated[[r]]$M, D = assay.aggregated[[r]]$D) %>%
      group_by(TPF, FPF) %>%
      dplyr::summarise(c = max(c, na.rm = TRUE)) %>%
      dplyr::mutate(data = ppi_prediction_result$assay[r])
  }

  annotations <- list()
  for(t in seq(n.assays)) {
    if(t == 1) {
      text <- paste0("AUC = ", round(DescTools::AUC(x = prob.aggregated.roc$FPF, y = prob.aggregated.roc$TPF), 2))
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
    rbind(purrr::list_rbind(assay.aggregated.roc)) %>%
    dplyr::mutate(data = factor(data, levels = c("probability", ppi_prediction_result$assay))) %>%
    ggplot(aes(x = FPF, y = TPF, color = data)) +
    geom_line() +
    scale_color_viridis(discrete = TRUE) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.1)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.1)) +
    annotations +
    theme_pubr() +
    ggplot2::theme(text = element_text(family = "Avenir"),
                   plot.title = element_text(size = 12),
                   plot.subtitle = element_text(size = 8),
                   axis.title = element_text(family = "Avenir Medium"),
                   legend.position = "right") +
    labs(color = "")
}
