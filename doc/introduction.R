## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.width = 6, 
  fig.height = 4,
  out.width = "100%"
)

## ----requiredPackages, echo = FALSE, warning=FALSE, results="hide"------------
suppressPackageStartupMessages({
library(tidyverse)
library(ggpubr)
library(viridis)
library(binaryPPIclassifier)
library(varhandle)
library(ggnewscale)
library(cowplot)
})


## ----loadData, echo=TRUE, message=FALSE, warning=FALSE------------------------
data("luthy_reference_sets")

## ----inputTable, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE----
knitr::kable(luthy_reference_sets %>% select(-plate, -PPI) %>% filter(interaction == "BAD + BCL2L1"),
             caption = 'Input table showing results for the interaction BAD + BCL2L1 from `luthy_reference_sets`.')

## ----luthyBretExample, echo=FALSE, message=FALSE, warning=FALSE, fig.cap="LuTHy-BRET results for one positive (BAD + BCL2L1) and one random (ARSA + DBN1) PPI"----
luthy_reference_sets %>% 
  filter(interaction == "BAD + BCL2L1" & data == "mean_cBRET" |
           interaction == "ARSA + DBN1" & data == "mean_cBRET" |
           interaction == "BAD + BCL2L1" & data == "mean_mCit" |
           interaction == "ARSA + DBN1" & data == "mean_mCit") %>%
  dplyr::mutate(interaction = paste0(complex, ": ", interaction)) %>%
  ggplot(aes(x = orientation, y = score, fill = score)) +
  geom_bar(stat = "identity", fill = "#6CA6C1") +
  theme_pubr() +
  facet_grid(data ~ interaction, scales = "free_y") +
  ggplot2::theme(text = element_text(family = "Avenir"),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8),
          axis.title = element_text(family = "Avenir Medium"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "LuTHy-BRET measurements", x = "tagging configuration")

## ----loadExamplePPIprediction, echo=TRUE, message=FALSE, warning=FALSE--------
data("example_ppi_prediction")

## ----MLresultPlot, fig.cap="Training results as (A) probability grid; (B) probability distribtuoin; (C) ROC curve plotting the training features and the probability; (D) recovery rates at 50%, 75% and 95% probability", fig.height=7.5, fig.width=10, message=FALSE, warning=FALSE, echo=TRUE----

#cowplot::plot_grid(
#    probGrid.plot(example_ppi_prediction, ylim = c(-0.1, 0.8)),
#    probDis.plot(example_ppi_prediction),
#    roc.plot(example_ppi_prediction),
#    recovery.plot(example_ppi_prediction),
#  ncol = 2, align = "hv", axis = "tlrb", labels = "AUTO"
#)


## ----exampleConfMat, echo=TRUE, message=FALSE, warning=FALSE------------------
example_confusion_matrix <- confusion.matrix(example_ppi_prediction)

## ----exampleConfMatDf, echo=FALSE, message=FALSE, warning=FALSE---------------
knitr::kable(example_confusion_matrix$confusionMatrixDf[,1:2], caption = "Tabulated confusion matrix.")

## ----exampleConfMatLis, echo=FALSE, message=FALSE, warning=FALSE--------------
example_confusion_matrix$confusionMatrixList$`95%`

## ----predExamplePPI, echo=FALSE, message=FALSE, warning=FALSE, fig.height=4, fig.width=5, fig.cap="(A) Predicted probabilities to be true-positive interactions and (B) mean cBRET values for the PRS-v2 interaction BAD + BCL2L1 and the RRS-v2 interaction ARSA + DBN1"----

cowplot::plot_grid(
  example_ppi_prediction$predDf %>%
  filter(interaction == "BAD + BCL2L1" |
           interaction == "ARSA + DBN1" |
           interaction == "BAD + BCL2L1" |
           interaction == "ARSA + DBN1") %>%
  ggplot(aes(x = orientation, y = interaction, fill = predMat)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), 
                       limits = c(0,1), na.value = "#FBF49C") +
  theme_pubr() +
  ggplot2::theme(text = element_text(family = "Avenir"),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8),
          axis.title = element_text(family = "Avenir Medium"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "right") +
  labs(fill = "probability"),
  
  example_ppi_prediction$predDf %>%
  filter(interaction == "BAD + BCL2L1" |
           interaction == "ARSA + DBN1" |
           interaction == "BAD + BCL2L1" |
           interaction == "ARSA + DBN1") %>%
  ggplot(aes(x = orientation, y = interaction, fill = mean_cBRET)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("#FBF49C", "#ff6b6b", "#f7fff7", "#48C2C5"), 
                       na.value = "#FBF49C") +
  theme_pubr() +
  ggplot2::theme(text = element_text(family = "Avenir"),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8),
          axis.title = element_text(family = "Avenir Medium"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "right") +
  labs(fill = "mean cBRET"), 
  ncol = 1, labels = "AUTO"
)


## ----exampleLearningCurve, message=FALSE, warning=FALSE, fig.cap="(A) Accuracy, (B) Hinge Loss, (C) Binary Cross-Entropy Loss", fig.width=7, fig.height=3----
example_learning_curve$learning_plot

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

