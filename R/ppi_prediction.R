usethis::use_package('dplyr')
usethis::use_package('tidyr')
usethis::use_package('tibble')
usethis::use_package('stringr')
usethis::use_package('e1071')
usethis::use_package('randomForest')

#' Function to use a support vector or random forest machine learning algorithm or to classify quantitative protein-protein interaction data.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import e1071
#' @import randomForest
#'
#' @param PPIdf: binary PPI data set containing interactions to be classified
#' @param referenceSet: reference PPI data set containing reference interactions used to train the svm models
#' @param seed: set seed
#' @param method.scaling: accepted scaling arguments are: "none", "standardize", "robust.scaler", "construct", "orientation"
#' @param iter.scaler: if TRUE and when using "robust.scaler" it iteratively performs robust scaler normalization until the IQR of each construct is within the IQR of all loaded data sets
#' @param range: IQR range used in "robust.scaler"
#' @param data.scaling: speficies which data provided under 'assay' is scaled. When "main" only the first assay is scaled. When "all", all the assays are scaled.
#' @param negative.reference: string in the column "complex" to specify the negative/random interactions
#' @param assay: assay parameters used for training
#' @param all.configurations: if TRUE all orientations for each interactions are used; if FALSE only the highest scoring orientation for each interaction is used
#' @param sampling: use "weighted" or "unweighted" sampling to generate the independent training sets
#' @param weightBy: assay parameter used for weighted sampling
#' @param weightHi: if TRUE weighted sampling uses preferentially higher values; if FALSE weighted sampling uses preferentially smaller values
#' @param model.type: the machine learning algorithm used. Support are support vector machines: "svm" and random fores: "randomForest
#' @param kernelType: the kernel used in training and predicting, see ?e1071::svm for details
#' @param svm.parameters: if TRUE, the best parameters (degree, gamma, coef0, cost) are calculated; if FALSE the parameters must be provided manually
#' @param C: cost of constraints violation (default: 1)—it is the ‘C’-constant of the regularization term in the Lagrange formulation; see ?e1071::svm for details
#' @param gamma: parameter needed for all kernels except linear (default: 1/(data dimension)); see ?e1071::svm for details
#' @param coef0: parameter needed for kernels of type polynomial and sigmoid (default: 0); see ?e1071::svm for details
#' @param degree: parameter needed for kernel of type polynomial (default: 3); see ?e1071::svm for details
#' @param ensembleSize: number of independent training sets assembled
#' @param top: number of highest scoring prs and rrs interactions to randomly sample from; if NULL, only interactions above the 'cs' will be used
#' @param inclusion: number of interactions to include during training of each model (>30); if NULL, number of interactions are calculated from the training sets assembled (ensembleSize) so that each interaction has been sampled with a 99.99% probability.
#' @param cs: sample interactions with quantitative scores above the specified 'cs' from the main assay
#' @param iter: number of iterations performed to reclassify the training set
#' @param verbose: give detailed information
#'
#' @return A list containing the results from the machine learning prediction classes and the parameters used. Input for the learning.curve function.
#' @export
#'
#' @examples
#' load('data/luthy_reference_sets.RData')
#' example <- ppi.prediction(PPIdf = luthy_reference_sets, referenceSet = luthy_reference_sets,
#'                           assay = c("mean_cBRET", "mean_mCit"), negative.reference = c("RRS"),
#'                           model.type = "svm", ensembleSize = 50,
#'                           sampling = 'weighted', weightBy = "mean_cBRET")
#'
ppi.prediction <- function(PPIdf = NULL, referenceSet = NULL, seed = 555,
                           method.scaling = "robust.scaler",
                           iter.scaler = TRUE, range = c(0.25, 0.75), data.scaling = "main",
                           negative.reference = c("RRS", "inter-complex"),
                           assay = c("mean_cBRET", "mean_mCit"), all.configurations = TRUE,
                           sampling = "weighted", weightBy = "mean_cBRET", weightHi = TRUE,
                           model.type = "svm", kernelType = "linear", svm.parameters = FALSE, C = 100, gamma = NULL, coef0 = 0, degree = 2,
                           ensembleSize = 50, top = NULL, inclusion = NULL, cs = "median", iter = 5, verbose = TRUE) {

  base::set.seed(seed)
  #define functions for multi-adaptive sampling
  multiAdaSampling <- function(train.mat, test.mat, label, model.type, kernelType = kernelType, iter = iter, cost = C, gamma = gamma, degree = degree, coef0 = coef0) {

    X <- train.mat
    Y <- label

    model <- NULL
    pred.mat <- NULL

    for (i in seq_len(iter)) {
      tmp <- X
      rownames(tmp) <- NULL

      if (is.null(gamma)) {
        gamma <- 0.1 / length(tmp)
      }

      if (model.type == "svm") {
        model <- e1071::svm(Y ~ ., data = data.frame(tmp, Y), type = "C-classification", kernel = kernelType, probability = TRUE, cost = C, gamma = gamma, degree = degree, coef0 = coef0)
        pred.train <- stats::predict(model, train.mat, decision.values = TRUE, probability = TRUE)
        pred.mat <- attr(pred.train, "probabilities")
      } else if (model.type == "randomforest") {
        model <- randomForest::randomForest(X, factor(Y))
        pred.train <- stats::predict(model, train.mat, type = "prob")
        pred.mat <- pred.train[,c(2,1)]
      } else {
        stop("Invalid model.type. Please choose 'svm' or 'randomforest'.")
      }

      X <- c()
      Y <- c()

      for (j in seq_len(base::ncol(pred.mat))) {
        voteClass <- pred.mat[label == base::colnames(pred.mat)[j], ]
        idx <- base::sample(seq_len(base::nrow(voteClass)),
                      size = base::nrow(voteClass), replace = TRUE,
                      prob = voteClass[, j])
        X <- base::rbind(X, base::as.matrix(train.mat[base::rownames(voteClass)[idx], ]))
        Y <- c(Y, label[base::rownames(voteClass)[idx]])
      }
      base::colnames(X) <- base::colnames(train.mat)
    }

    if(model.type == "svm") {
      pred.values <- stats::predict(model, newdata = test.mat, decision.values = TRUE, probability = TRUE)
      pred <- attr(pred.values, "probabilities")
    } else if (model.type == "randomforest") {
      pred.values <- stats::predict(model, newdata = test.mat, type = "prob")
      pred <- pred.values[,c(2,1)]
    } else {
      stop("Invalid model.type. Please choose 'svm' or 'randomforest'.")
    }

    return(list(pred, pred.values, model))
  }

  `%ni%` <- Negate(`%in%`)
  colnames_necessary <- c("Donor", "Donor_tag", "Donor_protein", "Acceptor", "Acceptor_tag","Acceptor_protein", "complex", "interaction", "sample", "orientation", "data", "score")

  #check reference set training set
  if(is.null(referenceSet)) {
    stop("User must provide a reference data set.")
    }
  else {
    base::message("User defined-reference set provided.")
    referenceSet <- referenceSet %>%
      dplyr::filter(data %in% assay) %>%
      dplyr::mutate(reference = base::ifelse(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|")), "RRS", "PRS"))
  }

  #check for luthy assay and additional parameters
  if(any(assay %ni% referenceSet$data)) {
    stop("the specified 'assay(s)' can not be found in the 'data' column of your 'referenceSet' ", call. = FALSE)
  }
  if(any(colnames_necessary %ni% colnames(referenceSet))) {
    missing_colnames <- colnames_necessary[base::which(colnames_necessary %ni% base::colnames(referenceSet))]
    stop(base::paste0("User provided reference set does not contain all necessary columns. The column(s) ", missing_colnames, " are missing.", "\n",
                "Necessary columns: ", base::paste(colnames_necessary, collapse = ", "), "\n"))
  }

  donor.n <- referenceSet %>%
    dplyr::filter(data == assay[1]) %>%
    dplyr::group_by(Donor) %>%
    dplyr::count()
  acceptor.n <- referenceSet %>%
    dplyr::filter(data == assay[1]) %>%
    dplyr::group_by(Acceptor) %>%
    dplyr::count()

  #check PPI test set
  if(is.null(PPIdf)) {
    stop("No user PPI test-dataset loaded: LuTHy reference set used as test set")
  }
  else {
    base::message("User-defined PPI test-dataset loaded.")
    PPIdf <- PPIdf %>%
      dplyr::filter(data %in% assay)
  }
  if(any(assay %ni% PPIdf$data)) {
    stop("the specified 'assay(s)' can not be found in the 'data' column of your 'PPIdf' ", call. = FALSE)
  }
  if(any(colnames_necessary %ni% base::colnames(PPIdf))) {
    missing_colnames <- colnames_necessary[base::which(colnames_necessary %ni% base::colnames(PPIdf))]
    stop(base::paste0("User provided PPI data set does not contain all necessary columns. The column(s) ", missing_colnames, " are missing.", "\n",
                "Necessary columns: ", base::paste(colnames_necessary, collapse = ", "), "\n"))
  }

  #check on which features data scaling should be applied
  if(data.scaling == "main") {
    assay.scaling <- assay[1]
  }
  if(data.scaling == "all") {
    assay.scaling <- assay
  }

  #check for the method to scale the data
  if(method.scaling == "none") {
    base::message("Data of the training and test sets are not scaled.")
  }
  if(method.scaling %ni% c("none", "standardize", "robust.scaler", "construct", "orientation")) {
    stop("Invalid 'method.scaling'. Choose between 'none', 'standardize', 'robust.scaler', 'construct' or 'orientation'.")
  }

  #z-score scale data to compare data measured at different plate readers
  if(method.scaling == "standardize") {
    if(verbose)
      base::message("Reference-set is z-score scaled.")
    for(i in assay.scaling) {
      referenceSet <- referenceSet %>%
        dplyr::filter(data %in% assay) %>%
        tidyr::pivot_wider(names_from = data, values_from = score) %>%
        dplyr::mutate(!!(rlang::sym(assay[1])) := scale(!!(rlang::sym(assay[1])))) %>%
        tidyr::pivot_longer(cols = assay, names_to = "data", values_to = "score")
    }
  }

  #construct normalization
  if(method.scaling == "construct") {
    #stop("This function is under construction")
    if(donor.n == FALSE | acceptor.n == FALSE) {
      base::message("'Donor' and 'Acceptor' constructs in the reference set are normalized for their median + 3SD.")
      for(i in assay.scaling) {
        scaler.global <- stats::median(referenceSet %>%
                                         tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
                                         dplyr::pull(i),
                                       na.rm = TRUE)
        mean.global <- base::mean(referenceSet %>%
                                    tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
                                    dplyr::pull(i),
                                  na.rm = TRUE)
        sd.global <- stats::sd(referenceSet %>%
                                 tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
                                 dplyr::pull(i),
                               na.rm = TRUE)
        donor.cf <- referenceSet %>%
          dplyr::filter(data %in% assay) %>%
          tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
          dplyr::group_by(Donor) %>%
          dplyr::summarise(median = stats::median(!!(rlang::sym(i)), na.rm = TRUE),
                           mean = base::mean(!!(rlang::sym(i)), na.rm = TRUE),
                           sd = stats::sd(!!(rlang::sym(i)), na.rm = TRUE)) %>%
          dplyr::mutate(cf.donor = median - scaler.global)
        referenceSet <- referenceSet %>%
          dplyr::filter(data %in% assay) %>%
          tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
          dplyr::left_join(donor.cf %>% dplyr::select(Donor, cf.donor), by = "Donor") %>%
          dplyr::mutate(!!(rlang::sym(i)) := !!(rlang::sym(i)) - cf.donor, .keep = "unused")
        acceptor.cf <- referenceSet %>%
          dplyr::group_by(Acceptor) %>%
          dplyr::summarise(median = stats::median(!!(rlang::sym(i)), na.rm = TRUE),
                           mean = base::mean(!!(rlang::sym(i)), na.rm = TRUE),
                           sd = stats::sd(!!(rlang::sym(i)), na.rm = TRUE)) %>%
          dplyr::mutate(cf.acceptor = median - scaler.global)
        referenceSet <- referenceSet %>%
          dplyr::left_join(acceptor.cf %>% dplyr::select(Acceptor, cf.acceptor), by = "Acceptor") %>%
          dplyr::mutate(!!(rlang::sym(i)) := !!(rlang::sym(i)) - cf.acceptor, .keep = "unused") %>%
          tidyr::pivot_longer(cols = assay, names_to = "data", values_to = "score")
      }
    } else(base::message("Not enough 'Donor' and 'Acceptor' constructs in the reference set to construct normalize."))

  }

  #orientation normalization
  if(method.scaling == "orientation") {
    base::message("median scaling the reference set by orientation.")
    for(i in assay.scaling) {
      scaler.global <- stats::median(referenceSet %>%
                                       tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
                                       dplyr::pull(i),
                                     na.rm = TRUE)
      mean.global <- base::mean(referenceSet %>%
                                  tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
                                  dplyr::pull(i),
                                na.rm = TRUE)
      sd.global <- stats::sd(referenceSet %>%
                               tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
                               dplyr::pull(i),
                             na.rm = TRUE)
      orientation.cf <- referenceSet %>%
        dplyr::filter(data %in% assay) %>%
        tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
        dplyr::group_by(orientation) %>%
        dplyr::summarise(median = stats::median(!!(rlang::sym(i)), na.rm = TRUE),
                         mean = base::mean(!!(rlang::sym(i)), na.rm = TRUE),
                         sd = stats::sd(!!(rlang::sym(i)), na.rm = TRUE)) %>%
        dplyr::mutate(cf = median - scaler.global)
      referenceSet <- referenceSet %>%
        dplyr::filter(data %in% assay) %>%
        tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
        dplyr::left_join(orientation.cf %>% dplyr::select(orientation, cf), by = "orientation") %>%
        dplyr::mutate(!!(rlang::sym(i)) := !!(rlang::sym(i)) - cf, .keep = "unused") %>%
        tidyr::pivot_longer(cols = assay, names_to = "data", values_to = "score")
    }
  }

  #robust scaler for each donor and acceptor construct: scaling to median and quantile
  IQR.custom <- function(x, na.rm = FALSE, type = 7, range = range)  {
    base::diff(stats::quantile(as.numeric(x), probs = range, na.rm = na.rm, names = FALSE,
                               type = type))
  }
  scaler.global <- list()
  IQR.global <- list()
  referenceSet.list <- c()

  for(j in assay.scaling) {
    scaler.global[[j]] <- stats::median(c(referenceSet %>% tidyr::pivot_wider(names_from = "data", values_from = "score") %>% dplyr::pull(j),
                                   PPIdf %>% tidyr::pivot_wider(names_from = "data", values_from = "score") %>% dplyr::pull(j)), na.rm = TRUE)
    IQR.global[[j]] <- IQR.custom(c(referenceSet %>% tidyr::pivot_wider(names_from = "data", values_from = "score") %>% dplyr::pull(j),
                                    PPIdf %>% tidyr::pivot_wider(names_from = "data", values_from = "score")  %>% dplyr::pull(j)), na.rm = TRUE, range = range)
  }
  if(method.scaling == "robust.scaler") {
    for(j in assay.scaling) {
      if(any(donor.n$n >= 20) == FALSE | any(acceptor.n$n >= 20) == FALSE) {
        base::message("Not enough 'Donor' and 'Acceptor' constructs in the provided 'reference.set' to perform robust scaler.")
        referenceSet.scaler <- base::data.frame(base::matrix(nrow = 0, ncol = base::ncol(referenceSet)))
        base::colnames(referenceSet.scaler) <- base::colnames(referenceSet)
      }

      if(any(donor.n$n >= 20) == TRUE | any(acceptor.n$n >= 20) == TRUE){

        donor.constructs <- donor.n %>% dplyr::filter(n >= 20) %>% dplyr::pull(Donor)
        acceptor.constructs <- acceptor.n %>% dplyr::filter(n >= 20) %>% dplyr::pull(Acceptor)

        referenceSet.scaler <- referenceSet %>%
          dplyr::filter(data %in% j) %>%
          dplyr::filter(Donor %in% donor.constructs | Acceptor %in% acceptor.constructs)

        scaler <- TRUE
        count <- 0
        while(scaler) {
          count <- count + 1
          donor.cf <- referenceSet.scaler %>%
            dplyr::filter(data %in% j & Donor %in% donor.constructs) %>%
            tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
            dplyr::filter(!is.na(!!(rlang::sym(j)))) %>%
            dplyr::group_by(Donor) %>%
            dplyr::summarise(median = stats::median(!!(rlang::sym(j)), na.rm = TRUE),
                             IQR = IQR.custom(!!(rlang::sym(j)), na.rm = TRUE, range = range)) %>%
            dplyr::mutate(cf.donor = median - scaler.global[[j]],
                          median = base::ifelse(is.na(median), 0, median),
                          IQR = base::ifelse(is.na(IQR), IQR.global[[j]], IQR),
                          cf.donor = base::ifelse(is.na(cf.donor), 0, cf.donor))
          referenceSet.scaler <- referenceSet.scaler %>%
            dplyr::filter(data %in% j) %>%
            tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
            dplyr::filter(!is.na(!!(rlang::sym(j)))) %>%
            dplyr::left_join(donor.cf %>% dplyr::select(Donor, IQR, cf.donor), by = "Donor") %>%
            dplyr::mutate(IQR = base::ifelse(is.na(IQR), 1, IQR),
                          cf.donor = base::ifelse(is.na(cf.donor), 0, cf.donor)) %>%
            dplyr::mutate(!!(rlang::sym(j)) := (!!(rlang::sym(j)) - cf.donor) / base::ifelse(IQR > IQR.global[[j]], (IQR / IQR.global[[j]]), 1), .keep = "unused")
          acceptor.cf <- referenceSet.scaler %>%
            dplyr::filter(Acceptor %in% acceptor.constructs) %>%
            dplyr::group_by(Acceptor) %>%
            dplyr::summarise(median = stats::median(!!(rlang::sym(j)), na.rm = TRUE),
                             IQR = IQR.custom(!!(rlang::sym(j)), na.rm = TRUE, range = range)) %>%
            dplyr::mutate(cf.acceptor = median - scaler.global[[j]],
                          median = base::ifelse(is.na(median), 0, median),
                          IQR = base::ifelse(is.na(IQR), IQR.global[[j]], IQR),
                          cf.acceptor = base::ifelse(is.na(cf.acceptor), 0, cf.acceptor))
          referenceSet.scaler <- referenceSet.scaler %>%
            dplyr::left_join(acceptor.cf %>% dplyr::select(Acceptor, IQR, cf.acceptor), by = "Acceptor") %>%
            dplyr::mutate(IQR = base::ifelse(is.na(IQR), 1, IQR),
                          cf.acceptor = base::ifelse(is.na(cf.acceptor), 0, cf.acceptor)) %>%
            dplyr::mutate(!!(rlang::sym(j)) := (!!(rlang::sym(j)) - cf.acceptor) / base::ifelse(IQR > IQR.global[[j]], (IQR / IQR.global[[j]]), 1), .keep = "unused") %>%
            tidyr::pivot_longer(cols = j, names_to = "data", values_to = "score")

          constructs <- referenceSet.scaler %>%
            dplyr::filter(data %in% j & Donor %in% donor.constructs) %>%
            tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
            dplyr::group_by(construct = Donor) %>%
            dplyr::summarise(median = stats::median(!!(rlang::sym(j)), na.rm = TRUE),
                             IQR = IQR.custom(!!(rlang::sym(j)), na.rm = TRUE, range = range)) %>%
            base::rbind(
              referenceSet.scaler %>%
                dplyr::filter(data %in% j & Acceptor %in% acceptor.constructs) %>%
                tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
                dplyr::group_by(construct = Acceptor) %>%
                dplyr::summarise(median = stats::median(!!(rlang::sym(j)), na.rm = TRUE),
                                 IQR = IQR.custom(!!(rlang::sym(j)), na.rm = TRUE, range = range))
            ) %>%
            mutate(scaler.IQR = base::ifelse(round(IQR, 4) > round(IQR.global[[j]],4), TRUE, FALSE),
                   scaler.median = base::ifelse(round(median, 4) > round(scaler.global[[j]],4), TRUE, FALSE))

          if(iter.scaler == TRUE) {
            if(any(constructs$scaler.IQR, na.rm = TRUE)) {
              scaler <- TRUE
            } else {
              scaler <- FALSE
            }
          } else {
            scaler <- FALSE
          }

        }
        base::message("Construct-specific median normalization for ", shQuote(j)," is done for ", length(base::which(donor.n$n >= 20)), " 'Donor construcst and ", length(base::which(acceptor.n$n >= 20)), " ' Acceptor' construcsts. \n",
                "Construct-specific IQR normalization for ", shQuote(j)," is done for ", length(base::which(donor.cf$IQR >= IQR.global[[j]])), " 'Donor' constructs and ", length(base::which(acceptor.cf$IQR >= IQR.global[[j]])), " 'Acceptor' constructs that showed a higher IQR than the global IQR of the 'reference.Set' in ", count, " iteration(s).")
      }

      if(any(donor.n$n < 20) == FALSE & any(acceptor.n$n < 20) == FALSE) {
        referenceSet.notScaler <- base::data.frame(base::matrix(nrow = 0, ncol = base::ncol(referenceSet)))
        base::colnames(referenceSet.notScaler) <- base::colnames(referenceSet)
      }

      if(any(donor.n$n < 20) == TRUE & any(acceptor.n$n < 20) == TRUE) {
        base::message("There are ", length(base::which(donor.n$n < 20)), " 'Donor' constructs with not enough interactions (<20) to perform 'Robust scaler' normalization. \n",
                "There are ", length(base::which(acceptor.n$n < 20)), " 'Acceptor' constructs with not enough interactions (<20) to perform 'Robust scaler' normalization. \n")

        donor.constructs <- donor.n %>% dplyr::filter(n < 20) %>% dplyr::pull(Donor)
        acceptor.constructs <- acceptor.n %>% dplyr::filter(n < 20) %>% dplyr::pull(Acceptor)

        if(length(donor.constructs) != 0 & length(acceptor.constructs) != 0) {
          referenceSet.notScaler <- referenceSet %>%
            dplyr::filter(data %in% j) %>%
            dplyr::filter(Donor %in% donor.constructs & Acceptor %in% acceptor.constructs)
        }
      }

      referenceSet.list[[j]] <- referenceSet.scaler %>% base::rbind(referenceSet.notScaler)

    }

    referenceSet.unscaled <- referenceSet %>%
      dplyr::filter(data %ni% assay.scaling)

    referenceSet <- dplyr::bind_rows(referenceSet.list) %>%
      base::rbind(referenceSet.unscaled)

  }

  #use all configurations/orientations used to measure "one" interactions
  if(all.configurations == FALSE) {
    if(verbose)
      base::message("Only the highest scoring orientation for each interaction is used for the reference set.")
    for(j in assay.scaling) {
      referenceSet <- referenceSet %>%
        dplyr::filter(data %in% j) %>%
        tidyr::pivot_wider(names_from = data, values_from = score) %>%
        dplyr::group_by(complex, interaction) %>%
        dplyr::summarise(!!(rlang::sym(j)) := max(!!(rlang::sym(j)), na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(referenceSet %>%
                    dplyr::filter(data == j) %>%
                    tidyr::pivot_wider(names_from = data, values_from = score),
                  by = c("complex", "interaction", j)) %>%
        tidyr::pivot_longer(cols = assay, names_to = "data", values_to = "score")
    }
  }

  donor.n <- PPIdf %>% dplyr::filter(data == assay[1]) %>% dplyr::group_by(Donor) %>% dplyr::count()
  acceptor.n <- PPIdf %>% dplyr::filter(data == assay[1]) %>% dplyr::group_by(Acceptor) %>% dplyr::count()

  #z-score scaling of test set
  if(method.scaling == "standardize") {
    if(verbose)
      base::message("Test-set is z-score scaled.")
    for(i in assay.scaling) {
      PPIdf <- PPIdf %>%
        dplyr::filter(data %in% assay) %>%
        tidyr::pivot_wider(names_from = data, values_from = score) %>%
        dplyr::mutate(!!(rlang::sym(i)) := scale(!!(rlang::sym(i)))) %>%
        tidyr::pivot_longer(cols = assay, names_to = "data", values_to = "score")
    }
  }

  #construct normalization
  if(method.scaling == "construct") {
    #stop("This function is under construction")
    if(donor.n == FALSE | acceptor.n == FALSE) {
      base::message("'Donor' and 'Acceptor' constructs in the reference set are normalized for their median + 3SD.")
      for(i in assay.scaling) {
        scaler.global <- stats::median(PPIdf %>% tidyr::pivot_wider(names_from = "data", values_from = "score") %>% dplyr::pull(i), na.rm = TRUE)
        mean.global <- base::mean(PPIdf %>% tidyr::pivot_wider(names_from = "data", values_from = "score") %>% dplyr::pull(i), na.rm = TRUE)
        sd.global <- stats::sd(PPIdf %>% tidyr::pivot_wider(names_from = "data", values_from = "score") %>% dplyr::pull(i), na.rm = TRUE)
        donor.cf <- PPIdf %>%
          dplyr::filter(data %in% assay) %>%
          tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
          dplyr::group_by(Donor) %>%
          dplyr::summarise(median = stats::median(!!(rlang::sym(i)), na.rm = TRUE),
                           mean = base::mean(!!(rlang::sym(i)), na.rm = TRUE),
                           sd = stats::sd(!!(rlang::sym(i)), na.rm = TRUE)) %>%
          dplyr::mutate(cf.donor = median - scaler.global)
        PPIdf <- PPIdf %>%
          dplyr::filter(data %in% assay) %>%
          tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
          dplyr::left_join(donor.cf %>% dplyr::select(Donor, cf.donor), by = "Donor") %>%
          dplyr::mutate(!!(rlang::sym(i)) := !!(rlang::sym(i)) - cf.donor, .keep = "unused")
        acceptor.cf <- PPIdf %>%
          dplyr::group_by(Acceptor) %>%
          dplyr::summarise(median = stats::median(!!(rlang::sym(i)), na.rm = TRUE),
                           mean = base::mean(!!(rlang::sym(i)), na.rm = TRUE),
                           sd = stats::sd(!!(rlang::sym(i)), na.rm = TRUE)) %>%
          dplyr::mutate(cf.acceptor = median - scaler.global)
        PPIdf <- PPIdf %>%
          dplyr::left_join(acceptor.cf %>% dplyr::select(Acceptor, cf.acceptor), by = "Acceptor") %>%
          dplyr::mutate(!!(rlang::sym(i)) := !!(rlang::sym(i)) - cf.acceptor, .keep = "unused") %>%
          tidyr::pivot_longer(cols = assay, names_to = "data", values_to = "score")
      }
    } else(base::message("Not enough 'Donor' and 'Acceptor' constructs in the reference set to construct normalize."))

  }

  #orientation normalization
  if(method.scaling == "orientation") {
    base::message("median scaling the test set by orientation.")
    for(i in assay.scaling) {
      scaler.global <- stats::median(PPIdf %>% tidyr::pivot_wider(names_from = "data", values_from = "score") %>% dplyr::pull(i), na.rm = TRUE)
      mean.global <- base::mean(PPIdf %>% tidyr::pivot_wider(names_from = "data", values_from = "score") %>% dplyr::pull(i), na.rm = TRUE)
      sd.global <- stats::sd(PPIdf %>% tidyr::pivot_wider(names_from = "data", values_from = "score") %>% dplyr::pull(i), na.rm = TRUE)
      orientation.cf <- PPIdf %>%
        dplyr::filter(data %in% assay) %>%
        tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
        dplyr::group_by(orientation) %>%
        dplyr::summarise(median = stats::median(!!(rlang::sym(i)), na.rm = TRUE),
                         mean = base::mean(!!(rlang::sym(i)), na.rm = TRUE),
                         sd = stats::sd(!!(rlang::sym(i)), na.rm = TRUE)) %>%
        dplyr::mutate(cf = median - scaler.global)
      PPIdf <- PPIdf %>%
        dplyr::filter(data %in% assay) %>%
        tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
        dplyr::left_join(orientation.cf %>% dplyr::select(orientation, cf), by = "orientation") %>%
        dplyr::mutate(!!(rlang::sym(i)) := !!(rlang::sym(i)) - cf, .keep = "unused") %>%
        tidyr::pivot_longer(cols = assay, names_to = "data", values_to = "score")
    }
  }

  #robust scaler for each donor and acceptor construct: scaling to median and quantile
  PPIdf.list <- c()

  if(method.scaling == "robust.scaler") {
    for(j in assay.scaling) {

      if(any(donor.n$n >= 20) == FALSE | any(acceptor.n$n >= 20) == FALSE) {
        base::message("Not enough 'Donor' and 'Acceptor' constructs in the provided 'PPIdf' to perform robust scaler.")
        PPIdf.scaler <- base::data.frame(base::matrix(nrow = 0, ncol = base::ncol(PPIdf)))
        base::colnames(PPIdf.scaler) <- base::colnames(PPIdf)
      }

      if(any(donor.n$n >= 20) == TRUE | any(acceptor.n$n >= 20) == TRUE){

        donor.constructs <- donor.n %>% dplyr::filter(n >= 20) %>% dplyr::pull(Donor)
        acceptor.constructs <- acceptor.n %>% dplyr::filter(n >= 20) %>% dplyr::pull(Acceptor)

        PPIdf.scaler <- PPIdf %>%
          dplyr::filter(data %in% j) %>%
          dplyr::filter(Donor %in% donor.constructs | Acceptor %in% acceptor.constructs)

        scaler <- TRUE
        count <- 0
        while(scaler) {
          count <- count + 1

          donor.cf <- PPIdf.scaler %>%
            dplyr::filter(data %in% j & Donor %in% donor.constructs) %>%
            tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
            dplyr::filter(!is.na(!!(rlang::sym(j)))) %>%
            dplyr::group_by(Donor) %>%
            dplyr::summarise(median = stats::median(!!(rlang::sym(j)), na.rm = TRUE),
                             IQR = IQR.custom(!!(rlang::sym(j)), na.rm = TRUE, range = range)) %>%
            dplyr::mutate(cf.donor = median - scaler.global[[j]],
                          median = base::ifelse(is.na(median), 0, median),
                          IQR = base::ifelse(is.na(IQR), IQR.global[[j]], IQR),
                          cf.donor = base::ifelse(is.na(cf.donor), 0, cf.donor))
          PPIdf.scaler <- PPIdf.scaler %>%
            dplyr::filter(data %in% j) %>%
            tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
            dplyr::filter(!is.na(!!(rlang::sym(j)))) %>%
            dplyr::left_join(donor.cf %>% dplyr::select(Donor, IQR, cf.donor), by = "Donor") %>%
            dplyr::mutate(IQR = base::ifelse(is.na(IQR), 1, IQR),
                          cf.donor = base::ifelse(is.na(cf.donor), 0, cf.donor)) %>%
            dplyr::mutate(!!(rlang::sym(j)) := (!!(rlang::sym(j)) - cf.donor) / base::ifelse(IQR > IQR.global[[j]], (IQR / IQR.global[[j]]), 1), .keep = "unused")
          acceptor.cf <- PPIdf.scaler %>%
            dplyr::filter(Acceptor %in% acceptor.constructs) %>%
            dplyr::group_by(Acceptor) %>%
            dplyr::summarise(median = stats::median(!!(rlang::sym(j)), na.rm = TRUE),
                             IQR = IQR.custom(!!(rlang::sym(j)), na.rm = TRUE, range = range)) %>%
            dplyr::mutate(cf.acceptor = median - scaler.global[[j]],
                          median = base::ifelse(is.na(median), 0, median),
                          IQR = base::ifelse(is.na(IQR), IQR.global[[j]], IQR),
                          cf.acceptor = base::ifelse(is.na(cf.acceptor), 0, cf.acceptor))
          PPIdf.scaler <- PPIdf.scaler %>%
            dplyr::left_join(acceptor.cf %>% dplyr::select(Acceptor, IQR, cf.acceptor), by = "Acceptor") %>%
            dplyr::mutate(IQR = base::ifelse(is.na(IQR), 1, IQR),
                          cf.acceptor = base::ifelse(is.na(cf.acceptor), 0, cf.acceptor)) %>%
            dplyr::mutate(!!(rlang::sym(j)) := (!!(rlang::sym(j)) - cf.acceptor) / base::ifelse(IQR > IQR.global[[j]], (IQR / IQR.global[[j]]), 1), .keep = "unused") %>%
            tidyr::pivot_longer(cols = j, names_to = "data", values_to = "score")


          constructs <- PPIdf.scaler %>%
            dplyr::filter(data %in% j & Donor %in% donor.constructs) %>%
            tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
            dplyr::group_by(construct = Donor) %>%
            dplyr::summarise(median = stats::median(!!(rlang::sym(j)), na.rm = TRUE),
                             IQR = IQR.custom(!!(rlang::sym(j)), na.rm = TRUE, range = range)) %>%
            base::rbind(
              PPIdf.scaler %>%
                dplyr::filter(data %in% j & Acceptor %in% acceptor.constructs) %>%
                tidyr::pivot_wider(names_from = "data", values_from = "score") %>%
                dplyr::group_by(construct = Acceptor) %>%
                dplyr::summarise(median = stats::median(!!(rlang::sym(j)), na.rm = TRUE),
                                 IQR = IQR.custom(!!(rlang::sym(j)), na.rm = TRUE, range = range))
            ) %>%
            mutate(scaler.IQR = base::ifelse(round(IQR, 4) > round(IQR.global[[j]],4), TRUE, FALSE),
                   scaler.median = base::ifelse(round(median, 4) > round(scaler.global[[j]],4), TRUE, FALSE))

          if(iter.scaler == TRUE) {
            if(any(constructs$scaler.IQR, na.rm = TRUE)) {
              scaler <- TRUE
            } else {
              scaler <- FALSE
            }
          } else {
            scaler <- FALSE
          }
        }
        base::message("Construct-specific median normalization for ", shQuote(j)," is done for ", length(base::which(donor.n$n >= 20)), " 'Donor construcst and ", length(base::which(acceptor.n$n >= 20)), " ' Acceptor' construcsts. \n",
                "Construct-specific IQR normalization for ", shQuote(j)," is done for ", length(base::which(donor.cf$IQR >= IQR.global[[j]])), " 'Donor' constructs and ", length(base::which(acceptor.cf$IQR >= IQR.global[[j]])), " 'Acceptor' constructs that showed a higher IQR than the global IQR of the 'PPIdf' in ", count, " iteration(s).")
      }

      if(any(donor.n$n < 20) == FALSE & any(acceptor.n$n < 20) == FALSE) {
        PPIdf.notScaler <- base::data.frame(base::matrix(nrow = 0, ncol = base::ncol(PPIdf)))
        base::colnames(PPIdf.notScaler) <- base::colnames(PPIdf)
      }

      if(any(donor.n$n < 20) == TRUE & any(acceptor.n$n < 20) == TRUE) {
        base::message("There are ", length(base::which(donor.n$n < 20)), " 'Donor' constructs with not enough interactions (<20) to perform 'Robust scaler' normalization. \n",
                "There are ", length(base::which(acceptor.n$n < 20)), " 'Acceptor' constructs with not enough interactions (<20) to perform 'Robust scaler' normalization. \n")

        donor.constructs <- donor.n %>% dplyr::filter(n < 20) %>% dplyr::pull(Donor)
        acceptor.constructs <- acceptor.n %>% dplyr::filter(n < 20) %>% dplyr::pull(Acceptor)

        if(length(donor.constructs) != 0 & length(acceptor.constructs) != 0) {
          PPIdf.notScaler <- PPIdf %>%
            dplyr::filter(data %in% j) %>%
            dplyr::filter(Donor %in% donor.constructs & Acceptor %in% acceptor.constructs)
        }

      }

      PPIdf.list[[j]] <- PPIdf.scaler %>% base::rbind(PPIdf.notScaler)

    }

    PPIdf.unscaled <- PPIdf %>%
      dplyr::filter(data %ni% assay.scaling)

    PPIdf <- dplyr::bind_rows(PPIdf.list) %>%
      base::rbind(PPIdf.unscaled)
  }

  #test all or highest scoring configurations per interaction
  if(all.configurations == FALSE) {
    if(verbose)
      base::message("Only the highest scoring orientation for each interaction is used for the test set.")

    for(j in assay.scaling) {
      PPIdf <- PPIdf %>%
        dplyr::filter(data %in% j) %>%
        tidyr::pivot_wider(names_from = data, values_from = score) %>%
        dplyr::group_by(complex, interaction) %>%
        dplyr::summarise(!!(rlang::sym(j)) := max(!!(rlang::sym(j)), na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(PPIdf %>%
                    dplyr::filter(data %in% j) %>%
                    tidyr::pivot_wider(names_from = data, values_from = score),
                  by = c("complex", "interaction", j)) %>%
        dplyr::filter(is.finite(!!(rlang::sym(assay[1])))) %>%
        tidyr::pivot_longer(cols = assay, names_to = "data", values_to = "score")
    }
  }

  #Assemble training and test matrix
  TrainTestMat <- referenceSet %>%
    base::unique() %>%
    dplyr::filter(data %in% assay) %>%
    tidyr::unite(complex, reference, interaction, sample, orientation, col = "sample", sep = ";", remove = FALSE) %>%
    tidyr::pivot_wider(names_from = data, values_from = score) %>%
    dplyr::filter(across(.cols = assay, ~!is.na(.x))) %>%
    tibble::column_to_rownames("sample") %>%
    dplyr::select(assay) %>%
    base::as.matrix()

  predTrainMat <- base::matrix(0, nrow = base::nrow(TrainTestMat), ncol = 1)
  base::rownames(predTrainMat) <- base::rownames(TrainTestMat)

  testMat <- PPIdf %>%
    base::unique() %>%
    dplyr::filter(data %in% assay) %>%
    tidyr::unite(complex, interaction, sample, orientation, col = "sample", sep = ";") %>%
    tidyr::pivot_wider(names_from = data, values_from = score) %>%
    dplyr::filter(across(.cols = assay, ~!is.na(.x))) %>%
    tibble::column_to_rownames("sample") %>%
    dplyr::select(assay) %>%
    base::as.matrix()

  predTrainMat <- list()
  predTrain.values.e <- list()
  predTrain.model.e <- list()
  pred.values.e <- list()
  model.e <- list()
  predMat <- list()
  TrainMat <- TrainTestMat
  TestMat <- testMat
  counter = 0

  #assemble training set by randomly sampling n ('top') prs and rrs interactions for the multi adaptive sampling and repeat according to ensemble Size
  if(svm.parameters == TRUE) {
    if(!is.null(C) | !is.null(gamma)) {
      base::message("User-provided parameters will be overwritten: \n")
    }
    base::message("Calculating best parameters for ", shQuote(kernelType)," kernel type. \n")
    tune.y <- base::factor(stringr::str_extract(base::rownames(TrainTestMat), "PRS|RRS"))
    names(tune.y) <- base::rownames(TrainTestMat)
    tune.out <- e1071::tune.svm(x = TrainTestMat, y = tune.y,
                                type = "C-classification", kernel = kernelType,
                                cost = 10^(-1:2),
                                gamma = base::ifelse(kernelType == "linear", 0, 10^(-4:1)),
                                coef0 = base::ifelse(stringr::str_detect(kernelType, "polynomial|sigmoid"), 0:3, 0),
                                degree = base::ifelse(kernelType == "polynomial", 1:5, 2))
    gamma <- tune.out$best.parameters$gamma
    C <- tune.out$best.parameters$cost
    coef0 <- tune.out$best.parameters$coef0
    degree <- tune.out$best.parameters$degree
    base::message("Type: ",  classificationType, "\n",
            "kernel: ", kernelType, "\n",
            "Cost: ", C, "\n",
            "gamma: ", gamma, "\n",
            "coef0: ", coef0, "\n",
            "degree: ", degree, "\n")
  }

  if(sampling != "unweighted" && sampling != "weighted") {
    stop("Invalid 'sampling'. Choose between 'unweighted' or 'weighted'.")
  }
  training.sets <- list()
  ensembleSize_total <- ensembleSize
  while(any(base::rowSums(!is.na(TrainMat[,-c(1:length(assay))]), na.rm = TRUE) == 0) | any(base::rowSums(!is.na(TestMat[,-c(1:length(assay))]), na.rm = TRUE) == 0)) {
    counter = counter+1
    if(counter > 1) {
      weightBy <- NULL
      sampling <- "unweighted"
      ensembleSize <- 10
      base::message("Not all interactions have been predicted. Sampling repeated ", ensembleSize, " times. \n")
      ensembleSize_total <- ensembleSize_total+ensembleSize
      probability <- base::rowMeans(TrainMat[,-c(1:length(assay))], na.rm = TRUE)
    }
    for (e in seq(from = ifelse(ensembleSize_total > ensembleSize, ensembleSize_total - ensembleSize, 1),
                  to = ifelse(ensembleSize_total > ensembleSize, ensembleSize_total, ensembleSize))) {

      if(ensembleSize > 1) {
        if(verbose)
          if(e %% 10 == 0) {
            message(e, ".")
          }
        if(cs == "median") {
          cs <- median(referenceSet %>% filter(data == assay[1]) %>% pull(score), na.rm = TRUE)
        }

        if(cs == "all") {
          cs <- min(referenceSet %>% filter(data == assay[1]) %>% pull(score), na.rm = TRUE)
        }

        n.ppis <- referenceSet %>%
          filter(data == assay[1] & score > cs) %>% nrow()

        if(is.null(top)) {
          top <- round(n.ppis, 0)
        }
        if(is.null(inclusion)) {
          inclusion <- round(log(1-0.9999) / log(1-ensembleSize/top),0)
          if(inclusion < 30) {
            cat(paste0("minimum number of interactions to include are 30", "\n"))
            inclusion <- 30
          }
        }
        if(is.null(weightBy)){
          weightBy <- assay[1]
        }
        if(sampling == "unweighted") {
          if(counter == 1) {
            prs.interactions <- referenceSet %>%
              filter(data == weightBy & !is.na(score) & reference == "PRS") %>%
              slice_max(n = top, order_by = score) %>%
              slice_sample(n = inclusion) %>%
              dplyr::select(Donor, Donor_tag, Donor_protein, Acceptor, Acceptor_protein, Acceptor_tag, orientation, complex, reference)
          }
          if(counter > 1) {
            prs.interactions <- data.frame(probability) %>%
              rownames_to_column("sample") %>%
              separate(col = "sample", into = c("complex", "reference", "interaction", "sample", "orientation"), sep = ";") %>%
              right_join(referenceSet, by = c("complex", "reference", "interaction", "sample", "orientation")) %>%
              filter(data == weightBy & !is.na(weightBy) & reference == "PRS" & !is.na(probability)) %>%
              slice_max(n = top, order_by = score) %>%
              slice_sample(n = inclusion) %>%
              dplyr::select(Donor, Donor_tag, Donor_protein, Acceptor, Acceptor_protein, Acceptor_tag, orientation, complex, reference)

          }
        } else if(sampling == "weighted") {
          if(weightHi == TRUE) {
            prs.interactions <- referenceSet %>%
              filter(data == weightBy & !is.na(score) & reference == "PRS") %>%
              slice_max(n = top, order_by = score) %>%
              slice_sample(n = inclusion, weight_by = score+abs(min(referenceSet %>% filter(data == weightBy) %>% pull(score), na.rm = TRUE))) %>%
              dplyr::select(Donor, Donor_tag, Donor_protein, Acceptor, Acceptor_protein, Acceptor_tag, orientation, complex, reference)
          }
          if(weightHi == FALSE) {
            prs.interactions <- referenceSet %>%
              filter(data == weightBy & !is.na(score) & reference == "PRS") %>%
              slice_min(n = top, order_by = score) %>%
              slice_sample(n = inclusion, weight_by = score+abs(min(referenceSet %>% filter(data == weightBy) %>% pull(score), na.rm = TRUE))) %>%
              dplyr::select(Donor, Donor_tag, Donor_protein, Acceptor, Acceptor_protein, Acceptor_tag, orientation, complex, reference)
          }
        }

        rrs.interactions <- referenceSet %>%
          filter(data == weightBy & !is.na(score) & reference == "RRS") %>%
          dplyr::select(Donor, Donor_tag, Donor_protein, Acceptor, Acceptor_protein, Acceptor_tag, orientation, complex, reference) %>%
          anti_join(prs.interactions, by = c("Donor", "Donor_tag", "Donor_protein", "Acceptor", "Acceptor_protein", "Acceptor_tag", "orientation", "complex", "reference")) %>%
          slice_sample(n = nrow(prs.interactions))
        }

      if(ensembleSize == 1) {
        prs.interactions <- referenceSet %>%
          filter(data == weightBy & !is.na(score) & reference == "PRS") %>%
          dplyr::select(Donor, Donor_tag, Donor_protein, Acceptor, Acceptor_protein, Acceptor_tag, orientation, complex, reference)

        rrs.interactions <- referenceSet %>%
          filter(data == weightBy & !is.na(score) & reference == "RRS") %>%
          dplyr::select(Donor, Donor_tag, Donor_protein, Acceptor, Acceptor_protein, Acceptor_tag, orientation, complex, reference) %>%
          anti_join(prs.interactions, by = c("Donor", "Donor_tag", "Donor_protein", "Acceptor", "Acceptor_protein", "Acceptor_tag", "orientation", "complex", "reference"))
        }

      positive.train <- referenceSet %>%
        filter(data %in% assay) %>%
        inner_join(prs.interactions, by = c("Donor", "Donor_tag", "Donor_protein", "Acceptor", "Acceptor_protein", "Acceptor_tag", "orientation", "complex", "reference")) %>%
        pivot_wider(names_from = "data", values_from = "score") %>%
        unite(complex, reference, interaction, sample, orientation, col = "sample", sep = ";") %>%
        column_to_rownames("sample") %>%
        dplyr::select(assay) %>%
        filter(across(.cols = assay, .fns = function(x) !is.na(x))) %>%
        as.matrix()
      positive.cls <- rep(2, nrow(positive.train))

      negative.train <- referenceSet %>%
        filter(data %in% assay) %>%
        inner_join(rrs.interactions, by = c("Donor", "Donor_tag", "Donor_protein", "Acceptor", "Acceptor_protein", "Acceptor_tag", "orientation", "complex", "reference")) %>%
        unite(complex, reference, interaction, sample, orientation, col = "sample", sep = ";") %>%
        pivot_wider(names_from = data, values_from = score) %>%
        column_to_rownames("sample") %>%
        dplyr::select(assay) %>%
        filter(across(.cols = assay, .fns = function(x) !is.na(x))) %>%
        as.matrix()
      negative.cls <- rep(1, nrow(negative.train))

      train.mat <- rbind(positive.train, negative.train)
      training.sets[[e]] <- train.mat
      cls <- as.factor(c(positive.cls, negative.cls))
      names(cls) <- rownames(train.mat)

      train.ppis <- as.data.frame(train.mat) %>% rownames_to_column("sample") %>%
        separate(col = "sample", into = c("complex", "reference", "interaction", "sample", "orientation"), sep = ";") %>%
        dplyr::select(-reference) %>%
        unite(complex, interaction, sample, orientation, col = "sample", sep = ";") %>%
        pull(sample)


      #predict positive interactions based on the training sets; also test performance on the reference set
      pred.reference <- multiAdaSampling(train.mat, test.mat = subset(TrainTestMat, rownames(TrainTestMat) %ni% rownames(train.mat)), model.type = model.type, label = cls, kernelType = kernelType, iter = iter,
                                         cost = C, gamma = gamma, degree = degree, coef0 = coef0)
      pred.test <- multiAdaSampling(train.mat, test.mat = subset(testMat, rownames(testMat) %ni% train.ppis), model.type = model.type, label = cls, kernelType = kernelType, iter = iter,
                                    cost = C, gamma = gamma, degree = degree, coef0 = coef0)

      predTrainMat[[e]] <- pred.reference[[1]][, 1]
      predTrain.values.e[[e]] <- pred.reference[[2]]
      predTrain.model.e[[e]] <- pred.reference[[3]]

      predMat[[e]] <- pred.test[[1]][, 1]
      pred.values.e[[e]] <- pred.test[[2]]
      model.e[[e]] <- pred.test[[3]]

      TrainMat <- as.data.frame(TrainMat) %>% rownames_to_column("id") %>%
        left_join(as.data.frame(predTrainMat[[e]]) %>% rownames_to_column("id"), by = "id") %>%
        column_to_rownames("id") %>%
        as.matrix()
      colnames(TrainMat)[[length(assay)+e]] <- paste0("ensembleSize_",e)

      TestMat <- as.data.frame(TestMat) %>% rownames_to_column("id") %>%
        left_join(as.data.frame(predMat[[e]]) %>% rownames_to_column("id"), by = "id") %>%
        column_to_rownames("id") %>%
        as.matrix()
      colnames(TestMat)[[length(assay)+e]] <- paste0("ensembleSize_",e)

    }

  }

  if(verbose)
    base::message(ensembleSize_total, " independent training sets assembled each consisting of ", inclusion, " randomly sampled interactions. Each training set was used to predict positive interactions in the test set through multi-adaptive sampling in ", iter, " iteration(s).")

  predTrainMat <- base::rowMeans(TrainMat[,-c(1:length(assay))], na.rm = TRUE)
  predMat <- base::rowMeans(TestMat[,-c(1:length(assay))], na.rm = TRUE)
  TrainMatCount <- cbind(TrainMat[,c(1:length(assay))], base::rowSums(!is.na(TrainMat[,-c(1:length(assay))]), na.rm = TRUE), predTrainMat)
  TestMatCount <- cbind(TestMat[,c(1:length(assay))], base::rowSums(!is.na(TestMat[,-c(1:length(assay))]), na.rm = TRUE), predMat)

  #save probability score for train and test data together with the respective input data
  predTrainDf <- base::data.frame(predTrainMat) %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::separate(col = "sample", into = c("complex", "reference", "interaction", "sample", "orientation"), sep = ";") %>%
    dplyr::right_join(referenceSet %>%
                        tidyr::pivot_wider(names_from = "data", values_from = "score"),
                      by = c("complex", "reference", "interaction", "sample", "orientation"))
  predDf <- base::data.frame(predMat) %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::separate(col = "sample", into = c("complex", "interaction", "sample", "orientation"), sep = ";") %>%
    dplyr::right_join(PPIdf %>%
                        tidyr::pivot_wider(names_from = "data", values_from = "score"),
                      by = c("complex", "interaction", "sample", "orientation"))

  #save data to file

  return(list(predTrainDf = predTrainDf,
              predDf = predDf,
              TrainMatCount = TrainMatCount,
              TestMatCount = TestMatCount,
              pred.reference = pred.reference,
              pred.test = pred.test,
              predTrain.values.e = predTrain.values.e,
              predTrain.model.e = predTrain.model.e,
              pred.values.e = pred.values.e,
              training.sets = training.sets,
              negative.reference = negative.reference,
              model.type = model.type,
              model.e = model.e,
              predMat = predMat,
              testMat = testMat,
              trainMat = train.mat,
              label = cls,
              inclusion = inclusion, ensembleSize = ensembleSize, sampling = sampling,
              kernelType = kernelType, iter = iter, C = C, gamma = gamma, coef0 = coef0, degree = degree,
              top = top,
              seed = seed,
              system.time = Sys.time()))
}
