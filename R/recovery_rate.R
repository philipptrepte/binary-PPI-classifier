usethis::use_package('dplyr')
usethis::use_package('tidyr')
usethis::use_package('tibble')
usethis::use_package('stringr')
usethis::use_package('rlang')
usethis::use_package('Rmisc')

#' Function to calculate the recovery rate from binary PPI assays.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import rlang
#' @import Rmisc
#'
#' @param referenceSet: reference PPI data set used to determine cutoffs to calculate recovery rates
#' @param PPIdf: optional binary PPI data set containing interactions classified using the cutoffs from the reference sets
#' @param assay: assay in the referenceSet and PPIdf to calculate recovery rates
#' @param construct.scaling: scaling of constructs by "median.normalization" or "robust.scaler"
#' @param range: probs in stats::quantile(x, ...): numeric vector of probabilities
#' @param iter.scaler: if TRUE and when using "robust.scaler" it iteratively performs robust scaler normalization until the IQR of each construct is within the IQR of all loaded data sets
#' @param cutoff: cutoffs are calculated based on "all" interactions, or for each "construct" or for each "orientation"
#' @param fns: function used to calculate the cutoffs: "median", "mean", "CI" or "max"
#' @param n.sd: number of sd used to calculate the cutoffs for functions "median" and "mean"
#' @param n.CI: confidence intervall used to calculate the cutoff for function "CI"
#' @param negative.reference: string in the column "complex" to specify the negative/random interactions
#' @param complex.parameters: additional grouping parameters from which the recovery rates are calculated
#' @param verbose: give detailed information
#'
#' @return
#' @export
#'
#' @examples
recovery.rate <- function(referenceSet = NULL,
                          PPIdf = NULL,
                          assay = "mean_cBRET",
                          construct.scaling = NULL,
                          range = c(0.25, 0.75), iter.scaler = TRUE,
                          cutoff = "all",
                          fns = "median",
                          n.sd = 1,
                          n.CI = 0.99,
                          negative.reference = "RRS",
                          complex.parameters = NULL,
                          verbose = TRUE) {
  `%ni%` <- Negate(`%in%`)

  #use provided reference set data if needed
  if(is.null(referenceSet)) {
    if(verbose)
      base::message("No user reference set data provided. Published reference data from Trepte et al. is used.")
    data("luthy_reference_set")

    referenceSet <- luthy_reference_sets %>%
      dplyr::filter(data %in% assay) %>%
      dplyr::mutate(reference = base::ifelse(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|")), "RRS", "PRS"))
  }

  if(!any(stringr::str_detect(base::colnames(referenceSet), "data"))) {
    stop("no column 'data' present in your 'referenceSet'. Must contain the 'assay' that you want to analyze.")
  }
  if(assay %ni% referenceSet$data) {
    stop("specify 'assay' that can be found in the 'data' column of your 'referenceSet' ", call. = FALSE)
  }
  if(cutoff %ni% c("construct", "orientation", "all")) {
    stop("'cutoff' is not specified. Choose between 'construct' and 'configuration'.")
  }

  if(!base::exists("negative.reference")){
    negative.reference = "."
  }
  if(is.null(negative.reference)){
    negative.reference = "."
  }
  if(!base::exists("data.scaling")){
    data.scaling <- "no"
  }
  if(is.null(data.scaling)) {
    data.scaling <- "no"
  }

  referenceSet <- referenceSet %>%
    dplyr::filter(data == assay[1]) %>%
    tidyr::pivot_wider(names_from = "data", values_from = "score")

  if(!is.null(PPIdf)) {
    PPIdf <- PPIdf %>%
      dplyr::filter(data == assay[1]) %>%
      tidyr::pivot_wider(names_from = "data", values_from = "score")
  }

  if(!is.null(complex.parameters)) {
    if(!any(stringr::str_detect(base::colnames(referenceSet), complex.parameters))) {
      stop("'complex.parameters': '", complex.parameters, "' was not found in the provided dataset.")
    }
  }

  if(is.null(construct.scaling)) {
    base::message("'Donor' and 'Acceptor' constructs are not normalized.")
    construct.scaling <- "none"
  }

  if(construct.scaling == "median.normalization") {
    base::message("'Donor' and 'Acceptor' constructs are median normalized.")
    median.global <- stats::median(referenceSet %>% dplyr::pull(assay[1]), na.rm = TRUE)

    donor.cf <- referenceSet %>%
      dplyr::group_by(Donor) %>%
      dplyr::summarise(median = rlang::exec(.fn = "median", !!(rlang::sym(assay[1])), na.rm = TRUE)) %>%
      dplyr::mutate(cf.donor = median - median.global)

    referenceSet <- referenceSet %>%
      dplyr::left_join(donor.cf %>% dplyr::select(-median), by = "Donor") %>%
      dplyr::mutate(!!(rlang::sym(assay[1])) := !!(rlang::sym(assay[1])) - cf.donor, .keep = "unused")

    acceptor.cf <- referenceSet %>%
      dplyr::group_by(Acceptor) %>%
      dplyr::summarise(median = rlang::exec(.fn = "median", !!(rlang::sym(assay[1])), na.rm = TRUE)) %>%
      dplyr::mutate(cf.acceptor = median - median.global)

    referenceSet <- referenceSet %>%
      dplyr::left_join(acceptor.cf %>% dplyr::select(-median), by = "Acceptor") %>%
      dplyr::mutate(!!(rlang::sym(assay[1])) := !!(rlang::sym(assay[1])) - cf.acceptor, .keep = "unused")

    if(!is.null(PPIdf)) {
      median.global <- stats::median(PPIdf %>% dplyr::pull(assay[1]), na.rm = TRUE)
      donor.cf <- PPIdf %>%
        dplyr::group_by(Donor) %>%
        dplyr::summarise(median = rlang::exec(.fn = "median", !!(rlang::sym(assay[1])), na.rm = TRUE)) %>%
        dplyr::mutate(cf.donor = median - median.global)

      PPIdf <- PPIdf %>%
        dplyr::left_join(donor.cf %>% dplyr::select(-median), by = "Donor") %>%
        dplyr::mutate(!!(rlang::sym(assay[1])) := !!(rlang::sym(assay[1])) - cf.donor, .keep = "unused")

      acceptor.cf <- PPIdf %>%
        dplyr::group_by(Acceptor) %>%
        dplyr::summarise(median = rlang::exec(.fn = "median", !!(rlang::sym(assay[1])), na.rm = TRUE)) %>%
        dplyr::mutate(cf.acceptor = median - median.global)

      PPIdf <- PPIdf %>%
        dplyr::left_join(acceptor.cf %>% dplyr::select(-median), by = "Acceptor") %>%
        dplyr::mutate(!!(rlang::sym(assay[1])) := !!(rlang::sym(assay[1])) - cf.acceptor, .keep = "unused")
    }
  }

  if(construct.scaling == "robust.scaler") {
    IQR.custom <- function(x, na.rm = FALSE, type = 7, range = range)  {
      base::diff(stats::quantile(as.numeric(x), probs = range, na.rm = na.rm, names = FALSE,
                    type = type))
    }

    scaler.global <- stats::median(referenceSet %>% dplyr::pull(assay[1]),na.rm = TRUE)
    IQR.global <- IQR.custom(referenceSet %>% dplyr::pull(assay[1]), na.rm = TRUE, range = range)

    if(!is.null(PPIdf)) {
      scaler.global <- stats::median(c(referenceSet %>% dplyr::pull(assay[1]), PPIdf %>% dplyr::pull(assay[1])), na.rm = TRUE)
      IQR.global <- IQR.custom(c(referenceSet %>% dplyr::pull(assay[1]), PPIdf %>% dplyr::pull(assay[1])), na.rm = TRUE, range = range)
    }

    donor.n.reference <- referenceSet %>% dplyr::group_by(Donor) %>% dplyr::count()
    acceptor.n.reference <- referenceSet %>% dplyr::group_by(Acceptor) %>% dplyr::count()

    if(any(donor.n.reference$n >= 20) == FALSE | any(acceptor.n.reference$n >= 20) == FALSE) {
      base::message("Not enough 'Donor' and 'Acceptor' constructs in the provided 'reference.set' to perform robust scaler.")
      referenceSet.scaler <- base::data.frame(base::matrix(nrow = 0, ncol = bae::ncol(referenceSet)))
      base::colnames(referenceSet.scaler) <- base::colnames(referenceSet)
    }

    if(any(donor.n.reference$n >= 20) == TRUE | any(acceptor.n.reference$n >= 20) == TRUE){

      donor.reference.constructs <- donor.n.reference %>% dplyr::filter(n >= 20) %>% dplyr::pull(Donor)
      acceptor.reference.constructs <- acceptor.n.reference %>% dplyr::filter(n >= 20) %>% dplyr::pull(Acceptor)

      referenceSet.scaler <- referenceSet %>%
        dplyr::filter(Donor %in% donor.reference.constructs & Acceptor %in% acceptor.reference.constructs)

      scaler <- TRUE
      count <- 0
      while(scaler) {
        count <- count + 1
        donor.cf <- referenceSet.scaler %>%
          dplyr::group_by(Donor) %>%
          dplyr::summarise(median = stats::median(!!(rlang::sym(assay[1])), na.rm = TRUE),
                           IQR = IQR.custom(!!(rlang::sym(assay[1])), na.rm = TRUE, range = range)) %>%
          dplyr::mutate(cf.donor = median - scaler.global)
        referenceSet.scaler <- referenceSet.scaler %>%
          dplyr::left_join(donor.cf %>% dplyr::select(Donor, IQR, cf.donor), by = "Donor") %>%
          dplyr::rowwise() %>%
          dplyr::mutate(!!(rlang::sym(assay[1])) := (!!(rlang::sym(assay[1])) - cf.donor) / base::ifelse(IQR > IQR.global, (IQR / IQR.global), 1), .keep = "unused")
        acceptor.cf <- referenceSet.scaler %>%
          dplyr::group_by(Acceptor) %>%
          dplyr::summarise(median = stats::median(!!(rlang::sym(assay[1])), na.rm = TRUE),
                           IQR = IQR.custom(!!(rlang::sym(assay[1])), na.rm = TRUE, range = range)) %>%
          dplyr::mutate(cf.acceptor = median - scaler.global)
        referenceSet.scaler <- referenceSet.scaler %>%
          dplyr::left_join(acceptor.cf %>% dplyr::select(Acceptor, IQR, cf.acceptor), by = "Acceptor") %>%
          dplyr::mutate(!!(rlang::sym(assay[1])) := (!!(rlang::sym(assay[1])) - cf.acceptor) / base::ifelse(IQR > IQR.global, (IQR / IQR.global), 1), .keep = "unused")

        constructs <- referenceSet.scaler %>%
          dplyr::group_by(construct = Donor) %>%
          dplyr::summarise(median = stats::median(!!(rlang::sym(assay[1])), na.rm = TRUE),
                           IQR = IQR.custom(!!(rlang::sym(assay[1])), na.rm = TRUE, range = range)) %>%
          base::rbind(
            referenceSet.scaler %>%
              dplyr::group_by(construct = Acceptor) %>%
              dplyr::summarise(median = stats::median(!!(rlang::sym(assay[1])), na.rm = TRUE),
                               IQR = IQR.custom(!!(rlang::sym(assay[1])), na.rm = TRUE, range = range))
          ) %>%
          dplyr::mutate(scaler.IQR = base::ifelse(round(IQR, 4) > round(IQR.global,4), TRUE, FALSE),
                 scaler.median = base::ifelse(round(median, 4) > round(scaler.global,4), TRUE, FALSE))

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
      base::message("Construct-specific median normalization was done for ", length(base::which(donor.n.reference$n >= 20)), " 'Donor construcst and ", length(base::which(acceptor.n.reference$n >= 20)), " ' Acceptor' construcsts. \n",
              "Construct-specific IQR normalization was done for ", length(base::which(donor.cf$IQR >= IQR.global)), " 'Donor' constructs and ", length(base::which(acceptor.cf$IQR >= IQR.global)), " 'Acceptor' constructs that showed a higher IQR than the global IQR of the 'reference.Set' in ", count, " iteration(s).")

    }

    if(any(donor.n.reference$n < 20) == FALSE & any(acceptor.n.reference$n < 20) == FALSE) {
      referenceSet.notScaler <- base::data.frame(base::matrix(nrow = 0, ncol = bae::ncol(referenceSet)))
      base::colnames(referenceSet.notScaler) <- base::colnames(referenceSet)
    }

    if(any(donor.n.reference$n < 20) == TRUE & any(acceptor.n.reference$n < 20) == TRUE) {
      base::message("There are ", length(base::which(donor.n.reference$n < 20)), " 'Donor' constructs with not enough interactions (<20) to perform 'Robust scaler' normalization. \n",
              "There are ", length(base::which(donor.n.reference$n < 20)), " 'Acceptor' constructs with not enough interactions (<20) to perform 'Robust scaler' normalization. \n")

      donor.reference.constructs <- donor.n.reference %>% dplyr::filter(n < 20) %>% dplyr::pull(Donor)
      acceptor.reference.constructs <- acceptor.n.reference %>% dplyr::filter(n < 20) %>% dplyr::pull(Acceptor)

      if(length(donor.reference.constructs) != 0 & length(acceptor.reference.constructs) != 0) {
        referenceSet.notScaler <- referenceSet %>%
          dplyr::filter(Donor %in% donor.reference.constructs & Acceptor %in% acceptor.reference.constructs)
      }
    }

    referenceSet <- referenceSet.scaler %>% base::rbind(referenceSet.notScaler) %>% dplyr::ungroup()


    if(!is.null(PPIdf)) {
      donor.n.PPIdf <- PPIdf %>% dplyr::group_by(Donor) %>% dplyr::count()
      acceptor.n.PPIdf <- PPIdf %>% dplyr::group_by(Acceptor) %>% dplyr::count()

      if(any(donor.n.PPIdf$n >= 20) == FALSE | any(acceptor.n.PPIdf$n >= 20) == FALSE) {
        base::message("Not enough 'Donor' and 'Acceptor' constructs in the provided 'PPIdf.set' to perform robust scaler.")
        PPIdf.scaler <- base::data.frame(base::matrix(nrow = 0, ncol = bae::ncol(PPIdf)))
        base::colnames(PPIdf.scaler) <- base::colnames(PPIdf)
      }

      if(any(donor.n.PPIdf$n >= 20) == TRUE | any(acceptor.n.PPIdf$n >= 20) == TRUE){

        donor.PPIdf.constructs <- donor.n.PPIdf %>% dplyr::filter(n >= 20) %>% dplyr::pull(Donor)
        acceptor.PPIdf.constructs <- acceptor.n.PPIdf %>% dplyr::filter(n >= 20) %>% dplyr::pull(Acceptor)

        PPIdf.scaler <- PPIdf %>%
          dplyr::filter(Donor %in% donor.PPIdf.constructs & Acceptor %in% acceptor.PPIdf.constructs) %>%
          dplyr::mutate(scaling = "median.scaled")

        scaler <- TRUE
        count <- 0
        while(scaler) {
          count <- count + 1

          donor.cf <- PPIdf.scaler %>%
            dplyr::group_by(Donor) %>%
            dplyr::summarise(median = stats::median(!!(rlang::sym(assay[1])), na.rm = TRUE),
                             IQR = IQR.custom(!!(rlang::sym(assay[1])), na.rm = TRUE, range = range)) %>%
            dplyr::mutate(cf.donor = median - scaler.global)
          PPIdf.scaler <- PPIdf.scaler %>%
            dplyr::left_join(donor.cf %>% dplyr::select(Donor, IQR, cf.donor), by = "Donor") %>%
            dplyr::rowwise() %>%
            dplyr::mutate(!!(rlang::sym(assay[1])) := (!!(rlang::sym(assay[1])) - cf.donor) / base::ifelse(IQR > IQR.global, (IQR / IQR.global), 1),
                          scaling = base::ifelse(IQR > IQR.global, "IQR.scaled", scaling), .keep = "unused")
          acceptor.cf <- PPIdf.scaler %>%
            dplyr::group_by(Acceptor) %>%
            dplyr::summarise(median = stats::median(!!(rlang::sym(assay[1])), na.rm = TRUE),
                             IQR = IQR.custom(!!(rlang::sym(assay[1])), na.rm = TRUE, range = range)) %>%
            dplyr::mutate(cf.acceptor = median - scaler.global)
          PPIdf.scaler <- PPIdf.scaler %>%
            dplyr::left_join(acceptor.cf %>% dplyr::select(Acceptor, IQR, cf.acceptor), by = "Acceptor") %>%
            dplyr::mutate(!!(rlang::sym(assay[1])) := (!!(rlang::sym(assay[1])) - cf.acceptor) / base::ifelse(IQR > IQR.global, (IQR / IQR.global), 1),
                          scaling = base::ifelse(IQR > IQR.global, "IQR.scaled", scaling), .keep = "unused")

          constructs <- PPIdf.scaler %>%
            dplyr::group_by(construct = Donor) %>%
            dplyr::summarise(median = stats::median(!!(rlang::sym(assay[1])), na.rm = TRUE),
                             IQR = IQR.custom(!!(rlang::sym(assay[1])), na.rm = TRUE, range = range)) %>%
            base::rbind(
              PPIdf.scaler %>%
                dplyr::group_by(construct = Acceptor) %>%
                dplyr::summarise(median = stats::median(!!(rlang::sym(assay[1])), na.rm = TRUE),
                                 IQR = IQR.custom(!!(rlang::sym(assay[1])), na.rm = TRUE, range = range))
            ) %>%
            dplyr::mutate(scaler.IQR = base::ifelse(round(IQR, 4) > round(IQR.global,4), TRUE, FALSE),
                   scaler.median = base::ifelse(round(median, 4) > round(scaler.global,4), TRUE, FALSE))

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
        base::message("Construct-specific median normalization is done for ", length(base::which(donor.n.PPIdf$n >= 20)), " 'Donor construcst and ", length(base::which(acceptor.n.PPIdf$n >= 20)), " ' Acceptor' construcsts. \n",
                "Construct-specific IQR normalization is done for ", length(base::which(donor.cf$IQR >= IQR.global)), " 'Donor' constructs and ", length(base::which(acceptor.cf$IQR >= IQR.global)), " 'Acceptor' constructs that showed a higher IQR than the global IQR of the 'PPIdf.Set' in ", count, " iteration(s).")

      }

      if(any(donor.n.PPIdf$n < 20) == FALSE & any(acceptor.n.PPIdf$n < 20) == FALSE) {
        PPIdf.notScaler <- base::data.frame(base::matrix(nrow = 0, ncol = bae::ncol(PPIdf)))
        base::colnames(PPIdf.notScaler) <- base::colnames(PPIdf)
      }

      if(any(donor.n.PPIdf$n < 20) == TRUE & any(acceptor.n.PPIdf$n < 20) == TRUE) {
        base::message("There are ", length(base::which(donor.n.PPIdf$n < 20)), " 'Donor' constructs with not enough interactions (<20) to perform 'Robust scaler' normalization. \n",
                "There are ", length(base::which(acceptor.n.PPIdf$n < 20)), " 'Acceptor' constructs with not enough interactions (<20) to perform 'Robust scaler' normalization. \n")

        donor.PPIdf.constructs <- donor.n.PPIdf %>% dplyr::filter(n < 20) %>% dplyr::pull(Donor)
        acceptor.PPIdf.constructs <- acceptor.n.PPIdf %>% dplyr::filter(n < 20) %>% dplyr::pull(Acceptor)

        if(length(donor.PPIdf.constructs) != 0 & length(acceptor.PPIdf.constructs) != 0) {
          PPIdf.notScaler <- PPIdf %>%
            dplyr::filter(Donor %in% donor.PPIdf.constructs & Acceptor %in% acceptor.PPIdf.constructs) %>%
            dplyr::mutate(scaling = "not.scaled")
        }
      }

      PPIdf <- PPIdf.scaler %>% base::rbind(PPIdf.notScaler) %>% dplyr::ungroup()

    }
  }

  if(cutoff == "all") {
    if(negative.reference == ".") {
      if(verbose) {
        base::message("Cut-offs are calculated based on all interactions in the 'referenceSet'.")
      }
    }
    else {
      if(verbose) {
        base::message("Cut-offs are calculated based on the provided 'negative.reference' in the 'referenceSet'.")
      }
    }

    if(fns %in% c("median", "mean")) {
      if(verbose) {
        base::message("Cutoffs are based on the ", fns, " + ", n.sd, " standard deviations.")
      }
      if(!base::exists("n.sd") | is.null(n.sd)) {
        stop("'n.sd' must be defined.")
      }
      cutoff.all <- rlang::exec(.fn = fns,
                         referenceSet %>%
                           dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|"))) %>%
                           dplyr::pull(rlang::sym(assay[1])), na.rm = TRUE) +
        n.sd*rlang::exec(.fn = sd, referenceSet %>%
                    dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|"))) %>%
                    dplyr::pull(rlang::sym(assay[1])), na.rm = TRUE)
      referenceSet <- referenceSet %>%
        dplyr::mutate(cutoff = cutoff.all,
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }

    if(fns == "CI") {
      if(verbose) {
        base::message("Cutoffs are based on the ", n.CI*100, "% confidence intervall.")
      }
      cutoff.all <- rlang::exec(.fn = fns, referenceSet %>%
                           dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|"))  & !is.na(!!(rlang::sym(assay[1])))) %>%
                           dplyr::pull(rlang::sym(assay[1])), ci = n.CI)[1]
      referenceSet <- referenceSet %>%
        dplyr::mutate(cutoff = cutoff.all,
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }

    if(fns == "max") {
      if(negative.reference == "."){
        stop("'negative.reference' is not specified. Must be provided to calculate cut-offs that are based on the interaction with the maximumn score in the respective 'negative.reference'.")
      }
      if(verbose) {
        base::message("Cutoffs are based on the highest interaction score per 'Donor' or 'Acceptor' construct in the 'negative.reference'." )
      }
      cutoff.all <- rlang::exec(.fn = max, referenceSet %>%
                           dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|"))  & !is.na(!!(rlang::sym(assay[1])))) %>%
                           dplyr::pull(rlang::sym(assay[1])))
      referenceSet <- referenceSet %>%
        dplyr::mutate(cutoff = cutoff.all,
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }

    if(!is.null(PPIdf)) {
      PPIdf <- PPIdf %>%
        dplyr::mutate(cutoff = cutoff.all,
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }
  }

  if(cutoff == "construct") {
    if(!is.null(PPIdf)) {
      donor.reference <- referenceSet %>% dplyr::pull(Donor) %>% base::unique()
      acceptor.reference <- referenceSet %>% dplyr::pull(Acceptor) %>% base::unique()
      if(!any(PPIdf$Donor %in% donor.reference)) {
        stop("The 'referenceSet' and 'PPIdf' do not share common 'Donor' constructs.")
      }
      if(!any(PPIdf$Acceptor %in% acceptor.reference)) {
        stop("The 'referenceSet' and 'PPIdf' do not share common 'Acceptor' constructs.")
      }
    }
    donor.n <- min(referenceSet %>% dplyr::group_by(Donor) %>% dplyr::count() %>% dplyr::pull(n), na.rm = TRUE) <= 20
    acceptor.n <- min(referenceSet %>% dplyr::group_by(Acceptor) %>% dplyr::count() %>% dplyr::pull(n), na.rm = TRUE) <= 20
    if(donor.n == TRUE | acceptor.n == TRUE) {
      stop("At least one 'Donor' or 'Acceptor' construct is not tested in enough interactions (>20) to calculate construct-specific cutoffs.")
    }

    if(negative.reference == ".") {
      base::message("Construct specific cut-offs are calculated based on all interactions for the respective construct.")
    }
    else {
      base::message("Construct specific cut-offs are calculated based on the provided 'negative.reference'.")
    }

    if(fns %in% c("median", "mean")) {
      if(verbose) {
        base::message("Construct-specific cutoffs are based on the ", fns, " + ", n.sd, " standard deviations.")
      }
      if(!base::exists("n.sd") | is.null(n.sd)) {
        stop("'n.sd' must be defined.")
      }
      donor.cutoff <- referenceSet %>%
        dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|"))) %>%
        dplyr::group_by(Donor) %>%
        dplyr::summarise(!!(rlang::sym(fns)) := rlang::exec(.fn = fns, !!(rlang::sym(assay[1])), na.rm = TRUE),
                         sd = stats::sd(!!(rlang::sym(assay[1])), na.rm = TRUE)) %>%
        dplyr::mutate(donor.cutoff = !!(rlang::sym(fns)) + n.sd*sd, .keep = "unused")
      acceptor.cutoff <- referenceSet %>%
        dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|"))) %>%
        dplyr::group_by(Acceptor) %>%
        dplyr::summarise(!!(rlang::sym(fns)) := rlang::exec(.fn = fns, !!(rlang::sym(assay[1])), na.rm = TRUE),
                         sd = stats::sd(!!(rlang::sym(assay[1])), na.rm = TRUE)) %>%
        dplyr::mutate(acceptor.cutoff = !!(rlang::sym(fns)) + n.sd*sd, .keep = "unused")
      referenceSet <- referenceSet %>%
        dplyr::left_join(donor.cutoff, by = "Donor") %>%
        dplyr::left_join(acceptor.cutoff, by = "Acceptor") %>%
        dplyr::mutate(donor.cutoff = base::ifelse(is.na(donor.cutoff), acceptor.cutoff, donor.cutoff),
                      acceptor.cutoff = base::ifelse(is.na(acceptor.cutoff), donor.cutoff, acceptor.cutoff),
                      cutoff = base::ifelse(is.na(donor.cutoff) & is.na(acceptor.cutoff), NA,
                                      base::ifelse(donor.cutoff > acceptor.cutoff, donor.cutoff, acceptor.cutoff)),
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }

    if(fns == "CI") {
      base::message("Construct-specific cutoffs are based on the ", n.CI*100, "% confidence intervall.")
      donor.cutoff <- referenceSet %>%
        dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|")) & !is.na(!!(rlang::sym(assay[1])))) %>%
        dplyr::group_by(Donor) %>%
        dplyr::summarise(donor.cutoff = Rmisc::CI(!!(rlang::sym(assay[1])), ci = n.CI)[1])
      acceptor.cutoff <- referenceSet %>%
        dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|")) & !is.na(!!(rlang::sym(assay[1])))) %>%
        dplyr::group_by(Acceptor) %>%
        dplyr::summarise(acceptor.cutoff = Rmisc::CI(!!(rlang::sym(assay[1])), ci = n.CI)[1])
      referenceSet <- referenceSet %>%
        dplyr::left_join(donor.cutoff, by = "Donor") %>%
        dplyr::left_join(acceptor.cutoff, by = "Acceptor") %>%
        dplyr::mutate(donor.cutoff = base::ifelse(is.na(donor.cutoff), acceptor.cutoff, donor.cutoff),
                      acceptor.cutoff = base::ifelse(is.na(acceptor.cutoff), donor.cutoff, acceptor.cutoff),
                      cutoff = base::ifelse(is.na(donor.cutoff) & is.na(acceptor.cutoff), NA,
                                      base::ifelse(donor.cutoff > acceptor.cutoff, donor.cutoff, acceptor.cutoff)),
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }

    if(fns == "max") {
      if(negative.reference == "."){
        stop("'negative.reference' is not specified. Must be provided to calculate cut-offs that are based on the interaction with the maximumn score in the respective 'negative.reference'.")
      }
      base::message("Construct-specific cutoffs are based on the highest interaction score per 'Donor' or 'Acceptor' construct in the 'negative.reference'." )
      donor.cutoff <- referenceSet %>%
        dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|"))) %>%
        dplyr::group_by(Donor) %>%
        dplyr::summarise(donor.cutoff = max(!!(rlang::sym(assay[1])), na.rm = TRUE))
      acceptor.cutoff <- referenceSet %>%
        dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|"))) %>%
        dplyr::group_by(Acceptor) %>%
        dplyr::summarise(acceptor.cutoff = max(!!(rlang::sym(assay[1])), na.rm = TRUE))
      referenceSet <- referenceSet %>%
        dplyr::left_join(donor.cutoff, by = "Donor") %>%
        dplyr::left_join(acceptor.cutoff, by = "Acceptor") %>%
        dplyr::mutate(donor.cutoff = base::ifelse(is.na(donor.cutoff), acceptor.cutoff, donor.cutoff),
                      acceptor.cutoff = base::ifelse(is.na(acceptor.cutoff), donor.cutoff, acceptor.cutoff),
                      cutoff = base::ifelse(is.na(donor.cutoff) & is.na(acceptor.cutoff), NA,
                                      base::ifelse(donor.cutoff > acceptor.cutoff, donor.cutoff, acceptor.cutoff)),
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }

    if(!is.null(PPIdf)) {
      PPIdf <- PPIdf %>%
        dplyr::left_join(donor.cutoff, by = "Donor") %>%
        dplyr::left_join(acceptor.cutoff, by = "Acceptor") %>%
        dplyr::mutate(donor.cutoff = base::ifelse(is.na(donor.cutoff), acceptor.cutoff, donor.cutoff),
                      acceptor.cutoff = base::ifelse(is.na(acceptor.cutoff), donor.cutoff, acceptor.cutoff),
                      cutoff = base::ifelse(is.na(donor.cutoff) & is.na(acceptor.cutoff), NA,
                                      base::ifelse(donor.cutoff > acceptor.cutoff, donor.cutoff, acceptor.cutoff)),
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }
  }

  if(cutoff == "orientation") {
    if(fns %in% c("median", "mean")) {
      if(verbose) {
        base::message("Orientation-specific cutoffs are based on the ", fns, " + ", n.sd, " standard deviations.")
      }
      if(!base::exists("n.sd")) {
        stop("'n.sd' must be defined.")
      }
      if(!is.numeric(n.sd)){
        stop("'n.sd' must be numeric")
      }
      orientation.cutoff <- referenceSet %>%
        dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|"))) %>%
        dplyr::group_by(orientation) %>%
        dplyr::summarise(!!(rlang::sym(fns)) := rlang::exec(.fn = fns, !!(rlang::sym(assay[1])), na.rm = TRUE),
                         sd = stats::sd(!!(rlang::sym(assay[1])), na.rm = TRUE)) %>%
        dplyr::mutate(orientation.cutoff = !!(rlang::sym(fns)) + n.sd*sd, .keep = "unused")
      referenceSet <- referenceSet %>%
        dplyr::left_join(orientation.cutoff, by = "orientation") %>%
        dplyr::mutate(cutoff = base::ifelse(is.na(orientation.cutoff), NA, orientation.cutoff),
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
      if(!is.null(PPIdf)) {
        PPIdf <- PPIdf %>%
          dplyr::left_join(orientation.cutoff, by = "orientation") %>%
          dplyr::mutate(cutoff = base::ifelse(is.na(orientation.cutoff), 0, orientation.cutoff),
                        pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                        pos = base::ifelse(is.finite(pos), pos, NA))
      }
    }

    if(fns == "CI") {
      if(verbose) {
        base::message("Orientation-specific cutoffs are based on the ", n.CI*100, "% confidence intervall.")
      }
      orientation.cutoff <- referenceSet %>%
        dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|")) & !is.na(!!(rlang::sym(assay[1])))) %>%
        dplyr::group_by(orientation) %>%
        dplyr::summarise(orientation.cutoff = Rmisc::CI(!!(rlang::sym(assay[1])), ci = n.CI)[1])
      referenceSet <- referenceSet %>%
        dplyr::left_join(orientation.cutoff, by = "orientation") %>%
        dplyr::mutate(cutoff = base::ifelse(is.na(orientation.cutoff), NA, orientation.cutoff),
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
      if(!is.null(PPIdf)) {
        PPIdf <- PPIdf %>%
          dplyr::left_join(orientation.cutoff, by = "orientation") %>%
          dplyr::mutate(cutoff = base::ifelse(is.na(orientation.cutoff), NA, orientation.cutoff),
                        pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                        pos = base::ifelse(is.finite(pos), pos, NA))
      }
    }

    if(fns == "max") {
      if(negative.reference == "."){
        stop("'negative.reference' is not specified. Must be provided to calculate cut-offs that are based on the interaction with the maximumn score in the respective 'negative.reference'.")
      }
      if(verbose) {
        base::message("Orientation-specific cutoffs are based on the highest interaction score per 'Donor' or 'Acceptor' construct in the 'negative.reference'." )
      }
      orientation.cutoff <- referenceSet %>%
        dplyr::filter(stringr::str_detect(complex, base::paste(negative.reference, collapse = "|")) & !is.na(!!(rlang::sym(assay[1])))) %>%
        dplyr::group_by(orientation) %>%
        dplyr::summarise(orientation.cutoff = max(!!(rlang::sym(assay[1])), na.rm = TRUE))
      referenceSet <- referenceSet %>%
        dplyr::left_join(orientation.cutoff, by = "orientation") %>%
        dplyr::mutate(cutoff = base::ifelse(is.na(orientation.cutoff), NA, orientation.cutoff),
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }

    if(!is.null(PPIdf)) {
      PPIdf <- PPIdf %>%
        dplyr::left_join(orientation.cutoff, by = "orientation") %>%
        dplyr::mutate(cutoff = base::ifelse(is.na(orientation.cutoff), NA, orientation.cutoff),
                      pos = base::ifelse(!!(rlang::sym(assay[1])) > cutoff, 1, 0),
                      pos = base::ifelse(is.finite(pos), pos, NA))
    }
  }

  if(base::exists("complex.parameters")) {
    grouping.parameters <- c("interaction", "complex", complex.parameters)
  }
  else {
    grouping.parameters <- base::paste(c("interaction", "complex"))
  }

  recovery.reference <- referenceSet %>%
    dplyr::group_by(across(.cols = grouping.parameters)) %>%
    dplyr::summarise(pos = max(pos, na.rm = TRUE)) %>%
    dplyr::group_by(across(.cols = grouping.parameters[-1])) %>%
    dplyr::mutate(pos = base::ifelse(is.finite(pos), pos, 0)) %>%
    dplyr::summarise(total = dplyr::n(),
                     pos = sum(pos, na.rm = TRUE),
                     recovery = round(pos/total*100, 1))

  if(!is.null(PPIdf)) {
    recovery.PPIdf <- PPIdf %>%
      dplyr::group_by(across(.cols = grouping.parameters)) %>%
      dplyr::summarise(pos = max(pos, na.rm = TRUE)) %>%
      dplyr::group_by(across(.cols = grouping.parameters[-1])) %>%
      dplyr::summarise(total = dplyr::n(),
                       pos = sum(pos, na.rm = TRUE),
                       recovery = round(pos/total*100, 1))
  }
  else {
    recovery.PPIdf <- c()
    PPIdf <- c()
  }


  return(list(recovery.reference = recovery.reference,
              recovery.PPIdf = recovery.PPIdf,
              referenceSet = referenceSet,
              PPIdf = PPIdf,
              fns = fns,
              cutoff = cutoff,
              negative.reference = negative.reference,
              system.time = Sys.time()
              ))
}
