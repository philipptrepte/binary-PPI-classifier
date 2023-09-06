
<!-- README.md is generated from README.Rmd. Please edit that file -->

# binaryPPIclassifier

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

The binaryPPIclassifier package uses machine learning algorithms to
prioritize protein-protein interactions (PPIs) by analyzing quantitative
data from binary PPI assays and AlphaFold-Multimer predictions.

## Installation

You can install the development version of binaryPPIclassifier from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("philipptrepte/binary-PPI-classifier")
```

## Requirements

``` r
library(binaryPPIclassifier)

#load example data
load(file = "data/luthy_reference_sets.RData")
```

The reference set used for training the machine learning algorithm and
the PPI data frame that it is applied to needs to have the following
format and columns names:

<img src="man/figures/README-inputRequirement-1.png" width="100%" />

## Usage

### ppi.prediction()

The function `ppi.prediction()` is used to train a classifier machine
learning algorithms on a set of reference interactions `referenceSet`
and apply it to a test set PPI data set `PPIdf`. Please see
`?ppi.prediction` for help on additional parameters that need to be
specified.

`YOUR_PREDICTION <- ppi.prediction(PPIdf = YOUR_TEST_SET, referenceSet = YOUR_TRAINING_SET)`

### learning.curve()

The function `learning.curve()` is used to calculate and plot the
accuracy, hinge loss and binary cross-entropy loss from the results of
your `ppi.prediction()` result. Please see `?learning.curve` for help on
additional parameters that can to be specified.

`YOUR_LEARNING_CURVE <- learning.curve(ppi_prediction_result = YOUR_PREDICTION)`

The results are saved in a list and the plot can be accessed using
`YOUR_LEARNING_CURVE$learning_plot`

### recovery.rate()

The function `recovery.rate()` is used to determine fixed cut-offs using
a set of reference interactions `referenceSet` that can be applied to a
PPI data set `PPIdf` (optional) and the recovery rates are calculated
for both data sets. Please see `?recovery.rate` for help on additional
parameters that can to be specified.

`YOUR_RECOVERY_RATE <- recovery.rate(PPIdf = YOUR_TEST_SET, referenceSet = YOUR_TRAINING_SET)`

## Vignette

Please see the Vignette (under development) for a detailed description.

## Application

Development, benchmarking and example applications have been published
here:

<u>**AI-guided pipeline for protein-protein interaction drug discovery
identifies a SARS-CoV-2 inhibitor**</u>

Philipp Trepte#, Christopher Secker#, Simona Kostova, Sibusiso
B. Maseko, Soon Gang Choi, Jeremy Blavier, Igor Minia, Eduardo,
Silva Ramos, Patricia Cassonnet, Sabrina Golusik, Martina Zenkner, Stephanie Beetz, Mara
J. Liebich, Nadine Scharek, Anja Schütz, MarcelSperling, Michael Lisurek, Yang Wang, Kerstin Spirohn, Tong Hao, Michael
A. Calderwood, David
E. Hill, Markus Landthaler, Julien Olivet, Jean-Claude Twizere, Marc Vidal, Erich
E. Wanker

bioRxiv 2023.06.14.544560; doi: <https://doi.org/10.1101/2023.06.14.544560>

## License

Distributed under the MIT License. See `License.md` for more
information.

## Contact

Philipp Trepte - <philipp.trepte@imba.oeaw.ac.at> -
[LinkedIn](https://www.linkedin.com/in/philipp-trepte/)

binaryPPIclassifier: <https://github.com/philipptrepte/binary-PPI-classifier>
