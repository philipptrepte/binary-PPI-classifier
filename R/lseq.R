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

