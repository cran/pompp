## usethis namespace: start
#' @useDynLib pompp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("pompp", libpath)
}
