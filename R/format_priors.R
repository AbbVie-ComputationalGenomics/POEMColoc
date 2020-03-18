#' format_priors
#'
#' @param dataset_list
#' list of datasets to be analyzed of length S
#' @param config
#' list of vectors of indices 1:S of length R
#' @param p1
#' vector of probability of length S or 1
#' @param p12
#' vector of probability of length R or 1
#' @return
#' list with elements
#' p1 vector of probability of length S
#' p12 vector of probability of length R

format_priors <- function(dataset_list, config, p1, p12) {
  if (!(length(p1)==1 | length(p1) == length(dataset_list))) {
    stop("p1 should be either the length of dataset_list or 1 if the same for all datasets")
  }
  if (!(length(p12)==1 | length(p12) == length(config))) {
    stop("p12 should either be length of config 1 if the same for all datasets")
  }
  if (length(p1)==1) {
    p1 <- rep(p1, length(dataset_list))
  }
  if (length(p12)==1) {
    p12 <- rep(p12, length(config))
  }
  return(list(p1=p1, p12=p12))
}
