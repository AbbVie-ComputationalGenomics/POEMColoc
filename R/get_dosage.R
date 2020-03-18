#' wrapper for seqGetData
#'
#' @param gds
#' gds file connection
#' @return
#' matrix with rows corresponding to samples and columns to positions
#'
#' @examples
get_dosage <- function(gds) {
  genotype = seqGetData(gds, "$dosage")
  return(genotype)
}
