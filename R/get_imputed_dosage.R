#' wrapper for seqGetData
#'
#' @param gds
#' gds file connection
#' @return
#' matrix with rows corresponding to samples and columns to positions
#'
#' @examples
get_imputed_dosage <- function(gds) {
  genotype = seqGetData(gds, "annotation/format/DS")$data
  return(genotype)
}
