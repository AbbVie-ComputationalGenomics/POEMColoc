#' Add some additional fields to datasets
#'
#' @param dataset
#' list of datasets formatted according to the requirements of coloc.abf with the following differences
#' pos, and chr are additionally required.
#' snp is ignored as it is determined by pos.
#' For top SNP datasets, R2 and MAF and gds_file and subset are optional, if supplied as arguments to the function.
#' For full summary statistic datasets, MAF is required.
#' @return
#' list of datasets with the following additions
#' snp is created from pos
#' top_pos is created as the position of the lowest p-value
#' imputation_class either all or top indicating if imputation is needed.
#' @examples
#' # p-value input
#' dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z")
#' dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
#' add_metadata(dataset_full)
#' add_metadata(dataset_top_SNP)
#' # beta-varbeta input
#' dataset_full_beta_varbeta = dataset_full
#' dataset_full_beta_varbeta$pvalues = NULL
#' dataset_full_beta_varbeta$beta = c(0.5, 0.5, 0.4, 0.5, 0.26)
#' dataset_full_beta_varbeta$varbeta = c(0.009, 0.008, 0.012, 0.01, 0.015)
#' dataset_full_beta_varbeta$sdY = 1
#' dataset_top_SNP_beta <- list(pos = 2, N= 10000, s =0.5, type="cc", beta = 2, varbeta = 0.05, chr = "Z")
#' add_metadata(dataset_full_beta_varbeta)
#' add_metadata(dataset = dataset_top_SNP_beta)
add_metadata <- function(dataset) {
  # add imputation class
  dataset$imputation_class <- detect_imputation_class(dataset)
  dataset$snp <- as.character(dataset$pos)
  # For those of imputation class all, create SNP from pos and add top SNP
  if (dataset$imputation_class == "all") {
    if (!is.null(dataset$pvalues)) {
      dataset$top_pos <- dataset$pos[which.min(dataset$pvalues)]
    } else {
      Zs <- dataset$beta / sqrt(dataset$varbeta)
      dataset$top_pos <- dataset$pos[which.max(abs(Zs))]
    }
  } else {
    dataset$top_pos <- dataset$pos
  }
  return(dataset)
}
