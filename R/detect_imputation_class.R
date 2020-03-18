#' Determine whether dataset has full or top SNP statistics
#'
#' @param dataset
#' POEMColoc formated dataset
#' @return
#' "top" or "all'
#' @examples
#' # p-value input
#' dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z")
#' dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
#' detect_imputation_class(dataset_full)
#' detect_imputation_class(dataset_top_SNP)
#' # beta-varbeta input
#' dataset_full_beta_varbeta = dataset_full
#' dataset_full_beta_varbeta$pvalues = NULL
#' dataset_full_beta_varbeta$beta = c(0.5, 0.5, 0.4, 0.5, 0.26)
#' dataset_full_beta_varbeta$varbeta = c(0.008, 0.008, 0.012, 0.01, 0.015)
#' dataset_full_beta_varbeta$sdY = 1
#' dataset_top_SNP_beta <- list(pos = 2, N= 10000, s =0.5, type="cc", beta = 2, varbeta = 0.05, chr = "Z")
#' detect_imputation_class(dataset_full_beta_varbeta)
#' detect_imputation_class(dataset = dataset_top_SNP_beta)

detect_imputation_class <- function(dataset) {
  nd <- names(dataset)
  if ("beta" %in% nd & "varbeta" %in% nd) {
    L <- length(dataset$beta)
    if (!length(dataset$varbeta)==L |!length(dataset$pos)==L) {
      stop("expect varbeta and pos to match length of beta")
    }
  } else {
    L <- length(dataset$pvalues)
    if (!length(dataset$pos)==L) {
      stop("expect pos to match length of pvalues")
    }
  }
  if (L==0) {
    stop("Need to supply nonzero length beta and varbeta or pvalues")
  }
  if (L==1) {
    impute_class <- "top"
  } else {
    impute_class <- "all"
  }
  return(impute_class)
}
