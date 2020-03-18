#' Remove positions from dataset
#'
#' @param dataset
#' POEM_output format dataset
#' Specifically this is a coloc format dataset with the following changes
#' top_pos giving the position of the top p-value
#' pos vector giving positions
#' chr giving chromosome
#' snp is mandatory
#' MAF is mandatory
#' imputation class that is either all or top
#' Unlike coloc, length of pvalues or beta L is permitted to be zero
#' @param keep
#' vector of length L with positions to keep
#' @return
#' POEM_output format dataset
#' @examples
#' # example with p-values
#' dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z", imputation_class = "all", top_pos = 1, snp = as.character(1:5))
#' dataset_filter(dataset = dataset_full, c(TRUE, FALSE, FALSE, TRUE, FALSE))
#' dataset_filter(dataset = dataset_full, rep(TRUE,5))
#' dataset_filter(dataset = dataset_full, rep(FALSE,5))
#' # example with beta
#' dataset_full_beta_varbeta = dataset_full
#' dataset_full_beta_varbeta$pvalues = NULL
#' dataset_full_beta_varbeta$beta = c(0.5, 0.5, 0.4, 0.5, 0.26)
#' dataset_full_beta_varbeta$varbeta = c(0.008, 0.008, 0.012, 0.01, 0.015)
#' dataset_full_beta_varbeta$sdY = 1
#' dataset_filter(dataset = dataset_full_beta_varbeta, c(TRUE, FALSE, FALSE, TRUE, FALSE))
#' dataset_filter(dataset = dataset_full_beta_varbeta, rep(TRUE,5))
#' dataset_filter(dataset = dataset_full_beta_varbeta, rep(FALSE,5))
#' # example with empty that should be TRUE
#' dataset_empty <- list(pos = integer(0), MAF = numeric(0), N=1000, type ="quant", pvalues = numeric(0), chr= "Z", imputation_class = "all", top_pos = 1, snp = character(0))
#' dataset_filter(dataset = dataset_empty, logical())

dataset_filter <- function(dataset, keep) {
  if (!empty_dataset_check(dataset)) {
    dataset_filter <- dataset
    if (!is.null(dataset_filter$pvalues)) {
      dataset_filter$pvalues <- dataset$pvalues[keep]
    }
    if (!is.null(dataset_filter$beta)) {
      dataset_filter$beta <- dataset$beta[keep]
    }
    if (!is.null(dataset_filter$varbeta)) {
      dataset_filter$varbeta <- dataset$varbeta[keep]
    }
    dataset_filter$snp <- dataset$snp[keep]
    dataset_filter$pos <- dataset$pos[keep]
    dataset_filter$MAF <- dataset$MAF[keep]
  } else {
    dataset_filter <- dataset
  }
  return(dataset_filter)
}
