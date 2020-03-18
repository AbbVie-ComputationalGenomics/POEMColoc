#' Perform some POEM specific checks on datasets
#'
#' @param dataset_list
#' list of datasets formatted according to the requirements of coloc.abf with the following differences
#' pos, and chr are additionally required.
#' snp is ignored as it is determined by pos.
#' For top SNP datasets, R2 and MAF and gds_file and subset are optional, if supplied as arguments to the function.
#' For full summary statistic datasets, MAF is required.
#' @return
#' NULL
#' @examples
#' dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z")
#' dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
#' dataset_top_SNP_no_chr <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9)
#' dataset_top_SNP_different_chr <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Y")
#' dataset_full_zero_p <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = rep(0,5), chr= "Z")
#' # Should pass
#' check_datasets(list(dataset_full, dataset_top_SNP))
#' # Should fail
#' try(check_datasets(list(dataset_full, dataset_top_SNP, dataset_top_SNP_no_chr)))
#' try(check_datasets(list(dataset_full, dataset_top_SNP, dataset_top_SNP_different_chr)))
#' try(check_datasets(list(dataset_full, dataset_top_SNP, dataset_full_zero_p)))

check_datasets <- function(dataset_list) {
  # Other sanity checks on input specific to POEM (pos, chr)
  sapply(dataset_list, single_dataset_check)
  if (!length(unique(sapply(dataset_list, function(x) x$chr)))==1) {
    stop("Expect all datasets to be from the same chromosome")
  }
  # zero p-values create problems always, especially for POEM
  sapply(dataset_list, check_zero_pvalues)
  # formatted input
}
