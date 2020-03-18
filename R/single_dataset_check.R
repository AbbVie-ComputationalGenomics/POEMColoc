#' Check a single dataset for POEM-specific requirements and fail if not met
#'
#' @param dataset
#' single POEMColoc dataset
#' @return
#' NULL
#' @examples
#' dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z")
#' dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
#' dataset_top_SNP_no_chr <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9)
#' dataset_top_SNP_no_pos <- list(N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
#' dataset_top_SNP_have_snp <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", snp = "abc")
#' dataset_top_SNP_too_many_chr <- list(pos = c(1,2), N= 10000, s =0.5, type="cc", pvalues = c(10^-8,10^-9), chr = c(1,2))
#' dataset_top_SNP_too_few_chr <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = c())
#' dataset_top_SNP_no_pvalues <- list(pos = 2, N= 10000, s =0.5, type="cc", chr = "Z")
#' # Should pass because correctly formatted
#' single_dataset_check(dataset_full)
#' single_dataset_check(dataset_top_SNP)
#' # Should pass because misformatting will be detected later
#' single_dataset_check(dataset_top_SNP_no_pvalues)
#' # Should warn
#' single_dataset_check(dataset_top_SNP_have_snp)
#' # Should fail
#' try(single_dataset_check(dataset_top_SNP_no_chr))
#' try(single_dataset_check(dataset_top_SNP_no_pos))
#' try(single_dataset_check(dataset_top_SNP_too_many_chr))
#' try(single_dataset_check(dataset_top_SNP_too_few_chr))
#'
single_dataset_check <- function(dataset) {
  if (!"pos" %in% names(dataset)) {
    stop("Expect snp position pos in dataset")
  }
  if (!"chr" %in% names(dataset)) {
    stop("Expect chromosome chr in dataset")
  }
  if ("snp" %in% names(dataset)) {
    warning("snp will be overriden by pos.  use pos to indicate same position in different datasets.")
  }
  if (!length(dataset$chr)==1) {
    stop("Expect chr to be length 1")
  }
}
