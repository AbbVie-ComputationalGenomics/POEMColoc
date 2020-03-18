#' Check zero p-values
#'
#' @param dataset
#' single POEMColoc dataset
#' @return
#' NULL
#' @examples
#' dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z")
#' dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
#' dataset_empty <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = c(), chr = "Z")
#' # Should pass
#' check_zero_pvalues(dataset_full)
#' check_zero_pvalues(dataset_top_SNP)
#' # Should fail
#' dataset_top_SNP_zero_p <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 0, chr = "Z")
#' dataset_full_zero_p <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = rep(0,5), chr= "Z")
#' try(check_zero_pvalues(dataset_top_SNP_zero_p))
#' try(check_zero_pvalues(dataset_full_zero_p))
#' try(check_zero_pvalues(dataset_empty))
#' # Should warn
#' dataset_full_zero_p2 <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(0, 0, 10^-5, 0, 0), chr= "Z")
#' check_zero_pvalues(dataset_full_zero_p2)

check_zero_pvalues <- function(dataset) {
  if ("pvalues" %in% names(dataset)) {
    if (all(dataset$pvalues==0)) {
      stop("No non-zero p-values.  Cannot run coloc.")
    }
    if (any(dataset$pvalues==0)) {
      warning("Zero p-values detected.  These will be removed by coloc")
    }
  }
}
