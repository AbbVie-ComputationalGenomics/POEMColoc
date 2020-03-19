#' Imputed coloc input
#'
#' @param dataset
#' POEM_input format dataset
#' Specifically, this is a coloc format dataset with the following additions
#' pos length L (where L is length of pvalues or beta)
#' chr length 1
#' snp length L equal to as.character(pos) is mandatory
#' imputation class either all or top
#' top_pos giving the position of the top p-value
#' if imputation class is top
#' -R2 which is either a named vector R2 with names being snp of length M or NA
#' -MAF which is either a named vector with names being snp of length P or NA
#' if imputation class is 'all'
#' -vector MAF of length L
#' @return
#' POEM_output format dataset
#' Specifically this is a coloc format dataset with the following changes
#' top_pos giving the position of the top p-value
#' pos vector giving positions
#' chr giving chromosome
#' snp is mandatory
#' MAF is mandatory
#' imputation class that is either all or top
#' Unlike coloc, length of pvalues or beta L is permitted to be zero
#' positions included will be those in R2 and MAF that are not NA
#' @examples
#' R2 = c(0.9, 1, 0.5, 0, 0.2)
#' names(R2) <- 1:5
#' MAF = c(0.1, 0.15, 0.3, 0.2, 0.4)
#' names(MAF) <- 1:5
#' # Top snp example
#' dataset_top_SNP <- list(pos = 2, MAF = MAF, R2=R2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=2, snp="2")
#' imputed_coloc_input(dataset = dataset_top_SNP)
#' # Full example
#' dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z", imputation_class = "all", top_pos = 1, snp=as.character(1:5))
#' imputed_coloc_input(dataset = dataset_full)
#' # Top snp example with beta
#' dataset_top_SNP_beta <- list(pos = 2, snp = "2", MAF = MAF, R2=R2, N= 10000, s =0.5, type="cc", beta = 2, varbeta = 0.05, chr = "Z", imputation_class = "top", top_pos=2)
#' imputed_coloc_input(dataset = dataset_top_SNP_beta)
#' # Full dataset example with beta
#' dataset_full_beta_varbeta = dataset_full
#' dataset_full_beta_varbeta$pvalues = NULL
#' dataset_full_beta_varbeta$beta = c(0.5, 0.5, 0.4, 0.5, 0.26)
#' dataset_full_beta_varbeta$varbeta = c(0.008, 0.008, 0.012, 0.01, 0.015)
#' dataset_full_beta_varbeta$sdY = 1
#' imputed_coloc_input(dataset = dataset_full_beta_varbeta)
#' # Top SNP dataset with NA R2
#' dataset_top_SNP_NA <- list(pos = 2, MAF = MAF, R2=NA, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=2)
#' imputed_coloc_input(dataset = dataset_top_SNP_NA)
#' # Top SNP dataset with missing MAF positions
#' dataset_top_SNP_beta2 <- list(pos = 2, snp = "2", MAF = MAF[-1], R2=R2, N= 10000, s =0.5, type="cc", beta = 2, varbeta = 0.05, chr = "Z", imputation_class = "top", top_pos=2)
#' imputed_coloc_input(dataset = dataset_top_SNP_beta2)
imputed_coloc_input <- function(dataset) {
  if (dataset$imputation_class == "all") {
    return(dataset)
  } else {
    if (all(is.na(dataset$R2))) {
      imputed_summary_stats_pval <- list(type=dataset$type, sdY=dataset$sdY, MAF=numeric(0),
                                         pos=integer(0), pvalues=numeric(0),
                                         N=dataset$N, s=dataset$s, snp=character(0), chr=dataset$chr,
                                         imputation_class = dataset$imputation_class, top_pos = dataset$top_pos)
    } else {
      if ("beta" %in% names(dataset) & "varbeta" %in% names(dataset)) {
        Zest <- dataset$beta / sqrt(dataset$varbeta)
      } else {
        Zest <- qt(dataset$pvalues / 2, dataset$N)
      }
      # Imputed summary stats
      point_est <- sqrt(dataset$R2) * Zest
      names(point_est) <- names(dataset$R2)
      # missing positions from SNPs not in reference panel
      keep <- !is.na(point_est)
      pos_keep <- (names(dataset$R2))[keep]
      # Ensure that MAF is also present for retained positions.
      pos_keep <- intersect(pos_keep, names(dataset$MAF[!is.na(dataset$MAF)]))
      imputed_summary_stats_pval <- list(type=dataset$type, sdY=dataset$sdY, MAF=dataset$MAF[pos_keep],
                                         pos=as.integer(pos_keep), pvalues=pt(-abs(point_est[pos_keep]), dataset$N) * 2,
                                         N=dataset$N, s=dataset$s, snp=pos_keep, chr=dataset$chr,
                                         imputation_class = dataset$imputation_class, top_pos = dataset$top_pos)
    }
    return(imputed_summary_stats_pval)
  }
}
