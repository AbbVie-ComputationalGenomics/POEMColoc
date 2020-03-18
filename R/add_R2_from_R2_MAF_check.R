#' Harmonize R2 and MAF inputs
#'
#' @param dataset_list
#' list of datasets formatted according to the requirements of coloc.abf with the following differences
#' pos, and chr are additionally required.
#' snp is ignored as it is determined by pos.
#' For top SNP datasets, R2 and MAF and gds_file and subset are optional, if supplied as arguments to the function.
#' For full summary statistic datasets, MAF is required.
#' imputation class either all or top
#' top_pos giving the position of the top p-value
#' @param gds_file
#' path to a gds file for the reference panel.  Can be NULL if specifying R2 and MAF, or if supplying within dataset.
#' @param R2
#' A matrix of R2 LD between SNP positions.  This should contain LD information for top SNPs in rows and LD for positions to be imputed in columns.  Can be NULL if specifying gds_file or if included in individual datasets.
#' Column and row names should be position.
#' @param MAF
#' optional input of minor allele frequencies.
#' only used for top SNP datasets
#' names should be position
#' @param subset
#' optional parameter that is the path to a file giving the names of samples to include as formatted in the gds file.  each line should be the name of a sample.
#' @param window_size
#' if NULL, window size is determined to cover all positions in dataset_list
#' if specified, covers all positions in dataset list plus window_size around top SNP dataset
#' @param get_dosage_fn
#' a function taking gds file connection and returning
#' a matrix with the genotype or imputed dosage at each position.  The function
#' should not perform any filtering or reordering of variants.  rows must correspond
#' to individuals and columns to positions.
#' @param dataset_R2_condition
#' vector of condition as output by R2_MAF_check_dataset.
#' If this output is inconsistent with the other things you supplied you
#' will probably get errors that will either break it or cause incorrect output
#' so don't use this independantly without a good reason.

add_R2_from_R2_MAF_check <- function(dataset_list, gds_file, R2, MAF, subset, window_size, get_dosage_fn, dataset_R2_condition) {
  # check how many distinct ones there are for top SNP datasets
  nfunction <- sum(grepl("function|none", dataset_R2_condition) & !grepl("MAF", dataset_R2_condition))
  if (nfunction > 0 & nfunction < sum(!grepl("MAF", dataset_R2_condition))) {
    warning("You have supplied R2 or gds_file at the dataset level for some datasets, but not all.  Function inputs will be used for those datasets for which it is not supplied.")
  }
  # Add function-level input to dataset
  for (i in seq_along(dataset_list)) {
    if (dataset_R2_condition[i] =="R2 from function") {
      dataset_list[[i]]$R2 <- R2
      dataset_list[[i]]$MAF <- MAF
    }
    if (dataset_R2_condition[i] == "gds from function") {
      dataset_list[[i]]$gds_file <- gds_file
      dataset_list[[i]]$subset <- subset
    }
  }
  # Format matrix R2 as vector
  for (i in which(grepl("R2", dataset_R2_condition))) {
    if (is.matrix(dataset_list[[i]]$R2)) {
      if (is.null(rownames(dataset_list[[i]]$R2))) {
        stop("Expect R2 matrix to have row names indicating position to impute from")
      }
      if (is.null(rownames(dataset_list[[i]]$R2))) {
        stop("Expect R2 matrix to have column names indicating position to impute to")
      }
      top_snp_pos <- as.character(dataset_list[[i]]$top_pos)
      if (top_snp_pos %in% rownames(dataset_list[[i]]$R2)) {
        dataset_list[[i]]$R2 <- dataset_list[[i]]$R2[top_snp_pos,]
      } else {
        dataset_list[[i]]$R2 <- NA
        warning("Top SNP not found in dataset R2, set to NA")
      }
    }
  }
  # Get the R2 from gds where needed
  dataset_list <- add_R2_from_gds(dataset_list, dataset_R2_condition, window_size, get_dosage_fn = get_dosage_fn)
  # Check and make sure same positions are in R2 and MAF
  return(dataset_list)
}
