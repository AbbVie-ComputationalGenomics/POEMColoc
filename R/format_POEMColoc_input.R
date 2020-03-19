#' Given arguments to POEMColoc, convert input to standardized format for imputation step.
#'
#' @param dataset_list
#' list of datasets formatted according to the requirements of coloc.abf with the following differences
#' pos, and chr are additionally required.
#' snp is ignored as it is determined by pos.
#' For top SNP datasets, R2 and MAF and gds_file and subset are optional, if supplied as arguments to the function.
#' For full summary statistic datasets, MAF is required.
#' @param gds_file
#' path to a gds file for the reference panel.  Can be NULL if specifying R2 and MAF, or if supplying within dataset.
#' @param subset
#' optional parameter that is the path to a file giving the names of samples to include as formatted in the gds file.  each line should be the name of a sample.
#' @param R2
#' A matrix of R2 LD between SNP positions.  This should contain LD information for top SNPs in rows and LD for positions to be imputed in columns.  Can be NULL if specifying gds_file or if included in individual datasets.
#' Column and row names should be position.
#' @param MAF
#' optional input of minor allele frequencies.
#' only used for top SNP datasets
#' names should be position
#' @param window_size
#' if NULL, window size is determined to cover all positions in dataset_list
#' if specified, covers all positions in dataset list plus window_size around top SNP dataset
#' @param get_dosage_fn
#' a function taking gds file connection and returning
#' a matrix with the genotype or imputed dosage at each position.  The function
#' should not perform any filtering or reordering of variants.  rows must correspond
#' to individuals and columns to positions.
#' @param flexible_input
#' if TRUE allows specifying gds_file and R2 by dataset.  if FALSE will only use function inputs.
#' the latter way to input has been better tested.  Further, it is not so clear if one ought to
#' attempt to colocalize things with divergent LD backgrounds.
#' @return
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
#' @examples

format_POEMColoc_input <- function(dataset_list, gds_file, R2, MAF, subset, window_size, get_dosage_fn = get_dosage_alt, flexible_input = TRUE) {
  check_datasets(dataset_list)
  dataset_list <- lapply(dataset_list, add_metadata)
  if (flexible_input) {
    dataset_list <- add_R2_to_dataset_list(dataset_list, gds_file = gds_file, R2 = R2, MAF = MAF,
                                           subset = subset, window_size = window_size, get_dosage_fn = get_dosage_fn)
  } else {
    dataset_list <- add_R2_to_dataset_list_function_input_only(dataset_list, gds_file = gds_file, R2 = R2, MAF = MAF,
                                                               subset = subset, window_size = window_size, get_dosage_fn = get_dosage_fn)
  }
  return(dataset_list)
}
