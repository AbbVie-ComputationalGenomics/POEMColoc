#' Perform coloc analysis using imputed p-values
#'
#' @param dataset_list
#' list of datasets formatted according to the requirements of coloc.abf with the following differences
#' pos, and chr are additionally required.
#' snp is ignored as it is determined by pos.
#' For top SNP datasets, R2 and MAF and gds_file and subset are optional, if supplied as arguments to the function.
#' For full summary statistic datasets, MAF is required.
#' Need to set flexible_input=TRUE if using this
#' @param config
#' a list of vectors of indices indicating which elements of dataset_list are to be run in a single colocalization analysis.
#' currently, they must be length 2 as we do not support multitrait colocalization
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
#' It is required if there are two top SNP datasets.
#' Default NULL will impute all positions present in the full datasets.
#' It should not be NULL if all datasets are top SNP datasets as this will impute everything between the two top SNPs only.
#' Otherwise, for each pairwise colocalization, only SNPs within window size of one of their top SNPs will be included, including all SNPs between them if disjoint.
#' @param min_MAF
#' Remove any site with minor allele frequency less than min_MAF.  Setting to zero means no filtering.
#' @param p1
#' prior probability a SNP is associated with trait t, where t is the index in dataset list, default 1e-4.  Either a single number or a vector of length same as dataset_list
#' @param p12
#' prior probability a SNP is associated with both traits, default 1e-5.  Either a single value or of the same length as config, with element i corresponding to the two sets in config[[i]]
#' @param get_dosage_fn
#' a function taking gds file connection and returning
#' a matrix with the genotype or imputed dosage at each position.  The function
#' should not perform any filtering or reordering of variants.  rows must correspond
#' to individuals and columns to positions.
#' @param flexible_input
#' if TRUE allows specifying gds_file and R2 by dataset.  if FALSE will only use function inputs.
#' the latter way to input has been better tested.  Further, it is not so clear if one ought to
#' attempt to colocalize things with divergent LD backgrounds.
#' @details
#' dataset and (elements of) dataset_list follow the format for coloc.abf in the coloc R package, except snp is not optional and additional elements pos and chr giving genome position and chromosome id are required.
#' Additionally, MAF need not be supplied for datasets to be imputed because this will be obtained elsewhere.
#' One or both of the two datasets can be in top SNP format.  snp will be overridden with pos so it is expected that all sources will use the same coordinate system.
#' Multiple alternative allele per site are collapsed.
#' If you have a vcf file first convert to gds using SeqArray:::seqVCF2GDS(your_vcf, your_gds)
#' prior parameters p1, p2, p12 are supplied to the coloc.abf function
#' POEMColoc creates imputed datasets for input to the coloc.abf function, refer to the coloc.abf documentation from the coloc R package for more information on this method and its inputs
#' pos and chr are additional requirements of POEMColoc in order to recognize each position in R2 and MAF, or in the gds file in the panel.
#' There are different ways to run POEMColoc.
#' POEMColoc can be run with R2 and MAF supplied as function arguments.  Note in this case R2 must be a matrix with both rownames including any top SNP positions to be imputed from and column names being positions of any SNPs to impute, including the top SNP itself.
#' R2 and MAF can also be supplied individually by dataset, in which case both can be named vectors with names being the positions to impute.  It is assume that you will impute from the top SNP position in the dataset.
#' Datasets consisting of more than one SNP will not be imputed.
#' Alternatively, R2 and MAF can be computed from a gds file, which can be created from a vcf file (see above).
#' The gds file can either be supplied as a function argument or separately by dataset.
#' Note that user-supplied R2 and MAF will take priority over gds file, and items specified by dataset will override options specified as function arguments.
#' @return
#' list of coloc.abf results.

POEMColoc_from_config <- function(dataset_list, config, gds_file = NULL, subset = NULL, R2 = NULL, MAF=NULL, window_size=10^6, min_MAF = 0, p1 = 10^-4, p12 = 10^-5, get_dosage_fn = get_dosage_alt, flexible_input = TRUE) {
  dataset_list <- format_POEMColoc_input(dataset_list = dataset_list, gds_file = gds_file, R2 = R2, MAF = MAF, subset = subset, window_size = window_size, get_dosage_fn = get_dosage_fn, flexible_input = flexible_input)
  imputed_dataset <- lapply(dataset_list, imputed_coloc_input)
  imputed_dataset <- lapply(imputed_dataset, MAF_filter, min_MAF = min_MAF)
  priors <- format_priors(dataset_list, config, p1, p12)
  coloc_res <- vector("list", length(config))
  for (i in seq_along(config)) {
    dataset_pair <- position_filter(imputed_dataset[config[[i]]], window_size)
    is_empty <- sapply(dataset_pair, empty_dataset_check)
    # additional check for window size being NULL with two top SNP
    if (!any(is_empty)) {
      coloc_string = capture.output(res <- coloc.abf(dataset1 = dataset_pair[[1]], dataset2 = dataset_pair[[2]],
                                                     p1 = priors$p1[config[[i]][1]], p2 = priors$p1[config[[i]][2]],
                                                     p12 = priors$p12[i]))
      coloc_res[[i]] <- res
    } else {
      coloc_res[[i]] <- NA
    }
  }
  return(coloc_res)
}

