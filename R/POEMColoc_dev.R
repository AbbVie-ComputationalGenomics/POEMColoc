#' Perform coloc analysis using imputed p-values - development version supporting separate imputation panel for each dataset

#' @param dataset
#' single dataset formatted according to the requirements of coloc.abf with the following differences
#' pos, and chr are additionally required.
#' snp is ignored as it is determined by pos.
#' For top SNP datasets, R2 and MAF and gds_file and subset are optional, if supplied as arguments to the function.
#' For full summary statistic datasets, MAF is required.
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
#' It is required if there are two top SNP datasets.
#' Default NULL will impute all positions present in the full datasets.
#' It should not be NULL if all datasets are top SNP datasets as this will impute everything between the two top SNPs only.
#' Otherwise, for each pairwise colocalization, only SNPs within window size of one of their top SNPs will be included, including all SNPs between them if disjoint.
#' @param min_MAF
#' Remove any site with minor allele frequency less than min_MAF.  Setting to zero means no filtering.
#' @param p1
#' prior probability a SNP is associated with trait in dataset, default 1e-4.
#' @param p2
#' prior probability a SNP is associated with trait t, where t is the index in dataset list, default 1e-4.  Either a single number or a vector of length same as dataset_list.
#' @param p12
#' prior probability a SNP is associated with both traits, default 1e-5.  Either a single value or of the same length as dataset_list.
#'@param dosage_method
#' Method to compute dosage if gds_file input mode is used.  one of 'dosage', 'dosage_alt' or 'dosage_imputed', see details.
#' @details
#' dataset and (elements of) dataset_list follow the format for coloc.abf in the coloc R package, except snp is not optional and additional elements pos and chr giving genome position and chromosome id are required.
#' Additionally, MAF need not be supplied for datasets to be imputed because this will be obtained elsewhere.
#' One or both of the two datasets can be in top SNP format.  snp will be overridden with pos so it is expected that all sources will use the same coordinate system.
#' Multiple alternative allele per site are collapsed.
#' If you have a vcf file first convert to gds using SeqArray:::seqVCF2GDS(your_vcf, your_gds)
#' prior parameters p1, p2, p12 are supplied to the coloc.abf function
#' POEMColoc_dev creates imputed datasets for input to the coloc.abf function, refer to the coloc.abf documentation from the coloc R package for more information on this method and its inputs
#' pos and chr are additional requirements of POEMColoc_dev in order to recognize each position in R2 and MAF, or in the gds file in the panel.
#' There are different ways to run POEMColoc_dev.
#' POEMColoc_dev can be run with R2 and MAF supplied as function arguments.  Note in this case R2 must be a matrix with both rownames including any top SNP positions to be imputed from and column names being positions of any SNPs to impute, including the top SNP itself.
#' R2 and MAF can also be supplied individually by dataset, in which case both can be named vectors with names being the positions to impute.  It is assume that you will impute from the top SNP position in the dataset.
#' Datasets consisting of more than one SNP will not be imputed.
#' Alternatively, R2 and MAF can be computed from a gds file, which can be created from a vcf file (see above).
#' The gds file can either be supplied as a function argument or separately by dataset.
#' Note that user-supplied R2 and MAF will take priority over gds file, and items specified by dataset will override options specified as function arguments.
#' subset will follow the same convention as gds_file such that if the gds_file is specified at the function level, subset will also be taken from the function
#' dosage_method is designed to support different gds input possibilities.
#' 'dosage' will compute reference allele dosage, 'dosage_alt' will compute dosages of alternative allele, imputed_dosage will use 'annotation/format/DS' to get the dosage.
#' note that not all gds files will have all nodes available, so choose an option that is compatible with the gds_file supplied.
#' @return
#' list of coloc.abf results.
#' @export
#'
#' @examples
#' # Example using R2 input mode
#' dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z")
#' dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
#' R2 <- matrix(data = c(0.94, 1, 0.65, 0, 0.51, 0.32, 0.64, 0.65, 1, 0, 0.21, 0.12), nrow=2, byrow=TRUE)
#' colnames(R2) <- 1:6
#' rownames(R2) <- 2:3
#' MAF = c(0.14, 0.15, 0.25, 0.2, 0.4, 0.25)
#' names(MAF) <- 1:6
#' POEMColoc_dev(dataset_top_SNP, dataset_full, gds_file = NULL, R2 = R2, MAF = MAF, subset=NULL)
#' # Example using R2 input mode and dataset_full is a list
#' dataset_full2 <- list(pos = c(1, 3, 4, 6), MAF = c(0.14, 0.25, 0.2, 0.25), type = "quant", pvalues = c(3 * 10^-4, 2 * 10^-8, 0.8, 0.03), N =350, chr= "Z")
#' dataset_full3 <- list(pos = c(1,4,5), MAF = c(0.14, 0.2, 0.4), type = "quant", pvalues = c(0.81, 10^-6, 0.97), N =350, chr= "Z")
#' POEMColoc_dev(dataset_top_SNP, list(dataset_full, dataset_full2, dataset_full3), gds_file = NULL, R2 = R2, subset=NULL, MAF = MAF)
#' # Example with two top SNP
#' dataset_top_SNP2 <- list(pos = 3, type= "quant", pvalues =2*10^-8, N = 350, chr = "Z")
#' POEMColoc_dev(dataset_top_SNP, dataset_top_SNP2, gds_file = NULL, R2 = R2, subset=NULL, MAF = MAF)
#' # Example with a mix of top SNP and SNP
#' POEMColoc_dev(dataset_top_SNP, list(dataset_full, dataset_top_SNP2, dataset_full3), gds_file = NULL, R2 = R2, subset=NULL, MAF = MAF)
#' # Example using min_MAF: should use 3 SNPs
#' POEMColoc_dev(dataset_top_SNP, dataset_full, gds_file = NULL, R2 = R2, MAF = MAF, subset=NULL, min_MAF=0.2)
#' # Example using min_MAF: should return NA because no SNPs meet filter
#' POEMColoc_dev(dataset_top_SNP, dataset_full, gds_file = NULL, R2 = R2, MAF = MAF, subset=NULL, min_MAF=0.5)
#' # Example using beta and varbeta
#' dataset_full_beta_varbeta = dataset_full
#' dataset_full_beta_varbeta$pvalues = NULL
#' dataset_full_beta_varbeta$beta = c(0.5, 0.5, 0.4, 0.5, 0.26)
#' dataset_full_beta_varbeta$varbeta = c(0.008, 0.008, 0.012, 0.01, 0.015)
#' dataset_full_beta_varbeta$sdY = 1
#' POEMColoc_dev(dataset_top_SNP, dataset_full_beta_varbeta, gds_file = NULL, R2 = R2, MAF = MAF, subset=NULL)
#' # Example that should reduce to full coloc assuming pos is snp.
#' POEMColoc_dev(dataset_full, dataset_full)
#' # Window size should also apply to datasets without top SNP
#' POEMColoc_dev(dataset_full, dataset_full, window_size = 1)
#' # Example of changing prior from coloc default
#' POEMColoc_dev(dataset_top_SNP, dataset_full, gds_file = NULL, R2 = R2, MAF = MAF, subset=NULL, p1 = 10^-6, p2 =10^-6, p12 = 10^-7)
#' # Example with R2 specified at the dataset level
#' POEMColoc_dev(dataset =c(dataset_top_SNP, list(R2=R2, MAF = MAF)), dataset_list = list(c(dataset_top_SNP2, list(R2 = R2["3",], MAF = MAF)), dataset_full, dataset_full2), window_size =1)
#' # Same thing with R2 specified at function level that should be same (in this case we used the same R2 for both sets, but we could have set them different)
#' POEMColoc_dev(dataset = dataset_top_SNP, dataset_list = list(dataset_top_SNP2, dataset_full, dataset_full2), window_size =1, R2=R2, MAF=MAF)
#' # Example with full dataset colocalized against top SNP
#' POEMColoc_dev(dataset = dataset_full_beta_varbeta, list(dataset_top_SNP, dataset_top_SNP2), R2=R2, MAF=MAF)
#' # Example using gds_file input mode
#' gds_file <- system.file("extdata", "example.gds", package = "POEMColoc")
#' gds_file2 <- system.file("extdata", "example2.gds", package = "POEMColoc")
#' subset <- system.file("extdata", "subset.ped", package = "POEMColoc")
#' subset2 <- system.file("extdata", "subset2.ped", package = "POEMColoc")
#' POEMColoc_dev(dataset_top_SNP, dataset_full, gds_file = gds_file)
#' # Example using gds_file input mode with subset
#' POEMColoc_dev(dataset_top_SNP, dataset_full, gds_file = gds_file, subset = subset)
#' # Example using gds_file input mode and dataset_full is a list
#' POEMColoc_dev(dataset_top_SNP, list(dataset_full, dataset_full2, dataset_full3), gds_file = gds_file)
#' # Example that should return NA because the top SNP is not in the file
#' dataset_top_SNP_missing <- list(pos = 8, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
#' POEMColoc_dev(dataset_top_SNP_missing, dataset_full, gds_file = gds_file)
#' # Example that should return NA because none of the SNPs are in the file
#' dataset_full_missing <- list(pos = c(7,8,9), MAF = c(0.14, 0.2, 0.4), type = "quant", pvalues = c(0.81, 10^-6, 0.97), N =350, chr= "Z")
#' POEMColoc_dev(dataset_top_SNP_missing, dataset_full_missing, gds_file = gds_file)
#' # Example that should still work because some eQTL SNPs are in the file and the top SNP is in the file
#' dataset_partial_missing <- list(pos = c(1,2,9), MAF = c(0.14, 0.2, 0.4), type = "quant", pvalues = c(0.81, 10^-6, 0.97), N =350, chr= "Z")
#' POEMColoc_dev(dataset_top_SNP, dataset_partial_missing, gds_file = gds_file)
#' # Example with two top SNP
#' POEMColoc_dev(dataset_top_SNP, dataset_top_SNP2, gds_file = gds_file, window_size = 2)
#' POEMColoc_dev(dataset_top_SNP, dataset_top_SNP2, gds_file = gds_file, window_size = NULL)
#' # Example with a mix of top SNP and SNP
#' POEMColoc_dev(dataset_top_SNP, list(dataset_full, dataset_top_SNP2, dataset_full3), gds_file = gds_file, window_size = 3)
#' # Same example reducing window size
#' POEMColoc_dev(dataset_top_SNP, list(dataset_full, dataset_top_SNP2, dataset_full3), gds_file = gds_file, window_size = 1)
#' # Example using min_MAF
#' POEMColoc_dev(dataset_top_SNP, dataset_full, gds_file = gds_file, min_MAF = 0.3)
#' # Example using min_MAF that should return NA
#' POEMColoc_dev(dataset_top_SNP, dataset_full, gds_file = gds_file, min_MAF = 0.5)
#' # Example with dataset-specific gds file input
#' POEMColoc_dev(dataset_full, list(c(dataset_top_SNP, list(gds_file = gds_file, subset = subset)),
#' c(dataset_top_SNP2, list(gds_file = gds_file, subset = subset2)),
#' c(dataset_top_SNP, list(gds_file = gds_file2, subset = subset))), window_size =1)


POEMColoc_dev <- function(dataset, dataset_list, gds_file = NULL, R2 = NULL, MAF=NULL, subset=NULL, window_size=NULL, min_MAF = 0, p1 = 10^-4, p2 =10^-4, p12 = 10^-5, dosage_method = 'dosage_alt') {
  if (!detect_list_condition(dataset_list)) {
    dataset_list <- list(dataset_list)
  }
  if (detect_list_condition(dataset)) {
    stop("Only argument 2 dataset_list may be a list of datasets")
  }
  if (length(p2) == 1) {
    p2 <- rep(p2, length(dataset_list))
  } else {
    if (!length(p2) == length(dataset_list)) {
      stop("p12 should either be length 1 or length equal to the number of elements of dataset list")
    }
  }
  if (length(p12) == 1) {
    p12 <- rep(p12, length(dataset_list))
  } else {
    if (!length(p12) == length(dataset_list)) {
      stop("p12 should either be length 1 or length equal to the number of elements of dataset list")
    }
  }
  dosage_fn <- switch(dosage_method, "dosage_alt"= get_dosage_alt, "dosage"= get_dosage,
                      "dosage_imputed" = get_imputed_dosage)
  if (is.null(dosage_fn)) {
    stop("Unrecognized dosage method")
  }
  dataset_list <- c(list(dataset), dataset_list)
  config <- lapply(as.list(2:length(dataset_list)), function(x) c(1,x))
  POEMColoc_from_config(dataset_list, config, gds_file = gds_file, subset = subset,
                        R2 = R2, MAF = MAF, window_size = window_size, min_MAF =min_MAF, p1 = c(p1, p2), p12 = p12,
                        get_dosage_fn = dosage_fn, flexible_input=TRUE)
}
