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
#' @examples
#' # R2 and MAF from function
#' dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=2, snp="2")
#' dataset_full <- list(pos = c(1, 2, 4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 0.01), chr= "Z", snp = c("1","2","4"), imputation_class = "all", MAF = c(0.2,0.4, 0.3))
#' R2 <- matrix(data = c(0.94, 1, 0.65, 0, 0.51, 0.32, 0.64, 0.65, 1, 0, 0.21, 0.12), nrow=2, byrow=TRUE)
#' colnames(R2) <- 1:6
#' rownames(R2) <- 2:3
#' MAF = c(0.14, 0.15, 0.25, 0.2, 0.4, 0.25)
#' names(MAF) <- 1:6
#' add_R2_to_dataset_list(list(dataset_full, dataset_top_SNP), R2 = R2, MAF = MAF, subset = NULL, window_size= 2, gds_file=NULL)
#' # R2 and MAF from function some positions missing from to impute
#' add_R2_to_dataset_list(list(dataset_full, dataset_top_SNP), R2 = R2[,c("1","2")], MAF = MAF[c("1","2")], subset = NULL, window_size= 2, gds_file=NULL)
#' add_R2_to_dataset_list(list(dataset_full, dataset_top_SNP), R2 = R2[,c("1","2")], MAF = MAF[c("2")], subset = NULL, window_size= 2, gds_file=NULL)
#' # R2 and MAF from function top position missing
#' add_R2_to_dataset_list(list(dataset_full, dataset_top_SNP), R2 = R2[c("3"),,drop=FALSE], MAF = MAF, subset = NULL, window_size= 2, gds_file=NULL)
#' # R2 and MAF from input
#' dataset_top_SNP2 <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=2, snp="2", R2=R2, MAF=MAF)
#' add_R2_to_dataset_list(list(dataset_full, dataset_top_SNP2), gds_file = NULL, subset = NULL, R2 = NULL, MAF = NULL, window_size= 2)
#' dataset_top_SNP3 <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=2, snp="2", R2=R2["2",], MAF=MAF)
#' add_R2_to_dataset_list(list(dataset_full, dataset_top_SNP3), gds_file = NULL, subset = NULL, R2 = NULL, MAF = NULL, window_size= 2)
#' # R2 and MAF from input some positions missing
#' dataset_top_SNP4 <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=2, snp="2", R2=R2["2",c("1","2")], MAF=MAF[c("1","2")])
#' add_R2_to_dataset_list(list(dataset_full, dataset_top_SNP4), gds_file = NULL, subset = NULL, R2 = NULL, MAF = NULL, window_size= 2)
#' gds_file <- system.file("extdata", "example.gds", package = "POEMColoc")
#' gds_file2 <- system.file("extdata", "example2.gds", package = "POEMColoc")
#' subset <- system.file("extdata", "subset.ped", package = "POEMColoc")
#' subset2 <- system.file("extdata", "subset2.ped", package = "POEMColoc")
#' # R2 and MAF from gds function
#' dataset_top_SNP5 <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=2, snp="2")
#' add_R2_to_dataset_list(list(dataset_full, dataset_top_SNP5), gds_file = gds_file, subset = subset2, R2 = NULL, MAF = NULL, window_size= 2)
#' add_R2_to_dataset_list(list(dataset_full, dataset_top_SNP5), gds_file = gds_file, subset = subset2, R2 = NULL, MAF = NULL, window_size= 3)
#' # R2 and MAF from gds function some positions missing
#' dataset_full2 <- list(pos = c(1, 2, 7), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 0.01), chr= "Z", snp = c("1","2","7"), imputation_class = "all", MAF = c(0.2,0.4, 0.3))
#' add_R2_to_dataset_list(list(dataset_full2, dataset_top_SNP5), gds_file = gds_file, subset = subset2, R2 = NULL, MAF = NULL, window_size= 1)
#' # R2 and MAF from gds dataset
#' dataset_top_SNP6 <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=2, snp="2", gds_file = gds_file)
#' add_R2_to_dataset_list(list(dataset_full2, dataset_top_SNP6), gds_file = gds_file, subset = subset2, R2 = NULL, MAF = NULL, window_size= 1)
#' # MAF from function should fail for full summary statistics - not supported
#' dataset_full3 <- list(pos = c(1, 2, 4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 0.01), chr= "Z", snp = c("1","2","4"), imputation_class = "all")
#' try(add_R2_to_dataset_list(list(dataset_full3, dataset_top_SNP6), gds_file = gds_file, subset = subset2, R2 = NULL, MAF = MAF, window_size= 1))
#' # Should fail due to missing MAF
#' dataset_full4 <- list(pos = c(1, 2, 7), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 0.01), chr= "Z", snp = c("1","2","7"), imputation_class = "all")
#' try(add_R2_to_dataset_list(list(dataset_full4, dataset_top_SNP6), gds_file = gds_file, subset = subset2, R2 = NULL, MAF = NULL, window_size= 1))
#' # Same datasets in a list
#' add_R2_to_dataset_list(list(dataset_full, dataset_full2, dataset_top_SNP, dataset_top_SNP2, dataset_top_SNP3, dataset_top_SNP4, dataset_top_SNP5, dataset_top_SNP6), gds_file = gds_file, subset = subset2, R2 = R2, MAF = MAF, window_size= 3)
#' # complex input with different gds
#' add_R2_to_dataset_list(list(c(dataset_top_SNP, list(gds_file = gds_file, subset = subset)),
#' c(dataset_top_SNP, list(gds_file = gds_file, subset = subset2)),
#' c(dataset_top_SNP, list(gds_file = gds_file2, subset = subset)), dataset_full, dataset_full2), R2 = NULL, MAF = NULL, gds_file = NULL, subset=NULL, window_size =1)

add_R2_to_dataset_list <- function(dataset_list, gds_file, R2, MAF, subset, window_size, get_dosage_fn = get_dosage_alt) {
  # Check what will be the source of LD information for each dataset
  dataset_R2_condition <- sapply(dataset_list, function(x) R2_MAF_check_dataset(x,  gds_file, R2, MAF))
  # add that information at the dataset level
  add_R2_from_R2_MAF_check(dataset_list = dataset_list, gds_file = gds_file, R2 = R2,
                           MAF = MAF, subset = subset, window_size = window_size,
                           get_dosage_fn = get_dosage_fn, dataset_R2_condition = dataset_R2_condition)
}
