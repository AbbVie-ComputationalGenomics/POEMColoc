#' Compute R2 from gds input
#'
#' @param dataset_list
#' list of datasets formatted according to the requirements of coloc.abf with the following differences
#' pos, and chr are additionally required.
#' snp is ignored as it is determined by pos.
#' For top SNP datasets, R2 and MAF and gds_file and subset are optional, if supplied as arguments to the function.
#' For full summary statistic datasets, MAF is required.
#' imputation class either all or top
#' top_pos giving the position of the top p-value
#' @param dataset_R2_condition
#' indicates whether gds is to be used to get the MAF and R2 for the dataset.  A string containing gds means it is.
#' @param window_size
#' if NULL, window size is determined to cover all positions in dataset_list
#' if specified, covers all positions in dataset list plus window_size around top SNP dataset
#' @param get_dosage_fn
#' a function taking gds file connection and returning
#' a matrix with the genotype or imputed dosage at each position.  The function
#' should not perform any filtering or reordering of variants.  rows must correspond
#' to individuals and columns to positions.
#' @return
#' a list of datasets, in which any dataset with dataset_R2_condition containing gds string now has R2 in addition to gds_file
#' @examples
#' \dontrun{
#' dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=2, snp="2")
#' dataset_top_SNP2 <- list(pos = 4, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", imputation_class = "top", top_pos=4, snp="4")
#' dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z", top_pos = 1, imputation_class = "all", snp = as.character(1:5))
#' # Example in which gds file is not used
#' add_R2_from_gds(list(dataset_top_SNP, dataset_top_SNP, dataset_full), c('R2 from function', 'R2 from function', 'none'), 1)
#' # Example in which there is a single gds file and no subsets
#' add_R2_from_gds(list(c(dataset_top_SNP, list(gds_file = "data/example.gds")), c(dataset_top_SNP2, list(gds_file = "data/example.gds")), dataset_full),
#' c('gds from dataset', 'gds from function', 'none'), 1)
#' # Example with multiple gds files
#' add_R2_from_gds(list(c(dataset_top_SNP, list(gds_file = "data/example.gds")), c(dataset_top_SNP2, list(gds_file = "data/example2.gds")), dataset_full),
#' c('gds from dataset', 'gds from function', 'none'), 1)
#' # Example with one gds file and multiple subsets.  Note it only selects the variants once.
#' add_R2_from_gds(list(c(dataset_top_SNP, list(gds_file = "data/example.gds", subset = "data/subset.ped")),
#' c(dataset_top_SNP2, list(gds_file = "data/example.gds", subset = 'data/subset2.ped')), dataset_full),
#' c('gds from dataset', 'gds from function', 'none'), 1)
#' # Example with single subset and single gds file.  Note in this case it only selects the subset from the gds.
#' add_R2_from_gds(list(c(dataset_top_SNP, list(gds_file = "data/example.gds", subset = "data/subset.ped")),
#' c(dataset_top_SNP2, list(gds_file = "data/example.gds", subset = 'data/subset.ped')), dataset_full),
#' c('gds from dataset', 'gds from function', 'none'), 1)
#' #Example in which the gds file does not have one of the needed positions
#' dataset_full2 <- list(pos = 1:7, MAF = c(0.14, 0.15, 0.25, 0.2, 0.4, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03, 0.1, 0.2), chr= "Z", top_pos = 1, imputation_class = "all", snp = as.character(1:7))
#' add_R2_from_gds(list(c(dataset_top_SNP, list(gds_file = "data/example.gds", subset = "data/subset.ped")),
#' c(dataset_top_SNP2, list(gds_file = "data/example.gds", subset = 'data/subset.ped')),
#' c(dataset_top_SNP2, list(gds_file = "data/example.gds", subset = 'data/subset2.ped')), dataset_full, dataset_full2),
#' c('gds from dataset', 'gds from function', 'gds from function', 'none', 'none'), 1)
#' }
add_R2_from_gds <- function(dataset_list, dataset_R2_condition, window_size, get_dosage_fn = get_dosage_alt) {
  pos <- get_start_and_end_pos(dataset_list, window_size)
  has_gds <- grepl("gds", dataset_R2_condition)
  if (any(has_gds)) {
    dataset_list_by_gds <- split(dataset_list[has_gds], sapply(dataset_list[has_gds], function(x) x$gds_file))
    R2_from_gds <- vector("list", length(dataset_list_by_gds))
    names(R2_from_gds) <- names(dataset_list_by_gds)
    for (i in seq_along(dataset_list_by_gds)) {
      # if single or no subset use that
      null_subset <- sapply(dataset_list_by_gds[[i]], function(x) is.null(x$subset))
      unique_subset <- unique(unlist(lapply(dataset_list_by_gds[[i]], function(x) x$subset)))
      if ((!any(null_subset)) & (length(unique_subset) == 1)) {
        subset <- unique_subset
      } else {
        subset <- NULL
      }
      current_geno_matrix <- getSNP(gds_file = names(dataset_list_by_gds)[i], chr = dataset_list_by_gds[[i]][[1]]$chr,
                                   subset = subset, start = pos$start, end = pos$end, get_dosage_fn = get_dosage_fn)
      subset_top_snp <- sapply(dataset_list_by_gds[[i]], function(x) ifelse(is.null(x$subset), NA, x$subset))
      pos_top_snp <- sapply(dataset_list_by_gds[[i]], function(x) x$pos)
      combos <- unique(data.frame(subset = subset_top_snp, pos = pos_top_snp, stringsAsFactors = FALSE))
      R2_from_gds[[i]] <- list(R2 = vector("list", length(combos)), index = combos)
      for (j in seq_len(nrow(combos))) {
        if (!is.na(combos$subset[j])) {
          subset <- scan(combos$subset[j], what = "character", sep="\n")
          geno_matrix <- current_geno_matrix$genotype[current_geno_matrix$id %in% subset,,drop=FALSE]
        } else {
          geno_matrix <- current_geno_matrix$genotype
        }
        R2_from_gds[[i]]$R2[[j]] <- compute_R2_MAF_from_geno_df(list(genotype = geno_matrix, pos = current_geno_matrix$pos),
                                                                list(pos = combos$pos[j]))
      }
    }
  }
  # Note: this could be sped up further by looking for same top position and subset and only computing once.
  for (i in which(has_gds)) {
    fileR2 = R2_from_gds[[dataset_list[[i]]$gds_file]]
    if (is.null(dataset_list[[i]]$subset)) {
      index <- which(is.na(fileR2$index$subset) & fileR2$index$pos == dataset_list[[i]]$pos)
    } else {
      index <- which(fileR2$index$subset == dataset_list[[i]]$subset & fileR2$index$pos == dataset_list[[i]]$pos)
    }
    R2_MAF = fileR2$R2[[index]]
    dataset_list[[i]]$R2 <- R2_MAF$R2
    dataset_list[[i]]$MAF <- R2_MAF$MAF
  }
  return(dataset_list)
}
