#' Determine what the source of R2 and MAF wil lbe
#'
#' @param dataset
#' single dataset formatted according to the requirements of coloc.abf with the following differences
#' pos, and chr are additionally required.
#' snp is ignored as it is determined by pos.
#' For top SNP datasets, R2 and MAF and gds_file and subset are optional, if supplied as arguments to the function.
#' For full summary statistic datasets, MAF is required in the full coloc format, i.e. the ordering of MAF is supposed to
#' correspond to the order of positions elsewhere in the dataset.
#' @param gds_file
#' path to a gds file for the reference panel.  Can be NULL if specifying R2 and MAF, or if supplying within dataset.
#' @param R2
#' A matrix of R2 LD between SNP positions.  This should contain LD information for top SNPs in rows and LD for positions to be imputed in columns.  Can be NULL if specifying gds_file or if included in individual datasets.
#' Column and row names should be position.
#' @param MAF
#' optional input of minor allele frequencies.
#' names should be position
#' @return
#' a string giving the source of R2 and MAF
#' @examples
#' \dontrun{
#' R2_matrix <- matrix(data = c(1,0.5,0.5,1), nrow=2, ncol=2)
#' rownames(R2_matrix) <- 1:2
#' colnames(R2_matrix) <- 1:2
#' R2_vec <- c(1,0.5)
#' names(R2_vec) <- 1:2
#' MAF <- c(0.1, 0.2)
#' names(MAF) <- 1:2
#' gds_file = "data/example.gds"
#' dataset_full <- list(pos = c(1, 2), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8), chr= "Z", snp = c("1","2"), imputation_class = "all")
#' dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", snp = c("2"), imputation_class = "top")
#' par_mat <- expand.grid(type = c("top", "all"), R2 = c("dataset", "both", "function", "none"), MAF = c("dataset", "both", "function", "none"), gds = c("dataset", "both", "function", "none"), stringsAsFactors=FALSE)
#' result <- vector("list", nrow(par_mat))
#' for (i in 1:nrow(par_mat)) {
#'   if (par_mat$R2[i] == "both" | par_mat$R2[i] == "function") {
#'     R2_fn = R2_matrix
#'   } else {
#'     R2_fn <- NULL
#'   }
#'   if (par_mat$R2[i] == "both" | par_mat$R2[i] == "dataset") {
#'     R2_ds = R2_vec
#'   } else {
#'     R2_ds <- NULL
#'   }
#'   if (par_mat$MAF[i] == "both" | par_mat$MAF[i] == "function") {
#'     MAF_fn = MAF
#'   } else {
#'     MAF_fn <- NULL
#'   }
#'   if (par_mat$MAF[i] == "both" | par_mat$MAF[i] == "dataset") {
#'     MAF_ds = MAF
#'   } else {
#'     MAF_ds <- NULL
#'   }
#'   if (par_mat$gds[i] == "both" | par_mat$gds[i] == "function") {
#'     gds_fn = gds_file
#'   } else {
#'     gds_fn <- NULL
#'   }
#'   if (par_mat$gds[i] == "both" | par_mat$gds[i] == "dataset") {
#'     gds_ds = gds_file
#'   } else {
#'     gds_ds <- NULL
#'   }
#'   if (par_mat$type[i] == "top") {
#'    ds_tmp <- dataset_top_SNP
#'   } else {
#'    ds_tmp <- dataset_full
#'   }
#'   ds_tmp$gds_file <- gds_ds
#'   ds_tmp$MAF <- MAF_ds
#'   ds_tmp$R2 <- R2_ds
#'   result[[i]] <- try(R2_MAF_check_dataset(ds_tmp, gds_file = gds_fn, R2 = R2_fn, MAF = MAF_fn))
#' }
#' result_summary <- sapply(result, function(x) x[1])
#' is_error <- grepl("Error", result_summary)
#' cbind(par_mat[is_error,], err = result_summary[is_error])
#' cbind(par_mat[!is_error,], res = result_summary[!is_error])
#' }
R2_MAF_check_dataset <- function(dataset, gds_file, R2, MAF) {
  MAF_condition <- check_input_condition(dataset$MAF, MAF)
  if (dataset$imputation_class == "all") {
    if (MAF_condition == "dataset" | MAF_condition == "both") {
      class <- "MAF from dataset"
    }
    if (MAF_condition == "none" | MAF_condition == "function") {
      stop("need to supply MAF within the dataset for full summary statistic datasets")
    }
  } else {
    R2_condition <- check_input_condition(dataset$R2, R2)
    if (!R2_condition == MAF_condition) {
      stop("Please supply R2 and MAF consistently either through a dataset or through function arguments when using top SNP datasets")
    }
    if (R2_condition == "both") {
      warning("R2 supplied from both dataset and function, using dataset")
    }
    MAF_condition[MAF_condition == "both"] <- "dataset"
    R2_condition[R2_condition == "both"] <- "dataset"
    if (R2_condition == "dataset") {
      class <- "R2 from dataset"
    } else {
      gds_condition <- check_input_condition(dataset$gds_file, gds_file)
      if (gds_condition == "both") {
        warning("gds from dataset overrides R2, MAF, and gds from function")
        gds_condition <- "dataset"
      }
      if (R2_condition == "none") {
        if (gds_condition == "none") {
          stop("Need to supply R2 and MAF or gds either as function arguments or in the dataset for top SNP datasets")
        } else {
          class <- paste("gds from", gds_condition)
        }
      }
      if (R2_condition == "function") {
        if (gds_condition == "dataset") {
          class <- "gds from dataset"
        } else{
          class <- "R2 from function"
        }
      }
    }
  }
  return(class)
}
