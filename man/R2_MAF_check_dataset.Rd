% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R2_MAF_check_dataset.R
\name{R2_MAF_check_dataset}
\alias{R2_MAF_check_dataset}
\title{Determine what the source of R2 and MAF wil lbe}
\usage{
R2_MAF_check_dataset(dataset, gds_file, R2, MAF)
}
\arguments{
\item{dataset}{single dataset formatted according to the requirements of coloc.abf with the following differences
pos, and chr are additionally required.
snp is ignored as it is determined by pos.
For top SNP datasets, R2 and MAF and gds_file and subset are optional, if supplied as arguments to the function.
For full summary statistic datasets, MAF is required in the full coloc format, i.e. the ordering of MAF is supposed to
correspond to the order of positions elsewhere in the dataset.}

\item{gds_file}{path to a gds file for the reference panel.  Can be NULL if specifying R2 and MAF, or if supplying within dataset.}

\item{R2}{A matrix of R2 LD between SNP positions.  This should contain LD information for top SNPs in rows and LD for positions to be imputed in columns.  Can be NULL if specifying gds_file or if included in individual datasets.
Column and row names should be position.}

\item{MAF}{optional input of minor allele frequencies.
names should be position}
}
\value{
a string giving the source of R2 and MAF
}
\description{
Determine what the source of R2 and MAF wil lbe
}
\examples{
R2_matrix <- matrix(data = c(1,0.5,0.5,1), nrow=2, ncol=2)
rownames(R2_matrix) <- 1:2
colnames(R2_matrix) <- 1:2
R2_vec <- c(1,0.5)
names(R2_vec) <- 1:2
MAF <- c(0.1, 0.2)
names(MAF) <- 1:2
gds_file <- system.file("extdata", "example.gds", package = "POEMColoc")
dataset_full <- list(pos = c(1, 2), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8), chr= "Z", snp = c("1","2"), imputation_class = "all")
dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", snp = c("2"), imputation_class = "top")
par_mat <- expand.grid(type = c("top", "all"), R2 = c("dataset", "both", "function", "none"), MAF = c("dataset", "both", "function", "none"), gds = c("dataset", "both", "function", "none"), stringsAsFactors=FALSE)
result <- vector("list", nrow(par_mat))
for (i in 1:nrow(par_mat)) {
  if (par_mat$R2[i] == "both" | par_mat$R2[i] == "function") {
    R2_fn = R2_matrix
  } else {
    R2_fn <- NULL
  }
  if (par_mat$R2[i] == "both" | par_mat$R2[i] == "dataset") {
    R2_ds = R2_vec
  } else {
    R2_ds <- NULL
  }
  if (par_mat$MAF[i] == "both" | par_mat$MAF[i] == "function") {
    MAF_fn = MAF
  } else {
    MAF_fn <- NULL
  }
  if (par_mat$MAF[i] == "both" | par_mat$MAF[i] == "dataset") {
    MAF_ds = MAF
  } else {
    MAF_ds <- NULL
  }
  if (par_mat$gds[i] == "both" | par_mat$gds[i] == "function") {
    gds_fn = gds_file
  } else {
    gds_fn <- NULL
  }
  if (par_mat$gds[i] == "both" | par_mat$gds[i] == "dataset") {
    gds_ds = gds_file
  } else {
    gds_ds <- NULL
  }
  if (par_mat$type[i] == "top") {
   ds_tmp <- dataset_top_SNP
  } else {
   ds_tmp <- dataset_full
  }
  ds_tmp$gds_file <- gds_ds
  ds_tmp$MAF <- MAF_ds
  ds_tmp$R2 <- R2_ds
  result[[i]] <- try(R2_MAF_check_dataset(ds_tmp, gds_file = gds_fn, R2 = R2_fn, MAF = MAF_fn))
}
result_summary <- sapply(result, function(x) x[1])
is_error <- grepl("Error", result_summary)
cbind(par_mat[is_error,], err = result_summary[is_error])
cbind(par_mat[!is_error,], res = result_summary[!is_error])
}
