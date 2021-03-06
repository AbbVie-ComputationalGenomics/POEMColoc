% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/single_dataset_check.R
\name{single_dataset_check}
\alias{single_dataset_check}
\title{Check a single dataset for POEM-specific requirements and fail if not met}
\usage{
single_dataset_check(dataset)
}
\arguments{
\item{dataset}{single POEMColoc dataset}
}
\description{
Check a single dataset for POEM-specific requirements and fail if not met
}
\examples{
dataset_full <- list(pos = c(1, 2, 3, 4, 5), MAF = c(0.14, 0.15, 0.25, 0.2, 0.4), N=1000, type ="quant", pvalues = c(2 * 10^-8, 4 * 10^-8, 2 * 10^-4, 0.6, 0.03), chr= "Z")
dataset_top_SNP <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
dataset_top_SNP_no_chr <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9)
dataset_top_SNP_no_pos <- list(N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z")
dataset_top_SNP_have_snp <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = "Z", snp = "abc")
dataset_top_SNP_too_many_chr <- list(pos = c(1,2), N= 10000, s =0.5, type="cc", pvalues = c(10^-8,10^-9), chr = c(1,2))
dataset_top_SNP_too_few_chr <- list(pos = 2, N= 10000, s =0.5, type="cc", pvalues = 10^-9, chr = c())
dataset_top_SNP_no_pvalues <- list(pos = 2, N= 10000, s =0.5, type="cc", chr = "Z")
# Should pass because correctly formatted
single_dataset_check(dataset_full)
single_dataset_check(dataset_top_SNP)
# Should pass because misformatting will be detected later
single_dataset_check(dataset_top_SNP_no_pvalues)
# Should warn
single_dataset_check(dataset_top_SNP_have_snp)
# Should fail
try(single_dataset_check(dataset_top_SNP_no_chr))
try(single_dataset_check(dataset_top_SNP_no_pos))
try(single_dataset_check(dataset_top_SNP_too_many_chr))
try(single_dataset_check(dataset_top_SNP_too_few_chr))

}
