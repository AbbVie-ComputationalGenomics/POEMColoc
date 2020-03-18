#' Get genotypes in region
#'
#' @param gds_file
#' Path to gds file
#' @param chr
#' Chromesome id
#' @param start
#' start position of region on chromosome
#' @param end
#' end position of region on chromosome
#' @param subset
#' optional samples file for input to bcftools
#' @param get_dosage_fn
#' a function taking gds file connection and returning
#' a matrix with the genotype or imputed dosage at each position.  The function
#' should not perform any filtering or reordering of variants.  rows must correspond
#' to individuals and columns to positions.
#' @return
#' A data.frame of alleles
#' @examples
#'
getSNP <- function(gds_file, chr, start, end, subset=NULL, get_dosage_fn = get_dosage_alt) {
  gds = seqOpen(gds_file)
  on.exit(seqClose(gds))
  seqSetFilterChrom(gds, chr, from.bp = start, to.bp = end)
  if (!is.null(subset)) {
    subset <- scan(subset, what = "character", sep="\n")
    seqSetFilter(gds, sample.id = subset)
  }
  position = seqGetData(gds, "position")
  sample_id = seqGetData(gds, "sample.id")
  genotype = get_dosage_fn(gds)
  return(list(genotype = genotype, pos = position, id = sample_id))
}
