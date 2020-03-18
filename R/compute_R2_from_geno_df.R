compute_R2_MAF_from_geno_df <- function(gdsout, dataset) {
  if (is.null(gdsout$genotype)) {
    R2 <- NA
  } else {
    if (any(duplicated(gdsout$pos))) {
      warning("Dropping duplicated positions")
      keep <- !duplicated(gdsout$pos)
    } else {
      keep <- TRUE
    }
    site_total <- colSums(gdsout$genotype)
    variable_sites = (site_total > 0) & (site_total < 2 * nrow(gdsout$genotype))
    keep <- keep & variable_sites
    geno_matrix <- gdsout$genotype[,keep,drop=FALSE]
    pos <- gdsout$pos[keep]
    top_SNP_pos <- dataset$pos
    top_SNP_index <- which(pos == top_SNP_pos)
    if (length(top_SNP_index) == 1) {
      R2 <- as.numeric(cor(as.numeric(geno_matrix[,top_SNP_index]), (geno_matrix))^2)
      names(R2) <- pos
    } else {
      R2 <- NA
    }
    MAF <- colMeans(geno_matrix) / 2
    MAF <- pmin(MAF, 1 - MAF)
    names(MAF) <- pos
  }
  return(list(R2 = R2, MAF = MAF))
}
