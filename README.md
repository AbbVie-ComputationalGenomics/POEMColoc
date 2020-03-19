# POEMColoc-rpackage
Implements POEMColoc colocalization method described in King et al. 2020.  This method is based on the coloc R package, but takes as input datasets that may have full summary statistics, or may have summary statistics for a single top SNP.

# Install

Install dependencies
```
install.packages("coloc")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SeqArray")
```

Install POEMColoc
```
git clone https://pig.abbvienet.com/kingea/POEMColoc-rpackage
cd POEMColoc
R CMD build POEMColoc
R CMD INSTALL POEMColoc_1.0.0.tar.gz
```

# Vignette
### Datasets
We have supplied example GWAS Catalog and GTEx eQTL datasets in this R package.  The GWAS Catalog dataset comes from one row of the GWAS Catalog reporting a pleiotropic locus on chromosome 19 for multiple chronic inflammatory diseases from Ellighaus et al. Nature Genetics.  The eQTL dataset is from GTEx whole blood eQTL for genes with summary statistics overlapping this locus.  We have provided LD and minor allele frequency for these samples, and then later show how to do the same analysis using a reference panel.
```
library(POEMColoc)
data(gwas_stats)
data(eQTL_stats)
data(R2)
data(MAF)
```

### Colocalization for two datasets
Compute colocalization between the GWAS association and PDE4A whole blood eQTL.
```
PDE4A <- eQTL_stats[["ENSG00000065989.11"]]
PDE4A_coloc <- POEMColoc(gwas_stats, PDE4A, R2= R2, MAF=MAF)
```
Notice that we get a list of colocalization results available with one element.  We can look at the hypothesis posterior probabilities

```
PDE4A_coloc[[1]]$summary
```

### Colocalization for one GWAS dataset with multiple overlapping QTL
POEMColoc can also accept a list of datasets as the second input.  Compute colocalization between the GWAS association and all overlapping QTL.
```
all_coloc <- POEMColoc(gwas_stats, eQTL_stats, R2= R2, MAF=MAF)
```
We get back a list of colocalization results of the same length as eQTL_stats.  This approach does not improve computational time when we supply R2 and MAF to the function, but it's advantage when using a reference panel will be discussed in the following sections.  We can extract colocalization probabilties by gene like this

```
sapply(all_coloc, function(x) x$summary["PP.H4.abf"])
```

### Colocalization using a reference panel
Download 1000 genomes reference panel for chromosome 19.  Also download sample info.
```
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr19.1kg.phase3.v5a.vcf.gz
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples.20130502.ALL.ped
```

Convert the vcf to a gds.  This will take a few minutes.
```
library(SeqArray)
seqVCF2GDS('chr19.1kg.phase3.v5a.vcf.gz', 'chr19.1kg.phase3.v5a.gds')
```

Create European subset file in R
```
sample_info <- read.delim('integrated_call_samples.20130502.ALL.ped', stringsAsFactors=FALSE)
eur <- sample_info$Population %in% c('CEU','TSI','FIN','GBR','IBS')
unrelateds <-  (sample_info$Paternal.ID == 0) & (sample_info$Maternal.ID == 0) & (sample_info$Second.Order == 0) & (sample_info$Siblings == 0)
eur_unrelateds <- sample_info$Individual.ID[eur & unrelateds]
write.table(eur_unrelateds, file = "eur_unrelateds.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
```

Colocalization of PDE4A using full 1000 genomes, and the gds file we just created to compute R2 and MAF.
```
library(POEMColoc)
data(gwas_stats)
data(eQTL_stats)
PDE4A <- eQTL_stats[["ENSG00000065989.11"]]
PDE4A_coloc <- POEMColoc(gwas_stats, PDE4A, gds_file = 'chr19.1kg.phase3.v5a.gds')[[1]]$summary
```

Colocalization of PDE4A using European-only reference panel.
```
PDE4A_coloc_eur <- POEMColoc(gwas_stats, PDE4A, gds_file = 'chr19.1kg.phase3.v5a.gds', subset = 'eur_unrelateds.txt')[[1]]$summary
```

Colocalization of all overlapping whole blood eQTL using European reference panel.
```
POEMColoc(gwas_stats, eQTL_stats, gds_file = 'chr19.1kg.phase3.v5a.gds', subset = 'eur_unrelateds.txt')
```
Notice that we only extract data using SeqArray once in this case.  This saves time as can be seen by the following

```
pt <- proc.time()
coloc1 <- POEMColoc(gwas_stats, eQTL_stats, gds_file = 'chr19.1kg.phase3.v5a.gds', subset = 'eur_unrelateds.txt')
time1 <- proc.time() - pt

coloc2 <- vector("list", length(eQTL_stats))
pt <- proc.time()
for (i in 1:length(eQTL_stats)) {
  coloc2[[i]] <- POEMColoc(gwas_stats, eQTL_stats[[i]], gds_file = 'chr19.1kg.phase3.v5a.gds', subset = 'eur_unrelateds.txt')[[1]]
}
time2 <- proc.time() - pt
```

In general, as long as memory is not a problem, it is advantageous to compute colocalization for multiple molecular eQTL with the same GWAS locus with a single function call.

### Colocalization when neither dataset has full summary statistics

To illustrate this feature we will artificially reduce the PDE4A eQTL dataset to a single position.  We specify a window size around either top SNP to impute.

```
top_pos_PDE4A <- which.max(abs(PDE4A$beta / sqrt(PDE4A$varbeta)))
PDE4A_top <- list(beta = PDE4A$beta[top_pos_PDE4A], varbeta = PDE4A$varbeta[top_pos_PDE4A], pos = PDE4A$pos[top_pos_PDE4A], sdY = PDE4A$sdY, type =PDE4A$type, chr = PDE4A$chr, N = PDE4A$N)
PDE4A_POEMColoc2 <- POEMColoc(gwas_stats, PDE4A_top, gds_file = 'chr19.1kg.phase3.v5a.gds', subset = 'eur_unrelateds.txt', window_size = 10^4)[[1]]$summary
```
In this case, the estimate is close to what we obtained before.
We can use multiple top SNP only datasets as well.
```
eQTL_top_only <- lapply(eQTL_stats, function(x) 
  {ii <- which.max(abs(x$beta / sqrt(x$varbeta))); list(beta = x$beta[ii], varbeta = x$varbeta[ii], 
                                                        pos = x$pos[ii], chr=x$chr, sdY=x$sdY, type=x$type, N=x$N)})
POEMColoc2_all <- POEMColoc(gwas_stats, eQTL_top_only, gds_file = 'chr19.1kg.phase3.v5a.gds', subset = 'eur_unrelateds.txt', window_size = 10^4)
```

Note that some of the colocalization results are NA because the top SNP was not found.  We can extract colocalization probabiltiies by gene like this

```
sapply(POEMColoc2_all, function(x) ifelse(isTRUE(all.equal(x, NA)), NA, x$summary["PP.H4.abf"]))
```

Using the output from earlier in the vignette, we could compare what we get using one and two top-SNP only datasets (POEMColoc-1 and POEMColoc-2)

```
summary(sapply(all_coloc, function(x) x$summary["PP.H4.abf"]) - sapply(POEMColoc2_all, function(x) ifelse(isTRUE(all.equal(x, NA)), NA, x$summary["PP.H4.abf"])))
```
