# POEMColoc-rpackage
Implements POEMColoc colocalization method described in King et al. 2020.  This method is based on the coloc R package, but takes as input datasets that may have full summary statistics, or may have summary statistics for a single top SNP.

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
POEMColoc(gwas_stats, PDE4A, R2= R2, MAF=MAF)
```

### Colocalization for one GWAS dataset with multiple overlapping QTL
POEMColoc can also accept a list of datasets as the second input.  Compute colocalization between the GWAS association and all overlapping QTL.
```
POEMColoc(gwas_stats, eQTL_stats, R2= R2, MAF=MAF)
```
This approach does not improve computational time when we supply R2 and MAF to the function, but it's advantage when using a reference panel will be discussed in the following sections.

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
POEMColoc(gwas_stats, PDE4A, gds_file = 'chr19.1kg.phase3.v5a.gds')
```


Colocalization of PDE4A using European-only reference panel.
```
POEMColoc(gwas_stats, PDE4A, gds_file = 'chr19.1kg.phase3.v5a.gds', subset = 'eur_unrelateds.txt')
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

