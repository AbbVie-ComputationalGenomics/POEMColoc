# POEMColoc-rpackage
Implements POEMColoc colocalization method described in King et al. 2020.  This method is based on the coloc R package, but takes as input datasets that may have full summary statistics, or may have summary statistics for a single top SNP.

# Vignette
### Datasets
We have supplied example GWAS Catalog and GTEx eQTL datasets in this R package.  The GWAS Catalog dataset comes from one row of the GWAS Catalog reporting a pleiotropic locus on chromosome 19 for multiple chronic inflammatory diseases from Ellighaus et al. Nature Genetics.  The eQTL dataset is from GTEx whole blood eQTL for genes with summary statistics overlapping this locus.  We have provided LD and minor allele frequency for these samples, and then later show how to do the same analysis using a reference panel.
### Colocalization for two datasets
### Colocalization for one GWAS dataset with multiple overlapping QTL
### Colocalization using a reference panel
Download 1000 genomes reference panel for chromosome 19.  Also download sample info.
```
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr19.1kg.phase3.v5a.vcf.gz
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples.20130502.ALL.ped
```

Create European subset file in R
```
sample_info <- read.delim('integrated_call_samples.20130502.ALL.ped', stringsAsFactors=FALSE)
eur <- sample_info$Population %in% c('CEU','TSI','FIN','GBR','IBS')
unrelateds <-  (sample_info$Paternal.ID == 0) & (sample_info$Maternal.ID == 0) & (sample_info$Second.Order == 0) & (sample_info$Siblings == 0)
eur_unrelateds <- sample_info$Individual.ID[eur & unrelateds]
write.table(eur_unrelateds, file = "eur_unrelateds.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
```
### Colocalization when neither dataset has full summary statistics

