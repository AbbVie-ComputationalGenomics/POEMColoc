\name{gwas_stats}
\alias{gwas_stats}
\docType{data}
\title{
GWAS summary statistics POEMColoc-formatted dataset
}
\description{
A POEMColoc-formmated dataset containing a GWAS association.  Here, we do not have full summary statistics available, so we include the top SNP position and chromosone (in hg19) and the p-value.  USe of hg19 is not required, but the build used should be consistent between all datasets and the reference panel.  We also supply the case fraction s and the sample size N, and specify this is a case-control study.
}
\usage{data("gwas_stats")}
\format{
  The format is:
List of 6
 $ type   : "cc" or "quant" (here "cc")
 $ pvalues: top SNP p-value
 $ chr    : chromosome
 $ pos    : position
 $ N      : sample size
 $ s      : case fraction
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Buniello A, MacArthur JAL, Cerezo M, Harris LW, Hayhurst J, Malangone C, McMahon A, Morales J, Mountjoy E, Sollis E, Suveges D, Vrousgou O, Whetzel PL, Amode R, Guillen JA, Riat HS, Trevanion SJ, Hall P, Junkins H, Flicek P, Burdett T, Hindorff LA, Cunningham F and Parkinson H.
The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019.
Nucleic Acids Research, 2019, Vol. 47 (Database issue): D1005-D1012.

Original data source: Ellinghaus D, et al, 2016. www.ncbi.nlm.nih.gov/pubmed/26974007
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(gwas_stats)
}
\keyword{datasets}
