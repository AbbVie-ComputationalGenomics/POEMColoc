\name{eQTL_stats}
\alias{eQTL_stats}
\docType{data}
\title{
List of whole blood eQTL summary statistics overlapping a single GWAS locus
}
\description{
GTEx whole blood eQTL with summary statistics overlapping chr19:10600418. 
A list of POEMColoc input datasets with the names of the list giving the gene.
Because full summary statistics were available, we used coefficients and variannces.
We supply sdY = 1 because the values have been standardized and hg19 positions.
}
\usage{data("eQTL_stats")}
\format{
  The format is:
List of 72 with the following format for each element
 $ ensembl_id:List of 7
  ..$ type   : "cc" or "quant", here "quant"
  ..$ sdY    : 1
  ..$ MAF    : vector of minor allele frequency at pos
  ..$ beta   : vector of beta estimate at pos
  ..$ varbeta: vector of variance of beta estimate at pos
  ..$ pos    : position
  ..$ chr    : chromosome
 
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
GTEx v7.  GTEx  Consortiumet al.,  “The  genotype-tissue  expression  (GTEx)  pilot  analysis:   Multitissue  generegulation in humans,”Science, vol. 348, no. 6235, pp. 648–660, 2015.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(eQTL_stats)
}
\keyword{datasets}
