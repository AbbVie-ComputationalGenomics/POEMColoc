\name{R2}
\alias{R2}
\docType{data}
\title{
LD matrix
}
\description{
LD matrix to use in vignette, computed from 1000 genomes Phase 3 European ancestry individuals.
Row names are top SNP positions, column names are other positions to be imputed.
}
\usage{data("R2")}
\format{
  The format is:
 num [1:62, 1:22381] 1.05e-04 1.02e-04 8.40e-04 3.35e-05 4.14e-03 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:62] "10600418" "10542807" "10449392" "10541384" ...
  ..$ : chr [1:22381] "8649879" "8649943" "8650263" "8650651" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(R2)
## maybe str(R2) ; plot(R2) ...
}
\keyword{datasets}
