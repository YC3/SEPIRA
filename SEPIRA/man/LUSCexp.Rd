\name{LUSCexp}
\alias{LUSCexp}
\docType{data}
\title{
LUSC RNA-seq dataset
}
\description{
\code{LUSCmeth} is a subset of TCGA LUSC RNA dataset. It contains 391 genes and 518 samples with 45 normal sample and 473 cancer samples.
}
\usage{data("LUSCexp")}
\format{
  The format is:
num [1:391, 1:518] -0.0842 0.7222 1.1949 -1.1637 0.2717 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:391] "100009676" "100129842" "100130417" "100131726" ...
  ..$ : chr [1:518] "TCGA-22-5481-11" "TCGA-56-7222-11" "TCGA-77-7338-11" "TCGA-77-7138-11" ...
}
\source{
https://cancergenome.nih.gov/
}
\examples{
data(LUSCexp)
}
\keyword{datasets}
