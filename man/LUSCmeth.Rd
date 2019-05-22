\name{LUSCmeth}
\alias{LUSCmeth}
\docType{data}
\title{
LUSC DNA methylation dataset
}
\description{
\code{LUSCmeth} is a subset of TCGA LUSC DNA methylation dataset. It contains 316 genes and 370 samples with 41 normal sample and 275 cancer samples.
}
\usage{data("LUSCmeth")}
\format{
  The format is:
 num [1:333, 1:316] 0.154246 0.197384 0.000585 -0.033967 -0.000737 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:333] "100130417" "100287284" "10105" "1016" ...
  ..$ : chr [1:316] "TCGA-43-3394-11A" "TCGA-18-3417-11A" "TCGA-18-4721-11A" "TCGA-22-4599-11A" ...
}
\source{
https://cancergenome.nih.gov/
}
\examples{
data(LUSCmeth)
}
\keyword{datasets}
