\name{VALexp}
\alias{VALexp}
\docType{data}
\title{
RNA-seq data matrix from ProteinAtlas project
}
\description{
\code{VALexp} is a gene expression data matrix with rows (369 genes) and 35 columns (samples from 7 tissue types). This data matrix was obtained from EBI arrayexpress and is part of ProteinAtlas project. It is a very big dataset, so we only used a small subset as an example dataset here.
}
\usage{data("VALexp")}
\format{
  The format is:
 num [1:467, 1:200] 0 2.379 0 0.379 5.392 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:467] "9496" "22797" "7080" "6943" ...
  ..$ : chr [1:200] "saliva-secreting gland" "lung" "liver" "heart" ...
}
\source{
https://www.ebi.ac.uk/arrayexpress/
}
\examples{
data(VALexp)
## view the first 5 rows and columns
VALexp[1:5,1:5]
}
\keyword{datasets}
