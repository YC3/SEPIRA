\name{GeneExp}
\alias{GeneExp}
\docType{data}
\title{
RNA-seq data matrix from GTEx
}
\description{
\code{GeneExp} is a gene expression data matrix with rows (467 genes) and 520 columns. This data matrix is a subset of GTEx RNA-seq data which is a very big dataset. In this subset we only randomly chose 20 samples from each tissue type with more than 50 samples in the original GTEx data set.
}
\usage{data("GeneExp")}
\format{
  The format is:
num [1:467, 1:520] 0.0361 0.758 0.5484 2.3569 6.1973 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:467] "9496" "22797" "7080" "6943" ...
  ..$ : chr [1:520] "Adipose Tissue" "Adipose Tissue" "Adipose Tissue" "Adipose Tissue" ...
}
\source{
https://www.gtexportal.org/home/
}
\examples{
data(GeneExp)
## view the first 5 rows and columns
GeneExp[1:5,1:5]
}
\keyword{datasets}
