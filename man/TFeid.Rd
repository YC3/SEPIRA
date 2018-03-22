\name{TFeid}
\alias{TFeid}
\docType{data}
\title{
Transcription factor genes annotated using Entrez gene ID
}
\description{
\code{TFeid} is a vector of transcription factor Entrez IDs. This list of TFs is from MSigDB with 1385 TFs in it. It is of different length as the one we used in our paper because of different "MSigDB" versions.
}
\usage{data("TFeid")}
\format{
  A vecto with 1385 observations.
  chr [1:1385] "9684" "59269" "2965" "6935" "9247" "3661" "6239" "11151" ...
}
\examples{
data(TFeid)
## maybe str(TFeid) ; plot(TFeid) ...
}
\keyword{datasets}
