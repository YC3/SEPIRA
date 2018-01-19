#' @title Infer TF activity from gene expression/ DNA methylation profile
#'
#' @description \code{sepiraRegAct} calculates TF activity scores in user input data set. It could be a gene expression dataset or a DNA methylation dataset
#'
#' @param data A gene expression or DNA methylation data matrix, with rows referring to genes and columns to samples.
#' @param type A character, "mRNA" for gene expression data; "DNAm" for DNA methylation data.
#' @param regnet The regulatory network inferred from \code{sepiraInfNet} function.
#' @param norm The method used to normalize your input data set, "c" for "centering"; "z" for "z-score normalization".
#' @param ncores The number of cores to use. See \code{\link[parallel]{mclapply}}.
#'
#' @return A matrix of TF activity score with rows referring to TFs, columns to samples.
#'
#' @details \code{sepiraRegAct} is one of the two main functions in \code{SEPIRA} package. It takes the output regulatory network from \code{sepiraInfNet} as input, and computes the activity of all TFs in this network from user provided \code{data}.
#'
#' The \code{data} matrix could be gene expression data or DAN methylation data, with rows are genes and columns are samples. Duplicated row names are not allowed, so you should average the these rows before running \code{sepiraRegAct}.
#'
#' Note that it's very important that you use the same gene identifier through out the whole analysis.
#'
#' @import parallel
#' @importFrom stats cor
#' @importFrom stats lm
#' @importFrom stats model.matrix
#' @importFrom stats pbinom
#' @importFrom stats pnorm
#' @importFrom stats sd
#'
#' @export sepiraRegAct
#'
#' @examples
#' # gene expression dataset
#' data("GeneExp")
#' # TFs
#' data("TFeid")
#' # run the function
#' cf <- "Blood"
#' coln <- colnames(GeneExp)
#' degth <- c(0.3,0.3) # 'degth = c(0.05, 0.05)' is recommended
#' net.o <- sepiraInfNet(GeneExp,coln,"Lung",cf,TFeid,sigth=0.05,degth=degth,minNtgts=5,ncores=1)
#' # normalized LSCC DNAm data set from TCGA
#' data("LUSCmeth")
#' # estimate TF activity
#' TFact.lusc <- sepiraRegAct(LUSCmeth,type="DNAm",regnet=net.o$netTOI,norm="z",ncores=1)
#' TFact.gtex <- sepiraRegAct(GeneExp,type="exp",regnet=net.o$netTOI,norm="z",ncores=1)

sepiraRegAct <- function(data, type = c("mRNA", "DNAm"), regnet, norm = c("c", "z"), ncores = 4) {

  if (type == "DNAm") {
    regnet <- -regnet
  }

  if (is.vector(data)) {
    inter <- intersect(names(data), rownames(regnet))
    data <- data[inter]
    regnet <- regnet[inter, ]
    actTF <- InferTFact(data, regnet)
    names(actTF) <- colnames(regnet)
  }
  else if (is.matrix(data)) {
    inter <- intersect(rownames(data), rownames(regnet))
    data <- data[inter, ]
    regnet <- regnet[inter, ]
    ndata <- data - rowMeans(data)
    if (norm == "z") {
      ndata <- (data - rowMeans(data)) / apply(data, 1, sd)
    }
    idx.l <- as.list(1:ncol(data))
    prl.o <- mclapply(idx.l, InferTFactPRL, ndata, regnet, mc.cores = ncores)
    actTF <- matrix(unlist(prl.o), nrow = ncol(regnet), ncol=  length(prl.o), byrow = FALSE)
    rownames(actTF) <- colnames(regnet)
    colnames(actTF) <- colnames(data)
  }

  return(actTF)
}
