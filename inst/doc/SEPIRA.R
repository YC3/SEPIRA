## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(fig.pos="h", out.extra='', fig.align="center")
library(SEPIRA)
data("GeneExp")
data.m <- GeneExp
dataVAL.m <- GeneExp
data("TFeid")

## ----warning=F-----------------------------------------------------------
net.o <- sepiraInfNet(data=data.m, tissue=colnames(data.m), toi = "Lung", cft = "Blood",
         TFs = TFeid, sdth = 0.25, sigth = 0.05, pcorth = 0.2, degth = c(0.05, 0.05),
         lfcth = c(log2(1.5), 0), minNtgts = 5, ncores = 1)
## Note: `data.m` should be a normalized gene expression data set.
## Parameters used here are not recommended. See "?sepiraInfNet" for more info.

## ------------------------------------------------------------------------
act <- sepiraRegAct(data = data.m, type = "DNAm", regnet = net.o$netTOI, norm = "z", ncores = 1)

## ----warning=F-----------------------------------------------------------
net.o <- sepiraInfNet(data=data.m, tissue=colnames(data.m), toi = "Lung", cft = "Blood",
         TFs = TFeid, sdth = 0.25, sigth = 0.05, pcorth = 0.2, degth = c(0.05, 0.05),
         lfcth = c(log2(1.5), 0), minNtgts = 3, ncores = 1)

## ----out.width = 400, fig.retina = NULL, echo=F--------------------------
knitr::include_graphics("Figures/LungNet_last_version.png")

## ------------------------------------------------------------------------
# estimating transcription factor activity in data.m
TFact_gtex.m <- sepiraRegAct(data.m,type="exp",regnet=net.o$netTOI,norm="z",ncores=1)

## ----out.width = 650, fig.retina = NULL, echo=F--------------------------
knitr::include_graphics("Figures/LungNet_box_GTEx.png")

## ------------------------------------------------------------------------
# estimating transcription factor activity in a validation dataset
TFact_val.m <- sepiraRegAct(dataVAL.m,type="exp",regnet=net.o$netTOI,norm="z",ncores=1)

## ----out.width = 650, fig.retina = NULL, echo=F--------------------------
knitr::include_graphics("Figures/LungNet_box_ProAtl.png")

## ----out.width = 300, fig.retina = NULL, echo=F--------------------------
knitr::include_graphics("Figures/SCM2.png")

## ----sessionInfo---------------------------------------------------------
sessionInfo()

