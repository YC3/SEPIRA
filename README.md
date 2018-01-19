# SEPIRA-package
Systens EPigenomics Inference of Regulatory Activity

`SEPIRA` is a novel algorithm which estimates transcription factor activity in any given sample from its genome-wide mRNA expression or DNA methylation profile[@Chen2017]. It encompasses two main steps:

1. Construction of a tissue-specific transcription factor regulatory network, consisting of transcription factors that are more highly expressed in the user-specified tissue type (the 'tissue type of interest') compared to other tissue types, plus an associated set of high-confidence downstream targets.
2. Estimation of transcription factor activity in this network, in any given dataset consisting of gene expression or promoter DNA methylation profiles.

To infer the network, we use the large [GTEx](https://www.gtexportal.org/home/) RNA-seq data set encompassing 8555 samples from about 30 tissue types. With this inferred network, we can then estimate transcription factor activity in samples from other data sets. Due to the large size of the GTEx data set we did not include it in our package. However, in order to gain an appreciation for the SEPIRA algorithm we show the results obtained by applying it to this full dataset. Then we validate it using an RNA-seq data set from Protein Atlas and a DNA methylation data set from the Stem-cell matrix compendium-2 (SCM2).

## Usage
```{r echo=FALSE}
knitr::opts_chunk$set(fig.pos="h", out.extra='', fig.align="center")
```

#### Inferring tissue-specific network

```{r eval=FALSE}
net.o <- sepiraInfNet(data=data.m, tissue=colnames(data.m), toi = "Lung", cft = "Blood",
         TFs = TFeid, sdth = 0.25, sigth = 0.05, pcorth = 0.2, degth = c(0.05, 0.05),
         lfcth = c(log2(1.5), 0), minNtgts = 3, ncores = 1)
```
** Note: `data.m` should be a normalized gene expression data set.

#### Estimating transcription factor activity
```{r eval=FALSE}
sepiraRegAct(data = data.m, type = "DNAm", regnet = net.o$netTOI, norm = "z", ncores = 1)
```

## Installation

An easy way to install SEPIRA is by facilitating the devtools R package.

#install.packages("devtools")
library(devtools)
install_github("YC3/SEPIRA", build_vignettes=TRUE)

Alternatively, the package can also be cloned or downloaded from this github-rep, built via R CMD build and installed via the R CMD INSTALL command.

## Getting started
The SEPIRA package contains a tutorial showing people how to implement SEPIRA in their work. The tutorial can be found in the package-vignette:

library(SEPIRA)
vignette("SEPIRA")

## Acknowledgements

Thanks to my supervisor Andrew Teschendorff for reviewing and commenting on the package and providing the code.
