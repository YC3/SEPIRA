# SEPIRA-package
Systems EPigenomics Inference of Regulatory Activity

`SEPIRA` is a novel algorithm which estimates transcription factor activity in any given sample from its genome-wide mRNA expression or DNA methylation profile[1]. It encompasses two main steps:

1. Construction of a tissue-specific transcription factor regulatory network, consisting of transcription factors that are more highly expressed in the user-specified tissue type (the 'tissue type of interest') compared to other tissue types, plus an associated set of high-confidence downstream targets.
2. Estimation of transcription factor activity in this network, in any given dataset consisting of gene expression or promoter DNA methylation profiles.

## Usage
#### Inferring tissue-specific network

```R
net.o <- sepiraInfNet(data=data.m, tissue=colnames(data.m), toi = "Lung", cft = "Blood",
         TFs = TFeid, sdth = 0.25, sigth = 0.05, pcorth = 0.2, degth = c(0.05, 0.05),
         lfcth = c(log2(1.5), 0), minNtgts = 3, ncores = 1)
```
** Note: `data.m` should be a normalized gene expression data set.

#### Estimating transcription factor activity
```R
sepiraRegAct(data = data.m, type = "DNAm", regnet = net.o$netTOI, norm = "z", ncores = 1)
```

## Installation

An easy way to install SEPIRA is by facilitating the devtools R package.

```R
#install.packages("devtools")
devtools::install_github("YC3/SEPIRA", build_vignettes=TRUE)
```
Alternatively, the package can also be cloned or downloaded from this github-rep, built via R CMD build and installed via the R CMD INSTALL command.

## Getting started
The SEPIRA package contains a tutorial showing people how to implement SEPIRA in their work. The tutorial can be found in the package-vignette [here](https://bioconductor.org/packages/release/bioc/vignettes/SEPIRA/inst/doc/SEPIRA.html) or generate one with:

```R
vignette("SEPIRA")
```

## References

Chen Y, Widschwendter M, and Teschendorff AE. 2017. Systems-Epigenomics Inference of Transcription Factor Activity Implicates Aryl-Hydrocarbon-Receptor Inactivation as a Key Event in Lung Cancer Development. Genome Biol 18:236.
