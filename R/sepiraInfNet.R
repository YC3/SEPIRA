#' @title Infer tissue-specific regulatory network from gene expression data
#'
#' @description \code{sepiraInfNet()} is one of the two main functions in package \code{SEPIRA}. Using it you can estimate tissue-specific regulatory networks in any tissue type of interest.
#'
#' @param data The normalized gene expression data matrix, with rows referring to unique genes and columns to samples from different tissue types.
#' @param tissue A phenotype vector, indicating the tissue types of samples. It should have the same order as the columns of the matrix.
#' @param toi The tissue type of interest, a character telling the function the tissue type for which a user wants to estimate the network.
#' @param cft A vector of tissue types to be used to adjust for confounding by immune or stromal cells infiltration in \code{toi}. It can be blood and/or spleen, which we found using \code{ESTIMATE} package that they contain extremely high proportion of immune and stromal cells.
#' @param TFs A vector of TFs. Note that one should use the same annotation in different data sets throughout the analysis.
#' @param sdth The standard deviation threshold used to remove genes with little or zero standard deviation of its expression levels.
#' @param sigth The unadjusted p-value threshold used to call significant interactions after calculating the correlation coefficients between TFs and target genes. This threshold is used to binarize the correlation coefficient matrix. If this value is not specified by user, the function will do Bonferroni correction and then use 0.05 as the threshold.
#' @param pcorth The partial correlation threshold, in the range between 0 and 1, used to remove indirect interactions between TFs and their target genes.
#' @param degth A vector of length three, thresholds of adjusted p-value to call significant TFs in 1) comparison between \code{toi} and all other tissue types; 2) & 3) comparison between \code{toi} and blood/spleen in \code{cft}.
#' @param lfcth vector of length three, thresholds of log2(fold-change) to call significant TFs in 1) comparison between \code{toi} and all other tissue types; 2) & 3) comparison between \code{toi} and blood/spleen in \code{cft}.
#' @param minNtgts An integer used to filter out TFs with few targets. Only TFs with more than `minNtgts` target genes can be kept in the network.
#' @param ncores The number of cores to use when computing partial correlation. See \code{\link[parallel]{mclapply}}.
#'
#' @return The output is a list with three entries:
#'
#' \code{$netTOI} the tissue specific network, rows refer to TF target genes, while columns refer to TFs.
#'
#' \code{$sumnet} a matrix summarizing the number of TF target genes and the number of positively/negatively regulated target genes for each TF in the inferred network.
#'
#' \code{$top} a list, entries are the tables summarizing the results of differential expression analyses. The first is the table from comparison between \code{toi} and 1) all other tissue; The rest tables are resulted from comparison to 2) the blood or/ and 3) spleen.
#'
#' @details \code{sepiraInfNet} generates tissue specific TF regulatory networks from gene expression data across multi-tissue types.
#'
#' The gene expression data set \code{data} should be normalized by user before inputting to \code{sepiraInfNet}, with rows referring to genes and columns to samples from different tissue types. Duplicated gene names/IDs should be averaged before normalization.
#'
#' The user needs to input the tissue type of samples (\code{tissue}) in the data set as well as the tissue type of interest (\code{toi}). Please make sure the \code{toi} is in the \code{tissue} and spelled correctly.
#'
#' Using differential gene expression analysis, we detect TFs that are highly active in \code{toi} and less active in other tissue types. When doing such analyses, the results could be confounded due to cell-type heterogeneity. \code{sepiraInfNet} provides a way to adjust for immune/ stromal cell contamination by doing additional comparisons between \code{toi} and 1) blood; 2)spleen as long as expression data for any one/ both of the tissue types are available in \code{data}.
#'
#' \code{TFs} is a vector containing the identifiers of all TFs (regulators). In our paper we used the 1313 TFs annotated as "transcription factors" in MSigDB. You could input your own list of TFs to \code{sepiraInfNet}.
#'
#' \code{sdth} is a standard deviation threshold that is used to remove genes in user provided data set which are with small or close to zero standard deviation. By default the threshold is 0.25.
#'
#' From the gene expression data matrix \code{sepiraInfNet} estimates Pearson correlation coefficient between every TF-gene pair as well as corresponding p-value. The p-value threshold \code{sigth} binarizes the network into "regulation" (1) /"no regulation" (0). This binarized network is used to determine the covariants when estimating the partial correlation between target genes and their regulators (TFs).
#'
#' \code{pcorth} is the partial correlation coefficient threshold for calling significant direct TF-gene interactions. By default \code{pcorth} equals 0.2.
#'
#' \code{degth} and \code{lfcth} are vectors each contains the 3 thresholds for adjusted p-value/log2 fold-change to call significant TFs in comparisons between \code{toi} and 1) all other tissue types; 2) the 1st tissue type (blood) in \code{cft}; 3) the 2nd tissue type (spleen) in \code{cft}. These differential expression analyses are done to find tissue-specific TFs that are only highly activated in tissue type of interest.
#'
#' When having detected tissue-specific TFs, we could get a network with only these TFs and their target genes. However \code{sepiraInfNet} further remove TFs with less than \code{minNtgts} target genes. By default the minimal number of TF targets in the final network is 10.
#'
#' The step of calculating partial correlation coefficients is done by in parallel, by default \code{sepiraInfNet} splits the work into 4 sub-processes. User could use more cores by specifying parameter \code{ncores}.
#'
#' @import parallel
#' @import corpcor
#' @import limma
#'
#' @export sepiraInfNet
#'
#' @examples
#' # gene expression data set (a subset of GTEx data set)
#' data("GeneExp")
#' # TFs
#' data("TFeid")
#' # run the function
#' cf <- "Blood"
#' coln <- colnames(GeneExp)
#' degth <- c(0.3,0.3) # 'degth = c(0.05, 0.05)' is recommended
#' # The resulted network is small due to the limited size of the 'GeneExp' data set
#' net.o <- sepiraInfNet(GeneExp,coln,"Lung",cf,TFeid,sigth=0.05,degth=degth,minNtgts=5,ncores=1)

sepiraInfNet <- function(data, tissue, toi, cft = NULL, TFs, sdth = 0.25, sigth = NULL, pcorth = 0.2, degth = c(0.05, 0.05), lfcth = c(1, log2(1.5)), minNtgts = 10, ncores = 4) {

  tt.v <- levels(factor(tissue))
  if (length(intersect(toi,tt.v)) == 0) stop("Your tissue of interest is not in 'tissue', or you have mispelled toi")
  if (FALSE %in% (cft %in% tt.v)) stop(paste0("'",cft[!cft %in% tt.v],"'"," is not in 'tissue', or you have mispelled it."))

  ### remove genes with no or little variance
  sd.v <- apply(data, 1, sd)
  selG.idx <- which(sd.v > sdth)
  exp.m <- data[selG.idx, ]
  ### find representation of regulators in data, and define regulatees/targets
  tfEID.v <- TFs
  repTF.v <- intersect(tfEID.v, rownames(exp.m))
  tgtsEID.v <- setdiff(rownames(exp.m), tfEID.v)
  mapTF.idx <- match(repTF.v, rownames(exp.m))
  mapTGTS.idx <- match(tgtsEID.v, rownames(exp.m))
  ### compute correlations and estimate P-values
  corNET.m <- cor(t(exp.m[mapTGTS.idx, ]), t(exp.m[mapTF.idx, ]))
  zNET.m <- 0.5 * log( (1 + corNET.m) / (1 - corNET.m))
  stdev <- 1 / sqrt(length(tt.v) - 3) ### this is not the number of independent samples but the number of independent tissues. the latter is used because of the strong dependence between samples from the same tissue, but also because it leads to a more stringent significance threshold
  pvNET.m <- 2 * pnorm(abs(zNET.m), 0, stdev, lower.tail = FALSE)
  ### for each gene, now identify the TFs which are correlated univariately- for these then run multivariate regression
  if (is.null(sigth)) {
    sigth <- 0.05 / prod(dim(pvNET.m))
  }
  binNET.m <- pvNET.m
  binNET.m[pvNET.m < sigth] <- 1
  binNET.m[pvNET.m >= sigth] <- 0
  if(sum(binNET.m) > 0.01*prod(dim(pvNET.m))){### if more than 1% cap at 1%
    topE <- floor(0.01*prod(dim(pvNET.m)));

    zNET.v <- as.vector(abs(zNET.m));
    tmp.s <- sort(zNET.v,decreasing=TRUE,index.return=TRUE);
    out.v <- rep(0,length(zNET.v));
    out.v[tmp.s$ix[1:topE]] <- 1;
    binNET.m <- matrix(out.v,nrow=nrow(zNET.m));
    rownames(binNET.m) <- rownames(zNET.m);
    colnames(binNET.m) <- colnames(zNET.m);
   }

  ### number of targets per tf
  ntgTF.v <- apply(binNET.m, 2, sum)
  ### number of regulators per gene
  nregG.v <- apply(binNET.m, 1, sum)

  ### select TFs with at least minNtgts targets
  selTF.idx <- which(ntgTF.v >= minNtgts)
  selbinNET.m <- binNET.m[, selTF.idx]
  selpvNET.m <- pvNET.m[, selTF.idx]
  selzNET.m <- zNET.m[, selTF.idx]
  selcorNET.m <- corNET.m[, selTF.idx]

  mapTG.idx <- match(rownames(selbinNET.m), rownames(exp.m))
  mapTF.idx <- match(colnames(selbinNET.m), rownames(exp.m))

  idx.l <- as.list(1:nrow(selbinNET.m))
  #print("Computing Partial Correlations")
  pcor.l <- mclapply(idx.l, ComputePCOR, mapTG.idx, mapTF.idx, selbinNET.m, exp.m, mc.cores = ncores)

  pcorNET.m <- matrix(0, nrow = nrow(selbinNET.m), ncol = ncol(selbinNET.m))
  rownames(pcorNET.m) <- rownames(selbinNET.m)
  colnames(pcorNET.m) <- colnames(selbinNET.m)
  for (g in 1:length(pcor.l)) {
    reg.idx <- which(selbinNET.m[g,] == 1)
    if (length(reg.idx) >= 2) {
      pcorNET.m[g, reg.idx] <- pcor.l[[g]][1, -1]
    }
    else if (length(reg.idx) == 1) {
      pcorNET.m[g, reg.idx] <- selcorNET.m[g, reg.idx]
    }
  }

  gtexNET.m <- sign(pcorNET.m)
  gtexNET.m[abs(pcorNET.m) < pcorth] <- 0

  sumnetTF.m <- matrix(nrow = 4, ncol = ncol(gtexNET.m))
  rownames(sumnetTF.m) <- c("nTG", "nUP", "nDN", "P")
  colnames(sumnetTF.m) <- colnames(gtexNET.m)
  sumnetTF.m[1, ] <- apply(abs(gtexNET.m), 2, sum)

  for (tf in 1:ncol(gtexNET.m)) {
    sumnetTF.m[2, tf] <- length(which(gtexNET.m[, tf] == 1))
    sumnetTF.m[3, tf] <- length(which(gtexNET.m[, tf] == -1))
  }
  pv.v <- pbinom(apply(sumnetTF.m[2:3, ], 2, max), size=sumnetTF.m[1, ], prob = 0.5, lower.tail = FALSE)
  sumnetTF.m[4, ] <- pv.v

  gtexNETf.m <- gtexNET.m[, which(sumnetTF.m[1, ] >= minNtgts)]
  ### now which TFs are overexpressed in toi?
  topTOI.lm <- list()
  ### compare toi to all other tissues
  pheno.v <- rep(0, ncol(exp.m))
  pheno.v[which(tissue == toi)] <- 1
  lim.o <- LimmaFn(pheno.v, exp.m)
  topTOI.lm[[1]] <- lim.o$top[[1]]
  ### now compare toi to other tissues (in order to avoid confounding by immune or stromal cell infiltrates)
  if (!is.null(cft)) {
    ti <- 2;
    for(t in cft){
     sel.idx <- which(tissue %in% c(toi,t))
     tmp.v <- tissue[sel.idx]
     tmpPH.v <- rep(0,length(sel.idx))
     tmpPH.v[which(tmp.v==toi)] <- 1
     if(length(unique(tmp.v))<=1) stop(paste0("In comparison between ",toi," and ",t, ", all data are from ",unique(tmp.v),"."))
     lim.o <- LimmaFn(tmpPH.v,exp.m[,sel.idx])
     topTOI.lm[[ti]] <- lim.o$top[[1]]
     ti <- ti+1
    }
  }

  ### now find tissue-specific TFs
  toiTF.lv <- list()
  for (i in 1:length(topTOI.lm)) {
    statTF.m <- topTOI.lm[[i]][match(colnames(gtexNETf.m), rownames(topTOI.lm[[i]])), c(1, 3, 4, 5)]
    toiTF.idx <- intersect(which(statTF.m[, 4] < degth[i]), which(statTF.m[, 1] > lfcth[i]))
    toiTF.lv[[i]] <- rownames(statTF.m[toiTF.idx, ])
  }

  toiTF.v <- toiTF.lv[[1]]
  if (length(toiTF.lv) > 1) {
    for (i in 2:length(toiTF.lv)) {
      toiTF.v <- intersect(toiTF.v, toiTF.lv[[i]])
    }
  }
  map.idx <- match(toiTF.v, colnames(gtexNETf.m))

  if (length(map.idx) <= 1) stop(paste0("Only ", length(map.idx), "TFs in your network. Could try releasing the threshold 'lfcth'/'degth'/'minNtgts' to keep more."))

  #### tissue-specific regulatory network is:
  netTOI.m <- gtexNETf.m[, map.idx]

  distNet.m <- matrix(nrow = ncol(netTOI.m), ncol = 3)
  rownames(distNet.m) <- colnames(netTOI.m)
  colnames(distNet.m) <- c("nTGTS", "Act", "Rep")
  distNet.m[, 1] <- apply(abs(netTOI.m), 2, sum)
  for (c in 1:ncol(netTOI.m)) {
    act.idx <- which(netTOI.m[, c] == 1)
    inact.idx <- which(netTOI.m[, c] == -1)
    distNet.m[c, 2:3] <- c(length(act.idx), length(inact.idx))
  }

  return(list(netTOI = netTOI.m, sumnet = distNet.m, top = topTOI.lm))
}  ### end of function sepiraInfNet
